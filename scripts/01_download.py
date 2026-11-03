#!/usr/bin/env python3
"""01_download.py v2 - resilient FinnGen R11 download, multi-connection edition.

Priority: aria2c (16 connections per file, automatic resume via control file)
          -> fallback: resilient curl loop (single connection, -C - resume).
Both honor http_proxy / https_proxy environment variables.

Idempotent: files whose size matches config/endpoints.tsv are skipped.
Partial files from ANY previous attempt (curl or aria2) are resumed,
never restarted from zero.
Records SHA256 into results/metadata/input_manifest.json.
"""
import os, sys, subprocess, json, time, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

BASE_URL = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_{ep}.gz"


def have_aria2():
    return shutil.which("aria2c") is not None


def fetch_aria2(url, out, expect_bytes):
    """One aria2c invocation. -c resumes from existing partial (curl or aria2).
    aria2c retries internally forever (--max-tries=0); we only get control back
    on success or a hard error, in which case the caller re-invokes (resume)."""
    r = subprocess.run([
        "aria2c", "--continue=true",
        "-x", "16", "-s", "16", "-k", "4M",
        "--max-connection-per-server=16",
        "--max-tries=0", "--retry-wait=5",
        "--timeout=30", "--connect-timeout=20",
        "--file-allocation=none", "--auto-file-renaming=false",
        "--console-log-level=warn", "--summary-interval=30",
        "-d", os.path.dirname(out), "-o", os.path.basename(out), url])
    return (r.returncode == 0 and os.path.exists(out)
            and os.path.getsize(out) == expect_bytes)


def fetch_curl(url, out, expect_bytes):
    attempt = 0
    while True:
        if os.path.exists(out) and os.path.getsize(out) == expect_bytes:
            return True
        attempt += 1
        have = os.path.getsize(out) if os.path.exists(out) else 0
        common.log(os.path.basename(out),
                   f"curl attempt {attempt}, resuming from {have/1e6:.0f} MB ...")
        r = subprocess.run(["curl", "-L", "-C", "-", "--retry", "2", "--retry-delay", "5",
                            "--speed-limit", "30720", "--speed-time", "60",
                            "-o", out, url])
        if (r.returncode == 0 and os.path.exists(out)
                and os.path.getsize(out) == expect_bytes):
            return True
        common.log(os.path.basename(out),
                   f"curl code {r.returncode}; retrying in 15 s (progress is kept)")
        time.sleep(15)


def fetch(endpoint, expect_bytes):
    out = common.p("data", "raw", f"finngen_R11_{endpoint}.gz")
    url = BASE_URL.format(ep=endpoint)
    if os.path.exists(out) and os.path.getsize(out) == expect_bytes:
        common.log(endpoint, "already complete", f"({expect_bytes/1e9:.2f} GB)")
        return out
    if have_aria2():
        round_ = 0
        while True:
            round_ += 1
            have = os.path.getsize(out) if os.path.exists(out) else 0
            common.log(endpoint, f"aria2 round {round_}, have {have/1e6:.0f} MB / "
                                 f"{expect_bytes/1e9:.2f} GB ...")
            if fetch_aria2(url, out, expect_bytes):
                return out
            common.log(endpoint, "aria2 exited early; resuming in 10 s")
            time.sleep(10)
    else:
        common.log(endpoint, "aria2c not found; single-connection curl fallback "
                             "(install aria2 for multi-connection download)")
        if fetch_curl(url, out, expect_bytes):
            return out


def main():
    man_path = common.p("results", "metadata", "input_manifest.json")
    man = json.load(open(man_path)) if os.path.exists(man_path) else {}
    for ep in common.load_endpoints():
        fp = fetch(ep["endpoint"], ep["remote_bytes"])
        key = os.path.basename(fp)
        if man.get(key, {}).get("bytes") != ep["remote_bytes"]:
            common.log(ep["endpoint"], "sha256 ...")
            man[key] = {"bytes": ep["remote_bytes"], "sha256": common.sha256(fp),
                        "source": BASE_URL.format(ep=ep["endpoint"])}
            common.write_json(man, "results", "metadata", "input_manifest.json")
    common.log("all 11 raw files present and size-verified")


if __name__ == "__main__":
    main()
