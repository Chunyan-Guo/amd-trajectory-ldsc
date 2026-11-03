#!/usr/bin/env python3
"""00_check_environment.py - environment + input sanity check.
Writes results/metadata/environment_check.json. Hard-aborts on sentinel mismatch."""
import os, sys, platform, gzip, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

def main():
    info = {"python": sys.version, "platform": platform.platform(),
            "cwd": common.ROOT, "packages": {}, "checks": {}}
    for mod in ("pandas", "numpy", "scipy"):
        try:
            m = __import__(mod)
            info["packages"][mod] = m.__version__
        except ImportError:
            info["packages"][mod] = "MISSING"
    free = shutil.disk_usage(common.ROOT).free
    info["disk_free_GB"] = round(free / 1e9, 1)

    # sentinel 1: LD score files present, header exact
    ld_dir = common.p("resources", "ldscores")
    lds = sorted(f for f in os.listdir(ld_dir) if f.endswith(".l2.ldscore.gz")) if os.path.isdir(ld_dir) else []
    info["checks"]["ld_files"] = len(lds)
    if len(lds) != 22:
        sys.exit("FATAL: resources/ldscores must contain 22 ldN.l2.ldscore.gz files "
                 f"(found {len(lds)}). Did you unzip BOTH LD resource zips into the same folder?")
    with gzip.open(os.path.join(ld_dir, "ld1.l2.ldscore.gz"), "rt") as f:
        hdr = f.readline().rstrip("\n")
    if hdr != "CHR\tSNP\tBP\tL2":
        sys.exit(f"FATAL: BUILD_MISMATCH sentinel - ld1 header is {hdr!r}, expected 'CHR\\tSNP\\tBP\\tL2'")
    info["checks"]["ld_header"] = "OK"

    # sentinel 2: snp list count
    snpgz = common.p("resources", "ldsnps.txt.gz")
    n = 0
    with gzip.open(snpgz, "rt") as f:
        for _ in f:
            n += 1
    if n != 9_997_231:
        sys.exit(f"FATAL: SNP list sentinel - expected 9,997,231 lines, got {n}")
    info["checks"]["ldsnps_lines"] = n

    # raw data presence (not fatal here; 01_download fills these)
    raw = []
    for ep in common.load_endpoints():
        fp = common.p("data", "raw", f"finngen_R11_{ep['endpoint']}.gz")
        ok = os.path.exists(fp) and os.path.getsize(fp) == ep["remote_bytes"]
        raw.append({"endpoint": ep["endpoint"], "ready": bool(ok)})
    info["checks"]["raw_ready"] = sum(r["ready"] for r in raw)
    info["checks"]["raw_detail"] = raw

    common.write_json(info, "results", "metadata", "environment_check.json")
    common.log("environment check written:",
               "pkgs", info["packages"],
               "| disk free", info["disk_free_GB"], "GB",
               "| raw ready", info["checks"]["raw_ready"], "/ 11")
    if any(v == "MISSING" for v in info["packages"].values()):
        sys.exit("FATAL: missing python packages - run bash install_ubuntu.sh first")
    if info["disk_free_GB"] < 25:
        common.log("WARNING: less than 25 GB free disk; downloads need ~9 GB, outputs ~3 GB")

if __name__ == "__main__":
    main()
