#!/usr/bin/env python3
"""Shared paths/config helpers. All paths resolve from the package root
(the directory containing this script's parent), never hard-coded."""
import os, sys, csv, json, hashlib, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def p(*parts):
    return os.path.join(ROOT, *parts)

def log(*a):
    m = " ".join(str(x) for x in a)
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {m}", flush=True)

def load_endpoints():
    """config/endpoints.tsv -> list of dicts."""
    out = []
    with open(p("config", "endpoints.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["N_total"] = int(row["N_total"])
            row["remote_bytes"] = int(row["remote_bytes"])
            out.append(row)
    if len(out) != 11:
        sys.exit(f"FATAL: endpoints.tsv must have 11 rows, got {len(out)}")
    return out

def sha256(path, _buf=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(_buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def write_json(obj, *parts):
    path = p(*parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1)
    os.replace(tmp, path)
