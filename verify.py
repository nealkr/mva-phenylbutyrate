#!/usr/bin/env python3
"""Verify release-file hashes and inspect text and Office XML."""
import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def excluded():
    provenance = json.loads((ROOT / "provenance.json").read_text())
    fetched = {
        item["path"]
        for item in provenance["inputs"] + provenance["derived_inputs"]
        if item["distribution"] == "fetch"
    }
    return fetched | set(provenance.get("generated_fetch_artifacts", []))


def inventory():
    omitted = excluded()
    rows = []
    for path in sorted(ROOT.rglob("*")):
        name = path.relative_to(ROOT).as_posix()
        if name.startswith(("results/", ".venv/", ".git/")) or "__pycache__" in path.parts:
            continue
        if path.is_symlink():
            raise ValueError("Symlinks are not included in the release")
        if not path.is_file() or name == "manifest.json" or name in omitted:
            continue
        if name.endswith((".pyc", ".part", ".zip")):
            continue
        content = path.read_bytes()
        rows.append({
            "path": name,
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        })
    return rows


def text_parts(path):
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for name in archive.namelist():
                if name.endswith((".xml", ".rels")):
                    yield archive.read(name).decode("utf-8", errors="replace")
    elif path.suffix.lower() not in {".xls", ".png", ".pdf"}:
        yield path.read_text(encoding="utf-8", errors="replace")


def privacy(rows):
    patterns = [
        r"/(?:home|Users)/[^/\s]+/",
        r"/mnt/[a-z]/Users/",
        r"[A-Z]:[/\\]Users[/\\]",
        r"(?:wsl\$|wsl\.localhost)[/\\]",
        r"-----BEGIN [A-Z ]*PRIVATE KEY-----",
        r"\b(?:hf_|ghp_|github_pat_|sk-)[A-Za-z0-9_\-]{20,}",
    ]
    found = []
    for row in rows:
        if any(
            re.search(pattern, text, re.I)
            for text in text_parts(ROOT / row["path"])
            for pattern in patterns
        ):
            found.append(row["path"])
    if found:
        raise ValueError("Review local paths or credential-like text in: " + ", ".join(found))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seal", action="store_true")
    arguments = parser.parse_args()
    rows = inventory()
    privacy(rows)
    manifest = ROOT / "manifest.json"
    if arguments.seal:
        manifest.write_text(json.dumps({"schema_version": 1, "files": rows}, indent=2) + "\n")
    elif json.loads(manifest.read_text())["files"] != rows:
        raise ValueError("Release differs from manifest.json")
    print(json.dumps({
        "status": "verified",
        "files": len(rows),
        "bytes": sum(row["bytes"] for row in rows),
        "path_and_secret_scan": "passed",
    }))


if __name__ == "__main__":
    main()
