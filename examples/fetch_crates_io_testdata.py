#!/usr/bin/env python
#
# Create some simple test-data by downloading crates.io download-stats. These provide some nice
# weekly seasonality and trends as new versions are released.

import requests

crate_name = "hashbrown"

d = requests.get(f"https://crates.io/api/v1/crates/{crate_name}/versions").json()
version_by_id = {}
for c in d["versions"]:
    version_by_id[c["id"]] = c["num"]

d = requests.get(f"https://crates.io/api/v1/crates/{crate_name}/downloads").json()
dl_by_id = {}
for dl in d["version_downloads"]:
    if not dl["version"] in dl_by_id:
        dl_by_id[dl["version"]] = []
    dl_by_id[dl["version"]].append((dl["date"], float(dl.get("downloads", 0))))
    
for k, v in dl_by_id.items():
    sorted(v, key=lambda x: x[0])
    print(f"\n\n{crate_name} version {version_by_id[k]}\ndownloads: {', '.join([str(x[1]) for x in v])}")
