import os
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, unquote, parse_qs, urlparse


# Base URL for the UNC FounderProbs page
base_url = "https://csbio.unc.edu/CCstatus/index.py?run=FounderProbs"
root = "https://csbio.unc.edu/CCstatus/"

# Output folder
os.makedirs("/Users/pathumnawarathna/Documents/Research/Simulation/keele/Probability_matrices_from_UNC/CC0_downloads", exist_ok=True)

# Load the webpage
response = requests.get(base_url)
soup = BeautifulSoup(response.text, "html.parser")

# Find all buttons on the page
buttons = soup.find_all("input", attrs={"type": "button"})

download_count = 0

suffixes = ["Tau", "TauUnc", "GeniUnc", "Unc"]
patterns = []

for i in range(1, 85):
    prefix = f"CC{i:03d}"
    for suffix in suffixes:
        patterns.append(f"{prefix}/{suffix}{"b38V01"}")



for btn in buttons:
    onclick = btn.get("onclick")
    if not onclick or "location.href=" not in onclick:
        continue

    # Extract URL from JavaScript onclick string
    start = onclick.find("'") + 1
    end = onclick.rfind("'")
    relative_url = onclick[start:end]
    full_url = urljoin(root, relative_url)

    # Parse the 'download=' value from the URL
    parsed = parse_qs(urlparse(full_url).query)
    download_val = parsed.get("download", [""])[0]


    #print(download_val)
    
    if download_val in patterns:
        strain_name = download_val.split("/")[0]  # e.g., CC001
        build = parsed.get("v", [""])[0]  # 37 or 38
        version = "b" + build if build else "unknown"
        filename = f"{strain_name}_{version}.csv"
        filepath = os.path.join("/Users/pathumnawarathna/Documents/Research/Simulation/keele/Probability_matrices_from_UNC/CC0_downloads", filename)

        print(f"Downloading {strain_name} ({version}) from {full_url}...")
        r = requests.get(full_url)
        if r.status_code == 200:
            with open(filepath, 'wb') as f:
                f.write(r.content)
            download_count += 1
        else:
            print(f"Failed to download: {full_url}")

print(f"Done. {download_count} files downloaded to 'CC0_downloads' folder.")

