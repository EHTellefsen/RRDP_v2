# Base URL for L3C CryoSat-2 data
$baseUrl = "https://dap.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness/L3C/cryosat2/v3.0/SH"

# Local directory to save downloaded files
$outputDir = "D:\RRDP_validation\RRDPp_v2\RawData\Antarctic\CryoSat_L3C"
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

# Retrieve the list of year directories
$yearHtml = Invoke-WebRequest "$baseUrl/"
$years = ($yearHtml.Links | Where-Object { $_.href -match "^\d{4}/$" }).href

foreach ($year in $years) {
    $yearUrl = "$baseUrl/$year"
    $yearPath = Join-Path $outputDir $year.TrimEnd('/')
    New-Item -ItemType Directory -Force -Path $yearPath | Out-Null

    # Retrieve the list of .nc files for the year
    $fileHtml = Invoke-WebRequest "$yearUrl"
    $ncFiles = ($fileHtml.Links | Where-Object { $_.href -match "\.nc$" }).href

    foreach ($file in $ncFiles) {
        $fileUrl = "$yearUrl$file"
        $localFile = Join-Path $yearPath $file
        Write-Host "Downloading $fileUrl → $localFile"
        curl.exe -L -o $localFile $fileUrl
    }
}