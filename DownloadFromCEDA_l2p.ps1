# Base URL
$baseUrl = "https://data.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness/L3C/cryosat2/v3.0/NH"

# Local output folder
$outputDir = "D:\RRDP_validation\RRDPp_v2\RawData\Arctic\CryoSat_L3C"
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

# Get subdirectories (years)
$yearHtml = Invoke-WebRequest "$baseUrl/"
$years = ($yearHtml.Links | Where-Object { $_.href -match "^\d{4}/$" }).href

foreach ($year in $years) {
    $yearUrl = "$baseUrl/$year"
    $yearPath = Join-Path $outputDir $year.TrimEnd('/')
    New-Item -ItemType Directory -Force -Path $yearPath | Out-Null

    # Get months
    $monthHtml = Invoke-WebRequest "$yearUrl"
    $months = ($monthHtml.Links | Where-Object { $_.href -match "^\d{2}/$" }).href

    foreach ($month in $months) {
        $monthUrl = "$yearUrl$month"
        $monthPath = Join-Path $yearPath $month.TrimEnd('/')
        New-Item -ItemType Directory -Force -Path $monthPath | Out-Null

        # Get .nc files
        $fileHtml = Invoke-WebRequest "$monthUrl"
        $ncFiles = ($fileHtml.Links | Where-Object { $_.href -match "\.nc$" }).href

        foreach ($file in $ncFiles) {
            $fileUrl = "$monthUrl$file"
            $localFile = Join-Path $monthPath $file
            Write-Host "Downloading $fileUrl → $localFile"
            curl.exe -L -o $localFile $fileUrl
        }
    }
}