#!/bin/bash

# Define the list of country codes
#all_countries=("AD" "AF" "AL" "AM" "AO" "AR" "AT" "AU" "AZ" "BA" "BD" "BE" "BF" "BG" "BI" "BJ" "BO" "BR" "BT" "BW" "BY" "BZ" "CA" "CD" "CF" "CG" "CH" "CI" "CK" "CL" "CM" "CN" "CO" "CR" "CV" "CY" "CZ" "DE" "DJ" "DK" "DZ" "EC" "EE" "EG" "ER" "ES" "ET" "FI" "FJ" "FM" "FO" "FR" "GA" "GB" "GE" "GF" "GH" "GL" "GN" "GQ" "GR" "GT" "GW" "GY" "HN" "HR" "HU" "ID" "IE" "IN" "IQ" "IR" "IS" "IT" "JO" "JP" "KE" "KG" "KH" "KI" "KP" "KR" "KZ" "LA" "LB" "LI" "LK" "LR" "LS" "LT" "LU" "LV" "LY" "MA" "MC" "MD" "ME" "MG" "MH" "MH" "MK" "ML" "MM" "MN" "MR" "MT" "MU" "MV" "MW" "MX" "MY" "MZ" "NA" "NC" "NE" "NG" "NI" "NL" "NO" "NP" "NR" "NU" "NZ" "PE" "PG" "PH" "PK" "PL" "PR" "PS" "PT" "PW" "PY" "QA" "RO" "RS" "RU" "RW" "SB" "SC" "SD" "SE" "SH" "SI" "SK" "SL" "SN" "SO" "SR" "SS" "ST" "SV" "SY" "SZ" "TD" "TG" "TH" "TJ" "TL" "TM" "TN" "TO" "TR" "TV" "TW" "TZ" "UA" "UG" "US" "UY" "UZ" "VE" "VI" "VN" "VU" "WS" "YE" "ZA" "ZM" "ZW")
# AQ = antarctica, takes too long to complete without significant usefulness

countries=("AD" "AF" "AL" "AM" "AO" "AR" "AT" "AU" "AZ" "BA" "BD" "BE" "BF" "BG" "BI" "BJ" "BO" "BR" "BT" "BW" "BY" "BZ" "CA" "CD" "CF" "CG" "CH" "CI" "CK" "CL" "CM" "CN" "CO" "CR" "CV" "CY" "CZ" "DE" "DJ" "DK" "DZ" "EC" "EE" "EG" "ER" "ES" "ET" "FI" "FJ" "FM" "FO" "FR" "GA" "GB" "GE" "GF" "GH" "GL" "GN" "GQ" "GR" "GT" "GW" "GY" "HN" "HR" "HU" "ID" "IE" "IN" "IQ" "IR" "IS" "IT" "JO" "JP" "KE" "KG" "KH" "KI" "KP" "KR" "KZ" "LA" "LB" "LI" "LK" "LR" "LS" "LT" "LU" "LV" "LY" "MA" "MC" "MD" "ME" "MG" "MH" "MH" "MK" "ML" "MM" "MN" "MR" "MT" "MU" "MV" "MW" "MX" "MY" "MZ" "NA" "NC" "NE" "NG" "NI" "NL" "NO" "NP" "NR" "NU" "NZ" "PE" "PG" "PH" "PK" "PL" "PR" "PS" "PT" "PW" "PY" "QA" "RO" "RS" "RU" "RW" "SB" "SC" "SD" "SE" "SH" "SI" "SK" "SL" "SN" "SO" "SR" "SS" "ST" "SV" "SY" "SZ" "TD" "TG" "TH" "TJ" "TL" "TM" "TN" "TO" "TR" "TV" "TW" "TZ" "UA" "UG" "US" "UY" "UZ" "VE" "VI" "VN" "VU" "WS" "YE" "ZA" "ZM" "ZW")
# countries=("DE")
for c in "${countries[@]}"
do
    fldr="./data/outputs/countries_h7/$c"

    # Check if folder exists
    if [ ! -d "$fldr" ]; then
        echo "Running script for country: $c"
        python cellxm/main.py locate "$c" --out-folder "$fldr"
    else
        echo "Country $c already exists. Skipping..."
    fi
done
