# 2010 FIPS Codes for Counties and County Equivalent Entities
# Details from census.gov:
# http://www.census.gov/geo/reference/codes/cou.html
library(dplyr)
col_names <- c("state", "state_fips", "county_fips", "county_name", "fips_class")
col_classes <- rep.int("character", times=length(col_names))
county_fips <- read.csv("http://www.census.gov/geo/reference/codes/files/national_county.txt",
                        header=T,
                        col.names=col_names,
                        colClasses=col_classes)

# County-level data from the Census Bureau
# Not included for Puerto Rico
census_data <- read.csv('http://www.census.gov/popest/data/counties/totals/2013/files/CO-EST2013-Alldata.csv',
                        header=TRUE)
# The variable SUMLEV indicates "Geographic summary level"
# 040 = State and/or Statistical Equivalent
# 050 = County and /or Statistical Equivalent
# We retain only the county data
census_data <- subset(census_data, SUMLEV == 50)

# For now, we are focused on county-level population from the 2010 census
census_data <- subset(census_data, select=c("STATE", "COUNTY", "CENSUS2010POP"))
colnames(census_data) <- c("state_fips", "county_fips", "population")

# Recode state_fips and county_fips to join against census_data
# Pads FIPS codes with the appropriate amount of 0's
census_data$state_fips <- sprintf("%02d", census_data$state_fips)
census_data$county_fips <- sprintf("%03d", census_data$county_fips)

# Merges the FIPS and Census data
# The left join induces NA for demographic/population data when not FIPS code
# not included (e.g., Puerto Rico).
counties <- left_join(county_fips, census_data)


# Next, we add the U.S. Census-defined geographic areas to the counties
# data.frame

# Combined Statistical Areas (CSAs)
csa <- read.csv("http://www.census.gov/popest/data/metro/totals/2013/files/CSA-EST2013-alldata.csv",
                header=TRUE)
combined_areas <- subset(csa,
                         LSAD == "Combined Statistical Area",
                         select=c(CSA, NAME))
colnames(combined_areas) <- c("CSA", "name")
rownames(combined_areas) <- NULL

# Core-based Statistical Areas (CBSAs)
cbsa <- read.csv("http://www.census.gov/popest/data/metro/totals/2013/files/CBSA-EST2013-alldata.csv",
                 header=TRUE)

# 1-to-many mapping of CBSA to CSA
CSA_to_CBSA <- subset(csa, !is.na(CBSA), select=c(CSA,CBSA))
CSA_to_CBSA <- CSA_to_CBSA[!duplicated(CSA_to_CBSA), ]

# CBSAs are the union of Metropolitan and Micropolitan Statistical Areas (MSAs)
metropolitan_areas <- subset(cbsa,
                             LSAD == "Metropolitan Statistical Area",
                             select=c(CBSA, NAME))
micropolitan_areas <- subset(cbsa,
                             LSAD == "Micropolitan Statistical Area",
                             select=c(CBSA, NAME))

# Combine metropolitan_areas and micropolitan_areas into single data.frame
corebased_areas <- rbind(
    data.frame(metropolitan_areas, type="Metropolitan"),
    data.frame(micropolitan_areas, type="Micropolitan")
)
corebased_areas <- left_join(corebased_areas, CSA_to_CBSA)
colnames(corebased_areas) <- c("CBSA", "name", "type", "CSA")
corebased_areas <- corebased_areas[c("CBSA", "CSA",  "name", "type")]
rownames(corebased_areas) <- NULL

# The STCOU variable is a concatenation of the state and county FIPS codes
# The 0 from single-digit state FIPS codes have been truncated.
# We pad the FIPS codes before extracting both the state and county FIPS codes.
stat_areas <- subset(cbsa,
                     LSAD == "County or equivalent",
                     select=c("CBSA", "STCOU"))
stat_areas$STCOU <- sprintf("%05d", stat_areas$STCOU)
stat_areas$state_fips <- substr(stat_areas$STCOU, start=1, stop=2)
stat_areas$county_fips <- substr(stat_areas$STCOU, start=3, stop=5)
stat_areas <- stat_areas[c("CBSA", "state_fips", "county_fips")]
stat_areas <- left_join(stat_areas, CSA_to_CBSA)

counties <- left_join(counties, stat_areas)


# Post-processing of counties to ensure things are juuuust right
counties <- counties[c("county_name", "state", "state_fips", "county_fips",
                       "fips_class", "CSA", "CBSA", "population")]
counties$state <- factor(counties$state)
counties$state_fips <- factor(counties$state_fips)
counties$county_fips <- factor(counties$county_fips)
counties$fips_class <- factor(counties$fips_class)
counties$CSA <- factor(counties$CSA)
counties$CBSA <- factor(counties$CBSA)

save(counties, file="../../data/counties.RData", compress='xz')

# Updates the encoding where Unicode is in a city/county name
# See: http://stackoverflow.com/a/19480595/234233
Encoding(levels(combined_areas$name)) <- "latin1"
levels(combined_areas$name) <- iconv(
      levels(combined_areas$name),
      "latin1",
      "UTF-8"
    )

Encoding(levels(corebased_areas$name)) <- "latin1"
levels(corebased_areas$name) <- iconv(
      levels(corebased_areas$name),
      "latin1",
      "UTF-8"
    )


# Saves the lookup tables of CSAs and CBSAs
save(combined_areas, file="../../data/combined_areas.RData", compress='xz')
save(corebased_areas, file="../../data/corebased_areas.RData", compress='xz')
