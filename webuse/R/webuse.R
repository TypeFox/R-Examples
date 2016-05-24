webuselist <- list(
    auto = "1978 Automobile Data",
    auto2 = "1978 Automobile Data",
    autornd = "Subset of 1978 Automobile Data",
    bplong = "fictional blood pressure data",
    bpwide = "fictional blood pressure data",
    cancer = "Patient Survival in Drug Trial",
    census = "1980 Census data by state",
    citytemp = "City Temperature Data",
    citytemp4 = "City Temperature Data",
    educ99gdp = "Education and GDP",
    gnp96 = "U.S. GNP, 1967-2002",
    lifeexp = "Life expectancy, 1998",
    network1 = "fictional network diagram data",
    network1a = "fictional network diagram data",
    nlsw88 = "U.S. National Longitudinal Study of Young Women (NLSW, 1988 extract)",
    nlswide1 = "U.S. National Longitudinal Study of Young Women (NLSW, 1988 extract)",
    pop2000 = "U.S. Census, 2000, extract",
    sandstone = "Subsea elevation of Lamont sandstone in an area of Ohio",
    sp500 = "S&P 500",
    surface = "NOAA Sea Surface Temperature",
    tsline1 = "simulated time-series data",
    tsline2 = "fictional data on calories consumed",
    uslifeexp = "U.S. life expectancy, 1900-1999",
    uslifeexp2 = "U.S. life expectancy, 1900-1940",
    voter = "1992 presidential voter data",
    xtline1 = "fictional data on calories consumed"
)
 
webuse <- function(data, version = 14, envir = parent.frame()) {
    d <- read_dta(paste0("http://www.stata-press.com/data/r",version,"/",data,".dta"))
    assign(data, d, envir = envir)
    invisible(d)
}
