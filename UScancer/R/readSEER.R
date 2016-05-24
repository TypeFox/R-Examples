# Read SEER data (assume downloaded locally)
#  file = vector of files to read from, i.e. files = c("usCancer/SEER_1973_2010_TEXTDATA/incidence/yr1992_2010.sj_la_rg_ak/RESPIR.TXT","usCancer/SEER_1973_2010_TEXTDATA/incidence/yr2000_2010.ca_ky_lo_nj_ga/RESPIR.TXT")

readSEER <- function(file,year,state,cancer,site,fips) {
  
  if(missing(file)) {
    print("Please specify SEER TXT file(s) to read from. Be sure the specified file includes the cancer and state of interest.")
  }
      
  # the file is in a weird format, extract data from these indicies
  startCols = c(sex=24, age=25, year=39, siteRecode=208, fips=246)
  endCols   = c(sex=25, age=28, year=43, siteRecode=212, fips=251)
  
  widths = c(min(startCols)-1, diff(sort(c(startCols,endCols))))
  theNames = c("junk1","sex","junk2","age","junk3","year","junk4","site","junk5","fips")
  
  allCancer = NULL
  for(D in file) {
    usCancer = read.fwf(D,widths = widths,
                        colClasses=c("character","integer","character","integer","character","integer","character","character","character","character"))
    names(usCancer) = theNames
    usCancer = usCancer[,c("year","site","sex","age","fips")]
    usCancer$site <- gsub(" ","",usCancer$site, fixed=TRUE) #remove some whitespacing
    usCancer = cbind(usCancer, fips[usCancer$fips,c("state","county")])
    rownames(usCancer) <- NULL #remove rownames that came from fipslookup
    allCancer = rbind(allCancer, usCancer)
  }
  usCancer = allCancer
  
  # remove cases that could not be mapped to a county/state
  usCancer = usCancer[!is.na(usCancer$state),] 
  
  # select only those rows that match 'years'
  if (missing(year)) {
    print("Note: year not supplied, including all years")
  } else {
    usCancer = usCancer[(usCancer$year %in% year),]
  }
  
  # select only those rows that match 'state'
  if (missing(state)) {
    print("Note: state not supplied, including all states")
  } else {
    usCancer = usCancer[(tolower(usCancer$state) %in% tolower(state)),]
  }
  
  # select only those rows that match the regular expression(s) in 'site'
  if (missing(cancer) & missing(site)) {
    print("Note: cancer/site not supplied, including all sites")
  } else {
    if (missing(site)) {
      site <- unlist(lapply(cancer,siteLookup))
    }
    siteind <- numeric(0)
    for (i in 1:length(site)) {
      siteind <- c(siteind, grep(site[i], usCancer$site))
    }
    usCancer = usCancer[siteind,]
  }
  
  rownames(usCancer) <- NULL #remove rownames that came from fipslookup
  
  return(usCancer)
}

# Test code
## files = c("../../../CCO/usCancer/SEER_1973_2010_TEXTDATA/incidence/yr1992_2010.sj_la_rg_ak/RESPIR.TXT","../../../CCO/usCancer/SEER_1973_2010_TEXTDATA/incidence/yr2000_2010.ca_ky_lo_nj_ga/RESPIR.TXT","../../../CCO/usCancer/SEER_1973_2010_TEXTDATA/incidence/yr1973_2010.seer9/RESPIR.TXT")
## seer <- readSEER(files)
## 
## dim(seer)



