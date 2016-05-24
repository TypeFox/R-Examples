##                        getHF                       ##
##      This code is part of the rusda package        ##
##      F.-S. Krah 2015 (last update: 2015-07-11)     ##  

getHF <- function(x, process, spec_type = c("plant","fungus")){
  if(process == TRUE){ message(paste(as.character(x[1]), as.character(x[2])))}
  if(spec_type == "fungus")
  {
    url <- paste("http://nt.ars-grin.gov/fungaldatabases/new_allViewGenBank.cfm?thisName=",
                 as.character(x[1]), "%20",
                 as.character(x[2]),
                 "&organismtype=Fungus&CFID=28606i&",
                 "CFTOKEN=5f05e06d6d0caa92-993D3B37-A643-36F0-49FFCF751D1BE3B1", sep="")
  }
  if(spec_type == "plant")
  {
    url <- paste("http://nt.ars-grin.gov/fungaldatabases/new_allViewGenBank.cfm?thisName=",
                 as.character(x[1]), "%20",
                 as.character(x[2]),
                 "&organismtype=Host&CFID=2980587&CFTOKEN=9092a7639c080a51-",
                 "32974067-A5AE-EA66-FE9A1BDF640011B5", sep="")
  }
  pars <- rawToChar(GET(url)$content)
  pars <- htmlTreeParse(pars, useInternalNodes = TRUE)
  pars <- xpathApply(pars, "//p", xmlValue)
  pars <- unlist(pars)
}
