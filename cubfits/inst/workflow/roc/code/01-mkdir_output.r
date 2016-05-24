### This script makes directories for outputs. 

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

### Set environment and data.
source("00-set_env.r")

### Check or create case directories.
for(i.case in case.names){
  dn <- paste(prefix$output, i.case, sep = "")
  if(!file.exists(dn)){
    cat(dn, " is made.\n", sep = "")
    dir.create(dn)
  }
}
