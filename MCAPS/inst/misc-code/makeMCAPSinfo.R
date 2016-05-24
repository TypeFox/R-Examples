library(stashR)

db <- new("localDB", dir = "MCAPSinfo", name = "MCAPSinfo")

sites <- dget("http://www.biostat.jhsph.edu/MCAPS/siteList.R")

APWdata <- lapply(sites, function(sitename) {
    cat(sitename, "\n")       
    read.csv(paste("../APWdata/", sitename, ".csv", sep = ""),
             colClasses = c("Date", rep("numeric", 5)), nrow = 1500)
})
names(APWdata) <- sites

cl <- c("factor", "factor", "factor", "numeric", "numeric", "integer", rep("factor", 4))
        
estimates.subset <- read.csv("http://www.biostat.jhsph.edu/MCAPS/estimates-subset.csv",
                             colClasses = cl)

estimates.full <- read.csv("http://www.biostat.jhsph.edu/MCAPS/estimates-full.csv",
                           colClasses = cl)


dbInsert(db, "APWdata", APWdata)
dbInsert(db, "siteList", sites)
dbInsert(db, "estimates.subset", estimates.subset)
dbInsert(db, "estimates.full", estimates.full)
