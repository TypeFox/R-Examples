scrapeToCsv <- function(Stations, get=seq(from=1,to=100000) , directory = "EnvCanada"){
   if (!file.exists(directory)) dir.create(directory)
   
  All <- 1:nrow(Stations)
  common <- intersect(All,get)
  if (length(common) == 0) stop(" get must be a subset of the total number of stations")
  
  Fname <- Stations$Name
  Fname <- gsub(" ","",Fname)
  Fname <- gsub(".","",Fname,fixed = TRUE)
  Fname <- gsub("(","",Fname,fixed = TRUE)
  Fname <- gsub(")","",Fname,fixed = TRUE)
  Fname <- gsub("/","",Fname,fixed = TRUE)
  Fname <- gsub("'","",Fname,fixed = TRUE)
  
  for( station in common){
  
    stationurl  <- paste(BASE.URL,Stations$WebId[station],YEAR.URL,Stations$Year[station], FORMAT.URL, sep = "")
    Destination <- paste(Stations$Id[station],"_",Fname[station],".Env",".csv", sep = "")
    Destination <- file.path(directory,Destination,fsep = .Platform$file.sep)   
    download.file(stationurl,Destination,mode = "wb")
    print(station)
   
  }
     
}