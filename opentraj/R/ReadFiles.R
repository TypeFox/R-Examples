ReadFiles <- 
  function( working_dir, ID, dates, tz) 
{
  combine.file.name <- paste("Rcombined_", ID, ".txt", sep="")
  dump.file.name <- paste("tdump_", ID, "_", "*", sep="")
  
  # find tdump files
  files <- list.files(path = working_dir, pattern = paste("tdump_", ID, sep=""))
  
  output <- file(combine.file.name, 'w')
  
  if(length(dates) != length(files)){
    print(length(dates))
    print(length(files))
    # cleanWD(hy.path, ID)
    stop("Please, make sure that all required meteorological files are available. Also, delete 
         all folders that starts with \"process_\".")
  }
  
  # read through them all, ignoring 1st 7 lines
  for (i in files){
    input <- readLines(i)
    input <- input[-c(1:7)] # delete header
    writeLines(input, output)
  }
  close(output)
  
  # read the combined txt file
  
  traj <- read.table(file.path(working_dir, combine.file.name), 
                     header = FALSE)
  
  traj <- subset(traj, select = -c(2, 7, 8))
  
  traj <- rename(traj, c(V1 = "receptor", V3 = "year", V4 = "month", V5 = "day",
                         V6 = "hour", V9 = "hour.inc", V10 = "lat", V11 = "lon",
                         V12 = "height", V13 = "pressure"))
  
  # hysplit uses 2-digit years ...
  year <- traj$year[1]
  
  if (year < 50) {
    traj$year <- traj$year + 2000 
  } else { 
    traj$year <- traj$year + 1900
  }
  
  # Setup the Canada, ON timezone
  traj$date2 <- with(traj, ISOdatetime(year, month, day, hour, min = 0, 
                                       sec = 0, tz = tz)) 
  
  # arrival time
  traj$date <- traj$date2 - 3600 * traj$hour.inc
  traj
}