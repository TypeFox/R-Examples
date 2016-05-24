## ---- eval=FALSE, install------------------------------------------------
#  library(abbyyR)

## ---- eval=FALSE, setapp-------------------------------------------------
#  # setapp(c("factbook", "7YVBc8E6xMricoTwp0mF0aH"))

## ---- eval=FALSE, comments_listtask--------------------------------------
#  "
#  all_tasks <- listTasks()
#  for (i in 1:nrow(all_tasks)) deleteTask(all_tasks$id[i])
#  "

## ---- eval=FALSE, iterate------------------------------------------------
#  # Set path to directory with all the images
#  path_to_img_dir <- paste0(path.package("abbyyR"),"/inst/extdata/wisc_ads/")
#  total_files <- length(dir(path_to_img_dir))
#  
#  # Iterate through the files and submit all the images
#  
#  # Monitor progress via progress bar package
#  library(progress)
#  pb <- progress_bar$new(format = "  downloading [:bar] :percent\n",
#  					    total = total_files,
#  					    clear = FALSE, width= 60)
#  
#  # Abbyy Fine API doesn't keep the file name so we have to keep track of it locally
#  tracker <- data.frame(filename=NA, taskid=NA)
#  
#  # Loop
#  j <- 1
#  
#  for (i in dir(path_to_img_dir)){
#  	
#  	# Assuming only 1 dot in the file name
#  	tracker[j,] <- c(unlist(strsplit(basename(i), "[.]"))[1], submitImage(file_path=paste0(path_to_img_dir, i))$id)
#  	j <- j + 1
#  
#  	# Prg. bar
#  	pb$tick()
#  	Sys.sleep(1/100)
#  }
#  

## ---- eval=FALSE, process------------------------------------------------
#  
#  for (i in 1:nrow(tracker)) processDocument(tracker$taskid[i])
#  

## ---- eval=FALSE, checktasks---------------------------------------------
#  "
#  i <- 1
#  
#  while(i < total_files){
#  	i <- nrow(listFinishedTasks())
#  	if (i == total_files){
#  		print("All Done!")
#  		break;
#  	}
#  
#  	Sys.sleep(5)
#  	}
#  "

## ---- eval=FALSE, download-----------------------------------------------
#  setwd(paste0(path.package("abbyyR"),"/inst/extdata/wisc_out/"))
#  
#  finishedlist <- listFinishedTasks()
#  results      <- merge(tracker, finishedlist, by.x="taskid", by.y="id")
#  
#  library(curl)
#  
#  for(i in 1:nrow(results)){
#  	curl_download(results$resultUrl[i], destfile=results$filename[i])
#  }
#  

