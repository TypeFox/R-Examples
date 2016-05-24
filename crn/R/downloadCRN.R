downloadCRN <-function(url = CRN.DAILY.URL, directory = DAILY_DIR, 
                       years = seq(from = 2000,to = 2011, by = 1)) {
  
  if (!file.exists(directory)) dir.create(directory)
  
  for( thisYear in 1:length(years)){
    urlList <- getUrlsCRN(url = url,year = years[thisYear])
    for( thisFile in 1: nrow(urlList)){
     download.file(urlList[thisFile], 
                   file.path(directory,basename(urlList[thisFile]),fsep = .Platform$file.sep))
    }
  }
} 
  
    
    
   