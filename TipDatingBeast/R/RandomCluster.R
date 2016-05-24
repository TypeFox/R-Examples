RandomCluster <- function(name, reps = 20, loadCluster = T) {
  
  options (warn = -1)
  
  setHook(packageEvent("mclust", "onLoad"), function(...) 
   sink(file(tempfile(), "w"), type = "message"))
  
   setHook(packageEvent("mclust", "attach"), function(...) 
    sink(file = NULL, type = "message"), "append")
  
  #library(mclust)
  
  inFileName <- paste0(name, ".xml")
  
  inFile <- readLines(inFileName)
  
  matchLines <- grep(pattern = "date value=", x = inFile, value = T)
  nameLinesRaw <- grep(pattern = "<taxon id=", x = inFile, value = T)
  
  nameLines = substring(nameLinesRaw, 14, nchar(nameLinesRaw))
  nameLines = substring(nameLines, 1, nchar(nameLines)-2)  
  
  if (length(matchLines) == 0) {stop(
    "No dates found, check the integrity of input files and Beast version.")}
  
  min <- 1
  
  tips <- length(matchLines)
  
  if (loadCluster == F){

    dates <- as.numeric(regmatches(matchLines[], 
                                   gregexpr('\\(?[0-9,.]+', matchLines[])))  
    
    d_clust <- Mclust(as.matrix(dates), G = min:tips)
    
    n <- dim(d_clust$z)[2]
    
    if (n == 1) {stop(
      "One cluster found, input clusters manually or use RandomDates.")}
    
    cat("Number of clusters found:", n, "\n")
    
    clusterOut <- cbind(nameLines, d_clust$classification)
    
    nameCluster <- paste0("clusters.", name,".csv")
    
    write.table(clusterOut, file = nameCluster, 
                row.names=F, col.names=F, sep=",")
    
  } else {
    
    nameCluster <- paste0("clusters.", name, ".csv")
    
    clusterOut <- read.table(nameCluster, sep=",")}

  matchLinesPosition <- grep(pattern = "date value=", x = inFile, value = F)
  
  matchFileName <- grep(pattern = "fileName", x = inFile, value = T)
  
  matchFileNamePosition <- grep(pattern = "fileName", x = inFile, value = F)
  
  matchLinesCluster <- cbind(matchLines, clusterOut[,2])
  
  cluster <- split(matchLinesCluster[,1], matchLinesCluster[,2])
  
  nCluster <- length(cluster)
  
  clusterList <- matchLinesCluster[, 2]
  
  if (loadCluster == T){
    if (nCluster == 1) {stop("Error, only One cluster found, check input.")}
        
    cat("Number of clusters loaded:", nCluster, "\n")}

  # Loops
  
  for (i in 1 : reps){
    
    newFile <- inFile
    
    fileRep <- paste0(" fileName=\"Rep", i, ".")
    
    newFile [matchFileNamePosition] <- gsub(" fileName=\"", 
                                            fileRep, matchFileName)
    
    for (j in 1 : tips){
      
      clusterL <- clusterList
      
      clusterId <- clusterOut[j,2]
      
      if (clusterId != 0){

        clusterL <- clusterL[clusterL != clusterId]
        
        if ("0" %in% clusterL) {
          clusterL <- clusterL[clusterL != "0"]}
        
        if ("0" %in% clusterList && nCluster == 2){
          clusterL <- clusterList[clusterList != "0"]}

        as.numeric(clusterL)
        
        pool <- matchLinesCluster[matchLinesCluster[,2] %in% clusterL,1]

        if (length(pool) == 1) {stop(
          "Single element in cluster, Check input or use RandomDates.")}
              
      newFile[matchLinesPosition[j]] <- sample(pool, size=1)

      }
  }
    
    outName <- paste0(name, "_Rep", i) 
    
    out <- paste0(name, ".Rep", i, ".xml")
    
    cat (newFile, file = out, sep = "\n")
  }
  
  cat ("Replicates done:", i,"\n")
}

