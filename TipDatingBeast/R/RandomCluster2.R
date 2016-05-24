RandomCluster2 <- function(name, reps = 20, loadCluster = T) {
    
  options (warn = -1)
  
  setHook(packageEvent("mclust", "onLoad"), function(...) 
    sink(file(tempfile(), "w"), type = "message")) 
  
  setHook(packageEvent("mclust", "attach"), function(...) 
    sink(file = NULL, type = "message"), "append") 
  
  #library(mclust)
  
  inFileName <- paste0(name, ".xml")
  
  inFile <- readLines(inFileName)
  
  tips <- length(grep("<sequence id", inFile))
  
  matchTraitDateLine <- grep(pattern = "beast.evolution.tree.TraitSet", 
                             x = inFile)
  
  if (length(matchTraitDateLine) == 0) {stop(
    "No date info found, check the integrity of files and Beast version.")}  
  
  datePositions <- seq(matchTraitDateLine + 1, matchTraitDateLine + tips)
  
  dateLines <- inFile[datePositions]
  
  date <- unlist(strsplit(dateLines, "="))
  
  dateHap <- date[c(T, F)]
  
  dateHap <- dateHap[1: tips]
  
  dateValues <- date[c(F, T)]
  
  lastLine <- length(grep("<taxa", dateValues[tips]))
    
  dateValues <- na.omit(as.numeric (gsub("[^\\d]+", "",
                                         dateValues, perl = T)))
  
  if (loadCluster == F){
  
    min <- 1
  
    d_clust <- Mclust(as.matrix(dateValues), G = min:tips)
  
    n <- dim(d_clust$z)[2]
      
    if (n == 1) {stop(
      "Single cluster found, input clusters manually or use RandomDates.")}
    
    cat("Number of clusters found:", n, "\n")
  
    clusterOut <- cbind(dateHap, d_clust$classification)
  
    nameCluster <- paste0("clusters.", name,".csv")
  
    write.table(clusterOut, file = nameCluster,
                row.names = F, col.names = F, sep = ",")
    
  } else {
    
    nameCluster <- paste0("clusters.", name, ".csv")
    
    clusterOut <- read.table(nameCluster, sep = ",")}
  
    matchLinesCluster <- cbind(dateValues, clusterOut[,2])
  
    cluster <- split(matchLinesCluster[,1], matchLinesCluster[,2])
  
    n <- length(cluster)
  
    check <- lapply(cluster, length)
  
    clusterList <- matchLinesCluster[, 2]
  
  if (loadCluster == T){
    if (n == 1) {stop("Error, only One cluster found, check input.")}
    
    cat("Number of clusters loaded:", n, "\n")}
  
  dateValues <- matrix(0, tips)
  
  # Loop
    
  for(i in 1 :reps) {
    
    newFile <- inFile 
    
    for (j in 1 : tips){
      
      clusterL <- clusterList
            
      clusterId <- as.numeric(as.character(matchLinesCluster[j, 2]))
      
      if (clusterId != 0){
        
        clusterL <- clusterL[clusterL != clusterId]
        
        if ("0" %in% clusterL) {
          clusterL <- clusterL[clusterL != "0"]}
        
        if ("0" %in% clusterList && n == 2){
          clusterL <- clusterList[clusterList != "0"]}
        
        as.numeric(clusterL)
        
        pool <- matchLinesCluster[matchLinesCluster[, 2] %in% clusterL, 1]
        
        if (length(pool) == 1) {stop(
          "Single element in cluster, check input or use RandomDates.")}
        
        dateValues[j] <- sample(pool, size = 1)
        
      }  else {

      dateValues[j] <- matchLinesCluster[j, 1]}
      
    }


    newDate <- paste0("\t\t\t", dateHap, "=", dateValues)
    
    newFile[datePositions] <- paste0(newDate, ",")

    if(lastLine == 0){newFile[(datePositions[tips])] <- 
                        paste0 (newDate[tips])
    
    } else {
    
    newFile[(datePositions[tips])] <-
                        paste0 (newDate[tips], "\t\t\t\t<taxa id=",
                                  date[tips*2+1], "=", date[tips*2+2])}
    
    fileRep <- paste0(" fileName=\"Rep", i, ".")
    
    newFile <- gsub(" fileName=\"", fileRep, newFile)
    
    outName <- paste0(name, "_Rep", i) 
    
    out <- paste0(name, ".Rep", i, ".xml")
    
    cat (newFile, file = out, sep = "\n")
    
}
cat ("Replicates done:", i,"\n")
}
