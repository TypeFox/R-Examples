scanCRAN <- function(folder,cutoffvalue=NULL,minArea=NULL,cores=1, gray=TRUE, stand=c(0,0), fact=0.25){
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("EBImage needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  tes <- dir(folder) # IS A PICTURE OR FOLDER??
  if(length(tes) == 0){# is a picture
    listp2 <- folder # pictures to analize
    #setwd(getwd()) # setwd
    folder <- getwd() # set a folder
  }else{# is a folder
    setwd(folder)
    listp2 <- dir(folder, '*.JPG')
  }
  
  ###################################
  ## is important for pictures to be taken with references at the top since bwlabel() function
  ## counts the reference circles in the right direction --> --> --> i.e.
  ##
  ##   O  O  O  O  O  O
  ##  .  .   .   . .   .
  ##   .   .      .  .
  ##       .   .     . 
  ##    .     .  .  .
  ##
  ##   O  O  O  O  O  O
  ##
  ## where "O" are reference circles and "." are fruits
  ###################################
  if(length(which(dir() %in%  "Pictures_analyzed")) > 0){
    i=1
    while(length(which(dir() %in%  paste("Pictures_analyzed",i, sep=""))) > 0){
      i=i+1
    }
    sto <- paste(getwd(), paste("Pictures_analyzed",i, sep=""),sep="/")
    dir.create(paste(getwd(), paste("Pictures_analyzed",i, sep=""),sep="/"))
    
  }else{sto <- paste(getwd(), "Pictures_analyzed",sep="/");  dir.create(paste(getwd(), "Pictures_analyzed",sep="/"))}
  
  
  x <- numeric()# creating a null vector or matrix
  #######################################################################################################
  scan <- function(folder,picture, cutoffvalue=NULL, minArea=NULL, grays=TRUE, stands=c(0,0), fact=0.25, stored=getwd()){
   
    #source('C:/Users/zalapalab/Desktop/R/computeFeatures.R')
    #picture <- dir()[1]
    A.1 <- EBImage::readImage(paste(folder,'/',picture, sep=''))
    resi <- dim(A.1[,,1])
    ####################################
    ## ----------------------------------
    ## orientate correctly
    ####################################
    if(resi[1] < resi[2]){
      ooo <- matrix(NA,nrow=resi[2], ncol=resi[1])
      A <- array(data=c(ooo,ooo,ooo), dim=c(nrow(ooo),ncol(ooo),3)) # creates an array
      for(i in 1:3){
        A[,,i] <- t(A.1[,,i])
      }
      ##
    }else{A <- A.1}
    
    
    
    ####################################
    ## ---------------------------------
    ## RESIZE
    ####################################
    resi3 <- dim(A[,,1]) 
    if(fact != 0){
      resi2 <- resi3 *fact
    }else{
      resi2 <- resi3
    }
    
    if(is.null(minArea)){
     minArea <- (resi2[1]/40) * (resi2[2]/40)
    }
    
    A <- EBImage::resize(A, w=resi2[1], h=resi2[2]) 
    
    B <- EBImage::channel(A,"gray") # transform to gray
    if(grays == TRUE){
      D <- list(r=A[,,1], g=A[,,2], b=A[,,3], g=B)
    }else{D <- list(r=A[,,1], g=A[,,2], b=A[,,3])}#, g=B)}
    
    #####################################
    
    col.indicator <- median(B) # if greater than 0.5 is white
    if(col.indicator > 0.35){
      back <- "w"
    }else{back <- "b"}
    #####################################
    if(is.null(cutoffvalue) & back == "b"){
      cutoffvalue <- 0.08
    }
    if(is.null(cutoffvalue) & back == "w"){
      cutoffvalue <- 0.6
    }
    #####################################
    if(back == "b"){
      B <- (B > cutoffvalue) # binary matrix
    }else{B <- (B < cutoffvalue)} # binary matrix}
    # if (interactive()) display(B, title="Cell nuclei")
    # cutoffvalue=0.08
    
    B <- EBImage::thresh(B, (sqrt((minArea*fact)/3.1459)/2), (sqrt((minArea*fact)/3.1459)/2), 0.05)
    B <-  EBImage::opening(B, EBImage::makeBrush((round(sqrt((minArea*fact)/3.1459)/10)), shape='disc')) # the brush disc has size r/10 = sqrt(A/pi)/10
    # fill holes in the picture and get label each bloob
    binaryImage <- EBImage::fillHull(B)     
    labeledImage <- EBImage::bwlabel(binaryImage) # this function identifies in which pixels the bloob is found
    # if (interactive()) EBImage::display(binaryImage, title="Cell nuclei")
    #str(binaryImage)
    www <- which(binaryImage@.Data == 0)
    col.zero <- lapply(D, function(x,www){x@.Data[www] <- 0; return(x) }, www)
    #str(col.zero)
    arr <- array(data=c(col.zero[[1]]@ .Data,col.zero[[2]]@ .Data,col.zero[[3]]@ .Data), dim=c(nrow(col.zero[[1]]@ .Data),ncol(col.zero[[1]]@ .Data),3))
    png::writePNG(arr, paste(stored,paste("used",picture,sep="-"),sep="/") )
    #writePNG(binaryImage, paste(getwd(),"measures_res",paste("used",picture,sep="-"),sep="/") )
    # target=paste(stored,  paste(picture,"at",cutoffvalue,"cutoffvalue.PNG", sep="_") , sep="/")
    
    #ptm <- proc.time()
    blobMeasurements <- EBImage::computeFeatures.shape(labeledImage)
    blobMeasurements <- blobMeasurements[,c("s.area", "s.perimeter", "s.radius.min", "s.radius.max")]
    blobMeasurements <- blobMeasurements * (1/fact) # adjust increasing the area by the factor you reduced the picture 
    #proc.time() - ptm
    # still need extract colors so corelation matches
    intensityMeasurements <- lapply(D,  EBImage::computeFeatures.basic, x=labeledImage, properties=F, basic.quantiles=0.05)
    intensityMeasurements <- lapply(intensityMeasurements, function(x){y <- x[,c("b.mean", "b.sd")]})
        
    #int <- intensityMeasurements[[1]]
    int <- intensityMeasurements[[1]]
    if(length(intensityMeasurements) > 1){
      for(l in 2:length(intensityMeasurements)){
        int <- cbind(int,intensityMeasurements[[l]])}}
    
    # getting rid of bad bloobs
    f <- which(as.vector(blobMeasurements[,1]) <= minArea) # remove bloobs that are smaller than minimum area specified
    if(length(f) > 0){
      blobMeasurements <- blobMeasurements[-f,]
      int <- int[-f,]
    } # is else necesary?
    # number of objects found in the picture after brushing
    numberOfBlobs <- dim(blobMeasurements)[[1]]
    # bind shape measurements and color measurements
    tmp <- cbind(blobMeasurements,int)
    # add shape and volume to the dataframe
    sha <- (2 * as.vector(tmp[,4])) / (2 * as.vector(tmp[,3])) 
    vol <- (4/3)*pi*sqrt(tmp[,1]/pi)^3 
    tmp <- cbind(tmp, sha, vol)
    
    
    # identify the references with respect to any other things in the picture
    two.groups <- (kmeans(tmp, 2))$cluster
    group1 <- tmp[which(two.groups == 1),]
    group2 <- tmp[which(two.groups == 2),]
    
    # variances for colors in both groups (including grey)
    #vvv1 <- mean(apply(group1, 1, function(y){ var(c(abs(y[5] - y[7]), abs(y[5] - y[9]), abs(y[7] - y[9])))  })) + var(group1[,"sha"])
    #vvv2 <- mean(apply(group2, 1, function(y){ var(c(abs(y[5] - y[7]), abs(y[5] - y[9]), abs(y[7] - y[9])))  })) + var(group2[,"sha"])
    
    if(back == "b"){ # if background is black the references tend to 1, then 1 - 1 tend to 0
    vvv1 <- 1 - mean(apply(group1[,c(5,7,9,11)],2,mean))
    vvv2 <- 1 - mean(apply(group2[,c(5,7,9,11)],2,mean))
    }else{ # if background is black the references tend to 0
      vvv1 <- mean(apply(group1[,c(5,7,9,11)],2,mean))
      vvv2 <- mean(apply(group2[,c(5,7,9,11)],2,mean))
    }
    
    #vvv1 <- var((group1[,1] - mean(group1[,1], na.rm=TRUE)), na.rm=TRUE)
    #vvv2 <- var((group2[,1] - mean(group2[,1], na.rm=TRUE)), na.rm=TRUE)
    #var(abs(y[5] - y[7]) + abs(y[5] - y[9]) + abs(y[7] - y[9]))
    
    variances <- c(vvv1, vvv2)
    # references identified
    ref <- which(variances == min(variances))
    refsCircles <- which(two.groups == ref)
    fru<- which(variances == max(variances))
    fruits <- which(two.groups == fru)

    # get medians for all features in reference and standarize
    ref.meds <- apply(tmp[refsCircles,], 2, FUN=mean, na.rm=T)
    
    
    ## standarize with respect to the references but don't do it for color
    if((stands[1] - stands[2]) == 0){ # if a circle
      tmp2 <- t(apply(tmp[fruits,], 1, function(x,x2){x3 <- x[c(1:4)]/x2[1:4]; x[c(1:4)] <- x3;  return(x)}, ref.meds))
    }else{ # if NOT a circle
      ref.len <- stands[1]/(2*ref.meds[4]) #equivalent length of one pixel
      ref.wid <- stands[2]/(2*ref.meds[3]) #equivalent length of one pixel
      sq.px <- ref.wid * ref.len          #equivalent area in cm2 of one pixel
      tmp[,3:4] <- tmp[,3:4] * 2 # change radius in length and width of the object
      tmp[,"s.area"] <- tmp[,"s.area"] * sq.px
      tmp[,"vol"] <- tmp[,"vol"] * sq.px
      tmp[,2:4] <- tmp[,2:4] * (mean(c(ref.len,ref.wid), na.rm=T))
      tmp2 <- tmp[fruits,]
    }
    
    # names for the dataframe depending if gray channel is included or not
    tmpe <- as.data.frame(tmp2)
    if(grays == T){
      names(tmpe) <- c("Area","Perimeter","Width","Length","Red","sd_Red","Green","sd_Green","Blue","sd_Blue","Gray","sd_Gray","Shape","Volume")
    }else{names(tmpe) <- c("Area","Perimeter","Width","Length","Red","sd_Red","Green","sd_Green","Blue","sd_Blue","Shape","Volume")}
    print("Plant done")
    return(tmpe)
  }
  
  x.core <- parallel::detectCores()
  n.cores <-  cores   # detects cores and store a number, all - 2
  cat(paste("\nYou are using",cores,"core(s) out of", x.core,"available\n"))
  cat(paste("\nPlease monitor the progress of pictures analyzed at the folder:\n"))
  cat(paste("   ", sto,"\n\n"))
  cl <- parallel::makeCluster(n.cores)     # number of cores to be used
  doParallel::registerDoParallel(cl)             # stablishing connection
  
  i=NULL
  
  x <- foreach::foreach(i=1:length(listp2)) %dopar% scan(folder,paste(listp2[i]), cutoffvalue=cutoffvalue, minArea=minArea, grays=gray, stands=stand, fact=fact, stored=sto)  # parallel function
  
  names(x) <- listp2
  
  parallel::stopCluster(cl)                # stops connection
  #stopCluster(cl)                # stops connection
  
  ################# -------------------------------
  ## make a dataframe instead of a list
  x <- lapply(x, function(x){rownames(x) <- NULL; return(x)}) 
  non <- length(x)
  dimens <- unlist(lapply(x,function(x){dim(x)[1]}))
  names.pics <- list(NA)
  for(j in 1:length(dimens)){names.pics[[j]] <- rep(names(dimens)[j],dimens[j])};# names.pics <- unlist(names.pics)
  
  for(k in 1:length(x)){
    rownames(x[[k]]) <- paste(names.pics[[k]],rownames(x[[k]]), sep=".")
  }
  ## based on number of pictures make a new matrix
  if(non > 1){
    x3 <- x[[1]]
    for(j in 2:length(x)){
      x3 <- rbind(x3,x[[j]])
    }
  }else{x3 <- x[[1]]}
  
  return(x3)
}