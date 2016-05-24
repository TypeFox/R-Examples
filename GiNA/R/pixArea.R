pixArea <- function(folder,cutoffvalue=0.5,cores=1, square=10, fact=0.25){
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
    setwd(getwd()) # setwd
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
  #dir.create(paste(getwd(),"measures_res",sep="/"))
 
  #listp2 <- dir(folder, '*.JPG')
  #nref <- nrefs/2
  
  x <- numeric()# creating a null vector or matrix for parallelized result
  #######################################################################################################
  scan <- function(folder,picture, squares, facts=0.25, cutoffvalue){
    #print(i)
    # picture <- listp2[1]
    A.1 <- EBImage::readImage(paste(folder,'/',picture, sep=''))
    resi <- dim(A.1[,,1])
    ####################################
    ## ----------------------------------
    ## orientate correctly
    ####################################
    if(resi[1] < resi[2]){
      ooo <- matrix(NA,nrow=resi[2], ncol=resi[1])
      A <- array(data=c(ooo,ooo,ooo), dim=c(nrow(ooo),ncol(ooo),3)) # creates an array
      for(h in 1:3){
        A[,,h] <- t(A.1[,,h])
      }
      ##
    }else{A <- A.1}
    
    
    ####################################
    ## ---------------------------------
    ## RESIZE
    ####################################
    resi2 <- dim(A[,,1]) *facts
    A <- EBImage::resize(A, w=resi2[1], h=resi2[2]) 
    #if (interactive()) display(B, title='resize(x, 128)')
    
    B <- EBImage::channel(A,"gray") # transform to gray
    
    ### find what is the background color
    col.indicator <- median(B) # if greater than 0.5 is white
    if(col.indicator > 0.35){
      backg <- "w"
    }else{backg <- "b"}
    
    #####################################
    if(backg == "b"){
      B <- (B > cutoffvalue) # binary matrix
    }else{B <- (B < cutoffvalue)} # binary matrix}
    B <- EBImage::thresh(B, squares, squares, 0.05)
    B <-  EBImage::opening(B, EBImage::makeBrush(5, shape='disc'))
    # if (interactive()) display(B, title="Cell nuclei")
    
    # the display() funtion let us see the images we are creating
    binaryImage <- EBImage::fillHull(B)     
    labeledImage <- EBImage::bwlabel(binaryImage) # this function identifies in which pixels the bloob is found
    #if (interactive()) EBImage::display(labeledImage, title="Cell nuclei")
    #writePNG(binaryImage, paste(getwd(),"measures_res",paste("used",picture,sep="-"),sep="/") )
    blobMeasurements <- EBImage::computeFeatures.shape(labeledImage)#, binaryImage)
    ## -------------------------
    ## -------------------------
    # identify the references with respect to any other things in the picture
    two.groups <- (kmeans(blobMeasurements, 2))$cluster
    group1 <- blobMeasurements[which(two.groups == 1),]
    group2 <- blobMeasurements[which(two.groups == 2),]
    
    vvv1 <- var((group1[,1] - mean(group1[,1], na.rm=TRUE)), na.rm=TRUE)
    vvv2 <- var((group2[,1] - mean(group2[,1], na.rm=TRUE)), na.rm=TRUE)
    
    variances <- c(vvv1, vvv2)
    # references identified
    refs <- which(variances == min(variances))
    references <- blobMeasurements[which(two.groups == refs),]
    frus<- which(variances == max(variances))
    fruits <- blobMeasurements[which(two.groups == frus),]
    #intensityMeasurements <- computeFeatures.basic(labeledImage, A)
    #hist(blobMeasurements[,"s.area"], main="Pixel area distribution", col=2:20)
    #plot(1:length(blobMeasurements[,"s.area"]),sort(blobMeasurements[,"s.area"]), yaxt="n", col="blue", ylab="Pixel area", xlab="sorted object measure by area", main=paste(picture))
    #axis(side=2,at=blobMeasurements[,"s.area"],labels=blobMeasurements[,"s.area"], las=2, cex.axis=.5)
    
    #references$type <- "ref"
    #fruits$type <- "fru"
    #resolution <- data.frame(rbind(references,fruits)) 
    #y <- blobMeasurements[,c("s.area","type")]
    #y[,1] <- y[,1]*(1/fact)
    
    y <- blobMeasurements#[,"s.area"]*(1/fact)
    return(y)
  }
  #######################################################################################################
  x.core <- parallel::detectCores()
  n.cores <-  cores   # detects cores and store a number, all - 2
  cat(paste("\nYou are using",cores,"core(s) out of", x.core,"available"))
  
  cl <- parallel::makeCluster(n.cores)     # number of cores to be used
  doParallel::registerDoParallel(cl)             # stablishing connection
  
  #x <- foreach(j=1:length(listp2)) %dopar% scan(folder=folder,picture=paste(listp2[j]), squares=square, facts=fact)  # parallel function
  #x <- list()
  i=NULL
  
  x <- foreach::foreach(i=1:length(listp2))%dopar%scan(folder=folder,picture=paste(listp2[i]), squares=square, facts=fact, cutoffvalue=cutoffvalue)
  
  ##################################
  parallel::stopCluster(cl) # stop connection
  ##################################
  # prepare data frmae for plot
  y2 <- data.frame(resp=x[[1]], group=rep(paste("P1"),dim(x[[1]])[1]), xx=1:dim(x[[1]])[1])
  if(length(x) > 1){
    for(g in 2:length(x)){
      prov <- data.frame(resp=x[[g]], group=rep(paste("P",g,sep=""),dim(x[[g]])[1]), xx=1:dim(x[[g]])[1])
      y2 <- rbind(y2,prov)
    }
  }else{y2 <- y2}
  #
  transp <- function (col, alpha = 0.5) {
    res <- apply(grDevices::col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                             c[3]/255, alpha))
    return(res)
  }
  # get colors and plot
  yoyoyo <- round(mean(y2$resp.s.area)/4, 0)
  col.scheme <- rep(brewer.pal(9,"Set1"),length(x))
  g1 <- seq(1,length(y2$resp.s.area),by=3)
  plot(jitter(y2$xx, factor=2), jitter(y2$resp.s.area, factor=2), yaxt="n",col=transp(col.scheme[factor(y2$group)],.3), cex=3, pch=20,ylab="Pixel area", xlab="sorted object measure by area", main="Picture area compendium")
  #axis(side=2,at=sort(y2$resp)[g1],labels=sort(y2$resp)[g1], las=2, cex.axis=.5)
  axis(side=2,at=seq(0,900000,yoyoyo),labels=seq(0,900000,yoyoyo), las=2, cex.axis=.5)
  grid()
  legend("topright",legend=levels(factor(y2$group)),col=col.scheme,pch=20,cex=0.6, bty="n")
  
  cat(paste("\n\nMinimum values were found at",min(y2$resp.s.area), "pixels (see the plot). \n\nPlease use an approximate value as a minimum area when using the scancran function.\n\n"))
  
  return(y2)
}