cutoff <- function(folder,cores=1, prefs=NULL){
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("EBImage needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("parallel needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  tes <- dir(folder) # IS A PICTURE NAME OR FOLDER??
  if(length(tes) == 0){# is a picture
    listp2 <- folder # pictures to analize
    setwd(getwd()) # setwd
    folder <- getwd() # set a folder
  }else{# is a folder
    setwd(folder)
    listp2 <- dir(folder, '*.JPG')
  }
  ###################
    
  if(length(which(dir() %in%  "cutoff_vals")) > 0){
    i=1
    while(length(which(dir() %in%  paste("cutoff_vals",i, sep=""))) > 0){
      i=i+1
    }
    sto <- paste(getwd(), paste("cutoff_vals",i, sep=""),sep="/")
    dir.create(paste(getwd(), paste("cutoff_vals",i, sep=""),sep="/"))
    
  }else{sto <- paste(getwd(), "cutoff_vals",sep="/");  dir.create(paste(getwd(), "cutoff_vals",sep="/"))}
  #Progress combine function
  
  prefuse <- prefs
  #backs <- back
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
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  scan <- function(folder,picture, pref, stored=getwd()){
    print(i)
    
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
    resi2 <- dim(A[,,1]) *0.5
    A <- EBImage::resize(A, w=resi2[1], h=resi2[2]) 
    #if (interactive()) display(y, title='resize(x, 128)')
    
    B <- EBImage::channel(A,"gray") # transform to gray
    
    ### find what is the background color
    col.indicator <- median(B) # if greater than 0.5 is white
    if(col.indicator > 0.35){
      backg <- "w"
    }else{backg <- "b"}
    
    ## based on the background color and the argument at the beggining decide the preferences
    if(is.null(pref)){
      
      if(backg == "w"){
        pref <- seq(0.5,0.8,by=0.05) ########!!!!!!!!!!!!!!!!! depending the background
      } 
      if(backg == "b"){
        pref <- seq(0.01,.12,by=0.02)
      }     
    }
    ####################################
    #if (interactive()) display(B, title="Cell nuclei")
    #B <- A[,,3]
    for(i in 1:length(pref)){
      
      
      
      cutoffvalue <- pref[i]
      if(backg == "b"){
        B.2 <- (B > cutoffvalue) ########!!!!!!!!!!!!!!!!! depending the background
      } else{B.2 <- (B < cutoffvalue)}
      B.2 <- EBImage::thresh(B.2, 50, 50, 0.05)
      B.2 <-  EBImage::opening(B.2, EBImage::makeBrush(5, shape='disc'))
      binaryImage <- EBImage::fillHull(B.2)
      
      B.out <- B
      B.out[which( binaryImage == 0, arr.ind=TRUE)] <- 0
      png::writePNG(B.out, target=paste(stored,  paste(picture,"at",cutoffvalue,"cutoffvalue.PNG", sep="_") , sep="/"))
      
    }
    # binary matrix
    #return(y)
  }
  ####################################################################################################### scan(folder,paste(listp2[i]), pref=prefuse, backg=backs)  # parallel function
  #######################################################################################################
  #######################################################################################################
  #######################################################################################################
  
  
  x.core <- parallel::detectCores()
  n.cores <-  cores   # detects cores and store a number, all - 2
  cat(paste("\nYou are using",cores,"core(s) out of", x.core,"available"))
  cat(paste("\n\nPlease monitor the progress of the segmentation for the different values of color threshold at the folder:\n"))
  cat(paste(sto))
  if(!is.null(prefs)){
  cat(paste("\n\nYou selected to test the following values:\n\n",paste(prefs, collapse = ", "), "\nPlease check the segmentation pictures and decide the best value of segmentation\n"))
  }else{cat(paste("\n\nYou selected to test the following values:\n\n",paste(seq(0.5,0.8,by=0.05), collapse=", "),"if background is white, or\n", paste(seq(0.01,.12,by=0.02), collapse=", "),"if background is black.","\n\nPlease check the segmentation pictures and decide the best value of segmentation for posterior analysis\n\n"))
  }
  cl <- parallel::makeCluster(n.cores)     # number of cores to be used
  doParallel::registerDoParallel(cl)             # stablishing connection
    
  x <- foreach::foreach(i=1:length(listp2)) %dopar% scan(folder,paste(listp2[i]), pref=prefuse, stored=sto)
  #x <- foreach::foreach(i=1:length(listp2), .combine=scan) %dopar%  {folder,paste(listp2[i]), pref=prefuse, stored=sto}
  
  #foreach(i=1:nrow(m), .combine=rbind) %dopar%    (m[i,] / mean(m[i,]))
  
  cat(paste("\nPlease remember to use the optimum cutoff value we have found for the posterior functions \n\nSee the results in the folder: \n\n", sto, "\n\n"))
  ##################################
  parallel::stopCluster(cl) # stop connection
  ##################################
  # prepare data frmae for plot
  #return(y2)
}