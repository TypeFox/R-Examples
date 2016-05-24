mixCov <- function(..., use = "everything", method = "pearson")
{
  ## convert arguments to named list
  covLIST <- list(...)
  covNAMES <- as.list(sys.call())[-1]
  names(covLIST)[1:length(covNAMES)] <- covNAMES
  
  ## function for inflating covariance matrices
  makeCov <- function(x) {
    ## calculate final dimension
    LEN <- sapply(x, function(a) NCOL(a))  
    DIM <- sum(LEN, na.rm = TRUE)
    covMAT <- matrix(0, DIM, DIM)
    NAMES <- names(x)
                
    ## initialize counter and names vector
    counter <- 1
    nameVEC <- vector("character", length = length(NAMES))
       
    ## fill covariance matrix
    for (i in 1:length(x)) {
      ITEM <- x[[i]]
      POS <- counter:(counter + LEN[i] - 1)
      covMAT[POS, POS] <- ITEM
      if (LEN[i] == 1) nameVEC[POS] <- NAMES[i] else nameVEC[POS] <- colnames(ITEM)
      counter <- counter + LEN[i]      
    } 
    
    colnames(covMAT) <- rownames(covMAT) <- make.unique(nameVEC)      
    return(covMAT)
  }
  
  ## preformat data, depending on input
  preLIST <- vector("list", length = length(covLIST))
  
  for (i in 1:length(covLIST)) {
    tempDAT <- covLIST[[i]]   
    
    ## check if input has features of a covariance matrix
    if (!is.null(rownames(tempDAT))) {
      if ((all(rownames(tempDAT) == colnames(tempDAT))) && (nrow(tempDAT) == ncol(tempDAT))) 
      isCov <- TRUE
    } else isCov <- FALSE
    
    ## case 1: input is raw data => covariance matrix
    if (NROW(tempDAT) > 1 && NCOL(tempDAT) > 1) preLIST[[i]] <- cov(tempDAT, use = use, method = method) 
    
    ## case 2: input data is covariance matrix
    if (isCov) preLIST[[i]] <- tempDAT
    
    ## case 3: input data is c(mean, sd) => sd^2
    if (NROW(tempDAT) == 2 && NCOL(tempDAT) == 1) preLIST[[i]] <- tempDAT[2]^2 
    
    ## case 4: input data is single error value
    if (NROW(tempDAT) == 1 && NCOL(tempDAT) == 1) preLIST[[i]] <- tempDAT^2    
  }
  
  names(preLIST) <- names(covLIST)  
  
  ## run makeCov over preLIST
  OUT <- makeCov(preLIST)
  
  return(OUT)  
}


