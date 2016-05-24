stratamean <-
function(y, h, Nh, wh, level = 0.95, eae = FALSE)
{
### input treatment
  # data and strata
  if(!(is.numeric(y) && is.vector(y))){
	stop("Invalid input: ", sQuote("y"), " has to be a numeric vector of observations.")
  }
  if(!is.vector(h)){
	stop("Invalid input: ", sQuote("h"), " has to be a vector indicating the strata.")
  }else{
	M <- nlevels(as.factor(h)) # number of strata
  }
  if(length(y) != length(h) ){
	stop("Invalid input: ", sQuote("y"), " and ", sQuote("h"), " must have the same length.")
  }
  # define data set
  data <- data.frame(target=y, strata=h)
  # omit NA values in data
  if(any(is.na(data))){
    data <- na.omit(data)
    warning("Missing values will be ignored in calculations.")
  }
  # each strata need more than one observation
  if(any(table(data$strata) < 2)){
    stop("Invalid input: stratum with only one observation was entered.")
  }

  ## strata weights
  if(missing(Nh) == missing(wh)){  
    stop("Invalid input: only ", sQuote("Nh"), " or ", sQuote("wh"), " valid inputs." )
  }
  # define strata weights by wh
  if(missing(Nh)){
    if( M != length(wh) ){
      stop("Invalid input: the number of strata given in ", sQuote("wh"), " is not valid according the sample data.")
    }
    if( sum(wh) != 1){
      stop("Invalid input: the sum of ", sQuote("wh"), " must equal 1.")
    }
    Nh <- NA
    fpc <- FALSE
  }
  # define strata weights by Nh
  else{ #if(missing(wh))
    if( M != length(Nh) ){
      stop("Invalid input: the number of strata given in ", sQuote("Nh"), " is not correct according the sample data.")
    }
    N <- sum(Nh)
    wh <- Nh/N
    fpc <- TRUE
  }
  # level and quantile
  if( level < 0 || level > 1 ){
    stop("Invalid input: ", sQuote("level"), " has to be probability between 0 and 1.")
  }else{
    q <- qnorm((1 + level) / 2)
  }

  ### calculate strata means
  splitted <- split(data, data$strata)	
  CIuh <- CIoh <- Varh <- Meanh <- vector(length = M) 	
  for(i in 1:M){
    size <- ifelse(fpc, Nh[i], Inf)
    if (dim(splitted[[i]])[1] > size){
      stop("Invalid input: the population of a stratum can not be less than the number of observations of this stratum")
    }
    Smean.i <- Smean(y = as.numeric(as.vector(splitted[[i]]$target)), N=size)
    Meanh[i] <- Smean.i$mean
    Varh[i] <- Smean.i$se^2
    CIuh[i] <- Smean.i$ci[1]
    CIoh[i] <- Smean.i$ci[2]
  }
  eaes <- cbind(Mean = Meanh, SE = sqrt(Varh), CIu = CIuh, CIo = CIoh)
  rownames(eaes) <- sort(unique(h))

  ### compute overall mean
  Mean <- sum(Meanh*wh)
  Var <- sum(Varh*((wh)^2))
  CIo <- Mean + q*sqrt(Var)
  CIu <- Mean - q*sqrt(Var)
  overall <- c(Mean, sqrt(Var), CIu, CIo)

  ### output
  if(eae){
    res <- rbind(eaes,overall)
    return(res) 
  }else{
    ret <- list()
    ret$call <- list(y = y, h = h, Nh = Nh, wh = wh, fpc = fpc, level = level)
    ret$mean <- Mean
    ret$se <- sqrt(Var)
    ret$ci <- c(CIu, CIo)
    structure(ret, class = "stratamean")   
  }
}
