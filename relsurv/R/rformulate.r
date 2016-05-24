rformulate <- function (formula, data = parent.frame(), ratetable, na.action, 
                        int, centered, cause) 
{ 
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$ratetable <- m$int <- m$centered <- NULL
  if(!missing(cause)){  					#NEW: ce cause obstaja
    Terms <- if (missing(data)) 
      terms(formula, "ratetable", "cause")
    else terms(formula, "ratetable", "cause",  data = data)
  }
  else{								#NEW: ce cause ne obstaja 
    Terms <- if (missing(data)) 
      terms(formula, "ratetable")
    else terms(formula, "ratetable", data = data)
  }
  rate <- attr(Terms, "specials")$ratetable
  if (length(rate) > 1) 
    stop("Can have only 1 ratetable() call in a formula")
  if (length(rate) == 0) {
    xx <- function(x) formula(x)
    if (is.ratetable(ratetable)) 
      varlist <- attr(ratetable, "dimid")
    else stop("Invalid rate table")
    ftemp <- deparse(formula)
    ftemp <- paste(ftemp,collapse="")
    formula <- xx(paste(ftemp, "+ ratetable(", paste(varlist, 
                                                     "=", varlist, collapse = ","), ")"))
    Terms <- if (missing(data)) 
      terms(formula, "ratetable")
    else terms(formula, "ratetable", data = data)
    rate <- attr(Terms, "specials")$ratetable
  }
  
  strats <- NULL 						#NEW in 2.05
  if (missing(cause)){ 
    Terms <-if (missing(data))
      terms(formula,c("strata","ratetable"))
    else terms(formula,c("strata","ratetable"),data=data)
    rate <- attr(Terms, "specials")$ratetable
    strats <- attr(Terms,"specials")$strata
  }
  
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  n <- nrow(m)
  Y <- model.extract(m, "response")

  if (!is.Surv(Y)) 
    stop("Response must be a survival object")
  Y.surv <- Y
  if (attr(Y, "type") == "right") {
    type <- attr(Y, "type")
    status <- Y[, 2]
    Y <- Y[, 1]
    start <- rep(0, n)
    ncol0 <- 2
  }
  else if (attr(Y, "type") == "counting") {
    type <- attr(Y, "type")
    status <- Y[, 3]
    start <- Y[, 1]
    Y <- Y[, 2]
    ncol0 <- 3
  }
  else stop("Illegal response value")
  if (any(c(Y, start) < 0)) 
    stop("Negative follow up time")
  if(max(Y)<30)
    warning("The event times must be expressed in days! (Your max time in the data is less than 30 days) \n")
  
  if (is.ratetable(ratetable)) {
    israte <- TRUE
    rtemp <- match.ratetable(m[, rate], ratetable)
    if(is.null(attributes(ratetable)$factor))attributes(ratetable)$factor <- attributes(ratetable)$type==1
    rtorig <- attributes(ratetable)
    nrt <- length(rtorig$dimid)
    R <- rtemp$R
    if (!is.null(rtemp$call)) {
      ratetable <- eval(parse(text = rtemp$call))
    }
    #checking if the ratetable variables are given in days
    wh.age <- which(attributes(ratetable)$dimid=="age")
    wh.year <- which(attributes(ratetable)$dimid=="year")
    if(length(wh.age)>0){
      if(max(rtemp$R[,wh.age])<150& median(diff(attributes(ratetable)$cutpoints[[wh.age]]))>12)
        warning("Age in the ratetable part of the formula must be expressed in days! \n (Your max age is less than 150 days) \n")
    }
    if(length(wh.year)>0){
      if(min(rtemp$R[,wh.year])>1850 & max(rtemp$R[,wh.year])<2020&class(attributes(ratetable)$cutpoints[[wh.year]])=="date")
        warning("The calendar year must be expressed in days since 1.1.1960! \n (Your variable seems to be expressed in years) \n")
    }
    #checking if one of the continuous variables is fixed:
    if(nrt!=ncol(R)){
      nonex <- which(is.na(match(rtorig$dimid,attributes(ratetable)$dimid)))
      for(it in nonex){
        if(rtorig$type[it]!=1)warning(paste("Variable ",rtorig$dimid[it]," is held fixed even though it changes in time in the population tables. \n (You may wish to set a value for each individual and not just one value for all)",sep=""))
      }
    }
    
  }
  else stop("Invalid ratetable argument")
    
 #NEW in 2.05 (strata)
if (missing(cause)){
  if (length(strats)){
    temp_str <- untangle.specials(Terms,"strata",1)
    dropx <- temp_str$terms
    if (length(temp_str$vars) == 1){
      strata.keep <- m[[temp_str$vars]]
    }
    else strata.keep <- strata(m[,temp_str$vars],shortlabel=TRUE,sep=",")
  }
  else strata.keep <- factor(rep(1,n)) # zgoraj ze definirano n = nrow(m)
}
  len_term.labels <- length(attr(Terms,"term.labels"))
  # del za X
  if (( (len_term.labels == 1)&(rate == 2) ) | ( (len_term.labels == 2)&(length(strats) == 1) )) {# spremenjeno ob dodajanju za strato
    # 1. pogoj: samo ratetable
    # 2. pogoj: samo ratetable in strata
    X <- NULL
    mm <- 0
  }
  else { 
    if (length(rate == 1)) {
      formula1 <- formula	#NEW: create object formula1
      #formula1[[3]] <- eval(parse(text=paste("formula1[[3]]",paste(rep("[[2]]",length(strats)+1),collapse = ""),sep=""))) # delete the ratetable (and strata) part
      f_temp1_part2 <- as.character(formula1[[2]])
      if (type == "right")f_temp1 <- paste(f_temp1_part2[1],"(",f_temp1_part2[2],",",f_temp1_part2[3],") ~ ",sep="")
      else if(type == "counting")f_temp1 <- paste(f_temp1_part2[1],"(",f_temp1_part2[2],",",f_temp1_part2[3],",",f_temp1_part2[4],") ~ ",sep="")
      
      f_temp2_temp <- attr(Terms,"term.labels")[-c(rate-1,strats-1)]
      if (length(f_temp2_temp) == 1){
        f_temp2 <- f_temp2_temp
      }
      else {
        f_temp2 <- f_temp2_temp[1]
        for (l in 2:length(f_temp2_temp)){
          f_temp2 <- paste(f_temp2,f_temp2_temp[l],sep=" + ")
        }
      }
      f_temp <- paste(f_temp1,f_temp2,paste="")
      f_temp <- formula(f_temp)
    }        
    
    
    #X <- as.data.frame(model.matrix(formula1, data = data))[,-1, drop = FALSE]
    X <- as.data.frame(model.matrix(f_temp, data = data))[,-1, drop = FALSE]		#NEW in 2.05
    if(nrow(X)!=n)warning("You have missing values in demographic variables \n")
    mm <- ncol(X)
  }
  mvalue <- rep(0,mm)
  if (!missing(centered)) {
    if (mm != 0 & centered == TRUE) {
      mvalue <- apply(as.matrix(X),2,mean)
      X <- apply(as.matrix(X), 2, function(x) x - mean(x))
    }
  }
  offset <- attr(Terms, "offset")
  tt <- length(offset)
  offset <- if (tt == 0) 
    rep(0, n)
  else if (tt == 1) 
    m[[offset]]
  else {
    ff <- m[[offset[1]]]
    for (i in 2:tt) ff <- ff + m[[offset[i]]]
    ff
  }
  keep <- Y > start
  cause <- model.extract(m, "cause")
  if(is.null(cause)) cause <- rep(2,nrow(m))					#NEW: ce cause manjka
  #status[cause==0] <- 0
  if (!missing(int)) {
    int <- max(int)
    status[Y > int * 365.241] <- 0
    Y <- pmin(Y, int * 365.241)
    keep <- keep & (start < int * 365.241)
  }
  if (any(start > Y) | any(Y < 0)) 
    stop("Negative follow-up times")
  X <- X[keep, , drop = FALSE]
  Y <- Y[keep]
  start <- start[keep]
  status <- status[keep]
  R <- R[keep, ,drop=FALSE] 
  strata.keep <- strata.keep[keep] # dodano za strato  #NEW in 2.05
  offset <- offset[keep]
  Y.surv <- Y.surv[keep, , drop = FALSE]
  cause <- cause[keep]
  n <- sum(keep)
  data <- data.frame(start = start, Y = Y, stat = status, R)
  if (mm != 0) 
    data <- cbind(data, X)
  out <- list(data = data, R = R, status = status, start = start, 
              Y = Y, X = as.data.frame(X), m = mm, n = n, type = type, Y.surv = Y.surv, 
              Terms = Terms, ratetable = ratetable, offset = offset, formula=formula,
              cause = cause,mvalue=mvalue,strata.keep=strata.keep) # dodano za strato  #NEW in 2.05
  na.action <- attr(m, "na.action")
  if (length(na.action)) 
    out$na.action <- na.action
  out
}
