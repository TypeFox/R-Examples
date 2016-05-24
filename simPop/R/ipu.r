#' iterative proportional updating
#' 
#' adjust sampling weights to given totals based on household-level and/or
#' individual level constraints
#' 
#' @name ipu
#' @param inp a \code{data.frame} or \code{data.table} containing household ids
#' (optionally), counts for household and/or personal level attributes that
#' should be fitted.
#' @param con named list with each list element holding a constraint total with
#' list-names relating to column-names in \code{inp}.
#' @param hid character vector specifying the variable containing household-ids
#' within \code{inp} or NULL if such a variable does not exist.
#' @param eps number specifiying convergence limit
#' @param verbose if TRUE, ipu will print some progress information.
#' @export
#' @author Bernhard Meindl
#' @keywords method
#' @examples
#' # basic example
#' inp <- as.data.frame(matrix(0, nrow=8, ncol=6))
#' colnames(inp) <- c("hhid","hh1","hh2","p1","p2","p3")
#' inp$hhid <- 1:8
#' inp$hh1[1:3] <- 1
#' inp$hh2[4:8] <- 1
#' inp$p1 <- c(1,1,2,1,0,1,2,1)
#' inp$p2 <- c(1,0,1,0,2,1,1,1)
#' inp$p3 <- c(1,1,0,2,1,0,2,0)
#' con <- list(hh1=35, hh2=65, p1=91, p2=65, p3=104)
#' res <- ipu(inp=inp, hid="hhid", con=con, verbose=FALSE)
#' 
#' # more sophisticated
#' # load sample and population data
#' data(eusilcS)
#' data(eusilcP)
#' 
#' # variable generation and preparation
#' eusilcS$hsize <- factor(eusilcS$hsize)
#' 
#' # make sure, factor levels in sample and population match
#' eusilcP$region <- factor(eusilcP$region, levels = levels(eusilcS$db040))
#' eusilcP$gender <- factor(eusilcP$gender, levels = levels(eusilcS$rb090))
#' eusilcP$hsize  <- factor(eusilcP$hsize , levels = levels(eusilcS$hsize))
#' 
#' # generate input matrix
#' # we want to adjust to variable "db040" (region) as household variables and
#' # variable "rb090" (gender) as individual information
#' samp <- data.table(eusilcS)
#' pop <-  data.table(eusilcP)
#' setkeyv(samp, "db030")
#' hh <- samp[!duplicated(samp$db030),]
#' hhpop <- pop[!duplicated(pop$hid),]
#' 
#' # reg contains for each region the number of households
#' reg <- data.table(model.matrix(~db040 +0, data=hh))
#' # hsize contains for each household size the number of households
#' hsize <- data.table(model.matrix(~factor(hsize) +0, data=hh))
#' 
#' # aggregate persons-level characteristics per household
#' # gender contains for each household the number of males and females
#' gender <- data.table(model.matrix(~db030+rb090 +0, data=samp))
#' setkeyv(gender, "db030")
#' gender <- gender[, lapply(.SD, sum), by = key(gender)]
#' 
#' # bind together and use it as input
#' inp <- cbind(reg, hsize, gender)
#' 
#' # the totals we want to calibrate to
#' con <- c(
#'   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$region)),
#'   as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$hsize)),
#'   as.list(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$gender))
#' )
#' # we need to have the same names as in 'inp'
#' names(con) <- setdiff(names(inp), "db030")
#' 
#' # run ipu und check results
#' res <- ipu(inp=inp, hid="db030", con=con, verbose=TRUE)
#' 
#' is <- sapply(2:(ncol(res)-1), function(x) { 
#'   sum(res[,x]*res$weights)
#' }) 
#' data.frame(required=unlist(con), is=is)
ipu <- function(inp, con, hid=NULL, eps=1e-07, verbose=FALSE) {
  if ( class(inp)[1] %in% c("data.frame", "data.table") ) {
    cnames <- colnames(inp)
    inp <- as.matrix(inp)
    colnames(inp) <- cnames
  }

  if ( is.null(hid) ) {
    hhid <- 1:nrow(inp)
  } else {
    ii <- match(hid, colnames(inp))
    hhid <- inp[,ii]
    inp <- inp[,-c(ii)]
  }
  inp <- inp[,match(names(con), colnames(inp))]
  con <- con[match(colnames(inp), names(con))]
  if ( any(colnames(inp) != names(con)) ) {
    stop("constraints (con) do not match input (inp) -> check your input!\n")
  }
  w <- as.numeric(rep(1,nrow(inp)))
  con <- as.numeric(unlist(con))
  w <- .Call("simPop_ipu_work", inp, con, w, eps=eps, verbose=ifelse(verbose, 1L, 0L), package="simPop")
  out <- cbind(hhid, inp)
  out <- as.data.frame(out)
  out$weights <- w
  invisible(out)
}
ipu2 <- function(dat,conP,conH,hid=NULL,epsP=1e-2,epsH=1e-2,verbose=FALSE,
    w=NULL,bound=4,maxIter=200,meanHH=TRUE){

  wvst <- melt <- calibWeight <- wValue <- f <- value <- baseWeight <- NULL

  # dat sollte ein data.table sein
  # w ein Name eines Basisgewichts oder NULL
  ncp <- length(conP) # number of constraints on person level
  nch <- length(conH) # number of constraints on household level
  dimncp <- sapply(conP,function(x)prod(dim(x)))
  dimnch <- sapply(conH,function(x)prod(dim(x)))
  valueP <- paste0("valueP",seq_along(conP))###fixed target value, should not be changed in iterations
  valueH <- paste0("valueH",seq_along(conH))
  ###Housekeeping of the varNames used
  usedVarNames <- c(valueP,valueH,"value","baseWeight","wvst","wValue")
  renameVars <- NULL
  if(any(names(dat)%in%usedVarNames)){
    renameVars <- names(dat)[names(dat)%in%usedVarNames]
    setnames(dat,renameVars,paste0(renameVars,"_safekeeping"))
  }
  ### Treatment of HID, creating 0,1 var for being the first hh member
  delVars <- c()
  if(is.null(hid)){
    delVars <- c("hid")
    hid <- "hid"
    dat[,hid:=1:nrow(dat)]
    dat[,wvst:=1]
  }else{
    dat[,wvst:=as.numeric(!duplicated(dat[,hid,with=FALSE]))]
  }


  boundsFak <- function(g1,g0,f,bound=4){ # Berechnet die neuen Gewichte (innerhalb 4, .25 Veraenderungsraten)
    g1 <- g1 * f
    g1[(g1/g0)>bound] <- bound*g0[(g1/g0)>bound]
    g1[(g1/g0)<(1/bound)] <- (1/bound)*g0[(g1/g0)<(1/bound)]
    return(g1)
  }
  mconP <- lapply(conP,melt)##convert tables to long form
  mconH <- lapply(conH,melt)

  for(i in seq_along(conP)){
    dat <- merge(dat,mconP[[i]],by=colnames(mconP[[i]])[-ncol(mconP[[i]])])
    setnames(dat,"value",valueP[i])
  }
  for(i in seq_along(conH)){
    dat <- merge(dat,mconH[[i]],by=colnames(mconH[[i]])[-ncol(mconH[[i]])])
    setnames(dat,"value",valueH[i])
  }
  pCalVar <- paste0("pcal",1:ncp)
  hCalVar <- paste0("hcal",1:nch)
  if(is.null(w)){
    if(!is.null(bound)&&is.null(w))
      stop("Bounds are only reasonable if base weights are provided")
    dat[,calibWeight:=1]
    delVars <- c(delVars,"baseWeight")
  }else{
    dat[,calibWeight:=dat[,w,with=FALSE]]
    setnames(dat,w,"baseWeight")
  }
  ## Names of the calibration variables for Person and household dimension
  pColNames <- lapply(conP,function(x)names(dimnames(x)))
  hColNames <- lapply(conH,function(x)names(dimnames(x)))
  ###Calib
  error <- TRUE
  calIter <- 1
  while(error&&calIter<=maxIter){
    error <- FALSE

    ### Person calib
    for(i in seq_along(conP)){
      dat[,wValue:=sum(calibWeight),by=eval(pColNames[[i]])]
      setnames(dat,valueP[i],"value")
      dat[,f:=value/wValue,by=eval(pColNames[[i]])]
      if(!is.null(bound)){
        dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound),by=eval(pColNames[[i]])]
      }else{
        dat[,calibWeight:=f*calibWeight,by=eval(pColNames[[i]])]
      }
      setnames(dat,"value",valueP[i])
      curEps <- dat[,max(abs(f-1))]
      if(curEps>epsP){
        error <- TRUE
      }
      if(verbose&&curEps&&calIter%%10==0){
        cat(calIter, ":Not yet converged for P-Constraint",i,"\n")
      }
    }
    ### Household calib
    for(i in seq_along(conH)){
      dat[,wValue:=sum(calibWeight*wvst),by=eval(hColNames[[i]])]
      setnames(dat,valueH[i],"value")
      dat[,f:=value/wValue,by=eval(hColNames[[i]])]
      if(!is.null(bound)){
        dat[,calibWeight:=boundsFak(calibWeight,baseWeight,f,bound=bound),by=eval(hColNames[[i]])]
      }else{
        dat[,calibWeight:=f*calibWeight,by=eval(hColNames[[i]])]
      }
      setnames(dat,"value",valueH[i])
      curEps <- dat[,max(abs(f-1))]
      if(curEps>epsH){
        error <- TRUE
      }
      if(verbose&&curEps&&calIter%%10==0){
        cat(calIter, ":Not yet converged for H-Constraint",i,"\n")
      }
      if(meanHH){
        dat[,calibWeight:=mean(calibWeight),by=eval(hid)]
      }
    }
    calIter <- calIter + 1
    if(verbose&&!error){
      cat("Convergence reached in ",calIter," steps \n")
    }else if(verbose&&maxIter==calIter){
      cat("Not converged in",maxIter,"steps \n")
    }
  }
  setnames(dat,"baseWeight",w)
  ##Housekeeping
  ###Housekeeping of the varNames used
  delVars <- c(delVars,valueP,valueH,"wvst","wValue","f")
  dat[,eval(delVars):=NULL]
  if(!is.null(renameVars)){
    setnames(dat,paste0(renameVars,"_safekeeping"),renameVars)
  }
  if(any(names(dat)%in%usedVarNames)){
    renameVars <- names(dat)[names(dat)%in%usedVarNames]
    setnames(dat,renameVars,paste0(renameVars,"_safekeeping"))
  }
  invisible(dat)
}



#rm(list=ls())
#source("/Users/alex/git/simpopulation2/R/ipu.R")
#require(simPop);require(reshape2)
#data(eusilcS)
#eusilcS$agecut <- cut(eusilcS$age, 7)
#eusilcS$nat <- sample(LETTERS[1:5],nrow(eusilcS),replace=TRUE)
#eusilcS$emp <- sample(LETTERS[1:5],nrow(eusilcS),replace=TRUE)
#totals1 <- tableWt(eusilcS[, c("rb090","agecut","db040")], weights=eusilcS$rb050)
#totals2 <- tableWt(eusilcS[, c("rb090","emp","db040")], weights=eusilcS$rb050)
#totals3 <- tableWt(eusilcS[, c("nat","db040")], weights=eusilcS$rb050)
#totals1h <- tableWt(eusilcS[!duplicated(eusilcS$db030), c("hsize","db040")], weights=eusilcS$rb050[!duplicated(eusilcS$db030)])
#conP <- list(totals1,totals2,totals3)
#conH <- list(totals1h)
#dat <- data.table(eusilcS)
#cal1 <- ipu2(dat,conP=list(totals1,totals2,totals3),conH=list(totals1h),verbose=TRUE,
#    w="rb050",bound=4,maxIter=200,meanHH=TRUE,hid="db030")
#
#cal2 <- ipu2(dat,conP=list(totals1,totals2,totals3),conH=list(totals1h),verbose=TRUE,
#    w="rb050",bound=4,maxIter=200,meanHH=FALSE,hid="db030")
