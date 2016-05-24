#' EM-based replacement of rounded zeros in compositional data
#' 
#' Parametric replacement of rounded zeros for compositional data using
#' classical and robust methods based on ilr-transformations with special
#' choice of balances.
#' 
#' Statistical analysis of compositional data including zeros runs into
#' problems, because log-ratios cannot be applied.  Usually, rounded zeros are
#' considerer as missing not at random missing values.
#' 
#' The algorithm iteratively imputes parts with rounded zeros whereas in each
#' step (1) an specific ilr transformation is applied (2) tobit regression is
#' applied (3) the rounded zeros are replaced by the expected values (4) the
#' corresponding inverse ilr transformation is applied. After all parts are
#' imputed, the algorithm starts again until the imputations do not change.
#' 
#' @param x data.frame or matrix
#' @param maxit maximum number of iterations
#' @param eps convergency criteria
#' @param method either \dQuote{lm}, \dQuote{MM} or \dQuote{pls}
#' @param dl Detection limit for each variable. zero for variables with
#' variables that have no detection limit problems.
#' @param variation matrix is used to first select number of parts
#' @param nComp if determined, it fixes the number of pls components. If
#' \dQuote{boot}, the number of pls components are estimated using a
#' bootstraped cross validation approach.
#' @param bruteforce sets imputed values above the detection limit to the
#' detection limit. Replacement above the detection limit are only exeptionally
#' occur due to numerical instabilities. The default is FALSE!
#' @param noisemethod adding noise to imputed values. Experimental
#' @param noise TRUE to activate noise (experimental)
#' @param R number of bootstrap samples for the determination of pls
#' components. Only important for method \dQuote{pls}.
#' @param correction normal or density
#' @param verbose additional print output during calculations.
#' @return \item{x }{imputed data} \item{criteria }{change between last and
#' second last iteration} \item{iter }{number of iterations} \item{maxit
#' }{maximum number of iterations} \item{wind}{index of zeros}
#' \item{nComp}{number of components for method pls} \item{method}{chosen
#' method}
#' @author Matthias Templ and Peter Filzmoser
#' @seealso \code{\link{impRZalr}}
#' @keywords manip multivariate
#' @export
#' @importFrom sROC kCDF
#' @importFrom MASS rlm
#' @examples
#' 
#' data(arcticLake)
#' x <- arcticLake
#' ## generate rounded zeros artificially:
#' #x[x[,1] < 5, 1] <- 0
#' x[x[,2] < 44, 2] <- 0
#' xia <- impRZilr(x, dl=c(5,44,0), eps=0.01, method="lm")
#' xia$x
#' 
`impRZilr` <-
  function(x, maxit=10, eps=0.1, method="pls", 
           dl=rep(0.05, ncol(x)), variation=FALSE,	nComp = "boot", 
           bruteforce=FALSE,  noisemethod="residuals", 
           noise=FALSE, R=10, correction="normal",
           verbose=FALSE){

    if( is.vector(x) ) stop("x must be a matrix or data frame")
    ## check if only numeric variables are in x:
    cl <- lapply(x, class)
    if(!all(cl %in% "numeric")) stop("some of your variables are not of class numeric.")
    stopifnot((method %in% c("lm", "MM", "pls", "variation")))
    if( length(dl) < ncol(x)) stop(paste("dl has to be a vector of ", ncol(x)))
    if(method=="pls" & ncol(x)<5) stop("too less variables/parts for method pls")
    if(any(is.na(x))) stop("missing values are not allowed. \n Use impKNNa or impCoda to impute them first.")
    if(is.null(nComp)){
	  pre <- FALSE
	  nC <- NULL
    } else if(nComp=="boot"){
	  nC <- integer(ncol(x))
	  pre <- TRUE
	} else if(length(nComp) == ncol(x)){
	  nC <- nComp
	  pre <- FALSE
	} else  {
	  pre <- FALSE	
	}
  if(!(correction %in% c("normal","density"))) stop("correction method must be normal or density")
#     pre <- TRUE
#      if(length(nComp) != ncol(x) & nComp!="boot") stop("nComp must be NULL, boot or of length ncol(x)")
#    } else if(nComp == "boot"){#
#		pre <- TRUE
#	} else {
#		pre <- FALSE
#	}
    
    #################
    ## store rowSums
    rs <- rowSums(x)
    
    #################
    ## zeros to NA:
    # check if values are in (0, dl[i]):
    check <- logical(ncol(x))
    for(i in 1:ncol(x)){
#      check[i] <- any(x[,i] < dl[i] & x[,i] != 0)
      x[x[,i] < dl[i],i] <- 0
    }
#    if(any(check)){warning("values below detection limit have been set to zero and will be imputed")}
    check2 <- any(x < 0)
    if(check2){warning("values below 0 set have been set to zero and will be imputed")}
    x[x == 0] <- NA
    x[x < 0] <- NA
    indexFinalCheck <- is.na(x)

    ## check if rows consists of only zeros:
    checkRows <- unlist(apply(x, 1, function(x) all(is.na(x))))
    if(any(checkRows)){ 
      w <- which(checkRows)
      cat("\n--------\n")
      message("Rows with only zeros are not allowed")
      message("Remove this rows before running the algorithm")
      cat("\n--------\n")      
      stop(paste("Following rows with only zeros:", w))
    }  

  ## check if cols consists of only zeros:
  checkCols <- unlist(apply(x, 2, function(x) all(is.na(x))))
  if(any(checkCols)){ 
    w <- which(checkCols)
    cat("\n--------\n")
    message("Cols with only zeros are not allowed")
    message("Remove this columns before running the algorithm")
    cat("\n--------\n")      
    stop(paste("Following cols with only zeros:", colnames(x)[w]))
  }  
 
    
    ################
    ## sort variables of x based on 
    ## decreasing number of zeros in the variables
    cn <- colnames(x)
    wcol <- - abs(apply(x, 2, function(x) sum(is.na(x))))
    o <- order(wcol)
    x <- x[,o]
    if(verbose) cat("variables with decreasing number of missings:\n", colnames(x))
    ## --> now work in revised order of variables
    ## dl must also be in correct order
    dlordered <- dl[o]
    
    #################
    ## index of missings / non-missings
    w <- is.na(x)
    wn <- !is.na(x)
#    w2 <- apply(x, 1, function(x){ sum(is.na(x)) })
    #	indNA <- apply(x, 2, function(x){any(is.na(x))})
    
    #################
    ## sort the columns of the data according to the amount of missings in the variables
#    wcol <- apply(x, 2, function(x) length(which(is.na(x))))
#    indM <- sort(wcol, index.return=TRUE, decreasing=TRUE)$ix
    xcheck <- x
    w2 <- is.na(x)

  
    
    ################
    ## initialisation
    indNA <- apply(x, 2, function(x){any(is.na(x))})
    print(indNA)
    for(i in 1:length(dl)){
      ind <- is.na(x[,i])
      #		if(length(ind) > 0) x[ind,i] <- dl[i]*runif(sum(ind),1/3,2/3)
      if(length(ind) > 0) x[ind,i] <- dlordered[i] *2/3
    }
    xOrig <- x

    ################
    ## check if for any variable with zeros,
    ## the detection limit should be larger than 0:
    if(any(dlordered[indNA]==0)){
      w <- which(dlordered[indNA]==0)
      invalidCol <- colnames(x)[w]
      for(i in 1:length(invalidCol)){
        cat("-------\n")
         cat(paste("Error: variable/part", invalidCol[i], 
                           "has detection limit 0 but includes zeros"))
        cat("\n-------\n")
      }
      stop(paste("Set detection limits larger than 0 for variables/part \n including zeros"))
    }

    
    ################
    n <- nrow(x) 
    d <- ncol(x)
    ###  start the iteration
    if(verbose) cat("\n start the iteration:")
    it <- 1; criteria <- 99999999
    while(it <= maxit & criteria >= eps){
      if(verbose) cat("\n iteration", it, "; criteria =", criteria)	
      xold <- x  
      for(i in which(indNA)){
        if(verbose) cat("\n replacement on part", i)
        ## if based on variation matrix:
        if(variation == TRUE){
        orig <- x  
        rv <- variation(x, robust = FALSE)[1,]
        s <- sort(rv)[11]
        cols <- which(rv <= s)[1:11]
        x <- x[, cols]
        }
        ## detection limit in ilr-space
        forphi <- cbind(rep(dlordered[i], n), x[,-i,drop=FALSE])
        if(any(is.na(forphi))) break()
        phi <- isomLR(forphi)[,1] 
        #		part <- cbind(x[,i,drop=FALSE], x[,-i,drop=FALSE])
        x[x < 2*.Machine$double.eps] <- 2*.Machine$double.eps
        xilr <- data.frame(isomLR(cbind(x[,i,drop=FALSE], x[,-i,drop=FALSE])))
        c1 <- colnames(xilr)[1]					
        colnames(xilr)[1] <- "V1"	
        response <- as.matrix(xilr[,1,drop=FALSE])
        predictors <- as.matrix(xilr[,-1,drop=FALSE])
        if(method=="lm"){ 
          reg1 <- lm(response ~ predictors)
          yhat <- predict(reg1, new.data=data.frame(predictors))
        } else if(method=="MM"){
          reg1 <- MASS::rlm(response ~ predictors, method="MM",maxit = 100)#rlm(V1 ~ ., data=xilr2, method="MM",maxit = 100)
          yhat <- predict(reg1, new.data=data.frame(predictors))
        } else if(method=="pls"){
          if(it == 1 & pre){ ## evaluate ncomp.
            nC[i] <- bootnComp(xilr[,!(colnames(xilr) == "V1"),drop=FALSE],y=xilr[,"V1"], R, 
					      plotting=FALSE)$res #$res2
          }
          if(verbose) cat("   ;   ncomp:",nC[i])
          reg1 <- mvr(as.matrix(response) ~ as.matrix(predictors), ncomp=nC[i], method="simpls")
          yhat <- predict(reg1, new.data=data.frame(predictors), ncomp=nC[i])
        }
        
        #		s <- sqrt(sum(reg1$res^2)/abs(nrow(xilr)-ncol(xilr))) ## quick and dirty: abs()
        s <- sqrt(sum(reg1$res^2)/nrow(xilr)) 
        ex <- (phi - yhat)/s 
        if(correction=="normal"){
          yhat2sel <- ifelse(dnorm(ex[w[, i]]) > .Machine$double.eps,
                             yhat[w[, i]] - s*dnorm(ex[w[, i]])/pnorm(ex[w[, i]]),
                             yhat[w[, i]])
        } else if(correction=="density"){
          den <- density(ex[w[,i]])
          distr <- sROC::kCDF(ex[w,i])
        }
        if(any(is.na(yhat)) || any(yhat=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
        # check if we are under the DL:
        if(any(yhat2sel >= phi[w[, i]])){
          yhat2sel <- ifelse(yhat2sel > phi[w[, i]], phi[w[, i]], yhat2sel)
        }
        xilr[w[, i], 1] <- yhat2sel
        xinv <- isomLRinv(xilr)
        ## if variation:
        if(variation == TRUE){
          orig[, cols] <- xinv
          xinv <- orig
        }
        ## reordering of xOrig
        if(i %in% 2:(d-1)){
          xinv <- cbind(xinv[,2:i], xinv[,c(1,(i+1):d)])
        }
        if(i == d){
          xinv <- cbind(xinv[,2:d], xinv[,1])
        }
 #       browser()
        x <- adjustImputed(xinv, xOrig, w2)
        #		x <- adjust3(xinv, xOrig, w2) 
        #		## quick and dirty:
        #		x[!w] <- xOrig[!w]
      }
      
      it <- it + 1
      criteria <- sum( ((xold - x)/x)^2, na.rm=TRUE) ## DIRTY: (na.rm=TRUE)
      if(verbose & criteria != 0) cat("\n iteration", it, "; criteria =", criteria)
    }
    
    #### add random error ###
    if(noise){
      for(i in which(indNA)){
        if(verbose) cat("\n add noise on variable", i)
        
        # add error terms
        inderr <- w[,i]
        if(noisemethod == "residuals") {
          error <- sample(residuals( reg1 )[inderr], 
                          size=wcol[i], replace=TRUE)
          reg1$res[inderr] <- error
        } else {
          mu <- median(residuals( reg1 )[inderr])
          sigma <- mad(residuals( reg1 )[inderr])
          error <- rnorm(wcol[i], mean=mu, sd=sigma)
          reg1$res[inderr] <- error		   
        }
        # return realizations
        yhat[inderr] <- yhat[inderr] + error
        
        
        s <- sqrt(sum(reg1$res^2)/nrow(xilr)) ## quick and dirty: abs()
        ex <- (phi - yhat)/s 
        yhat2sel <- ifelse(dnorm(ex[w[, i]]) > .Machine$double.eps,
                           yhat[w[, i]] - s*dnorm(ex[w[, i]])/pnorm(ex[w[, i]]),
                           yhat[w[, i]])
        if(any(is.na(yhat)) || any(yhat=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
        # check if we are under the DL:
        if(any(yhat2sel >= phi[w[, i]])){
          yhat2sel <- ifelse(yhat2sel > phi[w[, i]], phi[w[, i]], yhat2sel)
        }
        xilr[w[, i], 1] <- yhat2sel
        xinv <- isomLRinv(xilr)
        ## reordering of xOrig
        if(i %in% 2:(d-1)){
          xinv <- cbind(xinv[,2:i], xinv[,c(1,(i+1):d)])
        }
        if(i == d){
          xinv <- cbind(xinv[,2:d], xinv[,1])
        }
        
        #	   x <- adjust2(xinv, xOrig, w) 
        #	   ## quick and dirty:
        #	   x[!w] <- xOrig[!w]
      }
    }
    ### end add random error ###
 #   x <- adjust3(x, xOrig, w)
 #   x[!w] <- xOrig[!w] 
    x <- x[,order(o)] ## checked: reordering is OK!
    colnames(x) <- cn
    ## check if all is fine:
    # check if values are in (0, dl[i]):
    checkDL <- function(x, dl, indexNA){
      check <- logical(ncol(x))
      for(i in 1:ncol(x)){
         check[i] <- any(x[indexNA[,i],i] > dl[i])
         if(check[i]){ 
           x[which(x[indexNA[,i],i] > dl[i]),i] <- dl[i]
         }
      }
      if(any(check)){
        message("few replaced values have been corrected")      
      }
      return(x)
    }
    x <- checkDL(x, dl, indexFinalCheck)
    
  
    
    res <- list(x=x, criteria=criteria, iter=it, 
                maxit=maxit, wind=w, nComp=nC, method=method, dl=dl)
    class(res) <- "replaced"
    invisible(res)
  }

# 
# 
# #' Bootstrap to find optimal number of components
# #' 
# #' Combined bootstrap and cross validation procedure to find optimal number of
# #' PLS components
# #' 
# #' Heavily used internally in function impRZilr.
# #' 
# #' @param X predictors as a matrix
# #' @param y response
# #' @param R number of bootstrap samples
# #' @param plotting if TRUE, a diagnostic plot is drawn for each bootstrap
# #' replicate
# #' @return Including other information in a list, the optimal number of
# #' components
# #' @author Matthias Templ
# #' @seealso \code{\link{impRZilr}}
# #' @keywords manip
# #' @examples
# #' 
# #' ## we refer to impRZilr()
# #' 
# bootnComp <- function(X,y, R=99, plotting=FALSE){
#   ind <- 1:nrow(X)
#   d <- matrix(, ncol=R, nrow=nrow(X))#nrow(X))
#   nc <- integer(R)
#   for(i in 1:R){
#     bootind <- sample(ind)
# #    XX <- X
# #    yy <- y
#     ds <- cbind(X[bootind,], as.numeric(y[bootind]))
#     colnames(ds)[ncol(ds)] <- "V1"
#     res1 <- mvr(V1~., data=data.frame(ds), method="simpls", 
#                   validation="CV")
#     d[1:res1$ncomp, i] <- res1$validation$PRESS
#     nc[i] <- which.min(res1$validation$PRESS)
# #    d[1:reg1$ncomp,i] <- as.numeric(apply(reg1$validation$pred, 3, 
# #                                  function(x) sum(((y - x)^2)) ) )
#   }
#   d <- na.omit(d)
#   sdev <- apply(d, 1, sd, na.rm=TRUE)
#   means <- apply(d, 1, mean, na.rm=TRUE)
#   mi <- which.min(means)
#   r <- round(ncol(X)/20)
#   mi2 <- which.min(means[r:length(means)])+r-1
#   w <- means < min(means) + sdev[mi]
#   threshold <- min(means) + sdev[mi]
#   sdev <- sdev
#   means <- means
#   mi <- mi
#   means2 <- means
#   means2[!w] <- 999999999999999
#   res <- which.min(means2)
#   mi3 <- which.max(w)
# #  minsd <- means - sdev > means[mi]
# #  check <- means
# #  check[!minsd] <- 99999999
#   if(plotting) plot(means, type="l")
# #  res <- which.min(check)
#   list(res3=mi3, res2=mi2, res=res, means=means)
# }
# 
# 
# bootnCompHD <- function(X,y, R=99, plotting=FALSE){
#   ind <- 1:nrow(X)
#   d <- matrix(, ncol=R, nrow=nrow(X))#nrow(X))
#   for(i in 1:R){
#     bootind <- sample(ind)
#     XX <- X
#     yy <- y
#     ds <- cbind(X[bootind,], as.numeric(y[bootind]))
#     colnames(ds)[ncol(ds)] <- "V1"
#     reg1 <- mvr(V1~., data=data.frame(ds), method="simpls", validation="CV")
#     d[1:reg1$ncomp,i] <- as.numeric(apply(reg1$validation$pred, 3, function(x) sum(((y - x)^2)) ) )
#   }
#   d <- na.omit(d)
#   sdev <- apply(d, 1, mad, na.rm=TRUE)
#   means <- apply(d, 1, median, na.rm=TRUE)
#   mi <- which.min(means)
#   if(plotting) plot(means, type="l", col="blue", ylab="squared total prediction error", xlab="number of components")
#   themean <- mean(means)
#   thesd <- sd(means)
#   abovethreshold <- themean - sdev > means
#   check <- means
#   check[!abovethreshold] <- 9999999999
#   res <- which.min(check)
#   #	minsd <- means - sdev > means[mi]
#   #	check <- means
#   #	check[!minsd] <- 99999999
#   #	res <- which.min(check)
#   if(plotting){
#     abline(v=res, lwd=3)
#     abline(h=mi, col="red")
#     #		abline(h=means-sdev, lty=3)
#   }
#   list(res=res, mean=means)
# }
# 

#checkIfValuesUnderDL <- function(x, dl, wind){
#	check <- logical(ncol(x))
#	for(i in 1:ncol(x)){
#		check[i] <- any(x[,i] > dl[i])
#	}	
#	return(check)
#}

## test adjust2:

adjust2 <- function (xImp, xOrig, wind){
  ## aim: do not change original values
  ## adapt imputations
  xneu = xImp
  s1 <- rowSums(xOrig, na.rm = TRUE)
  ## per row: consider rowsums of imputed data
  ## example: 
  ## wind: F F T F F
  ## ganz orig:  3 5 NA 8 10   (sum=26)
  ## orig(init): 3 5 6.5 8 10  (sum=32.5)
  ## imp:        3 5 7 8 10    (sum=33)
  ## s: 26
  ## s2: 7
  ## fac: 26/(26+7)
  ## s1: 32.5/(26/(26+7)) =  41.25
  for (i in 1:nrow(xImp)) {
    if(any(wind[i,])) s <- sum(xImp[i, !wind[i, ]]) else s <- 1
    if(any(wind[i,])) s2 <- sum(xImp[i, wind[i, ]]) else s2 <- 0
    # how much is rowsum increased by imputation:
    fac <- s/(s + s2)
    # decrese rowsums of orig.
    s1[i] <- s1[i]/fac
  }
  ## impS: 41.25/33
  impS <- s1/rowSums(xImp)
  for (i in 1:ncol(xImp)) {
    xneu[, i] <- xImp[, i] * impS
  }
  xImp <- xneu
  return(xImp)
}


adjust3 <- function(xImp, xOrig, wind){
  xOrigSum <- rowSums(xOrig)
  # sum imputed without former zeros:
  xImpSum <- numeric(ncol(xOrig))
  for(i in 1:nrow(xOrig)){
    xImpSum[i] <- sum(xImp[i,!wind[i,]])
    fac <- xOrigSum[i] / xImpSum[i]
    xImp[i,wind[i,]] <- xImp[i,wind[i,]]  * fac
  }
  xImp[!wind] <- xOrig[!wind]
  xImp
}





#
### switch function to automatically select methods
#getM <- function(xReg, ndata, type, index,mixedTF,mixedConstant,factors,step,robust,noise,noise.factor=1,force=FALSE, robMethod="MM") {
#	switch(type,
#			numeric = useLM(xReg, ndata, index, mixedTF,mixedConstant,factors,step,robust,noise,noise.factor,force,robMethod),
#			factor  = useMN(xReg, ndata, index,factors,step,robust),
#			bin     = useB(xReg, ndata, index,factors,step,robust),
#			count   = useGLMcount(xReg, ndata, index, factors, step, robust)
#	)
#}
#
#
#
#f <- function(){	
##	x <- constSum(x)
##	dl <- dl/sum(dl)
#	## initialisation:
##		x[is.na(x)] <- 0.001
#	for(i in 1:length(dl)){
#		ind <- is.na(x[,i])
#		#PF# if(length(ind) > 0) x[ind,i] <- dl[i]/3*2 
#		if(length(ind) > 0) x[ind,i] <- dl[i]*runif(sum(ind),1/3,2/3)
#	}
#	
##		x <- constSum(x)
#	
#	## parameters:
#	it=0
#	criteria <- 10000000
#	error <- rep(0, ncol(x))
#	nComp <- normalerror <- numeric(ncol(x))
#	if(noisemethod=="residuals") residuals <- matrix(,ncol=ncol(x), nrow=nrow(x))
#	nComp[nComp==0] <- NA
#	criteria <- 1e+07
#	sigma <- mu <- rep(0, ncol(x))
#	
#	###########################################
#	###  start the iteration
#	if(verbose) cat("\n start the iteration:")
#	while(it <= maxit & criteria >= eps){
#		xold <- x
#		it=it+1
#		for(i in which(indNA)){
#			if(verbose) cat("\n column", i)
#			## change the first column with that one with the highest amount of NAs
#			## in the step
#			if(wcol[indM[i]] > 0){
#				xNA=x[,indM[i]]
#				x1=x[,1]
#				x[,1]=xNA
#				x[,indM[i]]=x1
#				
#				if(method == "roundedZero"){
#					xilr <- isomLR(x)
#					phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1] # TODO: phi auserhalb der Schleife!
#					## --> x hat sich geaendert aber dl nicht.
#					xilr2 <- data.frame(xilr)
#					c1 <- colnames(xilr2)[1]
#					colnames(xilr2)[1] <- "V1"
#					reg1 = lm(V1 ~ ., data=xilr2)
#					yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
#					if(bruteforce){ 
#						xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
#					} else {
#						s <- sqrt(sum(reg1$res^2)/(nrow(xilr2)-ncol(xilr2)))
#						ex <- (phi - yhat2)/s 
#						yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#								yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]),
#								yhat2[w[, indM[i]]])
#						if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
#						# check if we are under the DL:
#						if(any(yhat2sel >= phi[w[, indM[i]]])){
#							yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#						}
#						xilr2[w[, indM[i]], 1] <- yhat2sel
#					}
#				}
#				if (method == "pls") {	
#					xilr = isomLR(x)
#					ttt <- cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE])
#					phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1]
##				if(verbose) cat("\n phi", phi, "in iteration", it)
#					xilr2 <- data.frame(xilr)	
#					c1 <- colnames(xilr2)[1]					
#					colnames(xilr2)[1] <- "V1"				
#					if(it == 1){ ## evaluate ncomp.
#						test <- xilr2[,!(colnames(xilr2) == "V1")]
#						testy <- xilr2[,"V1"]
#						X <- xilr2[,!(colnames(xilr2) == "V1")]
#						y <- xilr2[,"V1"]
#						nComp[i] <- bootnComp(xilr2[,!(colnames(xilr2) == "V1")],y=xilr2[,"V1"], R)
#					}
#					reg1 <- mvr(V1~.,ncomp=nComp[i], data = xilr2, method="simpls")
#					yhat2 <- as.numeric(predict(reg1, new.data=xilr2[,-i], ncomp=nComp[i]) )
##				if(noisemethod=="residuals") residuals[,i] <- reg1$residuals[,,nComp[i]]
##				if(noisemethod=="normal"){
##					mu[i] <- median(residuals(reg1)[,,nComp[i]])
##					sigma[i] <- mad(residuals(reg1)[,,nComp[i]]) * noiseeffect		  
##				}
#					colnames(xilr2)[1] <- c1					
##				fit=reg1$fitted.values[,,nComp[i]]	
#					s <- sqrt(sum(reg1$residuals[,,nComp[i]]^2)/abs(nrow(xilr2)-ncol(xilr2)))
#					ex <- as.numeric((phi - yhat2)/s )
#					yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#							yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]),
#							yhat2[w[, indM[i]]])
#					if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of NaN or NA estimates")
#					# check if we are under the DL:
#					if(any(yhat2sel >= phi[w[, indM[i]]])){
#						yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#					}
##				if(verbose) cat("\n yhat2sel at iteration", it, "on column", i ,"is", yhat2sel)
#					xilr2[w[, indM[i]], 1] <- yhat2sel
#					
#					
#				}	
#				if(method == "roundedZeroRobust"){
#					xilr <- isomLR(x)
#					x[x < .Machine$double.eps] <- 0.00000000001  ## TODO: better solution 
#					phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1] # TODO: phi auserhalb der Schleife!
#					xilr2 <- data.frame(xilr)
#					c1 <- colnames(xilr2)[1]
#					colnames(xilr2)[1] <- "V1"
#					reg1 = rlm(V1 ~ ., data=xilr2, method="MM",maxit = 100)
##	            reg1 = lmrob(V1 ~ ., data=xilr2)
#					yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
#					if(bruteforce){ 
#						xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
#					} else {
#						s <- reg1$s
#						#PF# s <- IQR(reg1$resid)/1.349
#						ex <- (phi - yhat2)/s 
#						yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#								yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]), 
#								yhat2[w[, indM[i]]])
#						if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
#						# check if we are under the DL:
#						if(any(yhat2sel >= phi[w[, indM[i]]])){
#							yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#						}
#						xilr2[w[, indM[i]], 1] <- yhat2sel
#					}
#				}
#				
#				xilr <- xilr2 
#				x <- isomLRinv(xilr)	
#				
#				## return the order of columns:
#				xNA=x[,1]
#				x1=x[,indM[i]]
#				x[,1]=x1
#				x[,indM[i]]=xNA
#				x <- adjust2(x, xcheck, w)
#				print(x[1:2,1:2])
#				if(verbose) cat("\n iteration", it, "column", i, "value", x[2,1])
#			}
#			
#			
#		}
#		
#		
#		criteria <- sum( ((xold - x)/x)^2, na.rm=TRUE) ## DIRTY: (na.rm=TRUE)
#		colnames(x) <- colnames(xcheck)
#		
#	}
#	
#	res <- list(xOrig=xcheck, xImp=x[,order(o)], criteria=criteria, iter=it,  nComp=nComp,
#			maxit=maxit, w=length(which(w)), wind=w)
#	class(res) <- "imp"
##	res <- adjust(res)
#	invisible(res)
#}
#
#
#
##################################################################################
#
#`impRZilr` <-
#function(x, maxit=10, eps=0.1, method="roundedZero", 
#		dl=rep(0.05, ncol(x)), bruteforce=FALSE,  noise = TRUE, noisemethod="residuals", noiseeffect=1, R=10,
#		verbose=FALSE){
#
#	
#	if( is.vector(x) ) stop("x must be a matrix or data frame")
#	stopifnot((method %in% c("ltsReg", "ltsReg2", "classical", "lm", "roundedZero","roundedZeroRobust","pls")))
#    if( length(dl) < ncol(x)) stop(paste("dl has to be a vector of ", ncol(x)))
#	
#	#################
#	## zeros to NA:
#	x[x==0] <- NA
#
#	#################
#	## index of missings / non-missings
#	w <- is.na(x)
#	wn <- !is.na(x)
#	w2 <- apply(x, 1, function(x){
#          length(which(is.na(x)))
#	})
#    indNA <- apply(x, 2, function(x){any(is.na(x))})
#
#    #################
#	## sort the columns of the data according to the amount of missings in the variables
#	wcol <- apply(x, 2, function(x) length(which(is.na(x))))
#	indM <- sort(wcol, index.return=TRUE, decreasing=TRUE)$ix
#	cn <- colnames(x)
#	xcheck <- x
#	
#	################
#	## detection limit in ilr-space
##	for(i in which(indNA)){
##	   	
##	}
#	
#	
##	x <- constSum(x)
##	dl <- dl/sum(dl)
#	## initialisation:
##		x[is.na(x)] <- 0.001
#	    for(i in 1:length(dl)){
#		   ind <- is.na(x[,i])
#		   #PF# if(length(ind) > 0) x[ind,i] <- dl[i]/3*2 
#		   if(length(ind) > 0) x[ind,i] <- dl[i]*runif(sum(ind),1/3,2/3)
#	    }
#		
##		x <- constSum(x)
#		
#    ## parameters:
#		it=0
#		criteria <- 10000000
#		error <- rep(0, ncol(x))
#		nComp <- normalerror <- numeric(ncol(x))
#		if(noisemethod=="residuals") residuals <- matrix(,ncol=ncol(x), nrow=nrow(x))
#		nComp[nComp==0] <- NA
#		criteria <- 1e+07
#		sigma <- mu <- rep(0, ncol(x))
#		
#	###########################################
#	###  start the iteration
#	if(verbose) cat("\n start the iteration:")
#	while(it <= maxit & criteria >= eps){
#  		xold <- x
#  		it=it+1
#  		for(i in which(indNA)){
#			if(verbose) cat("\n column", i)
#		    ## change the first column with that one with the highest amount of NAs
#		    ## in the step
#			if(wcol[indM[i]] > 0){
#		    xNA=x[,indM[i]]
#		    x1=x[,1]
#		    x[,1]=xNA
#		    x[,indM[i]]=x1
#			
#			if(method == "roundedZero"){
#				xilr <- isomLR(x)
#				phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1] # TODO: phi auserhalb der Schleife!
#				## --> x hat sich geaendert aber dl nicht.
#				xilr2 <- data.frame(xilr)
#				c1 <- colnames(xilr2)[1]
#				colnames(xilr2)[1] <- "V1"
#				reg1 = lm(V1 ~ ., data=xilr2)
#				yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
#				if(bruteforce){ 
#					xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
#				} else {
#					s <- sqrt(sum(reg1$res^2)/(nrow(xilr2)-ncol(xilr2)))
#					ex <- (phi - yhat2)/s 
#                                        yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#                                                           yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]),
#                                                           yhat2[w[, indM[i]]])
#                                        if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
#                                        # check if we are under the DL:
#                                        if(any(yhat2sel >= phi[w[, indM[i]]])){
#						yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#					}
#                                        xilr2[w[, indM[i]], 1] <- yhat2sel
#		        }
#			}
#			if (method == "pls") {	
#				xilr = isomLR(x)
#				ttt <- cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE])
#				phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1]
##				if(verbose) cat("\n phi", phi, "in iteration", it)
#				xilr2 <- data.frame(xilr)	
#				c1 <- colnames(xilr2)[1]					
#				colnames(xilr2)[1] <- "V1"				
#				if(it == 1){ ## evaluate ncomp.
#					test <- xilr2[,!(colnames(xilr2) == "V1")]
#					testy <- xilr2[,"V1"]
#					X <- xilr2[,!(colnames(xilr2) == "V1")]
#					y <- xilr2[,"V1"]
#					nComp[i] <- bootnComp(xilr2[,!(colnames(xilr2) == "V1")],y=xilr2[,"V1"], R)
#				}
#				reg1 <- mvr(V1~.,ncomp=nComp[i], data = xilr2, method="simpls")
#				yhat2 <- as.numeric(predict(reg1, new.data=xilr2[,-i], ncomp=nComp[i]) )
##				if(noisemethod=="residuals") residuals[,i] <- reg1$residuals[,,nComp[i]]
##				if(noisemethod=="normal"){
##					mu[i] <- median(residuals(reg1)[,,nComp[i]])
##					sigma[i] <- mad(residuals(reg1)[,,nComp[i]]) * noiseeffect		  
##				}
#				colnames(xilr2)[1] <- c1					
##				fit=reg1$fitted.values[,,nComp[i]]	
#				s <- sqrt(sum(reg1$residuals[,,nComp[i]]^2)/abs(nrow(xilr2)-ncol(xilr2)))
#				ex <- as.numeric((phi - yhat2)/s )
#				yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#						yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]),
#						yhat2[w[, indM[i]]])
#				if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of NaN or NA estimates")
#				# check if we are under the DL:
#				if(any(yhat2sel >= phi[w[, indM[i]]])){
#					yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#				}
##				if(verbose) cat("\n yhat2sel at iteration", it, "on column", i ,"is", yhat2sel)
#				xilr2[w[, indM[i]], 1] <- yhat2sel
#				
#				
#			}	
#			if(method == "roundedZeroRobust"){
#				xilr <- isomLR(x)
#				x[x < .Machine$double.eps] <- 0.00000000001  ## TODO: better solution 
#				phi <- isomLR(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]))[,1] # TODO: phi auserhalb der Schleife!
#				xilr2 <- data.frame(xilr)
#				c1 <- colnames(xilr2)[1]
#				colnames(xilr2)[1] <- "V1"
#				reg1 = rlm(V1 ~ ., data=xilr2, method="MM",maxit = 100)
##	            reg1 = lmrob(V1 ~ ., data=xilr2)
#				yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
#				if(bruteforce){ 
#					xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
#				} else {
#					s <- reg1$s
#					#PF# s <- IQR(reg1$resid)/1.349
#					ex <- (phi - yhat2)/s 
#					yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#					                   yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]), 
#					                   yhat2[w[, indM[i]]])
#					if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
#					# check if we are under the DL:
#                                        if(any(yhat2sel >= phi[w[, indM[i]]])){
#						yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
#					}
#					xilr2[w[, indM[i]], 1] <- yhat2sel
#				}
#			}
#			
#			xilr <- xilr2 
#			x <- isomLRinv(xilr)	
#
#			## return the order of columns:
#			xNA=x[,1]
#			x1=x[,indM[i]]
#			x[,1]=x1
#			x[,indM[i]]=xNA
#			x <- adjust2(x, xcheck, w)
#			print(x[1:2,1:2])
#			if(verbose) cat("\n iteration", it, "column", i, "value", x[2,1])
#			}
#
#
# 	   }
#
#
#  	  criteria <- sum( ((xold - x)/x)^2, na.rm=TRUE) ## DIRTY: (na.rm=TRUE)
#	  colnames(x) <- colnames(xcheck)
#
#	}
#		
#	res <- list(xOrig=xcheck, xImp=x, criteria=criteria, iter=it,  nComp=nComp,
#			    maxit=maxit, w=length(which(w)), wind=w)
#	class(res) <- "imp"
##	res <- adjust(res)
#	invisible(res)
#}


# `impRZilr` <-
# function(x, maxit=10, eps=0.1, method="roundedZero", 
# 		dl=rep(0.05, ncol(x)), bruteforce=FALSE){
# 
# 	`ilrM` <-
# 			function(x, info=TRUE){
# 		x.ilr=matrix(NA,nrow=nrow(x),ncol=ncol(x)-1)
# 		D=ncol(x)
# 		for (i in 1:ncol(x.ilr)){
# 			x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop=FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
# 		} 
# #		invisible(-x.ilr)
#         if(info)  {res <- list(xilr=-x.ilr,
# 			 xOrig=x)
# 	         class(res) <- "ilrTransform"
# 	    } else {
# 			res <- -x.ilr
# 		}
# 		res
# 	}
# 	`invilrM` <-
# 			function(x.ilr){
# 		if(class(x.ilr) =="ilrTransform" ){
# 			fac <- rowSums(x.ilr$xOrig)
# 			x.ilr <- x.ilr$xilr
# 			y <- matrix(0,nrow=nrow(x.ilr),ncol=ncol(x.ilr)+1)
# 			D=ncol(x.ilr)+1
# 			y[,1]=-sqrt((D-1)/D)*x.ilr[,1]
# 			for (i in 2:ncol(y)){
# 				for (j in 1:(i-1)){
# 					y[,i]=y[,i]+x.ilr[,j]/sqrt((D-j+1)*(D-j))
# 				}
# 			}
# 			for (i in 2:(ncol(y)-1)){
# 				y[,i]=y[,i]-sqrt((D-i)/(D-i+1))*x.ilr[,i]
# 			}
# 			yexp=exp(-y)
# 			x.back=yexp/apply(yexp,1,sum) * fac # * rowSums(derOriginaldaten)
# 			invisible(x.back)			
# 		} else {
# 			y=matrix(0,nrow=nrow(x.ilr),ncol=ncol(x.ilr)+1)
# 			D=ncol(x.ilr)+1
# 			y[,1]=-sqrt((D-1)/D)*x.ilr[,1]
# 			for (i in 2:ncol(y)){
# 				for (j in 1:(i-1)){
# 					y[,i]=y[,i]+x.ilr[,j]/sqrt((D-j+1)*(D-j))
# 				}
# 			}
# 			for (i in 2:(ncol(y)-1)){
# 				y[,i]=y[,i]-sqrt((D-i)/(D-i+1))*x.ilr[,i]
# 			}
# 			yexp=exp(-y)
# 			x.back=yexp/apply(yexp,1,sum) # * rowSums(derOriginaldaten)
# 			invisible(x.back)
#         #return(yexp)
# 		}
# 		x.back
# 	}
# 	
# 	
# 	
# 	if( is.vector(x) ) stop("x must be a matrix or data frame")
# 	stopifnot((method %in% c("ltsReg", "ltsReg2", "classical", "lm", "roundedZero","roundedZeroRobust")))
#     if( length(dl) < ncol(x)) stop(paste("dl has to be a vector of ", ncol(x)))
# 	
# 
# 	## zeros to NA:
# 	x[x==0] <- NA
# 
# 	##index of missings / non-missings
# 	w <- is.na(x)
# 	wn <- !is.na(x)
# 	w2 <- apply(x, 1, function(x){
#           length(which(is.na(x)))
# 	})
# 
# 
# 	##sort the columns of the data according to the amount of missings in the variables
# 	wcol <- apply(x,2,function(x) length(which(is.na(x))))
# 	indM <- sort(wcol, index.return=TRUE, decreasing=TRUE)$ix
# 	cn <- colnames(x)
# 	xcheck <- x
# 
# 	## initialisation:
# #		x[is.na(x)] <- 0.001
# 	    for(i in 1:length(dl)){
# 		   ind <- is.na(x[,i])
# 		   #PF# if(length(ind) > 0) x[ind,i] <- dl[i]/3*2 
# 		   if(length(ind) > 0) x[ind,i] <- dl[i]*runif(sum(ind),1/3,2/3)
# 	    }
# 		
# #		x <- constSum(x)
# 		
#     ## parameters:
# 		it=0
# 		criteria <- 10000000
# 		error <- rep(0, ncol(x))
# 		
# 	###########################################
# 	###  start the iteration
# 	while(it <= maxit & criteria >= eps){
#   		xold <- x
#   		it=it+1
#   		for(i in 1:ncol(x)){
# 		    ## change the first column with that one with the highest amount of NAs
# 		    ## in the step
# 			if(wcol[indM[i]] > 0){
# 		    xNA=x[,indM[i]]
# 		    x1=x[,1]
# 		    x[,1]=xNA
# 		    x[,indM[i]]=x1
# 			
# 			if(method == "roundedZero"){
# 				xilr <- ilrM(x)
# 				phi <- ilrM(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]), info=FALSE)[,1] # TODO: phi auserhalb der Schleife!
# 				## --> x hat sich geaendert aber dl nicht.
# 				xilr2 <- data.frame(xilr$xilr)
# 				c1 <- colnames(xilr2)[1]
# 				colnames(xilr2)[1] <- "V1"
# 				reg1 = lm(V1 ~ ., data=xilr2)
# 				yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
# 				if(bruteforce){ 
# 					xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
# 				} else {
# #					s <- sd(reg1$res, na.rm=TRUE)
# 					s <- sqrt(sum(reg1$res^2)/(nrow(xilr2)-ncol(xilr2)))
# 					ex <- (phi - yhat2)/s 
# #####################################################
# # CHANGED PF:
# #PF#                                    if(all(dnorm(ex) > 5 * .Machine$double.eps)) yhat2 <- yhat2 - s*dnorm(ex)/pnorm(ex)
#                                         yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
#                                                            yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]),
#                                                            yhat2[w[, indM[i]]])
#                                         if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
#                                         # check if we are under the DL:
#                                         if(any(yhat2sel >= phi[w[, indM[i]]])){
# 						yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
# 					}
#                                         xilr2[w[, indM[i]], 1] <- yhat2sel
# #####################################################
# 		        }
# 			}
# 			if(method == "roundedZeroRobust"){
# 				xilr <- ilrM(x)
# 				x[x < .Machine$double.eps] <- 0.00000000001  ## TODO: better solution 
# 				phi <- ilrM(cbind(rep(dl[indM[i]], nrow(x)), x[,-1,drop=FALSE]), info=FALSE)[,1] # TODO: phi auserhalb der Schleife!
# 				xilr2 <- data.frame(xilr$xilr)
# 				c1 <- colnames(xilr2)[1]
# 				colnames(xilr2)[1] <- "V1"
# 				reg1 = rlm(V1 ~ ., data=xilr2, method="MM",maxit = 100)
# #	            reg1 = lmrob(V1 ~ ., data=xilr2)
# 				yhat2 <- predict(reg1, new.data=xilr2[,-i]) 	
# 				if(bruteforce){ 
# 					xilr2[w[, indM[i]], 1] <- ifelse(yhat2[w[, indM[i]]] <= phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2[w[, indM[i]]] )
# 				} else {
# #						s <- mad(reg1$res, na.rm=TRUE)
# 					s <- reg1$s
# 					#PF# s <- IQR(reg1$resid)/1.349
# 					ex <- (phi - yhat2)/s 
# #####################################################
# # CHANGED PF:
# #PF#					if(all(dnorm(ex) > 5 * .Machine$double.eps)) yhat2 <- yhat2 - s*dnorm(ex)/pnorm(ex)
# 					yhat2sel <- ifelse(dnorm(ex[w[, indM[i]]]) > .Machine$double.eps,
# 					                   yhat2[w[, indM[i]]] - s*dnorm(ex[w[, indM[i]]])/pnorm(ex[w[, indM[i]]]), 
# 					                   yhat2[w[, indM[i]]])
# 					if(any(is.na(yhat2)) || any(yhat2=="Inf")) stop("Problems in ilr because of infinite or NA estimates")
# 					# check if we are under the DL:
#                                         if(any(yhat2sel >= phi[w[, indM[i]]])){
# 						yhat2sel <- ifelse(yhat2sel > phi[w[, indM[i]]], phi[w[, indM[i]]], yhat2sel)
# 					}
# 					xilr2[w[, indM[i]], 1] <- yhat2sel
# #####################################################
# 				}
# 			}
# 			
# 			xilr$xilr <- xilr2 
# 			x=invilrM(xilr)		
# 			## return the order of columns:
# 			xNA=x[,1]
# 			x1=x[,indM[i]]
# 			x[,1]=x1
# 			x[,indM[i]]=xNA
# 			}
# 
# 
#  	   }
# 
# 
#   	  criteria <- sum( ((xold - x)/x)^2, na.rm=TRUE) ## DIRTY: (na.rm=TRUE)
# 	  colnames(x) <- colnames(xcheck)
# 
# 	}
# 		
# 	res <- list(xOrig=xcheck, xImp=x, criteria=criteria, iter=it, 
# 			    maxit=maxit, w=length(which(w)), wind=w)
# 	class(res) <- "imp"
# 	invisible(res)
# }
# 
