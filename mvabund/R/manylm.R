############################################################################
# print.manylm prints the manylm object in a nice way                      #
# so far this is identical with the print.lm method                        #
# maybe we want to add some more specific things later on                  #
# 05-Jan-2010
###############################################################################
print.manylm <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2,
          quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    invisible(x)   # Invisibly return the manylm object
}

################################################################################
# manylm: Fit a multilinear model for high dimensional data                    #
# the (default) methods coef, residuals, fitted values can be used             #
################################################################################

manylm <- function (formula, data=NULL, subset=NULL, weights=NULL, na.action=options("na.action"), method = "qr", model = FALSE, x = TRUE, y = TRUE, qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset, test="LR", cor.type= "I", shrink.param=NULL, tol=1.0e-5, ...) {
ret.x <- x
ret.y <- y
ret.qr <- qr
cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m  <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
data <- mf <- eval(mf, parent.frame())    # Obtain the model.frame.
#mf <- model.frame(formula)

if (method == "model.frame") { 
    return(mf)
} else if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method), domain = NA)
mt <-  attr(mf, "terms")  # Obtain the model terms.

abundances <- as.matrix(model.response(mf, "numeric"))
if(any(is.na(abundances)) & is.null(na.action)) 
	stop("There are NA values in the response. An 'na.action' is necessary.")
if(any(is.na(abundances)) & any(as.character(na.action)=="na.pass")) 
	stop("There are NA values in the response. 'na.action=na.pass' cannot be used.")

w <- model.weights(mf)
w <- unabund(w)

offset <- model.offset(mf)
offset <- unabund(offset)
N <- NROW(abundances)     # eg number of sites
p <- NCOL(abundances)     # number of organism types
if(p>0) {
if(is.null(colnames(abundances))){
    if (p==1)
     colnames(abundances) <- deparse(attr( mt,"variables")[[2]],
      width.cutoff = 500) else {
       colnames(abundances) <- paste (deparse(attr( mt,"variables")[[2]],
        width.cutoff = 500), 1:p, sep = "")   }
    
} } else stop("A model without response cannot be fitted.")
labAbund<-colnames(abundances)

if (!is.null(offset)) {
    if (length(offset) == 1) 
        offset <- matrix( offset, NROW(abundances), NCOL(abundances))
    else if (NCOL(offset) != p | NROW(offset)!=N  )
        stop(gettextf("dimension of 'offset's is (%s), should equal (%s) (dimension of observations)",
        paste(NROW(offset),"x", NCOL(offset), sep=""), paste(N,"x", p, sep="")), domain = NA)
}

################################### BEGIN Estimation ###################################            
if (is.empty.model(mt)) {
   X <- NULL
   z <- list(coefficients = if (is.matrix(abundances)) matrix(,0, 3) 
                            else numeric(0), residuals = unabund(abundances),
            fitted.values = 0*abundances, weights = w, rank = 0, df.residual = N)
   if (!is.null(offset))
       z$fitted.values <- offset
   if(is.null(w))
       z$hat.X <- matrix(0, ncol=N, nrow=N) 
   else
       z$hat.X <- matrix(0, ncol=sum( w != 0), nrow=sum( w != 0))    
} 
else {
   X <- model.matrix(mt, mf, contrasts )    # Obtain the Designmatrix.
   if(any(is.na(X)) & is.null(na.action)) 
      stop("There are NA values in the independent variables. An 'na.action' is necessary.")
   if(any(is.na(X)) & any(as.character(na.action)=="na.pass")) 
      stop("There are NA values in the independent variables. 'na.action=na.pass' cannot be used.")

   if (N!=nrow(X))
      stop("dimensions of response and independent variables in 'formula' do not match")
   q <- ncol(X)
   labX<-colnames(X)   

    # new codes added here to deal with w
   if (is.null(w)){ 
       ### Fit the multivariate LM.
       z <- manylm.fit(x=X, y=abundances, offset = offset, singular.ok = singular.ok, ...)
       w   <- rep(1, times=N)
       ok <- w!=0
       wIsNull <- 1
   } 
   else {
       if (!is.numeric(w))  stop("'weights' must be a numeric vector")
       if (any(w < 0)) stop("negative 'weights' not allowed")
       ok <- w!=0
       wIsNull <- 0
       if(any(!ok)) {
          abundances  <- abundances[ok,, drop=FALSE]
          X           <- X[ok,, drop=FALSE]
       } 
       if (NCOL(w) == 1) { 
          z <- manylm.wfit(x=X, y=abundances, w=w, offset = offset, 
                           singular.ok = singular.ok, ... ) } 
       else { 
          z <- manylm.multiwfit(x=X, y=abundances, w=w, offset = offset, 
                                singular.ok = singular.ok, ... ) }
   }
   # New codes added for estimating ridge parameter 
   if (cor.type=="shrink") {
      if (is.null(shrink.param)) {
         tX <- matrix(1, NROW(abundances), 1)
         shrink.param <- ridgeParamEst(dat=z$residuals, X=tX, weights=w, only.ridge=TRUE, doPlot=FALSE, tol=tol)$ridgeParameter
      }	 
      # to simplify later computation
      if(shrink.param == 0) cor.type <- "I"
      else if(shrink.param == 1) cor.type <- "R"      
      else if (abs(shrink.param)>1)
          stop("the absolute 'shrink.param' should be between 0 and 1")
   }
   else if (cor.type == "I") shrink.param <- NULL 
   else if (cor.type == "R") {
      if (nrow(abundances)<=ncol(abundances))
          stop("An unstructured correlation matrix should only be used if N>>number of variables.") 
      if(nrow(abundances) < 2*ncol(abundances)) 
          warning("the calculated p-values might be unreliable as the number of cases 
                   is not much larger than the number of variables") 
      shrink.param <- NULL
   }

}
################################### END Estimation #############################
# parameters that are needed for tests in anova.manylm and summary.manylm.
z$test          <- test
z$cor.type      <- cor.type
z$shrink.param  <- shrink.param
z$call          <- cl
z$terms         <- mt
z$data          <- data

z$na.action <- attr(mf, "na.action")
z$offset <- offset
z$contrasts <- attr(X, "contrasts")
z$xlevels <- .getXlevels(mt, mf)
if (model) z$model <- mf
if (ret.x) z$x <- X
if (ret.y) z$y <- abundances
if (!ret.qr) z$qr <- NULL

class(z) <- # c( if(is.matrix(abundances)|is.mvabund(abundances))
c("manylm", "mlm", "lm" )

return(z)

}


############################################################################
## manylm.fit                                                              #
############################################################################
manylm.fit <- function(x, y, offset = NULL, tol=1.0e-010, singular.ok = TRUE, ...) 
{
    y <- as.matrix( unabund(y) )
    if (is.null(n <- nrow(x)))
        stop("'x' must be a matrix")
    if (n == 0)
        stop("0 (non-NA) cases in 'x'")
    p <- ncol(x)
    ny <- ncol(y)
    if (p == 0) {
        return(list(coefficients = numeric(0), residuals = y,
            fitted.values = 0 * y, rank = 0, df.residual = NROW(y)))
    }

    if (!is.null(offset)) y<- y - offset
    if (nrow(y) != n)
        stop("dimensions of 'x' and 'y' do not match")

    if (length(list(...))){
	  tmp <- paste(names(list(...)), sep = ", ")
        warning("extra arguments ", tmp ,
            " are just disregarded.")
    }
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    ############ BEGIN Fit multivariate linear model  #######
    results     <- list( )
    results$assign    <- attr(x, "assign")  
    qrx         <- qr(x)
    rank        <- qrx$rank
    p1          <- 1:rank

    if (rank < p & !singular.ok ) 
       stop("singular fit encountered") 
    x.ind    <- x[,qrx$pivot[p1], drop=FALSE]
    xnames   <- dimnames(x)[[2]]
    if (is.null(xnames))
       xnames <- paste("x", 1:p, sep = "")
    colnames(qrx$qr) <- xnames[qrx$pivot]

    coefficients <- matrix(nrow=p,ncol=ny)
    coefficients[qrx$pivot[p1],] <- chol2inv(qrx$qr[p1,p1, drop=FALSE])%*%t(x.ind)%*% y

    dimnames(coefficients) <- list( xnames,dimnames(y)[[2]] )
    hat.X <- x.ind %*% chol2inv(qrx$qr[p1,p1, drop=FALSE]) %*% t(x.ind)
    fitted.values <- x.ind %*% chol2inv(qrx$qr[p1,p1, drop=FALSE])%*% t(x.ind)%*% y
     # less efficient than x%*%coefficients, but more stable.
    residuals     <- y - fitted.values
    colnames(residuals) <- colnames(y)
    df.residual   <- n - rank

    ########### END Fit multivariate linear model ############
    results$coefficients        <- coefficients
    results$residuals           <- residuals
    results$fitted.values       <- fitted.values
    results$rank                <- rank
    results$qr                  <- qrx
    results$df.residual         <- df.residual
    results$hat.X               <- hat.X  
    results$txX                 <- t(x.ind)%*% x.ind
    results$data                <- data

    # To add: results$effects
    # (for not null fits) 'n' vector of orthogonal single-df effects.
    # The first 'rank' of them correspond to non-aliased
    # coeffcients, and are named accordingly.
    return(results)
}
                
                
################################################################################
## manylm.wfit                                                             #
################################################################################

manylm.wfit<-function(x, y, w, offset = NULL, tol=1.0e-010,  singular.ok = TRUE, ...)
{
   y <- as.matrix(unabund(y))
   if (is.null(n <- nrow(x)))
      stop("'x' must be a matrix")
   if (n == 0)
      stop("0 (non-NA) cases in 'x'")
   p <- ncol(x)
   if (p == 0) {
      return(list(coefficients = numeric(0), residuals = y,
          fitted.values = 0 * y, rank = 0, df.residual = NROW(y)))
    }
    ny <- ncol(y)

    if (!is.null(offset)) y <- y - offset
    if (nrow(y) != n)
        stop("dimensions of 'x' and 'y' do not match")

    if (length(list(...))){
	  tmp <- paste(names(list(...)), sep = ", ")
        warning("extra arguments ", tmp ,
            " are just disregarded.")
    }	
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    n.ok    <- n - sum(w == 0)
    if(n.ok==0)
        stop("No informative observations. All 'weights' are equal to Zero.")

    pivot <- qrx$pivot

    ############# BEGIN Fit multivariate linear model  ##################
    results     <- list()
    results$assign <- attr(x, "assign")
    xnames      <- colnames(x)
    if (is.null(xnames))
        xnames <- paste("x", 1:p, sep = "")

    ok    <- w != 0
    qrx   <- qr(x[ok,,drop=FALSE] * sqrt(w[ok]))  # Exclude cases with weights=0.
    colnames(qrx$qr) <- xnames[qrx$pivot]

    rank  <- qrx$rank
    p1    <- 1:rank

    if (rank < p & !singular.ok ) 
      stop("singular fit encountered") 

    x.ind <- x[ok,qrx$pivot[p1], drop=FALSE] 

    coefficients  <- matrix(nrow=p,ncol=ny)
    coefficients[qrx$pivot[p1],]    <- 
    chol2inv(qrx$qr[p1,p1, drop=FALSE])%*%t(x.ind)%*% (y*w)[ok,,drop=FALSE]
    dimnames(coefficients) <- list(xnames,dimnames(y)[[2]])
    hat.X  <- x.ind %*% chol2inv(qrx$qr[p1,p1, drop=FALSE])%*% t(x.ind*w[ok])

    fitted.values <- x %*% ( chol2inv(qrx$qr[p1,p1, drop=FALSE])%*% t(x.ind) %*%
    (y*w)[ok,,drop=FALSE] ) # less efficient than x%*%coefficients, but more stable.

    residuals      <- y-fitted.values
    # Ensure numerical stability for weights=0 cases.
    residuals[!ok] <- y[!ok]
    colnames(residuals) <- colnames(y)

    df.residual   <- n.ok-rank
    ############## END Fit multivariate linear model ################

    results$coefficients    <- coefficients
    results$residuals       <- residuals
    results$fitted.values   <- fitted.values
    results$weights         <- w
    results$rank            <- rank
    results$qr              <- qrx # QR decomposition of x * sqrt(w)if (rank < p & !singular.ok ) 
      stop("singular fit encountered") 
  
    x.ind <- x[,pivot[p1], drop=FALSE]

    coefficients  <- matrix(nrow=p,ncol=ny)
    txy      <-  t( x.ind ) %*% ( y*w )
    chol2invqrx  <-  list()
    hat.X  <- txX <- list()

for( i in 1:(ncol(w)) ) {
   qrx[[i]]         <- qr( x * sqrt(w[,i]) )
   chol2invqrx[[i]] <- chol2inv( qrx[[i]]$qr[p1,p1, drop=FALSE] )
   coefficients[ pivot[p1], ] <-  chol2invqrx[[i]] %*% txy
   txX[[i]] <- t(x.ind*w[,i]) %*% x.ind
}

dimnames(coefficients) <- list(xnames,dimnames(y)[[2]])


fitted.values <- x %*% coefficients

residuals     <- y - fitted.values
# Ensure the numerical stability for weights=0 cases.
residuals[ok] <- y[ok]
colnames(residuals) <- colnames(y)

df.residual   <- n.ok - rank

################## END Fit multivariate linear model ##################

results$coefficients        <- coefficients
results$residuals           <- residuals
results$fitted.values       <- fitted.values
results$weights             <- w
results$rank                <- rank
results$qr                  <- qrx # QR decomposition of x * sqrt(w)
results$df.residual         <- df.residual

return(results)
    results$df.residual     <- df.residual
    results$hat.X           <- hat.X # weighted, cases with weights=0 are not included
    results$txX             <- t(x.ind*w[ok])%*% x.ind  # weighted
    # To add: results$effects
    # (for not null fits) 'n' vector of orthogonal single-df effects.
    # The first 'rank' of them correspond to non-aliased
    # coeffcients, and are named accordingly.
    return(results)

}

################################################################################
## manylm.multiwfit                                                           ## 
################################################################################
manylm.multiwfit <- function(x, y, w, offset = NULL, tol=1.0e-010, singular.ok = TRUE, ...)
{
   y <- as.matrix(unabund(y))
   w <- as.matrix(w)
   if (is.null(n <- nrow(x)))
        stop("'x' must be a matrix")
   if (n == 0)
        stop("0 (non-NA) cases in 'x'")
   p <- ncol(x)
 
   if (p == 0) {
      return(list(coefficients = numeric(0), residuals = y,
          fitted.values = 0 * y, rank = 0, df.residual = NROW(y)))
    }
    
    ny <- ncol(y)

    if (!is.null(offset)) y <- y - offset
    if (nrow(y) != n)
        stop("dimensions of 'x' and 'y' do not match")

    if( ncol(y) != ncol(w) )  stop("dimensions of 'y' and 'w' do not match")

    if (length(list(...))){
	  tmp <- paste(names(list(...)), sep = ", ")
        warning("extra arguments ", tmp ,
            " are just disregarded.")
	}

    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    ok      <- w == 0
    n.ok    <- n - colSums( w == 0 )
    if( any(n.ok == 0) )
      stop("No informative observations in the weights for variable ",
    which(n.ok==0), ". All 'weights' are equal to Zero.")
    
    ################### BEGIN Fit multivariate linear model  ##################
    ####### does the same as model.fit, but I can return txX and hat.X ########

    results         <- list()
    results$assign  <- attr(x, "assign")
    xnames          <- colnames(x)
    if(is.null( xnames ))
         xnames <- paste("x", 1:p, sep = "")

    qrx <- list()
    qrx[[1]]  <- qr( x * sqrt(w[,1]) )

    pivot <- qrx[[1]]$pivot
    rank  <- qrx[[1]]$rank
    p1    <- 1:rank

    if (rank < p & !singular.ok ) 
    stop("singular fit encountered") 
	
    x.ind <- x[,pivot[p1], drop=FALSE]
	
    coefficients  <- matrix(nrow=p,ncol=ny)
    txy      <-  t( x.ind ) %*% ( y*w )
     chol2invqrx  <-  list()
	
	hat.X  <- txX <- list()
	
	for( i in 1:(ncol(w)) ) {
	qrx[[i]]         <- qr( x * sqrt(w[,i]) )
	chol2invqrx[[i]] <- chol2inv( qrx[[i]]$qr[p1,p1, drop=FALSE] )
	coefficients[ pivot[p1], ] <-  chol2invqrx[[i]] %*% txy
	txX[[i]] <- t(x.ind*w[,i]) %*% x.ind
	}
	
	dimnames(coefficients) <- list(xnames,dimnames(y)[[2]])
	
	
	fitted.values <- x %*% coefficients
	
	residuals     <- y - fitted.values
	# Ensure the numerical stability for weights=0 cases.
	residuals[ok] <- y[ok]
	colnames(residuals) <- colnames(y)
	
	df.residual   <- n.ok - rank
	
	################## END Fit multivariate linear model ##################
	
	results$coefficients        <- coefficients
	results$residuals           <- residuals
	results$fitted.values       <- fitted.values
	results$weights             <- w
	results$rank                <- rank
	results$qr                  <- qrx # QR decomposition of x * sqrt(w)
	results$df.residual         <- df.residual
	
	return(results)
}
