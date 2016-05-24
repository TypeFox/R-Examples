#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by as.data.frame to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param overid A user-specified degree of overidentification (\code{overid}). See Example 3 below. 
#' @param print.miivs A logical indicating whether or not to display MIIVs in output.
#' @param varcov Option for estimating conditional variance and coavariance paramaters. Default is \code{NULL}.
#' @param bootstrap.se Option \code{"pairs"} or \code{"residual"} for obtaining bootstrap standard errors, t-tests, and bootstrap P values. Default is \code{NULL}.
#' @param cov A named numeric matrix. Default is \code{NULL}.
#' @param means A sample mean vector. Default is \code{NULL}.
#' @param N Numeric. The number of observations in the sample. Default is \code{NULL}.
#'
#' @return model
#' @return dat
#' @return modeeqns
#' 
#' @details 
#' \itemize{
#' \item{\code{overid}} {If the user-specified degree of overidentification (\code{overid}) exceeds the number of available MIIVs for a given equation, the maximum number of MIIVs will be used instead.  In this case, the \code{df} column for any equations in which the degrees of freedom for the \code{Sargan} test are less than the \code{overid} value will be appended with an \code{*}. A note will also  be displayed to alert the user to the changes. In the example below, the \code{overid} parameter is set to 2, however the \code{y1} equation has only 1 additional MIIV avaialable.}
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.  To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation, and not for a specific endogenous variable.}
#' \item{\code{varcov}} {Currently, only \code{"ML"} and \code{"ULS"} fitting functions are supported through the \code{\link[lavaan]{lavaan}} package.}
#' \item{\code{bootstrap.se}} {Currently, only \code{"pairs"} and \code{"residual"} unrestricted bootstrap are implemented using 999 bootstrap repititions.}
#' \item {Equality Constraints} Users can specify equality constraints directly using labels in the model syntax.  See Example 4 below. 
#' }
#' 
#' @references 
#' 
#' Bollen, K. A. 1996.	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen,	K. A. 2001.	Two-stage	Least	Squares	and	Latent	Variable	Models:	
#' Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#'
#' @examples
#' 
#' # Example 1
#' 
#'  bollen1989a_model <- '
#'
#'    Eta1 =~ y1 + y2  + y3  + y4  
#'    Eta2 =~ y5 + y6  + y7  + y8    
#'    Xi1  =~ x1 + x2 + x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8 
#'  '
#'  
#'   miive(model = bollen1989a_model, data = bollen1989a)
#'  
#'  
#' # Example 2
#' 
#'   my_instruments <- ' 
#'    y1 ~ x2 + x3                            
#'    y5 ~ y2 + y3 + y4 + x2                
#'    y2 ~ y3 + y7 + y8 + x2           
#'    y3 ~ y2 + y4 + y6 + y8        
#'    y4 ~ y3 + y6           
#'    y6 ~ y3 + y4 + y7 + x2            
#'    y7 ~ y2 + y4 + y6 + y8       
#'    y8 ~ y2 + y3 + y7 + x2          
#'    x2 ~ y1 + y5 + y2 + y3 + y4 + y6
#'    x3 ~ y1 + y5 + y2 + y3 + y4 + y6
#'  '
#'  
#' miive(model = bollen1989a_model, data = bollen1989a, 
#'       instruments = my_instruments)
#'  
#'  
#' # Example 3
#'  
#'  miive(model = bollen1989a_model, data = bollen1989a, overid = 2)
#'  
#'  
#' # Example 4
#'  bollen1989a_model_r <- '
#'
#'    Eta1 =~ y1 + l2*y2  + l3*y3  + l4*y4  
#'    Eta2 =~ y5 + l2*y6  + l3*y7  + l4*y8    
#'    Xi1  =~ x1 + x2 + l1*x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8
#'    
#'    # Equality Constraints
#'    l1   == 0.5
#'    l2   == l2
#'    l3   == l3
#'    l4   == l4 
#'  '
#'  
#'  miive(model = bollen1989a_model_r, data = bollen1989a)
#'  
#'  
#' @export
miive <- function(model = model, data = NULL, overid = NULL, varcov = NULL, 
                  print.miivs = FALSE, bootstrap.se = NULL, instruments = NULL,
                  cov = NULL, means = NULL, N = NULL){

  if ( "miivs" == class(model) ){ mod <- model; d <- mod$eqns; } 
  if ( "miivs" != class(model) ){ mod <- miivs(model); d <- mod$eqns; } 

  if (is.null(cov)) { covariance = FALSE }
  if (!is.null(cov)){ covariance = TRUE  }

  if (covariance == TRUE & !is.null(bootstrap.se)){
    stop(paste("Bootstrap procedures not supported when using covariance matrix as input."))
  }
  
  if (covariance == TRUE){
    if (is.null(means)){
      stop(paste("Must supply a vector of means when using covariance matrix as input."))
    }
    if (is.null(N)){
      stop(paste("Must supply the number of observations when using covariance matrix as input."))
    }
    means  <- setNames(means, colnames(cov))
  }

  if (covariance == FALSE){
    means  <- colMeans(data) 
    N      <- nrow(data)
    data.c <- apply(data, 2, function(y) y - mean(y)) 
  }

  means <- matrix(means, nrow = 1, ncol = length(means), dimnames = list("", names(means)) )

  if ( length(mod$constr) == 0 ) { restrictions <- FALSE }
  if ( length(mod$constr)  > 0 ) { restrictions <- TRUE  }
  
  if ( !is.null(overid) && !is.null(instruments)){
    stop(paste("Cannot supply both instruments list and overid."))
  }
  
  if (!is.null(instruments)){
    table   <- lavParTable(instruments)
    table   <- table[table$op == "~", ]
    dv_list <- unique(table$lhs)
    iv_list <- list(DV_user = "", IV_user = "")
    iv_list <- replicate(length(dv_list), iv_list, simplify = FALSE)
  
    for (i in 1:length(dv_list)){
      iv_list[[i]]$DV_user <- dv_list[i]
      iv_list[[i]]$IV_user <- c(table[table$lhs == dv_list[i],]$rhs)
    }
  
    dv_user <- unlist(lapply(iv_list, "[", c("DV_user")), use.names=FALSE)
    iv_user <- unlist(lapply(iv_list, "[", c("IV_user")), use.names=FALSE)
  
    eqns_not_to_estimate <- c()
    
    for (i in 1:length(d)){
    
      dv <- d[[i]]$DVobs
      index <- which(dv_user == dv)
    
      # if the equation isn't listed remove it
      if (is.integer(index) && length(index) == 0L) {
        eqns_not_to_estimate <- c(eqns_not_to_estimate, i) 
      }
      
      if (is.integer(index) && length(index) != 0L) {
      iv_user  <- iv_list[[index]]$IV_user
      iv_miivs <- d[[i]]$IV
      n_pred   <- length(d[[i]]$IVobs)
    
      if (length(iv_user) < n_pred) {
        stop(paste("Need at least ", n_pred," instruments for ", dv))}
    
      check <- iv_user[which(!(iv_user%in%iv_miivs))]
    
      if ( !(is.character(check) && length(check) == 0L) ) {  
         stop(paste("Instruments for ", dv," are not valid."))} 
  
      d[[i]]$IV <- iv_user
    
      }
    }
    if (length(eqns_not_to_estimate) > 0){d <- d[-eqns_not_to_estimate]}
  } 


  optim <- function(d, overid, data){
    for (i in 1:length(d)){
      y  <- as.matrix( cbind(data[,d[[i]]$DVobs] ) )
      P  <- d[[i]]$IVobs
      k  <- overid + length(d[[i]]$IVobs) 
      
      if (k > length(d[[i]]$IV) ) {
        d[[i]]$NOTE <- paste("* Maximum number of MIIVs is less than requested degree of overidentification. See df for degree of overidentification.", sep="")
        d[[i]]$MSG <- "*"
      }

      else if (k <= length(d[[i]]$IV) ) {
        cd <- cor(data)
        cd <- cbind(cd[,P])
        cd[cd == 1] <- 0
        cd_m <- as.matrix(cd)
    
        cd_m <- cbind( apply(cd_m, 1, max) ) 
        ord_names <- rownames(cd_m)[order(cd_m, 
                                        decreasing=TRUE)][1:nrow(cd_m)]
    
        temp <- d[[i]]$IV
        temp <- temp[order(match(temp,ord_names))]
        d[[i]]$IV <- temp[1:k]
      }
    }
    return(d)
  } 

  if (!is.null(overid)){ d <- optim(d, overid, data)}
  
  if (restrictions == TRUE){

    for(i in 1:length(d)){
      if (covariance == TRUE){
        gl0 <-  paste(d[[i]]$DVobs, "_",d[[i]]$IVobs, sep="")
      }
      if (covariance == FALSE){
        gl0 <- c(paste(d[[i]]$DVobs, "_Int", sep=""), 
                 paste(d[[i]]$DVobs, "_",d[[i]]$IVobs, sep=""))
      }
      if(i==1)(gl <- gl0)
      if(i >1)(gl <- c(gl, gl0))
    }
    
    con <- mod$constr
    R <- matrix(0, nrow = length(con), ncol = length(gl))
    L <- matrix(0, nrow = length(con), ncol = 1)
    dimnames(R) <- list(unlist(lapply(con, "[[", c("NAME"))), gl)
    dimnames(L) <- list(unlist(lapply(con, "[[", c("NAME"))), c(" "))
    
    for (r in 1:length(con)){
     if (con[[r]]$FIX == 0 ){ 
       R[r, paste(con[[r]]$SET[1],"_",con[[r]]$DV[1], sep="")] <-  1
       R[r, paste(con[[r]]$SET[2],"_",con[[r]]$DV[2], sep="")] <- -1
       L[r] <- 0
     }
  
     if (con[[r]]$FIX != 0 ){
       R[r, paste(con[[r]]$SET,"_",con[[r]]$DV, sep="")] <- 1
       L[r] <- con[[r]]$FIX
     }
    }
  } 
  
  if (restrictions == FALSE){
    R <- NULL
    L <- NULL
  }
  
  if (covariance == FALSE){ 
    for (i in 1:length(d)){
      if (i == 1){
        y1     <- cbind(data[,d[[i]]$DVobs])
        X1     <- as.matrix(cbind(1,data[,d[[i]]$IV]))
        Y1     <- as.matrix(cbind(1,data[,d[[i]]$IVobs]))
      }
      if (i >= 2){
        y1    <- rbind(y1, cbind(data[,d[[i]]$DVobs]))
        X1    <- as.matrix(bdiag(X1, as.matrix(cbind(1,data[,d[[i]]$IV]))))
        Y1    <- as.matrix(bdiag(Y1, as.matrix(cbind(1,data[,d[[i]]$IVobs]))))
      }
      d[[i]]$IVobsInt <- c(paste(d[[i]]$DVobs,"_Int", sep= ""), d[[i]]$IVobs)
    }
    
    y1 <- as(Matrix(y1), "dgeMatrix")
    X1 <- as(Matrix(X1), "sparseMatrix")
    Y1 <- as(Matrix(Y1), "sparseMatrix")
    
    X1inv    <- solve(crossprod(X1)) 
    Y1hat    <- X1%*% (X1inv %*% (t(X1) %*% Y1))
    
    if (restrictions == TRUE){
      b <- solve(rbind(cbind(as.matrix(t(Y1hat) %*% Y1hat), t(R)), cbind(R, matrix(0, 
           ncol=nrow(R), nrow=nrow(R))))) %*% rbind(as.matrix(t(Y1hat) %*% y1), L)
      b <- cbind(b[1:ncol(R)])
      b.unr <- solve(t(Y1hat)%*%Y1) %*% t(Y1hat)%*%y1
      b.unr <- cbind(b.unr[1:ncol(R)])
    }
    if (restrictions == FALSE){
      b <- solve(t(Y1hat)%*%Y1) %*% t(Y1hat)%*%y1
    }
  } 
 
  if (covariance == TRUE){ 
    for (i in 1:length(d)){
      if (i == 1){
        sYX <- as.matrix(cov[d[[i]]$IVobs, d[[i]]$IV,  drop = FALSE])
        sXY <- as.matrix(cov[d[[i]]$IV, d[[i]]$IVobs, drop = FALSE])
        sXX <- as.matrix(cov[d[[i]]$IV, d[[i]]$IV,  drop = FALSE])
        sYY <- as.matrix(cov[d[[i]]$IVobs, d[[i]]$IVobs,  drop = FALSE])
        sXy <- as.matrix(cov[d[[i]]$IV, d[[i]]$DVobs, drop = FALSE])
        syY <- as.matrix(cov[d[[i]]$DVobs, d[[i]]$IVobs, drop = FALSE])
        sy1 <- as.matrix(cov[d[[i]]$DVobs, d[[i]]$DVobs, drop = FALSE])
        sy1X <- as.matrix(cov[d[[i]]$DVobs, d[[i]]$IV, drop = FALSE])
        sXy1 <- as.matrix(cov[d[[i]]$IV, d[[i]]$DVobs, drop = FALSE])
      }
      if (i >= 2){
        sYX <- as.matrix(bdiag(sYX ,as.matrix(cov[d[[i]]$IVobs, d[[i]]$IV,  drop = FALSE])))
        sXY <- as.matrix(bdiag(sXY ,as.matrix(cov[d[[i]]$IV, d[[i]]$IVobs, drop = FALSE])))
        sXX <- as.matrix(bdiag(sXX ,as.matrix(cov[d[[i]]$IV, d[[i]]$IV,  drop = FALSE])))
        sYY <- as.matrix(bdiag(sYY ,as.matrix(cov[d[[i]]$IVobs, d[[i]]$IVobs,  drop = FALSE])))
        sXy <- as.matrix(rbind(sXy ,as.matrix(cov[d[[i]]$IV, d[[i]]$DVobs, drop = FALSE])))
        syY <- as.matrix(bdiag(syY ,as.matrix(cov[d[[i]]$DVobs, d[[i]]$IVobs, drop = FALSE])))
        sy1 <- as.matrix(rbind(sy1, as.matrix(cov[d[[i]]$DVobs, d[[i]]$DVobs, drop = FALSE])))
        sy1X <- as.matrix(bdiag(sy1X ,as.matrix(cov[d[[i]]$DVobs, d[[i]]$IV, drop = FALSE])))
        sXy1 <- as.matrix(bdiag(sXy1 ,as.matrix(cov[d[[i]]$IV, d[[i]]$DVobs, drop = FALSE])))
      }
      
    }
  
    if (restrictions == TRUE){
      left  <- solve(rbind(cbind(sYX %*% solve(sXX) %*% sXY, t(R)), 
               cbind(R, matrix(0, ncol=nrow(R), nrow=nrow(R)))))
      right <-  rbind(sYX %*% solve(sXX) %*% sXy, L)
      b <- left %*% right
      b.unr <- solve(sYX %*% solve(sXX) %*% sXY) %*% sYX %*% solve(sXX) %*% sXy
    }
    
    if (restrictions != TRUE){
      b <- solve(sYX %*% solve(sXX) %*% sXY) %*% sYX %*% solve(sXX) %*% sXy
    }
  }

  dvs <- unlist(lapply(d, function(x) unlist(x$DVobs))) 
  ivs <- unlist(lapply(d, function(x) unlist(x$IVobs)))
  if (covariance == FALSE){
    de <- do.call("rbind",lapply(d, function(x) cbind(x$DVobs, c(paste("Int_",x$DVobs,sep="") ,x$IVobs))))
  }
  if (covariance == TRUE){
    de <- do.call("rbind",lapply(d, function(x) cbind(x$DVobs ,x$IVobs)))
  }
  b <- data.frame(cbind(b[1:nrow(de),], de))
  colnames(b) <- c("b","dv", "iv")
  b$eq <- with(b, ave(as.character(dv), FUN = function(x) cumsum(!duplicated(x))))
  b$b <- as.numeric(as.character(b$b))
  b$id <- 1:nrow(b)
  if (restrictions == TRUE){ 
    b.ur <- b 
    b.ur$b <- b.unr
  }
  dups <- which(as.numeric(ave(paste(b$dv), b$dv, FUN = seq_along))!=1)
  ints <- which(as.numeric(ave(paste(b$dv), b$dv, FUN = seq_along))==1)
  nod <- setdiff(ivs, dvs)
  nod <- nod[!grepl("_Int", nod)]
  all <- c(dvs, nod)
  mx  <- as.matrix(b[,c("iv","dv","b")])

  B <- outer(all, all, 
    function(x, y) {
      mapply(function(x.sub, y.sub) {
        val <- mx[mx[, 1] == x.sub & mx[, 2] == y.sub, 3]
        if(length(val) == 0L) 0 else as.numeric(as.character(val))*-1
      }, x, y)
    } 
  )
  B <- as(B,"Matrix") + Diagonal(nrow(B)); dimnames(B) <- list(all,all)
  
  if (restrictions == TRUE){ 
    mx.ur <- as.matrix(b.ur[,c("iv","dv","b")])
    B.ur <- outer(all, all, 
    function(x, y) {
      mapply(function(x.sub, y.sub) {
        val <- mx.ur[mx.ur[, 1] == x.sub & mx.ur[, 2] == y.sub, 3]
        if(length(val) == 0L) 0 else as.numeric(as.character(val))*-1
      }, x, y)
    } 
  )
  B.ur <- as(B.ur,"Matrix") + Diagonal(nrow(B)) 
  }
  
  if (covariance == FALSE){
    Y         <- as(as.matrix(data)[,all, drop=FALSE], "Matrix")
    Y.c       <- as(apply(Y, 2, function(y) y - mean(y)), "Matrix")
    E         <- Y.c %*% B; colnames(E) <- all
    S.EE      <- (t(B) %*% crossprod(Y.c) %*% B) / N
    V         <- E %*% solve(B); colnames(V) <- all
    V.r       <- V[,unlist(dvs), drop = FALSE]
    S.MI      <- crossprod(V.r) / N
    S.OB      <- as(cov(data[, unlist(dvs), drop = FALSE]), "Matrix")
    S.RS      <- S.OB - S.MI
    u.hat     <- bdiag(split(t(E[,dvs]),1:ncol(E[,dvs, drop = FALSE])))# for sarg.
    ##P.z       <- X1 %*% X1inv %*% t(X1)
    ##srg       <- diag(crossprod(u.hat,P.z) %*% u.hat / (crossprod(u.hat)/ N))
    srg       <- diag(t(u.hat) %*% (X1 %*% (X1inv %*% (t(X1) %*% u.hat))) / (crossprod(u.hat)/ N) )
    srg.df    <- unlist(lapply(d, function(x) length(x$IV) - length(x$IVobs)))
    srg       <- data.frame(srg.df, srg)
    srg$srg.p <- apply(srg, 1, function(x)  1 - pchisq(x[2] ,df = x[1]))
    srg$dv    <- dvs
    b         <- merge(b, srg, all.x=TRUE)
    b         <- b[order(b$id), ]
    b[dups, c("srg","srg.df","srg.p")] <- NA
    
    if (restrictions == TRUE){
      S.EEur <- (t(B.ur) %*% t(Y.c) %*% Y.c %*% B.ur) / N
      I <- Diagonal(N)
      omega <- solve(kronecker(diag(diag(S.EE[unlist(dvs),unlist(dvs)])), I) )
      top <- cbind(as.matrix(t(Y1hat) %*% omega %*% Y1hat), t(R))
      bot <- cbind(R, matrix(0, ncol=nrow(R), nrow=nrow(R)))
      S.EEr <- as.matrix(solve(rbind(top, bot)) )
      S.EEr <- (S.EEr[1:length(gl),1:length(gl)])
      rownames(S.EEr) <- gl
      b$se <- cbind(sqrt(diag(S.EEr)))
      # for test of restrictions only
      dimnames(S.EEur) <- dimnames(S.EE) ## added this temporarily
      omega.ur <- solve(kronecker(diag(diag(S.EEur[unlist(dvs),unlist(dvs)])), I))
      S.EEur.r <- as.matrix(solve(t(Y1hat) %*% omega.ur %*% Y1hat))
    }
    
    if (restrictions == FALSE){
      I     <- Diagonal(N)
      omega <- solve(kronecker(diag(diag(S.EE[unlist(dvs),unlist(dvs), drop = FALSE]),
                                    nrow = length(dvs), ncol = length(dvs)), I))
      S.EEr <- as.matrix(solve(t(Y1hat) %*% omega %*% Y1hat))
      rownames(S.EEr) <- b$iv
      b$se <- cbind(sqrt(diag(S.EEr)))
    }
    
  }

  if (covariance == TRUE){
    S.EE <- t(B) %*% cov[colnames(B), colnames(B)] %*% B
    S.OB <- cov[colnames(B), colnames(B)]
    S.MI <-  t(solve(B)) %*% S.EE %*% solve(B)
    S.RS <- S.OB - S.MI
    se <- cbind(diag(S.EE[dvs,dvs]/N)) ##
    colnames(se) <- "se"
    rownames(se) <- dvs
    b <- merge(b, se, by.x = "dv", by.y = "row.names")
    b <- b[order(b$id), ]
    if (restrictions == TRUE){
      R0 <- matrix(0, ncol=nrow(R), nrow=nrow(R))
      omega <- diag((solve(rbind(cbind(sYX %*% solve(sXX) %*% sXY, t(R)), 
               cbind(R, R0))))[1:length(gl),1:length(gl)])
    }
    if (restrictions == FALSE){
      omega <- diag((solve(sYX %*% solve(sXX) %*% sXY))) 
    }
    b$se <- sqrt(b$se*omega)
    ybar <- t(means); colnames(ybar) <- "ybar"; 
    xbar <- ybar; colnames(xbar) <- "xbar";
    int <- merge(b, xbar, by.x = "iv", by.y = "row.names") 
    int <- int[order(int$id), ];
    int$yb <- int$xbar*int$b
    int <- aggregate(yb~dv+eq,int,sum)
    int <- merge(int, ybar, by.x = "dv", by.y = "row.names") 
    int$b <- int$ybar - int$yb; int$iv <- "Int";
    int$ybar <- int$yb <- NULL; int$se <- NA; int$id <- NA;
    int[order(as.numeric(as.character(int$eq))), ]
    b$id <- 2; int$id <- 1; 
    b <- rbind(b,int)
    b <- b[with(b, order(as.numeric(as.character(eq)), id)), ]
    b$id <- 1:nrow(b)
    b$srg <- NA
    b$srg.df <- NA
    b$srg.p <- NA
  }
  
  restests   <- NULL
  if (restrictions == TRUE & covariance == FALSE){ 
    restlabels <- unlist(lapply(con, "[[", c("NAME")))
    
    #lrtest.est <- N * (log(det(S.EE[dvs,dvs,drop=FALSE])) - log(det(S.EEur[dvs,dvs,drop=FALSE])))
    #lrtest.df <- nrow(R)
    #lrtest.p <- pchisq(lrtest.est, lrtest.df, lower.tail = FALSE)
    #lrtest.lab <- "Likelihood Ratio Test: Asymptotic Chi-squared"
    #lrtest <- list(lrtest.est, lrtest.p, lrtest.df, lrtest.lab)
    lrtest <- NULL
    
    waldtest.est <- t(R %*% b.ur$b - L) %*% solve(R %*% S.EEur.r %*% t(R)) %*% (R %*% b.ur$b - L)
    waldtest.df <- nrow(R)
    waldtest.p <- pchisq(waldtest.est, waldtest.df, lower.tail = FALSE)
    waldtest.lab <- "Wald Test: Asymptotic Chi-squared"
    waldtest <- list(waldtest.est, waldtest.p, waldtest.df, waldtest.lab)
    
    
    restests <- list(lrtest,waldtest, restlabels)
  }
  
  # get p values for coefficients
  b$z <- apply(b[,c("b","se")], 1, function(x) x[1]/x[2])
  b$p <- apply(b[,c("b","se")], 1, function(x) 2*(pnorm(abs(x[1]/x[2]), lower.tail=FALSE)))
  b$iv <- as.character(b$iv)
  b$dv <- as.character(b$dv)
  if (covariance == FALSE) {b[ints, "iv"] <- "Int"} 
  
  if (restrictions == TRUE){
    for (i in 1:length(con)){
      if (con[[i]]$FIX !=0) {
        b[b$dv == con[[i]]$SET & b$iv == con[[i]]$DV, c("se","z","p")]<-NA
      }
    }
  }
  
  for (i in 1:length(d)){
    d[[i]]$EST <- b[b$dv == d[[i]]$DVobs, "b"]
  }
  
  for (i in 1:length(d)){ # i =1
    
    if (d[[i]]$NOTE != "") {
      olddf <- b[b$eq == i & !is.na(b$srg.df),]$srg.df  
      newdf <- paste(olddf,"*",sep="")
      b[b$eq == i & !is.na(b$srg.df),]$srg.df <- newdf
    }
    
    if (covariance == FALSE){
      if (b[b$eq == i & !is.na(b$srg.df),]$srg.df == 0) {
        b[b$eq == i & !is.na(b$srg.df),]$srg    <- NA
        b[b$eq == i & !is.na(b$srg.df),]$srg.p  <- NA
        b[b$eq == i & !is.na(b$srg.df),]$srg.df <- NA
      }
    }
    
    if (!any(is.na(d[[i]]$IVlat))) {
      dvobs <- d[[i]]$DVobs
      for(z in 1:length(d[[i]]$IVlat)){
        ivobs <- d[[i]]$IVobs[z]
        ivlat <- d[[i]]$IVlat[z]
        b[b$dv == dvobs & b$iv == ivobs,]$iv <- ivlat
      }
    }
    
    if (!is.na(d[[i]]$DVlat)) {
      b[b$dv == d[[i]]$DVobs, "dv"] <- d[[i]]$DVlat
    }
  }
  

  dat <- b[,c("dv", "iv", "b", "se", "z", "p","srg", "srg.df", "srg.p")]
  colnames(dat) <- c("DV", "EV", "Estimate", "StdErr", 
                     "Z", "P(|Z|)","Sargan", "df", "P(Chi)")

  modeqns <- mod$df

  lavsyntax <- NULL
  if (!is.null(varcov)){
   
    estimator <- varcov
    
    fit <- lavaan(model, 
                auto.fix.first = TRUE, 
                auto.var = TRUE, 
                auto.cov.lv.x = TRUE)
    
    lavsyntax <- lavExport(fit, export = FALSE)
    
    ls <- strsplit(lavsyntax, "\n") 

    for (i in 1:nrow(dat)){
      dv  <- dat[i,"DV"]
      iv  <- dat[i,"EV"]
      est <- dat[i,"Estimate"]

      for (m in 1:length(ls[[1]])){ 
        line <- ls[[1]][m] 
          if (grepl("=~", line)){
            line <- strsplit(line, "=~")
            
              if (iv == gsub(".*\\*","",line[[1]][1]) & 
                  dv == gsub(".*\\*","",line[[1]][2])){
                ls[[1]][m] <- paste(iv, " =~ ", est, "*", dv,sep="")
              }
          }
          if (grepl(" ~ ", line)){
            line <- strsplit(line, "~")
    
              if (iv == gsub(".*\\*","",line[[1]][1]) & 
                  dv == gsub(".*\\*","",line[[1]][2])){
                  ls[[1]][m] <- paste(iv, " ~ ", est, "*", dv,sep="")
              }
          }
          if (grepl("==", line)){
            ls[[1]][m] <- NA
          }
      }
    }
      
    ls[[1]]  <- ls[[1]][!is.na(ls[[1]])]
    lavsyntax <- paste(unlist(ls), "\n", collapse="")
    lavsyntax <- gsub("NA", "1", lavsyntax)
    if (covariance == FALSE){
      fitcov   <- lavaan(lavsyntax, data = data, 
                  auto.fix.first = TRUE, estimator = estimator,
                  auto.var=TRUE, auto.cov.lv.x = TRUE)
    }
    if (covariance == TRUE){
      fitcov   <- lavaan(lavsyntax, sample.cov = cov, 
                         sample.nobs = N, auto.fix.first = TRUE, 
                         estimator = estimator, auto.var=TRUE, 
                         auto.cov.lv.x = TRUE)
    }
    pe  <-  parameterEstimates(fitcov)
    cpe <-  pe[pe$op == "~~" & pe$lhs != pe$rhs,]
    if (dim(cpe)[1] != 0){
      cpe$x <- paste(cpe$lhs, cpe$op, cpe$rhs, sep=" ")
      cpe <- cpe[,c("x","est", "se", "z", "pvalue")]
      colnames(cpe) <- c( "","Estimate", 
                            "StdErr", "Z", "P(|Z|)")
    }
    if (dim(cpe)[1] == 0){ cpe <- NULL }
   
    vpe <- pe[pe$op == "~~" & pe$lhs == pe$rhs,]
    vpe <- vpe[,c("lhs","est", "se", "z", "pvalue")]
    if (dim(vpe)[1] != 0){
      colnames(vpe) <- c( "","Estimate", "StdErr", "Z", "P(|Z|)")
    }
    if (dim(vpe)[1] == 0){ vpe <- NULL }
      vcov <- list(cov = cpe, var = vpe)
    
  }
  
  if (!is.null(bootstrap.se)){
    dat <- miivboot(d = d, dat = dat, restrictions = restrictions, data = data, 
                    bootstrap.se = bootstrap.se, reps=1000, R = R, L = L, B = B)
  }

  ctrlopts <- list(print.miivs = print.miivs, 
                   restrictions = restrictions, 
                   varcov = varcov, 
                   bootstrap.se =  bootstrap.se,
                   covariance = covariance)
  
  res <- list(model = d, dat = dat, modeqns = modeqns, 
              ctrlopts = ctrlopts, 
              restests = restests, vcov = vcov, 
              lavsyntax = lavsyntax)
  
  class(res) <- "miive"
  res
}