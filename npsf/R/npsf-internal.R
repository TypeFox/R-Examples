.prepareYX <- function(formula, data, subset, rts = c("C", "NI", "V"),
                       base = c("output", "input"), ref = NULL,
                       data.ref = NULL, subset.ref = NULL, print.level = 1, 
                       type = "RM", winw = 50, rts.subst = NULL, ...)
{
 # get y and x matrices
 
 # mf0 <- match.call(expand.dots = FALSE, call = sys.call(which = 1))
 needed.frame <- sys.nframe() - 1
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 
 # if(print.level >= 1) print(mf0)
 
 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 # subsetsupplied <- !(match("subset", names(mf0), 0) == 0)
 # print(subsetsupplied)
 # cat("1\n")
 
 # if data are supplied
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing
  
  # first using data and subset to get x without NA
  mf <- mf0
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf$formula <- Formula( formula )
  # mf <- eval(mf, parent.frame(n = 1))
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  x <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(x))
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  y <- as.matrix(model.part(Formula(formula), data = data[esample,], lhs = 1))
  x <- as.matrix(model.matrix(Formula(formula), data = data[esample,], rhs = 1)[,-1])
  
  # print(y)
  # print(x)
 }
 # if data are not supplied
 else {
  mf <- mf0
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  y <- as.matrix(model.response(mf))
  x <- as.matrix(model.matrix(mt, mf)[,-1])
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(x)
  # print(table(esample))
  # print(y)
  # print(x)
 }
 
 if( !is.numeric(y) ) stop("Some of the outputs are not numeric.")
 if( !is.numeric(x) ) stop("Some of the inputs are not numeric.")
 
 rts <- rts[1]
 base <- base[1]
 
 rts1 <- tolower(substr(rts, 1,1 ))
 base1 <- tolower(substr(base, 1,1 ))
 
 if(is.null(rts.subst)){
  if(rts1 == "c" | rts == 1){
   myrts <- 3
   myrts1 <- "CRS"
  } else if (rts1 == "n" | rts == 2){
   myrts <- 2
   myrts1 <- "NIRS"
  } else if (rts1 == "v" | rts == 3){
   myrts <- 1
   myrts1 <- "VRS"
  } else {
   stop("invalid 'rts'; 'rts' must be 'CRS', 'NIRS', or 'VRS'")
  }
 }
 else {
  myrts1 <- rts.subst
 }
 
 if (base1 == "o" | base == 2){
  mybase <- 2
  mybase1 <- "output"
 } else if (base1 == "i" | base == 1){
  mybase <- 1
  mybase1 <- "input"
 } else {
  stop("invalid 'base'; 'base' must be 'output' or 'input'")
 }
 
 # for printing
 if(print.level >= 1 & winw > 50){
  mymesage <- paste("\n",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency under assumption of ",myrts1," technology are computed for the following data:\n", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  # 		cat("\n ",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency \n", sep = "")
  # 		cat(" under assumption of ",myrts1," technology are computed for the \n", sep = "")
  # 		cat(" following data:\n\n", sep = "")
  cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
  cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
  cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")
 }
 
 # get y_ref and x_ref matrices
 
 if(!is.null(ref)){
  mf <- mf0
  
  # check if it is a matrix
  datasupplied <- !(match("data.ref", names(mf), 0) == 0)
  
  # if data are supplied
  if(datasupplied){
   N_all_ref <- nrow(data.ref)
   if(N_all_ref == 0) warning("Provided data for reference set does not have a signle data point", call. = FALSE)
   
   # begin get a logical vector equal TRUE if !missing
   
   # first using data and subset to get x without NA
   mf <- mf0
   m <- match(c("ref", "data.ref", "subset.ref"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   # change the names for eval
   names(mf)[which(names(mf) == "ref")] <- "formula"
   names(mf)[which(names(mf) == "data.ref")] <- "data"
   names(mf)[which(names(mf) == "subset.ref")] <- "subset"
   mf[[1L]] <- as.name("model.frame")
   # mf$formula <- Formula( formula )
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   x_ref <- as.matrix(model.matrix(mt, mf))
   if(length(x_ref) == 0){
    warning(" Given 'subset' reference set contains zero data points", call. = FALSE)
    # warning(" Reference set will be based on observations,")
    # warning(" for which efficiency is to be computed")
    y_ref <- y
    x_ref <- x
    esample_ref <- NULL
   } else {
    # now get the names in the entire data
    esample_ref <- seq_len( N_all_ref ) %in% as.numeric(rownames(x_ref))
    # print(table(esample_ref))
    # end get a logical vector equal TRUE if !missing
    y_ref <- as.matrix(model.part(Formula(ref), data = data.ref[esample_ref,], lhs = 1))
    x_ref <- as.matrix(model.matrix(Formula(ref), data = data.ref[esample_ref,], rhs = 1)[,-1])
   }
   # print(y_ref)
   # print(x_ref)
  }
  # if data are not supplied
  else {
   mf <- mf0
   subsetsupplied <- !(match("subset.ref", names(mf), 0) == 0)
   if(subsetsupplied) stop("Subset for reference set with matrices is not allowed. \n   Or: \nOption 'data.ref' must be provided if option 'ref' is provided", call. = FALSE)
   m <- match(c("ref", "data.ref", "subset.ref"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   names(mf)[which(names(mf) == "ref")] <- "formula"
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   mt <- attr(mf, "terms")
   y_ref <- as.matrix(model.response(mf))
   x_ref <- as.matrix(model.matrix(mt, mf)[,-1])
   # print(y_ref)
   # print(x_ref)
   # get a logical vector equal TRUE if !missing
   with.na <- model.frame(mt, na.action = na.pass)
   esample_ref <- rownames(with.na) %in% rownames(x_ref)
   # print(table(esample_ref))
  }
  # if reference set not provided
  # for printing
  # 		if(print.level >= 1){
  # 			if(!is.null(esample_ref)){
  # 				cat("\n Reference set is composed as follows:\n\n", sep = "")
  # 				cat(" Number of observations (K) = ",nrow(y_ref),"\n", sep = "")
  # 				cat(" Number of outputs      (M) = ",ncol(y_ref),"\n", sep = "")
  # 				cat(" Number of inputs       (N) = ",ncol(x_ref),"\n\n", sep = "")
  # 			}
  # 		}
  if(!is.null(y_ref)){
   if(ncol(y_ref) != ncol(y)) stop("Number of outputs for data points and reference set must be the same.")
   if(ncol(x_ref) != ncol(x)) stop("Number of inputs for data points and reference set must be the same.")
  }		
 } else {
  y_ref <- y
  x_ref <- x
  esample_ref <- NULL
 }
 
 if(print.level >= 1  & winw > 50){
  if(is.null(esample_ref)){
   mymesage <- paste("\nData for reference set are not provided. Reference set is formed by ",nrow(y)," data ",ngettext(nrow(y), "point", "point(s)"), " for which measures of technical efficiency are computed", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   #     cat("\n Data for reference set are not provided. Reference set is formed by ",nrow(y),"\n", sep = "")
   #     cat(" data points, for which measures of technical efficiency are computed.\n\n", sep = "")
  }
  else {
   mymesage <- paste("\nReference set is formed by ",nrow(y_ref)," provided reference data ",ngettext(nrow(y_ref), "point", "point(s)"), "", sep = "")
   cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
   # cat("\n Reference set is formed by ",nrow(y_ref)," provided reference data points.\n\n", sep = "")
  }
 }
 
 tymch <- list(y = y, x = x, y.ref = y_ref, x.ref = x_ref, esample = esample, esample.ref = esample_ref, myrts = myrts, rts.string = myrts1, mybase = mybase, base.string = mybase1)
 class(tymch) <- "npsf"
 return(tymch)
}

.prepareYXnoRef <- function(formula, data, subset, 
                            rts = c("C", "NI", "V"), rts.subst = NULL,
                            base = c("output", "input"), print.level = 1, 
                            type = "RM", winw = 50, ...)
{
 # get y and x matrices
 
 # mf0 <- match.call(expand.dots = FALSE, call = sys.call(which = 1))
 needed.frame <- sys.nframe() - 1
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 
 # if(print.level >= 1) print(mf0)
 
 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 # subsetsupplied <- !(match("subset", names(mf0), 0) == 0)
 # print(subsetsupplied)
 # cat("1\n")
 
 # if data are supplied
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing
  
  # first using data and subset to get x without NA
  mf <- mf0
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf$formula <- Formula( formula )
  # mf <- eval(mf, parent.frame(n = 1))
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  x <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(x))
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  y <- as.matrix(model.part(Formula(formula), data = data[esample,], lhs = 1))
  x <- as.matrix(model.matrix(Formula(formula), data = data[esample,], rhs = 1)[,-1])
  
  # print(y)
  # print(x)
 }
 # if data are not supplied
 else {
  mf <- mf0
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  y <- as.matrix(model.response(mf))
  x <- as.matrix(model.matrix(mt, mf)[,-1])
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(x)
  # print(table(esample))
  # print(y)
  # print(x)
 }
 
 if( !is.numeric(y) ) stop("Some of the outputs are not numeric.")
 if( !is.numeric(x) ) stop("Some of the inputs are not numeric.")
 
 rts <- rts[1]
 base <- base[1]
 
 rts1 <- tolower(substr(rts, 1,1 ))
 base1 <- tolower(substr(base, 1,1 ))
 
 if(is.null(rts.subst)){
  if(rts1 == "c" | rts == 1){
   myrts <- 3
   myrts1 <- "CRS"
  } else if (rts1 == "n" | rts == 2){
   myrts <- 2
   myrts1 <- "NIRS"
  } else if (rts1 == "v" | rts == 3){
   myrts <- 1
   myrts1 <- "VRS"
  } else {
   stop("invalid 'rts'; 'rts' must be 'CRS', 'NIRS', or 'VRS'")
  }
 }
 else {
  myrts1 <- rts.subst
 }
 
 if (base1 == "o" | base == 2){
  mybase <- 2
  mybase1 <- "output"
 } else if (base1 == "i" | base == 1){
  mybase <- 1
  mybase1 <- "input"
 } else {
  stop("invalid 'base'; 'base' must be 'output' or 'input'")
 } 
 
 # for printing
 #  if(print.level >= 1){
 #   cat("\n Radial (Debreu-Farrell) ",mybase1,"-based measures of technical efficiency \n", sep = "")
 #   cat(" under assumption of CRS, NIRS, and VRS technology are computed for the\n", sep = "")
 #   cat(" following data:\n\n", sep = "")
 #   cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
 #   cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
 #   cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")
 #  }
 
 if(print.level >= 1 & winw > 50){
  mymesage <- paste("\n",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency under assumption of ",myrts1," technology are computed for the following data:\n", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
  # 		cat("\n ",ifelse(type == "RM", "Nonradial (Russell)", "Radial Debrue-Farrell")," ",mybase1,"-based measures of technical efficiency \n", sep = "")
  # 		cat(" under assumption of ",myrts1," technology are computed for the \n", sep = "")
  # 		cat(" following data:\n\n", sep = "")
  cat("  Number of data points (K) = ",nrow(y),"\n", sep = "")
  cat("  Number of outputs     (M) = ",ncol(y),"\n", sep = "")
  cat("  Number of inputs      (N) = ",ncol(x),"\n", sep = "")
  
  mymesage <- paste("\nReference set is formed by ",nrow(y)," data point(s), for which measures of technical efficiency are computed", sep = "")
  cat("",unlist(strsplit(mymesage, " ")),"", sep = " ", fill = winw-10 )
 }
 
 #  if(print.level >= 1){
 #    cat("\n Reference set is formed by ",nrow(y)," data points,\n", sep = "")
 #    cat(" for which measures of technical efficiency are computed.\n\n", sep = "")
 #  }
 
 tymch <- list(y = y, x = x, esample = esample, mybase = mybase, base.string = mybase1)
 if(is.null(rts.subst)){
  tymch$myrts <- myrts
  tymch$rts.string <- myrts1
 }
 class(tymch) <- "npsf"
 return(tymch)
}

.su <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = FALSE){
 
 xvec2 <- xvec1 <- FALSE
 
 if(is.matrix(x)){
  if(min(dim(x)) == 1){
   # xvec1 <- TRUE
   if(which(dim(x) == 1) == 2){
    mynames <- colnames(x)
   } else {
    mynames <- rownames(x)
    x <- t(x)
   }
   # x <- as.vector(x)	
  } else {
   if(!mat.var.in.col){
    x <- t(x)
   }
   mynames <- colnames(x)
  }
  # print(mynames)
  if(is.null(mynames)) mynames <- paste("Var", seq_len(ncol(x)), sep = "")
  x <- as.data.frame(x)
  # print(x)
  # mynames <- colnames(x)
 } # end if matrix
 
 if(is.vector(x)){
  xvec2 <- TRUE
  mynames <- deparse(substitute(x))
  x <- data.frame(Var1 = x)	
 } # end if vector
 
 # cat("nymanes", sep ="")
 # print(mynames)
 
 if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
  stop("Provide vector, matrix, or data.frame")
 } else {
  t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
  # print(t1)
  # print(mynames)
  # print(class(t1))
  # print(dim(t1))
  if(xvec2 & !xvec1) colnames(t1) <- mynames
  if(print) print(t(t1), digits = digits)		
 }
 tymch <- t(t1)
 class(tymch) <- "npsf"
 return(tymch)
}

.teRad <- function(Y,X,M,N,K,
                   Yr,Xr,Kref,
                   rts,base,ifqh,
                   print.level=0){
 .C("radial",
    as.double(Y), 
    as.double(X), 
    as.integer(M), 
    as.integer(N), 
    as.integer(K), 
    as.double(Yr), 
    as.double(Xr), 
    as.integer(Kref), 
    as.integer(rts), 
    as.integer(base),
    as.integer(ifqh),
    as.integer(print.level),
    te = double(K) )$te
}

.teNonrad <- function(Y,X,M,N,K,
                      Yr,Xr,Kref,
                      rts,base,ifqh,
                      print.level=0){
 .C("nonradial",
    as.double(Y), 
    as.double(X), 
    as.integer(M), 
    as.integer(N), 
    as.integer(K), 
    as.double(Yr), 
    as.double(Xr), 
    as.integer(Kref), 
    as.integer(rts), 
    as.integer(base),
    as.integer(ifqh),
    as.integer(print.level),
    te = double(K) )$te
}

.dots <- function(nrep, message = NULL, width = 50, character = "."){
 if (!is.numeric(width)) {
  stop("'width' should be numeric")
 }
 if (width != 50 & width != 60 & width != 70 & width != 80 & width != 90  & width != 100){
  stop("'width' should be 50, 60, 70, 80, 90, or 100")
 }
 if (nrep == 0){
  if (!is.null(message)){
   cat("",message,"\n", sep = "")
  }
  cat("____|___ 1 ___|___ 2 ___|___ 3 ___|___ 4 ___|___ 5", sep = "")
  if (width == 60) {
   cat(" ___|___ 6", sep = "")
  }
  if (width == 70) {
   cat(" ___|___ 6 ___|___ 7", sep = "")
  }
  if (width == 80) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8", sep = "")
  }
  if (width == 90) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9", sep = "")
  }
  if (width == 100) {
   cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9 ___|___ 10", sep = "")
  }  
  cat("\n")
 }
 else {
  if (nrep/width != floor(nrep/width)){
   cat("",character,"", sep = "")
  }
  else {
   cat("",character," ",nrep,"\n", sep = "")
  }
 }
}

.smplHom <- function(teRef,terfl,Kr,mybw,scVarHom,YorX,ba){
 # sample with replacement from Farrell measures
 # (including reflected ones)
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 bStar     <- terfl[newSample]
 # generate the sequence randomly
 epsStar   <- rnorm(Kr)
 tStar     <- bStar + mybw * epsStar
 tStar     <- ifelse(tStar >= 1, tStar, 2-tStar)
 # correct the sequence
 tStar     <- scVarHom * (tStar - mean(bStar)) + mean(bStar)
 # get new xobs if input-based, new yobs if output-based
 if (ba == 1){
  MatStar  <- YorX * ( teRef * tStar )
 }
 else {
  MatStar  <- YorX * ( teRef / tStar )
 }
 return(MatStar)
}

.smplHomTE <- function(terfl,Kr,mybw,scVarHom,ba){
 # sample with replacement from Farrell measures
 # (including reflected ones)
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 bStar     <- terfl[newSample]
 # generate the sequence randomly
 epsStar   <- rnorm(Kr)
 tStar     <- bStar + mybw * epsStar
 tStar     <- ifelse(tStar >= 1, tStar, 2-tStar)
 # correct the sequence
 tStar     <- scVarHom * (tStar - mean(bStar)) + mean(bStar)
 # inverse if input-based
 if (ba == 1){
  tStar    <- 1/tStar
 }
 return(tStar)
}

.smplHet1 <- function(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M){
 ZStar1 <- NULL
 # print(class(Zt))
 if (ba == 1){
  ymax <- apply(Zt[, 1:M, drop = FALSE], MARGIN = 2, FUN = max)
  # print(ymax)
 } 	else {
  xmin <- apply(Zt[, M:(nZ-1), drop = FALSE], MARGIN = 2, FUN = min)
  # print(xmin)
 }
 cons3 <- cons1[1,] * cons2[1,]
 for(ww in seq_len(Kr)){
  flag <- TRUE
  while( flag ){
   # step 5
   newSample <- sample( seq_len(2*Kr), 1, replace = TRUE)
   ZStar <- Zt[newSample,]
   # step 6
   epsStar <- rnorm(nZ)
   if(newSample <= Kr){
    epsStar <- epsStar %*% L1
   } else {
    epsStar <- epsStar %*% L2
   }
   # step 7
   ZStar <- ZStar + as.vector(epsStar) * cons3
   if(ba == 1){
    flag0 <- max( ZStar[1:M] - ymax ) > 0
   }
   else {
    flag0 <- min( ZStar[M:(nZ-1)] - xmin ) < 0
   }
   flag <- min(ZStar) < 0 || flag0
  }
  ZStar1 <- rbind(ZStar1, ZStar)
 }
 
 ZStar <- ZStar1
 
 #  # toss possible negative values
 #  nonNegV <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
 #  # nonNegV <- rep(TRUE, Kr)
 #  ZStar <- ZStar[nonNegV,]
 nonNegV <- rep(TRUE, Kr)
 
 # step 8
 ZStar[,nZ] <- ifelse(ZStar[,nZ] > 1, ZStar[,nZ], 2 - ZStar[,nZ])
 
 # step 9
 if (ba == 1){
  teb <- 1 / ZStar[,nZ]
  yrb <- ZStar[, 1:M]
  xrb <- cbind( 1, tan(ZStar[,(M+1):(nZ-1)]) )
 } 	else {
  teb <- ZStar[,nZ]
  yrb <- cbind( 1, tan(ZStar[,1:(M-1)]) )
  xrb <- ZStar[, M:(nZ-1)]
 }
 return( list(Yrb = yrb, Xrb = xrb, teb = teb, Krb = sum(nonNegV)) )
}

.smplHet <- function(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N){
 MinZ <- ifelse(M == 1, M, M - 1)
 if (ba == 1){
  ymax <- apply(Zt[, 1:M, drop = FALSE], MARGIN = 2, FUN = max)
 } 	else {
  xmin <- apply(Zt[, (MinZ+1):(nZ-1), drop = FALSE], MARGIN = 2, FUN = min)
 }
 # first try
 # step 5
 newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
 ZStar <- Zt[newSample,]
 zBarStar <- matrix( colMeans( ZStar ), nrow  = 1)
 
 # step 6
 epsStar <- matrix(rnorm(Kr * nZ), nrow = Kr, ncol = nZ)
 for(ww in seq_len(Kr)){
  if(newSample[ww] <= Kr){
   epsStar[ww,] <- epsStar[ww, , drop = F] %*% L1
  } else {
   epsStar[ww,] <- epsStar[ww, , drop = F] %*% L2
  }
 }
 # step 7
 ZStar <- kronecker (onesN, zBarStar) + cons1 * ( M1 %*% ZStar + cons2 * epsStar)
 
 # toss values outside the frontier (including possible negative values)
 nonNeg <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
 if(ba == 1){
  withinFr <- apply(ZStar[,1:M, drop = FALSE], MARGIN = 1, FUN = function(z) min(ymax-z) > 0)
 }
 else {
  withinFr <- apply(ZStar[,(MinZ+1):(nZ-1), drop = FALSE], MARGIN = 1, FUN = function(z) max(xmin-z) < 0)
 }
 Zgood <- withinFr & nonNeg
 ZStar0 <- ZStar[Zgood,]
 
 # complete Zstar if not full
 if( nrow(ZStar0) < Kr){
  while (nrow(ZStar0) < Kr){
   # step 5
   newSample <- sample( seq_len(2*Kr), Kr, replace = TRUE)
   ZStar <- Zt[newSample,]
   zBarStar <- matrix( colMeans( ZStar ), nrow  = 1)
   # step 6
   epsStar <- matrix(rnorm(Kr * nZ), nrow = Kr, ncol = nZ)
   for(ww in seq_len(Kr)){
    if(newSample[ww] <= Kr){
     epsStar[ww,] <- epsStar[ww, , drop = F] %*% L1
    } else {
     epsStar[ww,] <- epsStar[ww, , drop = F] %*% L2
    }
   }
   # step 7
   ZStar <- kronecker (onesN, zBarStar) + cons1 * ( M1 %*% ZStar + cons2 * epsStar)
   # toss values outside the frontier (including possible negative values)
   nonNeg <- apply(ZStar, MARGIN = 1, FUN = function(x) ifelse(min(x) < 0, FALSE, TRUE))
   if(ba == 1){
    withinFr <- apply(ZStar[,1:M, drop = FALSE], MARGIN = 1, FUN = function(z) min(ymax-z) > 0)
   }
   else {
    withinFr <- apply(ZStar[,(MinZ+1):(nZ-1), drop = FALSE], MARGIN = 1, FUN = function(z) max(xmin-z) < 0)
   }
   Zgood <- withinFr & nonNeg
   ZStar0 <- rbind(ZStar0, ZStar[Zgood,])
  }
 }
 ZStar <- ZStar0[seq_len(Kr),]
 # print(nrow(ZStar))
 
 # step 8
 ZStar[,nZ] <- ifelse(ZStar[,nZ] > 1, ZStar[,nZ], 2 - ZStar[,nZ])
 
 # step 9
 if (ba == 1){
  teb <- 1 / ZStar[,nZ]
  if ( N == 1){
   xrb <- ZStar[, (M+1), drop = FALSE]
  }
  else {
   xrb <- cbind( 1, tan(ZStar[,(M+1):(nZ-1)]) )
  }
  yrb <- ZStar[, 1:M, drop = FALSE]
 }
 else {
  teb <- ZStar[,nZ]
  if ( M == 1){
   yrb <- ZStar[, 1, drop = FALSE]
  }
  else {
   yrb <- cbind( 1, tan(ZStar[,1:MinZ]) )
  }
  xrb <- ZStar[, (MinZ+1):(nZ-1), drop = FALSE]
 }
 return( list(Yrb = yrb, Xrb = xrb, teb = teb, Krb = nrow(ZStar)) )
}

.biasAndCI <- function(te, teboot, msub  = 0, K, Kr, M, N, 
                       level, smoothed, forceLargerOne){
 
 if ( forceLargerOne ){
  teff   <- 1/te
  teboot <- 1/teboot
 }
 else {
  teff <- te
 }
 
 if (smoothed) {
  con1 <- 1
  con2 <- 1
 }
 else {
  con1 <- msub ^ ( 2 / (M + N + 1) )
  con2 <- Kr ^ (-2 / (M + N + 1) )
 }
 con3 <- con1 * con2
 
 # 	# bias-correction
 # 	tebm <- colMeans(teboot, na.rm = TRUE)
 # 	bias <- con3 * (tebm - teff)
 # 	tebc <- teff - bias
 # 	vari <- apply(teboot, 2, var, na.rm = TRUE)
 # 	BovV <- bias^2 / vari / 3
 # 	
 # 	# CI
 # 	teOverTEhat <- sweep(teboot, 2, FUN = "/", teff)
 # 	teOverTEhat <- (teOverTEhat - 1) * con1
 # 	quans <- apply(teOverTEhat, 2, quantile, probs = (100 + c(-level, level))/200, na.rm = T)
 # 	LL <- teff / (1 + quans[2] * con2)
 # 	UU <- teff / (1 + quans[1] * con2)
 # 	
 # 	# drop if inference is based on less than 100
 # 	reaB <- apply(teboot, 2, function(x) sum( !is.na(x) ) )
 # 	
 # 	bias <- ifelse(reaB < 100, NA, bias)
 # 	vari <- ifelse(reaB < 100, NA, vari)
 # 	tebc <- ifelse(reaB < 100, NA, tebc)
 # 	LL <- ifelse(reaB < 100, NA, LL)
 # 	UU <- ifelse(reaB < 100, NA, UU)
 
 bias <- vari <- reaB <- BovV <- tebc <- LL <- UU <- numeric(K)
 
 for (i in seq_len(K) ) {
  TEi <- teboot[,i]
  TEi <- TEi[!is.na(TEi)]
  reaB[i] <- length(TEi)
  if( reaB[i] < 100 ){
   bias[i] = NA
   vari[i] = NA
   BovV[i] = NA
   tebc[i] = NA
   LL[i]   = NA
   UU[i]   = NA		
  }
  else {
   bias[i] = con3 * ( mean(TEi) - teff[i] )
   vari[i] = var(TEi)
   BovV[i] = (bias[i])^2 / vari[i] / 3
   tebc[i] = teff[i] - bias[i]
   teOverTEhat = (TEi / teff[i] - 1) * con1
   quans = quantile(teOverTEhat, probs = (100 + c(-level, level))/200, na.rm = TRUE)
   LL[i] = teff[i] / (1 + quans[2] * con2)
   UU[i] = teff[i] / (1 + quans[1] * con2)
  }
 }
 
 return(list(reaB=reaB, bias=bias, vari=vari, BovV=BovV, tebc=tebc, LL=LL, UU=UU))
 
}


.pvalsTestOne <- function(seCrs, seCrsMean, te1boot, te2boot, ba){
 seCrsB <- te1boot / te2boot
 seCrsMeanB <- rowMeans(te1boot, na.rm = TRUE) / rowMeans(te2boot, na.rm = TRUE)
 seCrsMeanB <- na.omit( seCrsMeanB )
 if(ba ==1){
  # global CRS
  pgCRS <- mean(seCrsMeanB <= seCrsMean)
  # local CRS
  plCRS <- rowMeans( apply(seCrsB, MARGIN = 1, FUN = function(z) z <= seCrs), na.rm = TRUE )
 }
 else {
  pgCRS <- mean(seCrsMeanB >= seCrsMean)
  # local CRS
  plCRS <- rowMeans( apply(seCrsB, MARGIN = 1, FUN = function(z) z >= seCrs), na.rm = TRUE )
 }
 
 # count those that are not NA
 nonNa <- apply(seCrsB, MARGIN = 2, FUN = function(z) length(na.omit(z)))
 plCRS <- ifelse(nonNa >= 100, plCRS, NA)
 
 return(list(pgCRS = pgCRS, plCRS = plCRS, nonNa = nonNa, reaB = length(seCrsMeanB)))
}

.pvalsTestTwo <- function(seNrs, seNrsMean, te1boot, te2boot, ba, performGlobal, s.inefficient){
 seNrsB <- te1boot / te2boot
 # cat("ncol(te1boot) = ",ncol(te1boot),"\n")
 # if(ncol(te1boot)==1) seNrsB <- matrix(seNrsB, ncol = 1)
 if ( performGlobal ){
  seNrsMeanB <- rowMeans(te1boot, na.rm = TRUE) / rowMeans(te2boot, na.rm = TRUE)
  seNrsMeanB <- na.omit( seNrsMeanB )
  
  if(ba ==1){
   # global NRS
   pgNRS <- mean(seNrsMeanB <= seNrsMean)
   # local NRS
   plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs), na.rm = TRUE )
  }
  else {
   # global NRS
   pgNRS <- mean(seNrsMeanB >= seNrsMean)
   # local NRS
   plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs), na.rm = TRUE )
  }
 }
 else {
  seNrs1 <- seNrs[s.inefficient]
  if(ba ==1){
   # local NRS
   if(ncol(te1boot) == 1){
    plNRS <- mean( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs1), na.rm = TRUE )
   }
   else {
    plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z <= seNrs1), na.rm = TRUE )
   }
  }
  else {
   # local NRS
   if(ncol(te1boot) == 1){
    plNRS <- mean( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs1), na.rm = TRUE )
   }
   else {
    plNRS <- rowMeans( apply(seNrsB, MARGIN = 1, FUN = function(z) z >= seNrs1), na.rm = TRUE )
   }
  }
 }
 # count those that are not NA
 nonNa <- apply(seNrsB, MARGIN = 2, FUN = function(z) length(na.omit(z)))
 plNRS <- ifelse(nonNa >= 100, plNRS, NA)
 
 # write into appropriate cells
 pineffdrs <- rep(NA, length(s.inefficient))
 
 if (performGlobal){
  pineffdrs <- ifelse(s.inefficient, plNRS, pineffdrs)
 }
 else {
  pineffdrs[s.inefficient] <- plNRS
 }
 if ( performGlobal ){
  tymch <- list(pgNRS = ifelse(performGlobal, pgNRS, NA), pineffdrs = pineffdrs, nonNa = nonNa, reaB = length(seNrsMeanB))
 }
 else {
  tymch <- list(pgNRS = ifelse(performGlobal, pgNRS, NA), pineffdrs = pineffdrs, nonNa = nonNa)
 }
 return(tymch)
}

# empirical distribution function
.edf <- function(Y,y){
 if(length(y) == 1){
  f1 <- Y<=y
  # f1 <- y-Y > 1e-6
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  # f2 <- apply(f3, MARGIN = 1, FUN = function(z)  z[1]-z[-1]> 1e-6)
  f1 <- rowMeans(f2)==1
 }
 return(mean(f1))
}

# empirical joint distribution function 
.ejdf <- function(Y, X, y, x){
 # 1
 if(length(y) == 1){
  f1 <- Y<=y
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  f1 <- rowMeans(f2)==1
 }
 # 2
 if(length(x) == 1){
  g1 <- X<=x
 }
 else {
  g3 <- cbind(x, t(X))
  g2 <- apply(g3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  g1 <- rowMeans(g2)==1
 }
 sum(f1 & g1) / length(f1)
}

# empirical joint distribution function 
.ejdfedf1 <- function(Y, X, y, x){
 # 1
 if(length(y) == 1){
  f1 <- Y<=y
 }
 else {
  f3 <- cbind(y, t(Y))
  f2 <- apply(f3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  f1 <- rowMeans(f2)==1
 }
 # 2
 if(length(x) == 1){
  g1 <- X<=x
 }
 else {
  g3 <- cbind(x, t(X))
  g2 <- apply(g3, MARGIN = 1, FUN = function(z) z[-1] <= z[1])
  g1 <- rowMeans(g2)==1
 }
 # mean(f1 & g1) - mean(f1)*mean(g1)
 mean(f1)
}

.ejdfedf <- function(Y, X, y, x){
 # 1
 f1 <- Y <= y
 # 2
 g2 <- X
 for(i in seq_len(ncol(X))){
  g2[,i] <- X[,i] <= x[i]
 }
 g1 <- ( rowMeans(g2) ==1 )
 mean(f1 & g1) - mean(f1)*mean(g1)
 # mean(f1)
}

.t4n <- function(w,d, print=FALSE){
 n <- length(d)
 tt <- numeric(n)
 for(i in 1:n){
  tt[i] <- .ejdfedf(Y = d, y = d[i], X = w, x = w[i,]) # - .edf(Y = d, y = d[i]) * .edf(Y = w, y = w[i,]) #
 }
 if(print) print(cbind(tt,d))
 return(sum(tt^2))
}

# begin parallel computing
.run.boots.hom.rts <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 ba <- args$ba
 teRef <- args$teRef
 terfl <- args$terfl
 mybw <- args$mybw
 scVarHom <- args$scVarHom
 
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 # Prepare output structure
 te1boot <- rep(NA, K)
 te2boot <- rep(NA, K)
 
 for(ii in seq_len(K)){
  # begin homogeneous
  # print(1)
  if (ba == 2) {
   # step 1: sampling
   Yb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,Y,ba)
   # step 2: applying DEA
   # CRS or NiRS
   te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,args$rtsHo,ba,
                         0,print.level=0)
   # VRS
   te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yb),t(X),K,1,ba,
                         0,print.level=0)
  }
  else {
   # step 1: sampling
   Xb <- .smplHom(teRef,terfl,Kr=K,mybw,scVarHom,X,ba)
   # step 2: applying DEA
   # CRS or NiRS
   te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,args$rtsHo,ba,
                         0,print.level=0)
   # VRS
   te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Y),t(Xb),K,1,ba,
                         0,print.level=0)
  }
  # end homogeneous
  # end for ii
 }
 return (c(te1boot, te2boot, K, 0))
}

.run.boots.hom.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 ba <- args$ba
 teRef <- args$teRef
 terfl <- args$terfl
 mybw <- args$mybw
 scVarHom <- args$scVarHom
 
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 
 # begin homogeneous
 # print(1)
 # teB <- numeric(K)
 if (ba == 2) {
  # step 1: sampling
  Yrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Yr,ba)
  # step 2: applying DEA
  teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xr),Kr,rt,ba,
                0,print.level=0)
 }
 else {
  # step 1: sampling
  Xrb <- .smplHom(teRef,terfl,Kr,mybw,scVarHom,Xr,ba)
  # step 2: applying DEA
  teB <- .teRad(t(Y),t(X),M,N,K,t(Yr),t(Xrb),Kr,rt,ba,
                0,print.level=0)
 }
 return (c(teB, K, 1))
}

.run.boots.het.rts <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 Zt <- args$Zt
 L1 <- args$L1
 L2 <- args$L2
 M1 <- args$M1
 cons1 <- args$cons1
 cons2 <- args$cons2
 onesN <- args$onesN
 ba <- args$ba
 
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 nZ <- ncol(Zt)
 # Prepare output structure
 te1boot <- rep(NA, K)
 te2boot <- rep(NA, K)
 
 for(ii in seq_len(K)){
  # begin heterogeneous
  # print(2)
  # step 1: sampling
  smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N)
  # step 2: get efficiency under H0
  teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                    t(Yr),t(Xr),Kr,rts=3,ba,0,print.level=0)
  # step 3: get efficient xRef or yRef
  if (ba == 1) {
   Yrb <- smplHet$Yrb
   Xrb <- smplHet$Xrb * teRefB / smplHet$teb
  }
  else {
   Yrb <- smplHet$Yrb * teRefB / smplHet$teb
   Xrb <- smplHet$Xrb
  }
  # handle infeasible
  teRefB.good <- !is.na(teRefB)
  # modify the existing matrices by tossing the missing values
  Yrb <- Yrb[teRefB.good, , drop = FALSE]
  Xrb <- Xrb[teRefB.good, , drop = FALSE]
  # new number of obs in reference
  KrB <- sum(teRefB.good)
  cat("my sample is ",KrB,"\n", sep = "")
  # step 4: get the bootstrapped TE
  # CRS or NiRS
  te1boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,args$rtsHo,ba,
                        0,print.level=0)
  # VRS
  te2boot[ii] <- .teRad(Y[ii,],X[ii,],M,N,1,t(Yrb),t(Xrb),KrB,1,ba,
                        0,print.level=0)
  # end heterogeneous
  # end for ii
 }
 return (c(te1boot, te2boot, K, 0))
}

.run.boots.het.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 Zt <- args$Zt
 L1 <- args$L1
 L2 <- args$L2
 M1 <- args$M1
 cons1 <- args$cons1
 cons2 <- args$cons2
 onesN <- args$onesN
 ba <- args$ba
 
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 nZ <- ncol(Zt)
 
 # begin heterogeneous
 # print(2)
 # step 1: sampling
 smplHet <- .smplHet(Zt,L1,L2,M1,cons1,cons2,onesN,Kr,nZ,ba,M,N)
 # step 2: get efficiency under H0
 teRefB  <- .teRad(t(smplHet$Yrb),t(smplHet$Xrb),M,N,smplHet$Krb,
                   t(Yr),t(Xr),Kr,rt,ba, 0,print.level=0)
 # step 3: get efficient xRef or yRef
 if (ba == 1) {
  Yrb <- smplHet$Yrb
  Xrb <- smplHet$Xrb * teRefB / smplHet$teb
 }
 else {
  Yrb <- smplHet$Yrb * teRefB / smplHet$teb
  Xrb <- smplHet$Xrb
 }
 # handle infeasible
 teRefB.good <- !is.na(teRefB)
 # modify the existing matrices by tossing the missing values
 Yrb <- Yrb[teRefB.good, , drop = FALSE]
 Xrb <- Xrb[teRefB.good, , drop = FALSE]
 # new number of obs in reference
 KrB <- sum(teRefB.good)
 # step 4: get the bootstrapped TE
 teB <- .teRad(t(Y),t(X),M,N,K,t(Yrb),t(Xrb),KrB,rt,ba,
               0,print.level=0)
 # end heterogeneous
 return (c(teB, KrB, 1))
}

.run.dots <- function(xs, nrep, args){
 width <- args$width
 character <- "."
 #  mylen <- length(xs[[nrep]])
 #  if (mylen > 0){
 #   if(xs[[nrep]][ mylen ] == 1){
 #    # boot for BC
 #    K <- mylen - 2
 #   }
 #   else {
 #    # boot for RTS test
 #    K <- (mylen - 2 )/2
 #   }
 #   # Pre-Last value in each output is 'KrB'
 #   KrB <- xs[[nrep]][ mylen-1 ]
 #   # cat("krB = ",KrB," K = ", K, " \n", sep = "")
 #   if (KrB/K < 0.80){
 #    character <- "?"
 #   }
 #   else if (KrB/K < 0.95){
 #    character <- "@"
 #   }
 #   else if (KrB/K < 1){
 #    character <- "x"
 #   }
 #  }
 
 #  if (!is.numeric(width)) {
 #   stop("'width' should be numeric")
 #  }
 #  if (width != 50 & width != 60 & width != 70 & 
 #      width != 80 & width != 90  & width != 100){
 #   stop("'width' should be 50, 60, 70, 80, 90, or 100")
 #  }
 #  if (nrep == 0){
 #   if (!is.null(message)){
 #    cat("",message,"\n", sep = "")
 #   }
 #   cat("____|___ 1 ___|___ 2 ___|___ 3 ___|___ 4 ___|___ 5", sep = "")
 #   if (width == 60) {
 #    cat(" ___|___ 6", sep = "")
 #   }
 #   if (width == 70) {
 #    cat(" ___|___ 6 ___|___ 7", sep = "")
 #   }
 #   if (width == 80) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8", sep = "")
 #   }
 #   if (width == 90) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9", sep = "")
 #   }
 #   if (width == 100) {
 #    cat(" ___|___ 6 ___|___ 7 ___|___ 8 ___|___ 9 ___|___ 10", sep = "")
 #   }  
 #   cat("\n")
 #  }
 #  else {
 if (nrep/width != floor(nrep/width)){
  cat("",character,"", sep = "")
 }
 else {
  cat("",character," ",nrep,"\n", sep = "")
 }
 # }
}

.run.boots.subs.bc <- function(inp, ...){
 # Read the arguments
 args <- list(...)[[1]]
 Y <- args$Y
 X <- args$X
 Yr <- args$Yr
 Xr <- args$Xr
 rt <- args$rt
 ba <- args$ba
 msub <- args$msub
 
 M  <- ncol(Y)
 N  <- ncol(X)
 K  <- nrow(Y)
 Kr <- nrow(Yr)
 
 # begin subsampling
 # print(3)
 # step 1: sub-sampling
 newSample <- sample( seq_len(Kr), msub, replace = TRUE)
 ystar <- Yr[newSample, , drop = FALSE]
 xstar <- Xr[newSample, , drop = FALSE]
 # step 2: applying DEA
 teB <- .teRad(t(Y),t(X),M,N,K,t(ystar),t(xstar),msub,rt,ba,
               0,print.level=0)
 mychar <- "."
 # end subsampling
 return (c(teB, K, 1))
}

# end parallel computing

.prepareYXZ <- function(formula, uhet = NULL, vhet = NULL, tmean = NULL, data, subset, ...) {
 needed.frame <- sys.nframe() - 1
 mf0 <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
 if ( is.null(uhet) ){
  uhet <- ~ 1
  ku <- 1
 } else {
  ku <- 17
 }
 if ( is.null(vhet) ){
  vhet <- ~ 1
  kv <- 1
 } else {
  kv <- 17
 }
 if ( is.null(tmean) ){
  tmean <- ~ 1
  kdel <- 1
 } else {
  kdel <- 17
 }
 
 form1 <- Formula(as.formula(paste("",deparse(formula, width.cutoff = 500L)," | ",uhet[2]," | ",vhet[2]," | ",tmean[2],"", sep = "")))
 
 form <- Formula(as.formula(paste("",deparse(formula, width.cutoff = 500L)," + ",uhet[2]," + ",vhet[2]," + ",tmean[2],"", sep = "")))
 
 # check if it is a matrix
 datasupplied <- !(match("data", names(mf0), 0) == 0)
 
 
 if(datasupplied){
  # begin get a logical vector equal TRUE if !missing
  # first using data and subset to get x without NA
  mf <- mf0
  mf$formula <- formula( form )
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # now get the names in the entire data
  esample <- seq_len( nrow(data) ) %in% as.numeric(rownames(X))
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  # print(form1)
  dataesample <- model.frame(form1, data = data[esample,])
  # print(dataesample)
  Y <- as.matrix(model.part(form1, data = dataesample, lhs = 1, drop = FALSE))
  # Y <- as.matrix( model.matrix(formula(form1, lhs = 1, rhs = 0), data = data[esample,]))
  # print(Y)
  X <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 1), data = dataesample))
  # print(X)
  n <- nrow(Y)
  k <- ncol(X)
  if(ku == 1){
   Zu <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zu <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 2), data = dataesample))
   ku <- ncol(Zu)
  }
  if(kv == 1){
   Zv <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zv <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 3), data = dataesample))
   kv <- ncol(Zv)
  }
  if(kdel == 1){
   Zdel <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zdel <- as.matrix( model.matrix(formula(form1, lhs = 0, rhs = 4), data = dataesample))
   kdel <- ncol(Zdel)
  }
  #   print(head(Y))
  #   print(head(X))
  #   print(head(Zu))
  #   print(head(Zv))
  #   print(head(Zdel))
 }
 # if data are not supplied
 else {
  # begin get a logical vector equal TRUE if !missing
  
  # first using data and subset to get XZ without NA
  mf <- mf0
  mf$formula <- formula( form )
  subsetsupplied <- !(match("subset", names(mf), 0) == 0)
  if(subsetsupplied) stop("Subset with matrices is not allowed")
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  Y <- as.matrix(model.response(mf))
  n <- nrow(Y)
  XZ <- as.matrix(model.matrix(mt, mf))
  # get a logical vector equal TRUE if !missing
  with.na <- model.frame(mt, na.action = na.pass)
  esample <- rownames(with.na) %in% rownames(XZ)
  if(sum(esample)==0){
   stop("no valid data points")
  }
  # print(table(esample))
  # end get a logical vector equal TRUE if !missing
  
  # get the data
  #   print(head(Y))
  #   print(head(XZ))
  
  # get X
  mf <- mf0
  m <- match(c("formula"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
  mt <- attr(mf, "terms")
  X <- as.matrix(model.matrix(mt, mf))
  # print(head(X))
  k <- ncol(X)
  
  # get Zu
  if(ku == 1){
   Zu <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",uhet[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZu <- as.matrix(model.matrix(mt, mf))
   Zu <- XZu[, -(2:k), drop = FALSE]
   ku <- ncol(Zu)
   # print(head(Zu))
  }
  
  # get Zv
  if(kv == 1){
   Zv <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   mf <- mf0
   mf$formula <- formula( Formula(as.formula(paste("",deparse(formula)," + ",uhet[2]," + ",vhet[2],"", sep = ""))))
   m <- match(c("formula"), names(mf), 0L)
   mf <- mf[c(1L, m)]
   mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, sys.frame(sys.parent(n = needed.frame)))
   mt <- attr(mf, "terms")
   XZuZv <- as.matrix(model.matrix(mt, mf))
   Zv <- XZuZv[, -(2:(k+ku-1)), drop = FALSE]
   kv <- ncol(Zv)
   # print(head(Zv))
  }
  
  # get Zdel
  if(kdel == 1){
   Zdel <- matrix(1, nrow = sum(esample), ncol = 1)
  }
  else {
   Zdel <- XZ[, -(2:(k+ku-1+kv-1)), drop = FALSE]
   # print(head(Zdel))
  }
  
  #   print(head(Zu))
  #   print(head(Zv))
  #   print(head(Zdel))
 }
 colnames(X)[1] <- colnames(Zu)[1] <- colnames(Zv)[1] <- colnames(Zdel)[1] <- "Intercept"
 
 tymch <- list(Y = Y, X = X, Zu = Zu, Zv = Zv, Zdel = Zdel, n = n, k = k, ku = ku, kv = kv, kdel = kdel, esample = esample)
 class(tymch) <- "npsf"
 return(tymch)
 
}

# Half-normal model

# Log-likelihood
.ll.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 s <- ifelse(prod, -1, 1) 
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y-X%*%beta
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 
 # Log-likelihood  
 llf <- sum(log(2) - log(sqrt(2*pi)) - log(sig) + pnorm(s*e*lmd/sig, log.p = TRUE) - 0.5*(e/sig)^2)
 return(llf)
}

# Gradient
.g.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 if(prod == TRUE){ s = -1 } else {s = 1}
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y - X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(expv + expu); lmd <- sqrt(expu/expv)
 ls <- lmd/sig; g <- dnorm(s*e*ls)/pnorm(s*e*ls)
 gels <- g*e*ls; 
 
 # Gradient
 gb <- t(X)%*%(-s*g*ls + e/sig^2)
 ggv <- 0.5*t(Zv)%*%(expv/sig^2 * ((e/sig)^2 - 1) - s*gels*(1 + expv/sig^2))
 ggu <- 0.5*t(Zu)%*%(expu/sig^2 * ((e/sig)^2 - 1) - s*gels*(expu/sig^2 - 1))                                                                                                                         
 grad <- rbind(gb, ggv, ggu)
 return(grad)
}

# Hessian
.hess.hn <- function(theta, prod, k, kv, ku, kdel = NULL, y, Zv, Zu, Zdel = NULL, X) {
 s <- ifelse(prod, -1, 1) 
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[-c(1:k, (k+1):(k+kv))]
 e <- y - X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv)
 sig <- sqrt(expv + expu); lmd <- sqrt(expu/expv)
 ls <- lmd/sig; els <- e*ls; 
 g <- dnorm(s*e*ls)/pnorm(s*e*ls); gels <- g*els
 
 # Hessian
 Hb <- t(as.numeric(-s*g*(els + s*g)*ls^2 - sig^(-2))*X)%*%X 
 Hbgu <- t(as.numeric(0.5*g*ls*(1 - expu/sig^2)*((g + s*els)*els - s*1) - e*expu/sig^4)*X)%*%Zu
 Hbgv <- t(as.numeric(-s*0.5*g*ls*(1 + expv/sig^2)*((els + s*g)*els - 1) - e*expv/sig^4)*X)%*%Zv
 Hgvu <- t(0.5*as.numeric(-s*0.5*gels*(1 - expu/sig^2)*(1 + expv/sig^2)*(1 - s*(g + s*els)*els) - expv*expu*(2*(e/sig)^2 - s*gels - 1)/sig^4)*Zv)%*%Zu
 Hgu <- t(0.5*as.numeric((1 - expu/sig^2)*(expu/sig^2 * ((e/sig)^2 - s*gels - 1) + s*0.5*gels*(1 - expu/sig^2)*(1 - s*(g + s*els)*els)) - (expu*e/sig^3)^2)*Zu)%*%Zu
 Hgv <- t(0.5*as.numeric(expv/sig^2 * ((1 - expv/sig^2)*((e/sig)^2 - s*gels - 1) - e^2*expv/sig^4) - s*0.5*gels*(1 + expv/sig^2)^2 * ((els + s*g)*els - 1))*Zv)%*%Zv
 H <- cbind(rbind(Hb, t(Hbgv), t(Hbgu)), rbind(Hbgv, Hgv, t(Hgvu)), rbind(Hbgu, Hgvu, Hgu))
 return(H)
} 

# Truncated-normal model

#Log-likelihood
.ll.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1) 
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta
 
 # Log-likelihood  
 llf <- sum(-log(sqrt(2*pi)) - log(sig) - pnorm(mu/sqrt(exp(Zu%*%gu)), log.p = TRUE) + pnorm(mu/(sig*lmd) + s*e*lmd/sig, log.p = TRUE) - 0.5*(((e - s*mu)^2)/sig^2))
 return(llf)
}

#Gradient

.g.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv) 
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta; ls <- lmd/sig
 g1 <- dnorm(ls*(mu/lmd^2 + s*e))/pnorm(ls*(mu/lmd^2 + s*e))
 g2 <- dnorm(mu/sqrt(expu))/pnorm(mu/sqrt(expu)); ml = mu/lmd^2
 d1 <- 0.5*mu/(sig*lmd)*(1 - expv/sig^2)
 d2 <- -0.5*s*e*ls * (1 + expv/sig^2)
 d3 = -0.5*mu/(sig*lmd)*(1 + expu/sig^2)
 d4 = 0.5*s*e*ls * (1 - expu/sig^2)
 
 # Gradient
 gb <- t(X)%*%((e - s*mu)/sig^2 - s*g1*ls)
 gdel <- t(Zdel)%*%(-g2/sqrt(expu) + g1/(sig*lmd) + s*(e - s*mu)/sig^2)
 gv <- t(Zv)%*%(-0.5*expv/sig^2 + g1*(d1 + d2) + 0.5*expv*(e - s*mu)^2/sig^4)
 gu  <- t(Zu)%*%(-0.5*expu/sig^2 + 0.5*g2*mu/sqrt(expu) + g1*(d3 + d4) + 0.5*expu*(e - s*mu)^2/sig^4)
 grad <- rbind(gb, gv, gu, gdel)
 return(grad)
}

#Hessian
.hess.tn <- function(theta, prod, k, kv, ku, kdel, y, Zv, Zu, X, Zdel) {
 s <- ifelse(prod, -1, 1)
 beta <- theta[1:k]
 gv <- theta[(k+1):(k+kv)]
 gu <- theta[(k+kv+1):(k+kv+ku)]
 delta <- theta[(k+kv+ku+1):(k+kv+ku+kdel)]
 e <- y-X%*%beta
 expu <- exp(Zu%*%gu); expv <- exp(Zv%*%gv) 
 sig <- sqrt(exp(Zv%*%gv) + exp(Zu%*%gu))
 lmd <- sqrt(exp(Zu%*%gu)/exp(Zv%*%gv))
 mu <- Zdel%*%delta; ls <- lmd/sig
 g1 <- dnorm(ls*(mu/lmd^2 + s*e))/pnorm(ls*(mu/lmd^2 + s*e))
 gr1u <- ls*(0.5*(1 - expu/sig^2)*(mu/lmd^2 + s*e) - mu/lmd^2)*g1*(-s*ls*(e + s*mu/lmd^2) - g1)
 gr1v <- ls*(mu/lmd^2 - 0.5*(1 + expv/sig^2)*(mu/lmd^2 + s*e))*g1*(-s*ls*(e + s*mu/lmd^2) - g1)
 g2 <- dnorm(mu/sqrt(expu))/pnorm(mu/sqrt(expu)); ml = mu/lmd^2
 
 # Hessian
 Hb <- t(as.numeric(-1/sig^2 - s*g1*ls^2*(ls*(e + s*ml) + s*g1))*X)%*%X 
 Hbdel <- t(as.numeric(-s/sig^2 - s*g1*(ls/lmd)^2*(-s*ls*(e + s*ml) - g1))*X)%*%Zdel
 Hbgv <- t(as.numeric(-(e - s*mu)*expv/sig^4 + s*0.5*ls*(1 + expv/sig^2)*g1 - s*ls*gr1v)*X)%*%Zv
 Hbgu <- t(as.numeric(-(e - s*mu)*expu/sig^4 - s*0.5*ls*(1 - expu/sig^2)*g1 - s*ls*gr1u)*X)%*%Zu
 Hdel <- t(as.numeric(g2/expu*(mu/sqrt(expu) + g2) - g1*ls/(lmd^3*sig) * (ls*(ml + s*e) + g1) - sig^(-2))*Zdel)%*%Zdel
 Hdelgv <- t(as.numeric(1/(lmd*sig) * (0.5*(1 - expv/sig^2)*g1 + gr1v) - s*(e - s*mu)*expv/sig^4)*Zdel)%*%Zv
 Hdelgu <- t(as.numeric(-s*(e - s*mu)*expu/sig^4 + 1/(sig*lmd)*(gr1u - 0.5*g1*(1 + expu/sig^2)) - 0.5*g2/sqrt(expu)*(-1 + mu/sqrt(expu)*(mu/sqrt(expu) + g2)))*Zdel)%*%Zu
 Hgv <- t(as.numeric(0.5*expv/sig^2 * ((1 - expv/sig^2)*((e - s*mu)^2/sig^2 - 1) - (e - s*mu)^2*expv/sig^4) + (ml - 0.5*(1 + expv/sig^2)*(ml + s*e))*(gr1v*ls - 0.5*g1*ls*(1 + expv/sig^2)) + g1*ls*(ml - 0.5*(expv/sig^2 * (1 - expv/sig^2)*(ml + s*e) + (1 + expv/sig^2)*ml)))*Zv)%*%Zv
 Hgu <- t(as.numeric((0.5*(1 - expu/sig^2)*(ml + s*e) - ml)*ls*(gr1u + 0.5*g1*(1 - expu/sig^2)) + g1*ls*(0.5*(expu/sig^2 - 1)*(expu/sig^2*(ml + s*e) + ml) + ml) + 0.5*(expu/sig^2 * ((1 - expu/sig^2)*((e - s*mu)^2/sig^2 - 1) - expu*(e - s*mu)^2/sig^4) + 0.5*g2*mu/sqrt(expu)*(-1 + mu/sqrt(expu)*(mu/sqrt(expu) + g2))))*Zu)%*%Zu
 Hgvgu <- t(as.numeric(0.5*expv*expu/sig^4 * (1 - 2*((e - s*mu)/sig)^2) + (ml - 0.5*(1 + expv/sig^2)*(ml + s*e))*ls*(gr1u + 0.5*(1 - expu/sig^2)*g1) + g1*ls*(-ml + 0.5*(expv*expu/sig^4*(ml + s*e) + (1 + expv/sig^2)*ml)))*Zv)%*%Zu
 
 H <- cbind(rbind(Hb, t(Hbgv),  t(Hbgu),t(Hbdel)), rbind(Hbgv, Hgv, t(Hgvgu), Hdelgv),rbind(Hbgu, Hgvgu, Hgu, Hdelgu), rbind(Hbdel, t(Hdelgv), t(Hdelgu), Hdel ))
 return(H)
}



# Technical efficiencies and prediction intervals
.u2efftnm <- function( e, su, sv, mu, alpha = 0.05, prod) {
 if(prod == T){sn = -1} else {sn = 1}
 s  <- sqrt(su^2 + sv^2);  m1 <- (sn*su^2 * e + mu*sv^2)/s^2
 s1 <- su * sv / s;  z  <- m1 / s1
 point.est.mean <- m1 + s1 * dnorm(z) / pnorm(z) 
 point.est.mode <- ifelse( m1 >= 0, m1, 0 ) 
 te_jlms_mean <- exp( -point.est.mean)
 te_jlms_mode <- exp( -point.est.mode)
 zl    <- qnorm( 1 - alpha / 2 * pnorm(z) )
 zu    <- qnorm( 1 - ( 1 - alpha/2 ) * pnorm(z) )
 te_l  <- exp( -m1 - zl * s1 )
 te_u  <- exp( -m1 - zu * s1 )
 te_bc <- exp(-m1 + .5 * s1^2) * pnorm(-s1 + z) / pnorm(z)
 tymch <- data.frame(te_l, te_jlms_mean, te_jlms_mode,te_bc,te_u)
 colnames(tymch) <- c("Lower bound","JLMS", "Mode", "BC","Upper bound" )
 return(tymch)
}



# Marginal effects
.me = function(theta, Zu, Zdel, ku, kdel, n, dist = c("h", "t")){
 
 mat.equal <- function(x, y) is.matrix(x) && is.matrix(y) && ncol(x) == ncol(y) && all(colnames(x) == colnames(y))
 
 gu = theta[1:ku,1,drop = F]
 expu = exp(Zu%*%gu)
 if(dist == "t"){
  delta = theta[-c(1:ku),1,drop = F]
  mu <- Zdel%*%delta
  if(length(delta) == 1) delta = rep(0, max(ku, kdel))
 }
 if(length(gu) == 1) gu = rep(0, max(ku, kdel))
 meff = matrix(NA, ncol = max(ku, kdel) - 1, nrow = n)
 
 if(dist == "h"){
  if(ncol(Zu) == 1) {
   warning("Marginal effects are not returned: scale of half-normal distribution of inefficiency term is not expressed as function of any exogenous variables", call. = FALSE)
   meff = NULL
  } else {
   arg = sqrt(1/(2*pi)) * sqrt(expu)
   for(i in 2:ku){
    meff[,i-1] = arg*gu[i]  
    meff = round(meff, digits = 5)
   }
   colnames(meff) = rownames(gu)[-1]
  }} else if(ncol(Zu) == 1 & ncol(Zdel) == 1){
   warning("Marginal effects are not returned: mean or variance of (pre-)truncated normal distribution of inefficiency term is not expressed as function of any exogenous variables", call. = FALSE) 
   meff = NULL
  }else  if(dist == "t" & (mat.equal(Zu, Zdel)| all(delta == 0) | all(gu == 0))){
   arg = mu/sqrt(expu)
   g = dnorm(arg)/pnorm(arg)
   arg1 = (1 - arg*g - g^2); arg2 = sqrt(expu)/2 * ((1 + arg^2)*g + arg*g^2)
   for(i in 2:max(ku, kdel)){
    meff[,i-1] = arg1*delta[i] + arg2*gu[i]
    meff = round(meff, digits = 5)
   }
   if(all(gu == 0)){colnames(meff) = rownames(delta)[-1]} else {colnames(meff) = rownames(gu)[-1]}
  }  else {
   warning("Marginal effects are not returned: mean and variance of (pre-)truncated normal distribution of inefficiency term are expressed as functions not of the same exogenous variables", call. = FALSE)
   meff = NULL
  }
 return( meff)
}


.skewness <- function (x, na.rm = FALSE, type = 3) 
{
 if (any(ina <- is.na(x))) {
  if (na.rm) 
   x <- x[!ina]
  else return(NA)
 }
 if (!(type %in% (1:3))) 
  stop("Invalid 'type' argument.")
 n <- length(x)
 # x <- x - mean(x)
 y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
 if (type == 2) {
  if (n < 3) 
   stop("Need at least 3 complete observations.")
  y <- y * sqrt(n * (n - 1))/(n - 2)
 }
 else if (type == 3) 
  y <- y * ((1 - 1/n))^(3/2)
 y
}

# Print the estimation results

.printoutcs = function(x, digits, k, kv, ku, kdel, na.print = "NA", dist){
 
 Cf = cbind(ifelse(x[,1, drop = F]> 999, formatC(x[,1, drop = F], digits = 1, format = "e",width = 10), formatC(x[,1, drop = F], digits = digits, format = "f", width = 10)),
            ifelse(x[,2, drop = F]>999, formatC(x[,2, drop = F], digits = 1, format = "e", width = 10), formatC(x[,2, drop = F], digits = digits, format = "f", width = 10)),
            ifelse(x[,3, drop = F]>999, formatC(x[,3, drop = F], digits = 1, format = "e", width = 7), formatC(x[,3, drop = F], digits = 2, format = "f", width = 7)),
            ifelse(x[,4, drop = F]>999, formatC(x[,4, drop = F], digits = 1, format = "e", width = 10), formatC(x[,4, drop = F], digits = digits, format = "f", width = 10)))
 
 cat("               Coef.        SE       z       P>|z|\n", sep = "")
 dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
 cat("__________________________________________________", "\n", "Frontier", "\n")
 print.default(Cf[1:k,,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 cat("--------------------------------------------------", "\n", "Random noise component: log(sigma_v^2)", "\n")
 # dimnames(Cf)[[2]] <- rep("", dim(Cf)[[2]])
 print.default(Cf[(k+1):(k+kv),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 cat("--------------------------------------------------", "\n", "Inefficiency component: log(sigma_u^2)", "\n")
 print.default(Cf[(k+kv+1):(k+kv+ku),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 if(dist == "t"){
  cat("--------------------------------------------------", "\n", "Mu (location parameter)", "\n")
  print.default(Cf[(k+kv+ku+1):(k+kv+ku+kdel),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 }
 
 if(nrow(Cf[-c(1:(k+kv+ku+kdel)),,drop=F]) >= 1){
  cat("--------------------------------------------------", "\n", "Parameters of compound error distribution", "\n")
  print.default(Cf[-c(1:(k+kv+ku+kdel)),,drop=F], quote = FALSE, right = TRUE, na.print = na.print)
 }
 cat("__________________________________________________", "\n")
 invisible(x)
}



.timing <- function(x, wording = "Time elapsed is")
{
 if(x < 60){
  # cat("\n")
  cat("",wording,"",x," seconds\n", sep = "")
  # cat("\n")
 } else {
  if(x >= 60 & x < 60*60){
   minutes <- floor(x/60)
   seconds <- round(x - minutes * 60,1)
   # cat("\n")
   cat("",wording,"",minutes," minute(s) and ",seconds," second(s)\n", sep = "")
   # cat("\n")
  } else {
   if(x >= 60*60 & x < 60*60*24){
    hours   <- floor(x / 60 / 60)
    minutes <- round( (x - hours * 60 *60) / 60, 1)
    seconds <- floor(x - hours * 60 *60 - minutes * 60)
    # cat("\n")
    cat("",wording,"",hours," hour(s) and ",minutes," minute(s) \n", sep = "")
    # cat("\n")
   } else {
    if(x >= 60*60*24){
     days    <- floor(x / 60 / 60 / 24)
     hours   <- round( (x - days * 60 * 60 * 24) / 60 /60 ,1)
     minutes <- floor( (x - days * 60 * 60 * 24 - hours * 60 *60) / 60)
     seconds <- floor(x - days * 60 * 60 * 24 - hours * 60 *60 - minutes * 60)
     # cat("\n")
     cat("",wording,"",days," day(s) and ",hours," hour(s)\n", sep = "")
     # cat("\n")
    }
   }
  }
 }
}

# library(matrixcalc)
is.negative.definite <- function (x, tol = 1e-08) 
{
 # if (!is.square.matrix(x)) 
 # stop("argument x is not a square matrix")
 # if (!is.symmetric.matrix(x)) 
 # stop("argument x is not a symmetric matrix")
 # if (!is.numeric(x)) 
 # stop("argument x is not a numeric matrix")
 eigenvalues <- eigen(x, only.values = TRUE)$values
 n <- nrow(x)
 for (i in 1:n) {
  if (abs(eigenvalues[i]) < tol) {
   eigenvalues[i] <- 0
  }
 }
 if (any(eigenvalues >= 0)) {
  return(FALSE)
 }
 return(TRUE)
}

.mlmaximize <- function(theta0, ll, gr = NULL, hess = NULL, alternate = NULL, BHHH = F, level = 0.99, step.back = .Machine$double.eps, reltol = .Machine$double.eps, lmtol = sqrt(.Machine$double.eps), steptol = sqrt(.Machine$double.eps), digits = 4, when.backedup = sqrt(.Machine$double.eps), max.backedup = 17, print.level = 6, only.maximize = FALSE, maxit = 150, n = 100, ...){
 
 theta00 <- theta0
 
 k4 <- length(theta0)
 
 if(print.level >= 6){
  cat("\n=================")
  cat(" Initial values:\n\n", sep = "")
  print(theta0)
 }
 
#  if(print.level >= 2){
#   cat("\n=================")
#   cat(" Maximization:\n\n", sep = "")
#  }
 
 # step.back = 2^-217
 
 ll0 <- ll(theta0, ...)
 ltol <- reltol * (abs(ll0) + reltol)
 typf <- ll0
 theta1 <- theta0
 
 iter <- iter.total <- backedup <- backedups <- wasconcave <- wasconcaves <- 0
 
 if( is.na(ll0) | ll0 == -Inf ){
  if(print.level >= 2){
   cat("Could not compute ll at starting values: trying something else\n")
  }
  iter1 <- backedups
  repeat{
   iter1 <- iter1 + 1
   theta0 <- theta0 * runif(length(theta0), 0.96, 1.05) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(!is.na(ll0)) break
   if(iter1 == 55){
    stop("it's not gonna happen...")
   }
  }
  # backedups <- iter1
 }
 
 delta1 <- gHg <- s1 <- 1
 h1 <- tryCatch( 2, error = function(e) e )
 cant.invert.hess <- FALSE
 
 if(print.level >= 2){
  cat(paste("Iteration ",formatC(iter, width = 3)," (at starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
 }
 
 repeat{
  iter.total <- iter.total + 1
  # cat("backedup = ",backedup," backedups = ",backedups,"\n", sep = "")
  if(print.level >= 6){
   print(theta0)
  }
  
  # cumulate how many times did it backed-up in a row
  if(s1 < when.backedup){
   backedup <- backedup + 1
  } else {
   backedup <- 0
  }
  # print(s1)
  # cumulate how many times was concave
  if( inherits(h1, "error") ){
   wasconcave <- wasconcave + 1
  } else {
   wasconcave <- 0
  }
  
  # try different values if was concave more than @@@ times
  if(wasconcave == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Not concave ",max.backedup," times in a row: trying something else (not concave ",wasconcaves+1," times in total)\n")
   }
   iter <- wasconcave <- backedup <- 0
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.96, 1.05) # not sure what to do here
   ll0 <- ll( theta0, ... )
   # if(print.theta) print(theta0)
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3),"  (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  
  # try different values if backed-up more than @@@ times
  if(backedup == max.backedup){
   # start over
   if(print.level >= 2){
    cat("Backed-up ",max.backedup," times in a row: trying something else (backup-up ",backedups+1," times in total)\n", sep = "")
   }
   iter <- backedup <- wasconcave <- 0
   backedups <- backedups + 1
   wasconcaves <- wasconcaves + 1
   theta0 <- theta0*runif(length(theta0), 0.96, 1.04) # not sure what to do here
   ll0 <- ll( theta0, ... )
   if(print.level >= 6){
    print(theta0)
   }
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  
  # see if it calculated ll
  if( is.na(ll0) | ll0 == -Inf | ll0 == Inf | ll0 == 0 ){
   if(print.level >= 2){
    cat("Could not compute ll: trying something else\n")
   }
   iter1 <- backedups
   repeat{
    iter1 <- iter1 + 1
    # theta0 <- c( cons0, beta0, mu = 0, eta = 0, lnsv2 = -1*iter1/2, lnsu2 = -1*iter1/2)
    theta0 <- theta00*runif(length(theta0), 0.96, 1.04) # not sure what to do here
    ll0 <- ll( theta0, ... )
    if(!is.na(ll0) & ll0 != 0) break
    if(iter1 == 15){
     stop("it's not gonna happen... could not compute at initial and find feasible values")
    }
   }
   iter <- backedup <- 0
   backedups <- iter1
   if(print.level >= 2){
    cat(paste("Iteration ",formatC(iter, width = 3)," (at slightly perturbed starting values):       log likelihood = ",format(ll0, digits = 13),"\n\n", sep = ""), sep = "")
   }
  }
  if(backedups == 5){
   stop("it's not gonna happen... backuped 5 times")
  }
  if(wasconcaves == 5){
   stop("it's not gonna happen... not concave 5 times")
  }
  iter <- iter + 1
  delta3 <- 1
  # step 2: direction vector
  
  # previous theta
  
  if(iter.total > 1) theta1 <- theta0 - s1 * d0
  
  # BHHH (faster, but different):
  # The Hessian is approximated as the negative 
  # of the sum of the outer products of the gradients 
  # of individual observations, 
  # -t(gradient) %*% gradient = - crossprod( gradient )
  g1 <- gr(theta0, ...)
  # print(g1)
  if(!is.null(alternate)) BHHH <- floor(iter/alternate) %% 2 == 1
  h0 <- hess(theta0,  ...)
  # print(h0)
  # check if negative definite
  eigen1 <- eigen( h0 )
  # eigen.tol <- k4 * max(abs(eigen1$values)) * .Machine$double.eps # this is for positive definiteness
  eigen.val <- ifelse(eigen1$values < .Machine$double.eps^.1, 0, eigen1$values)
  # hess.pos.def <- sum(eigen1$values > eigen.tol) == k4
  hess.neg.def <- !any(eigen.val >= 0)
  # 1. replace negative with small ones
  # eigen.val <- ifelse(eigen1$values < 0, .0001, eigen.val)
  # 2. replace negative with absolut values
  eigen.val <- abs(eigen1$values)
  # print(hess.neg.def)
  # make it negative definite if it is not already
  if(!hess.neg.def){
   h0_ <- matrix(0, k4, k4)
   # eigen1 <- eigen( h0 )
   for( i in seq_len( k4 ) ){
    # h0_ <- h0_ - abs(eigen1$values[i]) * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
    h0_ <- h0_ - eigen.val[i] * eigen1$vectors[,i] %*% t(eigen1$vectors[,i])
   }
   h0.old <- h0
   h0 <- h0_
  }
  # print( is.negative.definite(h0) )
  # remember hessian and negative of its inverse from previous iter that could have been inverted
  if( !cant.invert.hess ){
   h0_previous <- h0
   h1_previous <- h1
  } 
  # easier to invert positive definite matrix
  h1 <- tryCatch( qr.solve(-h0), error = function(e) e )
  # check if it can be inverted
  cant.invert.hess <- FALSE
  cant.invert.hess <- inherits(h1, "error")
  if( cant.invert.hess ){
   # print(h1)
   if(print.level >= 2){
    cat(paste("cannot invert Hessian, using eigenvalues\n", sep = ""), sep = "")
   }
   # this was just to get the uninvertable hessian
   # return(list(hess = h0, grad = g1))
   # stop("this")
   # @14@ this
   eig1 <- eigen( -h0_previous )
   d0 <- rep(0, length(theta0))
   # eig2 <- ifelse(eig1$values < eps1, 1, eig1$values)
   for (i in 1:length(eig1$values)){
    d0 <- d0 + (g1%*%eig1$vectors[,i])*eig1$vectors[,i] / eig1$values[i]
   }
   # @14@ could be done easier
   # d0 <- qr.solve(-h0, g1, tol = 1e-10)
   gHg <- sum( g1 * d0)
   # in the part of the ortogonal subspace where the eigenvalues
   # are negative or small positive numbers, use steepest ascent
   # in other subspace use NR step
   # d0 <- ifelse(eigen(-h0, only.values = TRUE)$values < reltol, g1, d0)
   gg <- sqrt( crossprod(g1) )
   gHg <- gg
   # d0 <- g1
   # d0
  } else {
   d0 <- as.vector( h1 %*% g1 )
   gg <- sqrt( crossprod(g1) )
   # h1.old <- solve(-h0.old)
   gHg <- as.vector( t(g1) %*% h1 %*% g1 )
  }
  # gg_scaled <- gg * max( crossprod(theta0), crossprod(theta1) ) / max( abs(ll0), abs(typf))
  # theta_rel <- max( abs(theta0 - theta1) / apply( cbind( abs(theta0),abs(theta1) ), 1, max) )
  theta_rel <- max( abs(theta0 - theta1) / (abs(theta1)+1) )
  
  
  # begin stopping criteria calculated using new values of g1 and h1
  if(s1 > when.backedup*10^-100 & delta1 != 17.17){ # if(s1 > when.backedup*10^-100 & !cant.invert.hess){
   if(abs(gHg) < lmtol & iter.total > 1){
    if(print.level >= 2){
     cat("\nConvergence given g inv(H) g' = ",abs(gHg)," < lmtol\n", sep = "")
    }
    break
   }
   if(theta_rel < steptol & iter.total > 2){
    # print(theta_rel)
    if(print.level >= 2){
     cat("\nConvergence given relative change in parameters = ",theta_rel," < steptol\n", sep = "")
    }
    break
   }
  }
  # end stopping criteria
  # use steepest ascent when backed-up
  if(s1 < when.backedup*10^-3){
   # eig1 <- eigen( -h0 )
   d0 <- g1
   # d0 <- ifelse(eig1$values < reltol, g1, d0)
   # theta0 <- theta0 - 1 * d0
  }
  # print(d0)
  # step 3: new guess
  # a: s = 1
  # b: funct(theta0 + d0) > funct(theta0)
  s1 <- 1
  theta1 <- theta0 + s1 * d0
  # print(12)
  # print(theta1)
  ll1 <- ll( theta1, ... )
  # print(13)
  delta2 <- ll1 - ll0
  flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
  # begin Cases
  if( flag ){
   # begin Case 1: f(theta1) > f(theta0)
   ll.temp <- ll0
   # check if s1 = 2, 3, ... increases f even more
   while( flag ){
    if(print.level >= 6){
     cat(paste("\t\tCase 1: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
    }
    ll0 <- ll1
    s1 <- s1 + 1
    theta1 <- theta0 + s1 * d0
    ll1 <- ll( theta1, ... )
    delta2 <- ll1 - ll0
    flag <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
   }
   # overall delta
   delta1 <- ll0 - ll.temp
   delta_rel <- abs(delta1 / ll.temp)
   # print(delta_rel)
   s1 <- s1 - 1
   # overwrite the values
   theta0 <- theta0 + s1 * d0
   # end Case 1: f(theta1) > f(theta0)
  } else {
   # begin Case 2: f(theta1) < f(theta0)
   # check only if s1=1/2 increases f
   s1 <- 0.5
   theta1 <- theta0 + s1 * d0
   ll1 <- ll( theta1, ... )
   # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
   delta2 <- ll1 - ll0
   if(print.level >= 6){
    cat(paste("\t\tCase 2: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
   }
   flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
   # end Case 2: f(theta1) < f(theta0)
   if( flag2 ){
    # begin Case 2a: f(theta1) > f(theta0)
    ll.temp <- ll0
    # check if s1=1/2^2,1/2^3,... increases f even more
    while( flag2 ){
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2a: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     ll0 <- ll1
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     # cat(" ll1 = ",ll1,", ll0 = ",ll0,", s = ",s1,"\n", sep = "")
     delta2 <- ll1 - ll0
     flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > ltol)
    }
    # overall delta
    delta1 <- ll0 - ll.temp
    delta_rel <- abs(delta1 / ll.temp)
    # print(delta_rel)
    s1 <- 2 * s1
    # overwrite the values
    theta0 <- theta0 + s1 * d0
    # end Case 2a: f(theta1) > f(theta0)
   } else {
    # begin Case 2b: f(theta1) < f(theta0)
    ll.temp <- ll0
    # try s1=1/2^2,1/2^3,... so that f(theta1) > f(theta0)
    while ( !flag2 & s1 > step.back ){
     s1 <- 0.5 * s1
     theta1 <- theta0 + s1 * d0
     ll1 <- ll( theta1, ... )
     delta2 <- ll1 - ll0
     if(print.level >= 6){
      cat(paste("\t\t\tCase 2b: s = ",s1,", delta = ",delta2,"\n", sep = ""), sep = "")
     }
     flag2 <- (!is.na(delta2) & delta2 != -Inf & delta2 > 0)
    } 		
    if( !flag2 | s1 < step.back ){
     # stop("provide different starting values")
     delta1 <- 17.17
    } else {
     # overwrite the values
     delta1 <- delta2
     delta_rel <- abs(delta1 / ll.temp)
     ll0 <- ll1
     theta0 <- theta0 + s1 * d0
    }
    # end Case 2b: f(theta1) < f(theta0)
   }
  }
  
  
  if(print.level >= 2){
   if( cant.invert.hess ){
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (not concave)\n", sep = ""), sep = "")
   } else if (s1 <= when.backedup) {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13)," (backed up)\n", sep = ""), sep = "")
   } else {
    cat(paste("Iteration ",formatC(iter, width = 3)," (hessian is ",ifelse(BHHH, "BHHH", "analytical"),", ",formatC(iter.total, width = 3)," in total):   log likelihood = ",format(ll0, digits = 13),"\n", sep = ""), sep = "")
   }
  }
  
  # printing criteria
  if(print.level >= 5){
   if( cant.invert.hess ){
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; quasi-gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 5.5){
     print(theta0)
     cat("\n")
    } 
   } else {
    cat(paste(" (in iter ",formatC(iter, width = 3),": delta = ",format(delta1, digits = 6),"; s = ",format(s1, digits = 6),"; gHg = ",format(gHg, digits = 6),"; theta_rel_change = ",format(theta_rel, digits = 6),")\n\n", sep = ""), sep = "")
    if(print.level >= 5.5){
     print(theta0)
     cat("\n")
    }
   }
  }
  # print(s1)
  if(s1 > when.backedup^2 & !cant.invert.hess){ # if(s1 > when.backedup^2 & !cant.invert.hess){
   # ltol <- reltol * (abs(ll0) + reltol)
   # print(cant.invert.hess)
   if(delta1 > 0 & !is.na(delta_rel) & delta_rel < ltol & iter.total > 1){
    if(print.level >= 2){
     cat("\nConvergence given relative change in log likelihood = ",delta_rel," < ltol\n", sep = "")
    }
    break
   }
  }
  if(iter.total > maxit){
   cat("\n Maximum number of iterations (",maxit,") reached without convergence\n", sep = "")
   break
  }
 } # end repeat
 
 if( !only.maximize & !cant.invert.hess){		
  names(ll0) <- NULL
  colnames(h1) <- rownames(h1) <- names(g1) <- names(theta0)
  
  # sqrt(crossprod(g1))
  
  b0 <- theta0
  sd0 <- sqrt( diag( h1 ) )
  t0 <- b0 / sd0
  p0 <- pt(abs(t0), n-length(b0), lower.tail = FALSE) * 2
  t10 <- qt((1-0.01*level)/2, n-length(b0), lower.tail = FALSE)
  t17 <- cbind( b0, sd0, t0, p0, b0 - t10*sd0, b0 + t10*sd0)
  # t17 <- cbind( b0, sd0, t0, p0)
  colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", paste("",level,"_CI_LB", sep = ""), paste("",level,"_CI_UB", sep = ""))
  # colnames(t17) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  # t17
  
  if(print.level >= 2){
   cat(paste("\nFinal log likelihood = ",format(ll0, digits = 5),"\n\n", sep = ""), sep = "")
  }
  # cat(paste("Stoc. frontier normal/",distribution,"\n", sep = ""), sep = "")
  if(print.level >= 5.5){
   cat("\nCoefficients:\n\n", sep = "")
   printCoefmat(t17[,1:4], digits = digits)
  }
  
  return(list(par = theta0, table = t17, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel))
 } else {
  return(list(par = theta0, gradient = g1, vcov = h1, ll = ll0, gg = gg, gHg = gHg, theta_rel_ch = theta_rel))
 }
 
 
}
.su1 <- function(x, mat.var.in.col = TRUE, digits = 4, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), print = TRUE){
 
 xvec2 <- xvec1 <- FALSE
 
 if(is.matrix(x)){
  if(min(dim(x)) == 1){
   xvec1 <- TRUE
   # mynames <- deparse(substitute(x))
   x <- as.vector(x)  
  } else {
   if(!mat.var.in.col){
    x <- t(x)
   }
   mynames <- paste("Var", seq_len(ncol(x)), sep = "")
  }
  # print(x)
  # mynames <- colnames(x)
 } # end if matrix
 
 if(is.vector(x)){
  xvec2 <- TRUE
  mynames <- deparse(substitute(x))
  x <- data.frame(Var1 = x)  
 } # end if vector
 
 # cat("nymanes", sep ="")
 # print(mynames)
 
 if(!is.vector(x) & !is.matrix(x) & !is.data.frame(x)){
  stop("Provide vector, matrix, or data.frame")
 } else {
  t1 <- apply(x, 2, function(x) c(Obs = length(x), NAs = sum(is.na(x)), Mean = mean(x, na.rm = TRUE), StDev = sd(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), quantile( x, probs = probs, na.rm = TRUE ), Max = max(x, na.rm = TRUE)))
  # print(class(t1))
  # print(dim(t1))
  if(xvec2 & !xvec1) colnames(t1) <- mynames
  if(print) print(t(t1), digits = digits)    
 }
 return(t1)
}

