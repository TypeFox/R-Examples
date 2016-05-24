.onLoad  <- function(libname, pkgname) {
    library.dynam ("kknn", pkgname, libname)
}


simulation <- function(formula, data, runs = 10, train = TRUE, k = 11, ...)
{
   .Deprecated("train.kknn or cv.kknn", "kknn", old="simulation")
   mf <- model.frame(formula, data=data)
   y <- model.response(mf)
   MISCLASS<-numeric(runs)
   MEAN.ABS<-numeric(runs)
   MEAN.SQU<-numeric(runs)
   for(i in 1:runs)
     {
     set.seed(i)
     m<-dim(data)[1]
     val<-sample(1:m, size=round(m/3), replace=FALSE, prob=rep(1/m, m)) 
     learn<-data[-val,]
     valid<-data[val,]
     ytmp<-y[val]
     
     if(train){
     	fit = train.kknn(formula , learn, kmax = k, ...)
     	pred <- predict(fit, valid)}
     if(!train) pred <- predict(kknn(formula, learn, valid, k = k, ...))
     if(is.factor(y)) MISCLASS[i]<-sum(ytmp != pred)/dim(valid)[1]
     if(is.numeric(y) | is.ordered(y)) MEAN.ABS[i]<-sum(abs(as.numeric(ytmp) -
		as.numeric(pred)))/dim(valid)[1]
     if(is.numeric(y) | is.ordered(y)) MEAN.SQU[i]<-sum((as.numeric(ytmp) -
		as.numeric(pred))^2)/dim(valid)[1]
     }  
   if(is.numeric(y)){ 
   		result <- matrix(data=c(mean(MEAN.ABS), sd(MEAN.ABS),
		mean(MEAN.SQU), sd(MEAN.SQU)), nrow=2, ncol=2)
		colnames(result)<-c("absolute distance", "squared distance")
		rownames(result)<-c("mean", "sd")                         
   }
   if(is.ordered(y)){ 
   		result<-matrix(data=c(mean(MISCLASS), sd(MISCLASS), 
		mean(MEAN.ABS), sd(MEAN.ABS),
		mean(MEAN.SQU), sd(MEAN.SQU)), 
		nrow=2, ncol=3)
		colnames(result)<-c("misclassification","absolute distance", "squared distance")
		rownames(result)<-c("Mean", "sd")
		} 
   if(is.factor(y) & !is.ordered(y)){ 
		result<-matrix(data=c(mean(MISCLASS), sd(MISCLASS)),	nrow=2, ncol=1)   
   		colnames(result)<-"misclassification"
   		rownames(result)<-c("mean", "sd")
   		} 
   result
}


contr.dummy <- function (n, contrasts = TRUE) 
{
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("contrasts are not defined for 0 degrees of freedom")
    }
    else levels <- n
    	lenglev <- length(levels)
        cont <- array(0, c(lenglev, lenglev), list(levels, levels))
        cont[col(cont) == row(cont)] <- 1
    cont
}


contr.ordinal <- function (n, contrasts = TRUE) 
{
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("contrasts are not defined for 0 degrees of freedom")
    }
    else levels <- n
    	lenglev <- length(levels)
        cont <- array(0.5, c(lenglev, lenglev - 1), list(levels, NULL))
        cont[lower.tri(cont)] <- -0.5
    cont
}


contr.metric <- function(n, contrasts = TRUE) 
{	
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- 1:n
        else stop("contrasts are not defined for 0 degrees of freedom")
    }
    else levels <- n
    	lenglev <- length(levels)
        cont <- array((1:lenglev)-(1+lenglev)/2 , c(lenglev,1), list(levels,NULL)) 
    cont	
}



contr.int <- function (n, contrasts = TRUE) 
{
    if (length(n) <= 1) {
        if (is.numeric(n) && length(n) == 1 && n > 1) 
            levels <- as.integer(1:n)
        else stop("contrasts are not defined for 0 degrees of freedom")
    }
    else levels <- n
    lenglev <- length(levels)
    cont <- array(as.integer(1:lenglev), c(lenglev, 1), 
                  list(levels, NULL))
    cont
}


optKernel <- function(k, d=1){
   1/k*(1 + d/2 - d/(2*k^(2/d)) * ( (1:k)^(1+2/d) - (0:(k-1))^(1+2/d)  ))
}


kknn <-  function (formula = formula(train), train, test, na.action=na.omit(), k = 7, distance = 2, 
    kernel = "optimal", ykernel = NULL, scale=TRUE, contrasts=c('unordered'="contr.dummy",ordered="contr.ordinal")) 
{
    if(is.null(ykernel)) ykernel=0

    weight.y = function(l=1,diff = 0){
        k=diff+1
        result=matrix(0,l,l)
        diag(result)=k
        for(i in 1:(k-1)){
            for(j in 1:(l-i)){
                result[j,j+i]=k-i
                result[j+i,j]=k-i
            }
        }
        result  
    }

    kernel <- match.arg(kernel, c("rectangular", "triangular", "epanechnikov", "biweight", 
        "triweight", "cos", "inv", "gaussian", "rank", "optimal"), FALSE)
    
    ca <- match.call()
    response = NULL
    old.contrasts<-getOption('contrasts')
    options(contrasts=contrasts)

    formula = as.formula(formula)

    mf <- model.frame(formula, data = train)    
    mt <- attr(mf, "terms")
#reformulate(, intercept = FALSE
    
    mt2 <- delete.response(mt)

    cl <- model.response(mf)

    d <- sum(attr(mt, "order"))

    if(is.ordered(cl)) {
        response<-"ordinal"
        lev <- levels(cl)
	}
    if(is.numeric(cl)) response<-"continuous"
    if(is.factor(cl) & !is.ordered(cl)){
        response<-"nominal"
	    lev <- levels(cl)
	}
	
    if(distance<=0)stop('distance must >0')
    if(k<=0)stop('k must >0')

    learn <- model.matrix(mt, mf)
    valid <-model.matrix(mt2,test)

    m <- dim(learn)[1]
    p <- dim(valid)[1]
    q <- dim(learn)[2]

    ind <- attributes(learn)$assign 

    d.sd<-numeric(length(ind))+1
    we<-numeric(length(ind))+1
    
    d.sd = apply(learn, 2, stats::var)
    for (i in unique(ind)){
    	d.sd[ind==i] = sqrt(mean(d.sd[ind==i]))
    	we[ind==i] = 1/sum(ind==i)
    	}

    we[d.sd==0]=0
    d.sd[d.sd==0]=1

    if(scale){
# change 5.3.2013      
        learn <- sweep(learn, 2L, d.sd, "/", check.margin = FALSE) 
        valid <- sweep(valid, 2L, d.sd, "/", check.margin = FALSE) 
    } 
# ordering allows branch and bound in distance computation
    ord = order(we * apply(learn, 2, sd), decreasing=TRUE)

    we = we[ord]
    learn = learn[,ord, drop=FALSE]
    valid = valid[,ord, drop=FALSE]

    Euclid <- FALSE
    if(distance==2) Euclid <- TRUE
    if(Euclid) dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
        as.integer(m), as.integer(p), as.integer(q), 
        dm=double((k+1L) * p), cl=integer((k+1L) * p), k=as.integer(k+1), 
        as.double(distance),as.double(we), PACKAGE='kknn')

    else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
        as.integer(m), as.integer(p), as.integer(q), 
        dm=double((k+1L) * p), cl=integer((k+1L) * p), k=as.integer(k+1), 
        as.double(distance),as.double(we), PACKAGE='kknn')
    
    D <- matrix(dmtmp$dm, nrow = p, ncol = k + 1)
    C <- matrix(dmtmp$cl, nrow = p, ncol = k + 1)
    maxdist <- D[, k + 1]
    maxdist[maxdist < 1.0e-6] <- 1.0e-6
    D <- D[, 1:k]
    C <- C[, 1:k]+1
    CL <- matrix(cl[C], nrow = p, ncol = k)     
    
    if(response!="continuous"){
        l <- length(lev)
        weightClass <- matrix(0, p, l)
    }
    if(response=="continuous"){
        weightClass <- NULL
    }

    W <- D/maxdist
    W <- pmin(W,1-(1e-6))
    W <- pmax(W,1e-6)

#	
# Kernels
#	
    if (kernel=="rank") W <- (k+1)-t(apply(as.matrix(D),1,rank))
    if (kernel=="inv") W <- 1/W
    if (kernel=="rectangular") W <- matrix(1,nrow = p, ncol = k)
    if (kernel=="triangular") W <- 1-W	 	
    if (kernel=="epanechnikov") W <- 0.75*(1-W^2)
    if (kernel=="biweight") W <- dbeta((W+1)/2,3,3)	 	
    if (kernel=="triweight") W <- dbeta((W+1)/2,4,4)	 	
    if (kernel=="cos") W <- cos(W*pi/2)
    if (kernel=="triweights") W <- 1
    if (kernel=="gaussian"){
	    alpha=1/(2*(k+1))
	    qua=abs(qnorm(alpha))
	    W=W*qua
        W = dnorm(W, sd = 1)
    }
    if (kernel == "optimal") {
        W = rep(optKernel(k, d=d), each=p)
    } 
    W <- matrix(W, p, k)

if(response!="continuous"){
    for (i in 1:l) {
		weightClass[, i] <- rowSums(W * (CL == lev[i]))	
    }
    weightClass <- weightClass/rowSums(weightClass)	
	colnames(weightClass) <- lev
    }
    
if (response=="ordinal") {
   
    blub = length(lev)
    weightClass= weightClass%*%weight.y(blub,ykernel)
    weightClass <- weightClass/rowSums(weightClass)	
    weightClass <- t(apply(weightClass, 1, cumsum))
    colnames(weightClass) <- lev
       	
    fit <- numeric(p)
    for (i in 1:p) fit[i] <- min((1:l)[weightClass[i, ] >= 0.5])
    fit <- ordered(fit, levels = 1:l, labels = lev)
    }
	if(response=="nominal"){ 
		fit <- apply(weightClass, 1, order, decreasing = TRUE)[1,]
		fit <- factor(fit, levels = 1:l, labels = lev)

	if(kernel=="rectangular" && k>1){
		blub <- apply(weightClass, 1, rank, ties.method = "max")
		indices = (1:p)[colSums(blub==l)>1]
		blub=t(blub)
		nM = matrix(0,p,l) 
		colnames(nM)=lev
		for(i in 1:l) nM[,i] = apply((CL==lev[i]) %*% diag(1:k) ,1,max)

		nM = (blub==l)*nM  
		nM[nM==0] <- k+1
		fitv = numeric(p)
		for(i in indices) fitv[i]=which(nM[i,]==min(nM[i,]))
		fit[indices] <- factor(fitv[indices], levels = 1:l, labels = lev)
		}
	}
	if(response=="continuous") fit <- rowSums(W*CL)/pmax(rowSums(W), 1e-6) 
    #fit <- rowSums(W*CL)/sapply(rowSums(W),'max',1e-6) 
    options('contrasts'=old.contrasts)	
	
    result <- list(fitted.values=fit, CL=CL, W=W, D=D, C=C, prob=weightClass, response=response, distance=distance, call=ca, terms=mt)	
    class(result)='kknn'
    result
}


# valid=NULL fuer leave one out?
# include in kknn, train.kknn? 
kknn.dist <- function(learn, valid, k = 10, distance = 2) 
{
    m <- dim(learn)[1]
    p <- dim(valid)[1]
    q <- dim(learn)[2]
  
    we <- rep(1.0, q)
  
    ord = order(we * apply(learn, 2, sd), decreasing=TRUE)
  
    learn = learn[,ord, drop=FALSE]
    valid = valid[,ord, drop=FALSE]
    
    Euclid <- FALSE
    if(distance==2) Euclid <- TRUE
    if(Euclid) dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
        as.integer(m), as.integer(p), as.integer(q), 
        dm=double(k * p), cl=integer(k * p), k=as.integer(k), 
        as.double(distance),as.double(we), PACKAGE='kknn')
    else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
        as.integer(m), as.integer(p), as.integer(q), 
        dm=double(k * p), cl=integer(k * p), k=as.integer(k), 
        as.double(distance),as.double(we), PACKAGE='kknn')
    D <- matrix(dmtmp$dm, nrow = p, ncol = k)
    C <- matrix(dmtmp$cl, nrow = p, ncol = k) + 1L
    list(C, D)
}  



print.kknn <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
   	cat("Response: ",deparse(x$response),"\n",sep="")
}



summary.kknn <- function(object, ...)
{

	cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
	cat("Response: ",deparse(object$response),"\n",sep="")
	digits = max(3, getOption("digits") - 3)
	if(object$response!="continuous")print(data.frame(fit=object$fitted.value,prob=object$prob),digits)
	fit <- object$fit
}


predict.kknn <- function(object, type = 'raw', ...) 
{ 
    pred <- switch(type, 
        raw = object$fit,
        prob = object$prob,
        stop('invalid type for prediction'))
    pred
}


predict.train.kknn <- function (object, newdata, ...) 
{
    if (missing(newdata)) 
        return(predict(object, ...)) 
    #    return(object$fit)
    res <- kknn(formula(terms(object)), object$data, newdata, 
        k = object$best.parameters$k, kernel = object$best.parameters$kernel, 
        distance = object$distance)
    return(predict(res, ...))
}


#predict.kknn <- function(object,...)return(object$fit)


#predict.train.kknn <- function (object, newdata, ...) 
#{
#    if (missing(newdata)) 
#        return(object$fit)
#    res = kknn(formula(terms(object)), object$data, newdata, 
#        k = object$best.parameters$k, kernel = object$best.parameters$kernel, 
#        distance = object$distance)
#   return(predict(res))
#}


train.kknn = function (formula, data, kmax = 11, ks = NULL, distance = 2, kernel = "optimal", ykernel = NULL, scale=TRUE, 
    contrasts = c(unordered = "contr.dummy", ordered = "contr.ordinal"), ...) 
{
	if(is.null(ykernel)) ykernel=0
    weight.y = function(l = 1, diff = 0) {
        k = diff + 1
        result = matrix(0, l, l)
        diag(result) = k
        for (i in 1:(k - 1)) {
            for (j in 1:(l - i)) {
                result[j, j + i] = k - i
                result[j + i, j] = k - i
            }
        }
        result
    }
#    if (is.null(kernel)) 
#        kernel = "triangular"
    kernel <- match.arg(kernel, c("rectangular", "triangular", "epanechnikov", "biweight", 
        "triweight", "cos", "inv", "gaussian", "rank", "optimal"), TRUE)
    
    if (is.null(ks)) {
      ks <- 1:kmax
      nk <- kmax
    } else {
      ks <- sort(ks)
      nk <- length(ks)
      kmax <- max(ks)
    }
    
    call <- match.call()
    mf <- model.frame(formula, data = data)
    mt <- attr(mf, "terms")

    y <- model.response(mf)
    cl <- model.response(mf)
    old.contrasts <- getOption("contrasts")
    options(contrasts = contrasts)
    mm.data <- model.matrix(mt, mf)
  
    d <- sum(attr(mt, "order"))

    r <- length(kernel)
    MISCLASS <- matrix(nrow = nk, ncol = r, dimnames = list(ks, 
        kernel))
    MEAN.ABS <- matrix(nrow = nk, ncol = r, dimnames = list(ks, 
        kernel))
    MEAN.SQU <- matrix(nrow = nk, ncol = r, dimnames = list(ks, 
        kernel))
    P <- list(nk * r)
    m <- dim(mm.data)[1]
    q <- dim(mm.data)[2]
    p <- m

    ind <- attributes(mm.data)$assign
  
    d.sd <- numeric(length(ind)) + 1
    we <- numeric(length(ind)) + 1

#    for (i in 1:max(ind)) {
#        d.sd[ind == i] = sqrt(mean(diag(cov(as.matrix(mm.data[, ind == i])))))
#        we[ind == i] = 1/sum(ind == i)
#    }

    d.sd = apply(mm.data, 2, stats::var)
    for (i in unique(ind)){
    	d.sd[ind==i] = sqrt(mean(d.sd[ind==i]))
    	we[ind==i] = 1/sum(ind==i)
    	}

    we[d.sd == 0] = 0
    d.sd[d.sd == 0] = 1
#    mm.data <- t(t(mm.data)/d.sd) 
 
#    raus 5.3.2013  
#    if(scale) mm.data <- mm.data %*% diag(1/d.sd)         
  	if(scale) mm.data <- sweep(mm.data, 2L, d.sd, "/", check.margin = FALSE) 
    ord = order(we * apply(mm.data, 2, sd), decreasing=TRUE)
# ordering
    mm.data <- mm.data[, ord, drop=FALSE]  
    we <- we[ord]

    Euclid <- FALSE
    if(distance==2) Euclid <- TRUE
    if(Euclid) dmtmp <- .C("dmEuclid", as.double(mm.data), as.double(mm.data), 
        as.integer(m), as.integer(p), as.integer(q), dm = double((kmax + 2L) * p), 
        cl = integer((kmax + 2L) * p), k = as.integer(kmax + 2), as.double(distance), 
        as.double(we), PACKAGE = "kknn")
    else dmtmp <- .C("dm", as.double(mm.data), as.double(mm.data), 
        as.integer(m), as.integer(p), as.integer(q), dm = double((kmax + 2L) * p), 
        cl = integer((kmax + 2L) * p), k = as.integer(kmax + 2), as.double(distance), 
        as.double(we), PACKAGE = "kknn")
    D <- matrix(dmtmp$dm, nrow = p, ncol = kmax + 2)
    C <- matrix(dmtmp$cl, nrow = p, ncol = kmax + 2)
    C <- C + 1
    CL <- matrix(cl[C], nrow = p, ncol = kmax + 2)  # y statt cl
    D <- D[, -1]
    C <- C[, -1]
    CL <- CL[, -1]
    if (is.ordered(y)) {
        response <- "ordinal"
        lev <- levels(y)
        l <- length(lev)
        weightClass <- matrix(0, m, l)
    }
    if (is.numeric(y)) {
        response <- "continuous"
        weightClass <- NULL
    }
    if (is.factor(y) & !is.ordered(y)) {
        response <- "nominal"
        lev <- levels(y)
        l <- length(lev)
        weightClass <- matrix(0, m, l)
    }

    for (k_i in 1:nk) {
        j <- ks[k_i]
        maxdist <- D[, j + 1]
        maxdist[maxdist < 1.0e-06] = 1.0e-06
        V <- D[, 1:j]/ maxdist # sapply(maxdist, "max", 1e-06)
#        V <- D[, 1:j]/sapply(maxdist, "max", 1e-06)
        V <- pmin(V, 1 - (1e-06))
        V <- pmax(V, 1e-06)
        for (s in 1:r) {
            if (kernel[s] == "rank") 
                W <- (j + 1) - t(apply(as.matrix(V), 1, rank))
            if (kernel[s] == "inv") 
                W <- 1/V
            if (kernel[s] == "rectangular") 
                W <- matrix(1, nrow = m, ncol = j)
            if (kernel[s] == "triangular") 
                W <- 1 - V
            if (kernel[s] == "epanechnikov") 
                W <- 0.75 * (1 - V^2)
            if (kernel[s] == "biweight") 
                W <- dbeta((V + 1)/2, 3, 3)
            if (kernel[s] == "triweight") 
                W <- dbeta((V + 1)/2, 4, 4)
            if (kernel[s] == "cos") 
                W <- cos(V * pi/2)
            if (kernel[s] == "gaussian") {
            	v <- j + 1
         	    alpha = 1/(2 * v)
             	qua = abs(qnorm(alpha))
             	W = V*qua
            	W = apply(as.matrix(W), 2, dnorm)
            }
            if (kernel[s] == "optimal") {
                W = rep(optKernel(j,d), each=m)
            }
            W <- matrix(W, m, j)
            if (response != "continuous") {
                for (i in 1:l) {
                  weightClass[, i] <- rowSums(W * (matrix(CL[, 
                    1:j], m, j) == lev[i]))
                }
                weightClass <- weightClass/rowSums(weightClass)
                colnames(weightClass) <- lev
            }
            if (response == "ordinal") {
                blub = length(lev)
                weightClass = weightClass %*% weight.y(blub, ykernel)
                weightClass <- weightClass/rowSums(weightClass)
                weightClass <- t(apply(weightClass, 1, cumsum))
                colnames(weightClass) <- lev
                fit <- numeric(m)             
				fit <- ((l+1)-(weightClass >= 0.5)%*%(numeric(l)+1))
                fit <- ordered(fit, levels = 1:l, labels = lev)
            }
            if (response == "nominal") {
				lwc = length(weightClass)
                fit <- apply(weightClass, 1, order, decreasing = TRUE)[1, ]
                fit <- factor(fit, levels = 1:l, labels = lev)
            }
            if (response == "continuous") {
#                fit <- rowSums(W * (matrix(CL[, 1:j], m, j)))/sapply(rowSums(matrix(W, m, j)), "max", 1e-06)
                fit <- rowSums(W * (matrix(CL[, 1:j], m, j)))/pmax(rowSums(matrix(W, m, j)), 1e-06)
                weightClass = fit
            }
            attr(fit, "kernel") = kernel[s]
            attr(fit, "k") = j
            P[[k_i + (s - 1) * nk]] = fit

        }
    }

    for (k_i in 1:nk) {
        j <- ks[k_i]
        for (s in 1:r) {
            if (is.factor(y)) 
                MISCLASS[k_i, s] <- sum(y != P[[k_i + (s - 1) * nk]])/m
            if (is.numeric(y) | is.ordered(y)) 
                MEAN.ABS[k_i, s] <- sum(abs(as.numeric(y) - as.numeric(P[[k_i + 
                  (s - 1) * nk]])))/m
            if (is.numeric(y) | is.ordered(y)) 
                MEAN.SQU[k_i, s] <- sum((as.numeric(y) - as.numeric(P[[k_i + 
                  (s - 1) * nk]]))^2)/m
        }
    }
    if (response == "nominal") 
        best <- which(MISCLASS == min(MISCLASS), arr.ind = TRUE)
    if (response == "ordinal") 
        best <- which(MEAN.ABS == min(MEAN.ABS), arr.ind = TRUE)
    if (response == "continuous") 
        best <- which(MEAN.SQU == min(MEAN.SQU), arr.ind = TRUE)
    best.parameters = list(kernel = kernel[best[1, 2]], k = ks[best[1, 
        1]])

    options('contrasts'=old.contrasts)
    
    result = list(MISCLASS = MISCLASS, MEAN.ABS = MEAN.ABS, MEAN.SQU = MEAN.SQU, 
        fitted.values = P, best.parameters = best.parameters, 
        response = response, distance = distance, call = call, 
        terms = mt, data = data)
    class(result) = c("train.kknn", "kknn")
    result
}


print.train.kknn <- function(x, ...)
{
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
	cat("Type of response variable: ", x$response,"\n", sep = "")
	if(x$response=='continuous'){cat("minimal mean absolute error: ",min(x$MEAN.ABS),"\n", sep = "")
		cat("Minimal mean squared error: ",min(x$MEAN.SQU),"\n", sep = "")}
	if(x$response=='nominal'){
		cat("Minimal misclassification: ",min(x$MISCLASS),"\n", sep = "")}
	if(x$response=='ordinal'){cat("minimal mean absolute error: ",min(x$MEAN.ABS),"\n", sep = "")
		cat("Minimal mean squared error: ",min(x$MEAN.SQU),"\n", sep = "")
		cat("Minimal misclassification: ",min(x$MISCLASS),"\n", sep = "")
		}	
	cat("Best kernel: ", x$best$kernel,"\n", sep = "")
	cat("Best k: ", x$best$k,"\n", sep = "")	
}


summary.train.kknn <- function(object, ...)
{
	cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
	cat("Type of response variable: ", object$response,"\n", sep = "")
	if(object$response=='continuous'){cat("minimal mean absolute error: ",min(object$MEAN.ABS),"\n", sep = "")
		cat("Minimal mean squared error: ",min(object$MEAN.SQU),"\n", sep = "")}
	if(object$response=='nominal'){
		cat("Minimal misclassification: ",min(object$MISCLASS),"\n", sep = "")}
	if(object$response=='ordinal'){cat("minimal mean absolute error: ",min(object$MEAN.ABS),"\n", sep = "")
		cat("Minimal mean squared error: ",min(object$MEAN.SQU),"\n", sep = "")
		cat("Minimal misclassification: ",min(object$MISCLASS),"\n", sep = "")
		}	
	cat("Best kernel: ", object$best$kernel,"\n", sep = "")
	cat("Best k: ", object$best$k,"\n", sep = "")	
}


plot.train.kknn <-function(x,...){
	if(x$response=='continuous'){		
		legend.text = colnames(x$MEAN.ABS)
		m = 1:length(colnames(x$MEAN.ABS))
		matplot(x = as.integer(rownames(x$MEAN.SQU)),
		        y = x$MEAN.SQU, xlab="k", ylab="mean squared error",pch = m,...)
		xy <- par("usr")
		legend(xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = legend.text, xjust = 1, yjust = 1,col=m,pch=m)
		}
	if(x$response=='ordinal'){
		legend.text = colnames(x$MISCLASS)
		m = 1:length(colnames(x$MISCLASS))
		matplot(x = as.integer(rownames(x$MEAN.ABS)),
		        y = x$MEAN.ABS, xlab="k", ylab="mean absolute error",pch = m,...)
		xy <- par("usr")
		legend(xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = legend.text, xjust = 1, yjust = 1,col=m,pch=m)
		}
	if(x$response=='nominal'){
		legend.text = colnames(x$MISCLASS)
		m = 1:length(colnames(x$MISCLASS))
		matplot(x = as.integer(rownames(x$MISCLASS)),
		        y = x$MISCLASS, xlab="k", ylab="misclassification",pch = m,...)
		xy <- par("usr")
		legend(xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = legend.text, xjust = 1, yjust = 1,col=m,pch=m)
		}	
}


cv.kknn <- function(formula, data, kcv = 10, ...)
{
  mf <- model.frame(formula, data=data) 
  # terms(formula, data = data) keine kopie der Daten?
  y <- model.response(mf)                 
  l <- length(y)    # nrow(data)                  
  val<-sample(kcv, size=l, replace=TRUE) 
  yhat <- numeric(l)
  for(i in 1:kcv){
    m<-dim(data)[1]
    learn<-data[val!=i,]
    valid<-data[val==i,]
    fit = kknn(formula , learn, valid, ...)
    yhat[val==i] <- predict(fit)
  }  
  if(is.factor(y)) MISCLASS <- sum(y != yhat)/l
  if(is.numeric(y) | is.ordered(y)) MEAN.ABS <- sum(abs(as.numeric(y) - as.numeric(yhat)))/l
  if(is.numeric(y) | is.ordered(y)) MEAN.SQU <- sum((as.numeric(y) - as.numeric(yhat))^2)/l
  if(is.numeric(y)) result <- c(MEAN.ABS, MEAN.SQU)
  if(is.ordered(y)) result <- c(MISCLASS, MEAN.ABS, MEAN.SQU)
  if(is.factor(y) & !is.ordered(y)) result<-MISCLASS 
  list(cbind(y=y, yhat=yhat), result)
}


prepare.Discrete <- function(data){
  if(class(data)=="factor")return(as.matrix(unclass(data)))
  if(class(data)=="data.frame")return(as.matrix(data.frame(lapply(data,unclass))))
}


kknn.dist2 <- function(learn, valid, learnD=NULL, validD=NULL, k = 10, distance = 2) 
{
  m <- dim(learn)[1]
  p <- dim(valid)[1]
  q <- dim(learn)[2]
  
  we <- rep(1.0, q)
  ord = order(we * apply(learn, 2, sd), decreasing=TRUE)
  
  if(!is.null(validD)){
    q2 <- dim(learnD)[2]
    we2 <- rep(1.0, ncol(learnD))
  } 
  
  learn = learn[,ord, drop=FALSE]
  valid = valid[,ord, drop=FALSE]
  
  Euclid <- FALSE
  if(distance==2) Euclid <- TRUE
  if(Euclid){
      if(is.null(learnD))
      dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
                         as.integer(m), as.integer(p), as.integer(q), 
                         dm=double(k * p), cl=integer(k * p), k=as.integer(k), 
                         as.double(distance),as.double(we), PACKAGE='kknn')
      else dmtmp <- .C("dmEuclid2", as.double(learn), as.double(valid), 
          as.integer(learnD), as.integer(validD),             
          as.integer(m), as.integer(p), as.integer(q), as.integer(q2), 
          dm=double(k * p), cl=integer(k * p), k=as.integer(k), 
          as.double(distance), as.double(we), as.double(we2), PACKAGE='kknn')  
#dmEuclid2( int *n, int *m, int *p, int *p2, double *dm, int *cl, int *k, double *mink, double *weights, double *weights2)
  }
  else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
                   as.integer(m), as.integer(p), as.integer(q), 
                   dm=double(k * p), cl=integer(k * p), k=as.integer(k), 
                   as.double(distance),as.double(we), PACKAGE='kknn')
  D <- matrix(dmtmp$dm, nrow = p, ncol = k)
  C <- matrix(dmtmp$cl, nrow = p, ncol = k) + 1L
  list(C, D)
}  
