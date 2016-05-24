added.variable.plots<-function(f,intercept=TRUE,...)
{
	UseMethod("added.variable.plots")
}

#S3 method for class "lm"
added.variable.plots.lm <-function(f,intercept=TRUE,...)
{
# Draws added variable plots (also known as partial regression plots)
if(class(f)[1]=="glm") stop("Added variable plots not valid for glm objects")	
	X <- model.matrix(f$terms, model.frame(f)) 
	y<-y<-fitted.values(f)+residuals(f)
	X<-X[,-1]
	k<-dim(X)[2]
		for(i in (1:k)) {
		z1 <- lsfit(X[,  - i], X[, i], intercept = intercept)
		z2 <- lsfit(X[,  - i], y, intercept = intercept)
		plotname <- dimnames(X)[[2]][i]
		plot(z1$residuals, z2$residuals, xlab = plotname, ylab = 
			"Residuals", main = paste("Partial plot of", plotname))
		abline(h=0,lty=1)
		abline(coef(lm(z2$residuals~z1$residuals)),lty=2)
		}
	invisible()
}


#S3 method for class "formula"
added.variable.plots.formula <-
function(f,intercept=TRUE, data, subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...) 
{
# Draws added variable plots (also known as partial regression plots)

#Code from lm
    			cl <- match.call()
    			mf <- match.call(expand.dots = FALSE)
   			m <- match(c("f", "data", "subset", "weights", "na.action", 
    			    "offset"), names(mf), 0L)
    			mf <- mf[c(1L, m)]
    			mf$drop.unused.levels <- TRUE
    			mf[[1L]] <- as.name("model.frame")
    			mf <- eval(mf, parent.frame())
    			if (method == "model.frame") 
        			return(mf)
    			else if (method != "qr") 
        			warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            			method), domain = NA)
    			mt <- attr(mf, "terms")
    			y <- model.response(mf, "numeric")
		if(any(y<=0)) stop("Responses must be positive")
    			w <- as.vector(model.weights(mf))
    			if (!is.null(w) && !is.numeric(w)) 
        			stop("'weights' must be a numeric vector")
    			offset <- as.vector(model.offset(mf))
    			if (!is.null(offset)) {
        			if (length(offset) != NROW(y)) 
            			stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                			length(offset), NROW(y)), domain = NA)
    			}

    			if (is.empty.model(mt)) stop("Model is empty")


	
	X <- model.matrix(mt, mf, contrasts) 
	X<-X[,-1]
	k<-dim(X)[2]
		for(i in (1:k)) {
		z1 <- lsfit(X[,  - i], X[, i], intercept = intercept)
		z2 <- lsfit(X[,  - i], y, intercept = intercept)
		plotname <- dimnames(X)[[2]][i]
		plot(z1$residuals, z2$residuals, xlab = plotname, ylab = 
			"Residuals", main = paste("Partial plot of", plotname))
		abline(h=0,lty=1)
		abline(coef(lm(z2$residuals~z1$residuals)),lty=2)
		}
	invisible()
}


################################################################

allpossregs<-function(f,best=1,Cp.plot=TRUE,text.cex=0.8, dp=3, cv.rep=50, nvmax=20,...)
{
	UseMethod("allpossregs")
}


#S3 method for class "formula" 
allpossregs.formula=function(f, best=1, Cp.plot=TRUE, text.cex=0.8, dp=3, cv.rep=50, nvmax=20, data, subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...)
{
# calculates model goodness statistics for All Possible Regressions

#Code from lm
    			cl <- match.call()
    			mf <- match.call(expand.dots = FALSE)
   			m <- match(c("f", "data", "subset", "weights", "na.action", 
    			    "offset"), names(mf), 0L)
    			mf <- mf[c(1L, m)]
    			mf$drop.unused.levels <- TRUE
    			mf[[1L]] <- as.name("model.frame")
    			mf <- eval(mf, parent.frame())
    			if (method == "model.frame") 
        			return(mf)
    			else if (method != "qr") 
        			warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            			method), domain = NA)
    			mt <- attr(mf, "terms")
    			y <- model.response(mf, "numeric")
		if(any(y<=0)) stop("Responses must be positive")
    			w <- as.vector(model.weights(mf))
    			if (!is.null(w) && !is.numeric(w)) 
        			stop("'weights' must be a numeric vector")
    			offset <- as.vector(model.offset(mf))
    			if (!is.null(offset)) {
        			if (length(offset) != NROW(y)) 
            			stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                			length(offset), NROW(y)), domain = NA)
    			}

    			if (is.empty.model(mt)) stop("Model is empty")


	

X <- model.matrix(mt, mf, contrasts) 
y= model.response(model.frame(f, data))
n<-dim(X)[1]
        selection.stuff<-summary(regsubsets(f, data, nbest=best, nvmax=nvmax))      
        Rsq<-selection.stuff$rsq
        rssp<-selection.stuff$rss
        p<- apply(selection.stuff$which, 1, sum)
        sigma2<-rssp/(n-p)
        adjRsq<- selection.stuff$adjr2
        Cp<-selection.stuff$cp
        AIC<-rssp/sigma2[length(p)] + 2*p
        BIC<-rssp/sigma2[length(p)] + log(n)*p
nmod<-dim(selection.stuff$which)[1]
k<-dim(selection.stuff$which)[2]


pred.error<-numeric(10)
CV.mat<-matrix(0, nmod,cv.rep)
m<-n%/%10

for(k in 1:cv.rep){

# randomise order
rand.order<-order(runif(n))

yr<-y[rand.order]
Xr<-X[rand.order,]

for(j in 1:nmod){
sample<-1:m
for(i in 1:10){
use.cols = selection.stuff$which[j,]
use.cols[1]=FALSE
use.mat<-as.matrix(Xr[-sample,use.cols])
colnames(use.mat) = colnames(X)[use.cols]
test.mat<-as.matrix(Xr[sample,use.cols])
colnames(test.mat) = colnames(use.mat)
use.data<-data.frame(y=yr[-sample], use.mat)
new.data<-data.frame(test.mat)

fit<-lm(y~., data=use.data)
pred.error[i]<-sum((yr[sample]-predict(fit, new.data))^2)
sample<-sample + m
}
CV.mat[j,k]<-mean(pred.error)
}}

CV=apply(CV.mat, 1,mean)
        result<-round(cbind(rssp,sigma2,adjRsq,Cp,AIC,BIC,CV, selection.stuff$which[,-1]),dp)
        plot(p-1,Cp,xlab="Number of variables",ylab="Cp",main="Cp Plot")
        textvec<-character(length(p))
        for(i in 1:length(p))textvec[i]<-paste((1:(max(p)-1))[selection.stuff$which[i,-1]],collapse=",")
        text(p-1,Cp,textvec,pos=3, cex=text.cex)
        lines(p-1,p,lty=2)
result
}


allpossregs.lm=function(f, best=1,Cp.plot=TRUE,text.cex=0.8, dp=3, cv.rep=50, nvmax=20, ...)
{
# calculates model goodness statistics for All Possible Regressions
if(class(f)[1]!="lm") stop("All possible regressions only implemented for lm objects")
X<-model.matrix(f$terms,model.frame(f))
y= fitted.values(f) + residuals(f)
n<-dim(X)[1]
        selection.stuff<-summary(regsubsets(X[,-1], y, nbest=best, nvmax=nvmax))      
        Rsq<-selection.stuff$rsq
        rssp<-selection.stuff$rss
        p<- apply(selection.stuff$which, 1, sum)
        sigma2<-rssp/(n-p)
        adjRsq<- selection.stuff$adjr2
        Cp<-selection.stuff$cp
        AIC<-rssp/sigma2[length(p)] + 2*p
        BIC<-rssp/sigma2[length(p)] + log(n)*p
nmod<-dim(selection.stuff$which)[1]
k<-dim(selection.stuff$which)[2]


pred.error<-numeric(10)
CV.mat<-matrix(0, nmod,cv.rep)
m<-n%/%10

for(k in 1:cv.rep){

# randomise order
rand.order<-order(runif(n))

yr<-y[rand.order]
Xr<-X[rand.order,]

for(j in 1:nmod){
sample<-1:m
for(i in 1:10){
use.cols = selection.stuff$which[j,]
use.cols[1]=FALSE
use.mat<-as.matrix(Xr[-sample,use.cols])
colnames(use.mat) = colnames(X)[use.cols]
test.mat<-as.matrix(Xr[sample,use.cols])
colnames(test.mat) = colnames(use.mat)
use.data<-data.frame(y=yr[-sample], use.mat)
new.data<-data.frame(test.mat)

fit<-lm(y~., data=use.data)
pred.error[i]<-sum((yr[sample]-predict(fit, new.data))^2)
sample<-sample + m
}
CV.mat[j,k]<-mean(pred.error)
}}

CV=apply(CV.mat, 1,mean)
        result<-round(cbind(rssp,sigma2,adjRsq,Cp,AIC,BIC,CV, selection.stuff$which[,-1]),dp)
        plot(p-1,Cp,xlab="Number of variables",ylab="Cp",main="Cp Plot")
        textvec<-character(length(p))
        for(i in 1:length(p))textvec[i]<-paste((1:(max(p)-1))[selection.stuff$which[i,-1]],collapse=",")
        text(p-1,Cp,textvec,pos=3, cex=text.cex)
        lines(p-1,p,lty=2)
result
}

######################################################################

boxcoxplot<- function(f, p = seq(-2, 2, length = 20),...)
{
	UseMethod("boxcoxplot")
}

#S3 method for class 'formula'
boxcoxplot.formula <-  function (f, p = seq(-2, 2, length = 20), data,  subset, weights, na.action, method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...) 

	{
	# Draws a box-Cox plot for transforming to normality

		#Code from lm
    			cl <- match.call()
    			mf <- match.call(expand.dots = FALSE)
   			m <- match(c("f", "data", "subset", "weights", "na.action", 
    			    "offset"), names(mf), 0L)
    			mf <- mf[c(1L, m)]
    			mf$drop.unused.levels <- TRUE
    			mf[[1L]] <- as.name("model.frame")
    			mf <- eval(mf, parent.frame())
    			if (method == "model.frame") 
        			return(mf)
    			else if (method != "qr") 
        			warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            			method), domain = NA)
    			mt <- attr(mf, "terms")
    			y <- model.response(mf, "numeric")
		if(any(y<=0)) stop("Responses must be positive")
    			w <- as.vector(model.weights(mf))
    			if (!is.null(w) && !is.numeric(w)) 
        			stop("'weights' must be a numeric vector")
    			offset <- as.vector(model.offset(mf))
    			if (!is.null(offset)) {
        			if (length(offset) != NROW(y)) 
            			stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                			length(offset), NROW(y)), domain = NA)
    			}

    			if (is.empty.model(mt)) stop("Model is empty")
			

		l <- length(p)
		boxcox <- seq(l)
		n <- length(y)		
		sumlog <- sum(log(y))
                x <- model.matrix(mt, mf, contrasts)
		for (i in seq(l)) 
			{
			y.p<- if(p[i]==0)log(y) else (y^p[i] -1)/p[i]
			z.trans <- if (is.null(w)) 
            			lm.fit(x, y.p, offset = offset, singular.ok = singular.ok, ...)
        			else lm.wfit(x, y.p, w, offset = offset, singular.ok = singular.ok, ...)
			trans.res<-residuals(z.trans)
			ResSS.p<-sum(trans.res^2)
			boxcox[i] <- n * log(ResSS.p)/2 - (p[i] - 1) * sumlog
		}
		plot(p, boxcox, type = "l", ylab = "Profile likelihood",
		main = "Box-Cox plot", ...)
	}


#S3 method for class "lm"
boxcoxplot.lm <- function(f,p = seq(-2, 2, length = 20), ...){
if(class(f)[1]!="lm") stop("must be an lm object")

	l <- length(p)
		boxcox <- seq(l)
		y<-fitted.values(f)+residuals(f)
		X <- model.matrix(f$terms, model.frame(f))
		n <- length(y)
		sumlog <- sum(log(y))
		for (i in seq(l)) 
			{
			y.p<- if(p[i]==0)log(y) else (y^p[i] -1)/p[i]
			trans.res<-residuals(lm(y.p~-1+X))
			ResSS.p<-sum(trans.res^2)
			boxcox[i] <- n * log(ResSS.p)/2 - (p[i] - 1) * sumlog
		}
		plot(p, boxcox, type = "l", ylab = "Profile likelihood",
		main = "Box-Cox plot", ...)
	}

#############################################################

cross.val<-function(f, nfold=10, nrep=20,...)
{
	UseMethod("cross.val")
}



#S3 method for class "formula"
cross.val.formula=function (f, nfold=10, nrep=20, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
if (!(family$family=="binomial"||family$family=="gaussian" )) {
        print(family)
        stop("'family' must be gaussian or binomial")
    }

 
    if (missing(data))
        data <- environment(f)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

if(family$family=="binomial"){
n<-dim(X)[1]
m<-n%/%nfold
Spec<-matrix(0,nrep,nfold)
Sense<-matrix(0,nrep,nfold)
Tot<-matrix(0,nrep,nfold)

 
for(i in 1:nrep){
rand.perm = order(runif(n))
X = X[rand.perm,]
Y = Y[rand.perm]
mysample = 1:m
for(j in 1:nfold){

 

fit = glm.fit(X[-mysample,], Y[-mysample], weights = rep(1, n-length(mysample)),
        start = NULL, etastart = NULL, mustart = NULL,
        offset = rep(0, n-length(mysample)), family = binomial(),
        control = list(), intercept = TRUE)
 

y.pred<-(X[mysample,]%*%fit$coefficients>0)*1
y.act<-Y[mysample]


 
Sense[i,j]<-sum(y.act*y.pred)/sum(y.act)
Spec[i,j]<-sum((1-y.pred)*(1-y.act))/sum(1-y.act)
Tot[i,j]<-(sum((1-y.pred)*(1-y.act))+sum(y.act*y.pred))/m
if(j<nfold) mysample = mysample + m else mysample = max(mysample)+1:n
}
}
cat("Mean Specificity = ",mean(Spec, na.rm=TRUE), "\n")
cat("Mean Sensitivity = ",mean(Sense, na.rm=TRUE), "\n")
cat("Mean Correctly classified = ",mean(Tot, na.rm=TRUE), "\n")
invisible(list(Specificity=Spec, Sensitivity=Sense, Correct=Tot))
} else {



    n <- dim(X)[1]
    CV <- numeric(nfold)
    pred.error <- numeric(nfold)
    m <- n%/%nfold
    for (k in 1:nrep) {
        rand.order <- order(runif(n))
        yr <- Y[rand.order]
        Xr <- X[rand.order, ]
        sample <- 1:m
        for (i in 1:nfold) {
              use.mat <- as.matrix(Xr[-sample,])
              test.mat <- as.matrix(Xr[sample,])
              y.use = yr[-sample]
              new.data <- data.frame(test.mat)
              fit <- lm(y.use ~ -1+use.mat)
              my.predict = test.mat%*%coefficients(fit) 
              pred.error[i] <- sum((yr[sample] - my.predict)^2)/m
              sample <- if(i==nfold) (max(sample)+1):n else sample + m

            }
            CV[k] <- mean(pred.error)
        }
cat("Cross-validated estimate of root \nmean square prediction error = ",
    sqrt(mean(CV)),"\n")
}

}





#S3 method for class "glm"
cross.val.glm<-function(f, nfold=10, nrep=20,...){

# f: a glm object
# data: the data frame containing the data
	
	classvec=class(f)
	if(class(f)[1]!="glm") stop("must be glm object")
	if(f$family$family!="binomial") stop("must be logistic regression") 
n<-length(residuals(f))
m<-n%/%nfold
Spec<-matrix(0,nrep,nfold)
Sense<-matrix(0,nrep,nfold)
Tot<-matrix(0,nrep,nfold)

for(i in 1:nrep){
rand.order<-order(runif(n))
newdata<-f$data[rand.order,]
sample<-1:m
for(j in 1:nfold){
fit<-glm(f$formula, family =binomial, data = newdata[-sample,]) 
y.pred<-(predict(fit, newdata[sample,],type="response")>0.5)*1
y.act<-(f$y[rand.order])[sample]

Sense[i,j]<-sum(y.act*y.pred)/sum(y.act)
Spec[i,j]<-sum((1-y.pred)*(1-y.act))/sum(1-y.act)
Tot[i,j]<-(sum((1-y.pred)*(1-y.act))+sum(y.act*y.pred))/m
sample<-sample + m
}
}
cat("Mean Specificity = ",mean(Spec, na.rm=TRUE), "\n")
cat("Mean Sensitivity = ",mean(Sense, na.rm=TRUE), "\n")
cat("Mean Correctly classified = ",mean(Tot, na.rm=TRUE), "\n")
invisible(list(Specificity=Spec, Sensitivity=Sense, Correct=Tot))
}




# lm object
cross.val.lm <- function (f, nfold = 10, nrep = 20, ...){
    X <- model.matrix(f$terms, model.frame(f))
    y = fitted.values(f) + residuals(f)
    n <- dim(X)[1]
    CV <- numeric(nfold)
    pred.error <- numeric(nfold)
    m <- n%/%nfold
    for (k in 1:nrep) {
        rand.order <- order(runif(n))
        yr <- y[rand.order]
        Xr <- X[rand.order, ]
        sample <- 1:m
        for (i in 1:nfold) {
              use.mat <- as.matrix(Xr[-sample,])
              test.mat <- as.matrix(Xr[sample,])
              y.use = yr[-sample]
              new.data <- data.frame(test.mat)
              fit <- lm(y.use ~ -1+use.mat)
              my.predict = test.mat%*%coefficients(fit) 
              pred.error[i] <- sum((yr[sample] - my.predict)^2)/m
              sample <- if(i==nfold) (max(sample)+1):n else sample + m

            }
            CV[k] <- mean(pred.error)
        }
cat("Cross-validated estimate of root \nmean square prediction error = ",
    sqrt(mean(CV)),"\n")
}
  


###################################################################
err.boot<-function(f,B=500,...)
{
	UseMethod("err.boot")
}



#S3 method for class "formula"
err.boot.formula = function(f,  B=500, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...){ 

call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
if (!(family$family=="binomial"||family$family=="gaussian")) {
        print(family)
        stop("'family' must be gaussian or binomial")
    }

 
    if (missing(data))
        data <- environment(f)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

nobs=length(Y)

  fit = glm.fit(X, Y,weights = rep(1, nobs),
	start = NULL, etastart = NULL, mustart = NULL,
        offset = rep(0, nobs), family = family,
        control = list(), intercept = TRUE)
if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = f, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))

if(family$family=="binomial"){

  Q = (fit$y!=(predict(fit)>0))
# bootstrap sample
  n = dim(data)[1]
  w.boot = numeric(B)
  for(b in 1:B){
    ni = rmultinom(1, n ,prob=rep(1/n,n))
    w.boot[b] = sum(((1/n) - (ni/n))* Q[rep(1:n,ni)])
  }
Err.B = mean(Q) + mean(w.boot)
list(err=mean(Q), Err = Err.B)
}

else {
n <- length(Y)
    w.boot <- numeric(B)
    for (b in 1:B) {
        use <- sample(n,n, replace=TRUE)
	    yr <- Y[use]
        Xr <- X[use,]
        boot.lm <- lm(yr~-1+Xr)
        beta.boot <- coefficients(boot.lm)
        res <- Y - X%*%beta.boot
        res.boot <- residuals(boot.lm)
        w.boot[b] <- mean(res^2) - mean(res.boot^2) 
            }
    res.old <- residuals(lm(Y~X))
    RSS.n <- mean(res.old^2)
    Err.B <- RSS.n + mean(w.boot)
    list(err = RSS.n, Err = Err.B)
}
}





#S3 method for class "glm"
err.boot.glm = function(f, B=500, ...){ 
  classvec=class(f)
  if(class(f)[1]!="glm") stop("must be glm object")
  if(f$family$family!="binomial") stop("must be logistic regression")
  Q = (f$y!=(predict(f)>0))
# bootstrap sample
  n = dim(f$data)[1]
  w.boot = numeric(B)
  for(b in 1:B){
    ni = rmultinom(1, n ,prob=rep(1/n,n))
    w.boot[b] = sum(((1/n) - (ni/n))* Q[rep(1:n,ni)])
  }
Err.B = mean(Q) + mean(w.boot)
list(err=mean(Q), Err = Err.B)
}


# lm object
err.boot.lm <-function (f, B = 500, ...) 
{
    classvec <- class(f)
    if (class(f)[1] != "lm") 
        stop("must be lm object")
    X <- model.matrix(f$terms, model.frame(f))
    y <- fitted.values(f) + residuals(f)
    n <- length(y)
    w.boot <- numeric(B)
    for (b in 1:B) {
        use <- sample(n,n, replace=TRUE)
	    yr <- y[use]
        Xr <- X[use,]
        boot.lm <- lm(yr~-1+Xr)
        beta.boot <- coefficients(boot.lm)
        res <- y - X%*%beta.boot
        res.boot <- residuals(boot.lm)
        w.boot[b] <- mean(res^2) - mean(res.boot^2) 
            }
    res.old <- residuals(lm(y~X))
    RSS.n <- mean(res.old^2)
    Err.B <- RSS.n + mean(w.boot)
    list(err = RSS.n, Err = Err.B)
}

###########################################################

funnel<-function(f,...)
{
	UseMethod("funnel")
}


#S3 method for class "lm"
funnel.lm <-function(f, ...)
{
# diagnostic plots for detecting non-constant variance

        # plots log sds vs log means 
        # and then squared residuals versus fitted values
        # returns  an estimate of the variance function
	
	classvec=class(f)
	if(length(classvec)!=1) stop ("must be lm object")
	if(classvec[1]!="lm") stop ("must be lm object")
	oldpar = par(mfrow=c(1,2))
        pred<-fitted(f)
        res<-residuals(f)
        cut.points<-quantile(pred,c(0.,.2,.4,.6,.8,1.0)) + c(-0.01,0,0,0,0,1.01)
        group<-cut(pred,cut.points)
        log.means<-log(tapply(pred+res,group,mean))
        log.sds<-log(sqrt(tapply(pred+res,group,var)))
        plot(log.means,log.sds,xlab="Log means",ylab="Log std. errors")
        res.sq<-res^2
        plot(pred,res.sq,ylab="Squared residuals",xlab="Fitted values")
        low.stuff<-loess(res.sq~pred,span=1)
        lines(sort(pred),fitted(low.stuff)[order(pred)],lty=2)

        cat("Slope:",lm(log.sds~log.means)$coef[2],"\n")
        par(oldpar)
        invisible(abs(fitted(low.stuff)))
}


#S3 method for class "formula"
funnel.formula <-
function(f, data,  subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...)
{
# diagnostic plots for detecting non-constant variance

        # plots log sds vs log means 
        # and then squared residuals versus fitted values
        # returns  an estimate of the variance function
	

  
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
            3) else numeric(), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w != 
            0) else if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }
			
	
	oldpar = par(mfrow=c(1,2))
	
        pred<- z$fitted.values
        res=  z$residuals
        cut.points<-quantile(pred,c(0.,.2,.4,.6,.8,1.0)) + c(-0.01,0,0,0,0,1.01)
        group<-cut(pred,cut.points)
        log.means<-log(tapply(pred+res,group,mean))
        log.sds<-log(sqrt(tapply(pred+res,group,var)))
        plot(log.means,log.sds,xlab="Log means",ylab="Log std. errors")
        res.sq<-res^2
        plot(pred,res.sq,ylab="Squared residuals",xlab="Fitted values")
        low.stuff<-loess(res.sq~pred,span=1)
        lines(sort(pred),fitted(low.stuff)[order(pred)],lty=2)

        cat("Slope:",lm(log.sds~log.means)$coef[2],"\n")
        par(oldpar)
        invisible(abs(fitted(low.stuff)))
}



###########################################################

HLstat<-function(f,...)
{
	UseMethod("HLstat")
}

#S3 method for class "glm"
HLstat.glm<-function(f,...){
	classvec=class(f)
	if(class(f)[1]!="glm") stop("must be glm object")
	if(f$family$family!="binomial") stop("must be logistic regression")
y<-f$y*1
probs<-predict(f, type="response")
prob.order<-order(probs)
n<-length(prob.order)
nn<-rep(n%/%10,10)
if(n>sum(nn)) {addvec<-1:(n-sum(nn))
			nn[addvec]<-nn[addvec]+1
			}
group<-rep(1:10,nn)
yi<-tapply(y[prob.order], group,sum)
pi<-tapply(probs[prob.order], group,sum)
HL<-sum((yi-pi)^2/(pi*(1-pi/nn)))
cat("Value of HL statistic = ",round(HL,3),"\n")
cat("P-value = ",round(1-pchisq(HL,8),3),"\n")
invisible()
}



#S3 method for class "formula"
HLstat.formula<-function(f,family=binomial,data=data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...){

call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
if (family$family!="binomial") {
        print(family)
        stop("'family' must be binomial")
    }

 
    if (missing(data))
        data <- environment(f)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")


nobs=length(Y)

fit<-glm.fit(X, Y,weights = rep(1, nobs),
	start = NULL, etastart = NULL, mustart = NULL,
        offset = rep(0, nobs), family = binomial(),
        control = list(), intercept = TRUE)


	classvec=class(fit)
	
if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = f, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))

y<-fit$y*1
probs<-predict(fit, type="response")
prob.order<-order(probs)
n<-length(prob.order)
nn<-rep(n%/%10,10)
if(n>sum(nn)) {addvec<-1:(n-sum(nn))
			nn[addvec]<-nn[addvec]+1
			}
group<-rep(1:10,nn)
yi<-tapply(y[prob.order], group,sum)
pi<-tapply(probs[prob.order], group,sum)
HL<-sum((yi-pi)^2/(pi*(1-pi/nn)))
cat("Value of HL statistic = ",round(HL,3),"\n")
cat("P-value = ",round(1-pchisq(HL,8),3),"\n")
invisible()
}



#########################################################

ROC.curve<-function(f,...)
{
	UseMethod("ROC.curve")
}


#S3 method for class "formula"
ROC.curve.formula = function(f, family=gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...){


call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
if (family$family!="binomial") {
        print(family)
        stop("'family' must be binomial")
    }

 
    if (missing(data))
        data <- environment(f)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")


nobs=length(Y)

fit<-glm.fit(X, Y,weights = rep(1, nobs),
	start = NULL, etastart = NULL, mustart = NULL,
        offset = rep(0, nobs), family = binomial(),
        control = list(), intercept = TRUE)


	classvec=class(fit)
	
if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = f, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))

    pred = 1/(1+exp(-predict(fit)))
    p = seq(0,1,length=51)
    ROC = matrix(0, 500,2)
    temp.table = matrix(0, 2,2)
    ROC[1,] = c(1,1)
    for(i in 499:2){
        y.levels  = sort(unique(Y))+1
        pred.levels  = sort(unique((pred > p[i])*1))+1
        temp.table[y.levels,pred.levels] = table( Y,1*( pred > p[i]))
        ROC[i,] = c(temp.table[1,2]/sum(temp.table[1,]), temp.table[2,2]/sum(temp.table[2,]))
    }
    plot(ROC, xlab = "False positive rate", ylab = "True positive rate", type="s")
    temp.table = table(Y,pred >0.5)
    TPR = temp.table[2,2]/sum(temp.table[2,])
    FPR = temp.table[1,2]/sum(temp.table[1,])
    points(FPR, TPR, pch = 19, cex = 1.3, col = "red")
    lines(c(0, 1), c(0, 1), col = "blue", lty = 2)
    x.vals = ROC[50:1,1]
    y.vals = ROC[50:1,2]
    Area = sum(diff(x.vals) * y.vals[-length(y.vals)])
    cat("Area under ROC curve = ", round(Area, 4), "\n")
    text(0.9, 0.1, paste("AUC = ", round(Area, 4)), pos = 2)

}


#S3 method for class "glm"
ROC.curve.glm = function(f, ...){
    classvec = class(f)
    if (class(f)[1] != "glm") 
        stop("must be glm object")
    if (f$family$family != "binomial") 
        stop("must be logistic regression")
    pred = predict(f, type="response")
    p = seq(0,1,length=51)
    ROC = matrix(0, 500,2)
    temp.table = matrix(0, 2,2)
    ROC[1,] = c(1,1)
    for(i in 499:2){
         y.levels  = sort(unique(f$y))+1
         pred.levels  = sort(unique((pred > p[i])*1))+1
         temp.table[y.levels,pred.levels] = table( f$y,1*( pred > p[i]))
         ROC[i,] = c(temp.table[1,2]/sum(temp.table[1,]), temp.table[2,2]/sum(temp.table[2,]))
    }
    plot(ROC, xlab = "False positive rate", ylab = "True positive rate", type="s")
    temp.table = table(f$y,pred >0.5)
    TPR = temp.table[2,2]/sum(temp.table[2,])
    FPR = temp.table[1,2]/sum(temp.table[1,])
    points(FPR, TPR, pch = 19, cex = 1.3, col = "red")
    lines(c(0, 1), c(0, 1), col = "blue", lty = 2)
    x.vals = ROC[50:1,1]
    y.vals = ROC[50:1,2]
    Area = sum(diff(x.vals) * y.vals[-length(y.vals)])
    cat("Area under ROC curve = ", round(Area, 4), "\n")
    text(0.9, 0.1, paste("AUC = ", round(Area, 4)), pos = 2)
}



########################################################

test.lc<-function(f, cc, c, ...)
{
	UseMethod("test.lc")
}

#S3 method for class "lm"
test.lc.lm = function (f, cc, c, ...){

# tests hypothesis cc^Tb=c
# cc: vector of coefficients specifying the linear comb in the hypothesis
# c: hypothetical value of the linear comb
# lm.obj: the lm object storing the regression fit

# returns list consisting of the estimated linear comb, its stadnard error, 
# the residual df, the test statistic and a two-sided p-value


classvec=class(f)
if(class(f)[1]!="lm") stop("must be lm object")
est=sum(cc*coef(f))
cov.beta = summary(f)$cov.unscaled* (summary(f)$sigma^2)
df = f$df.residual
std.err = sqrt(sum(cc*cov.beta%*%cc))
t.stat = (est-c)/std.err
p.val= 2*pt(-abs(t.stat),df)
list(est=est, std.err=std.err, df=df, t.stat = t.stat, p.val=p.val)
}


#S3 method for class "formula"
test.lc.formula = function (f,  cc, c, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...) {

# tests hypothesis cc^Tb=c0
# cc: vector of coefficients specifying the linear comb in the hypothesis
# c0: hypothetical value of the linear comb
# f: the lm object storing the regression fit

# returns list consisting of theestimated linear comb, its stadnard error, 
# the residual df, the test statistic and a two-sided p-value



    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) 
        return(mf)
    if (!is.character(method) && !is.function(method)) 
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit")) 
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

    fit <- eval(call(if (is.function(method)) "method" else method, 
        x = X, y = Y, weights = weights, start = start, etastart = etastart, 
        mustart = mustart, offset = offset, family = family, 
        control = control, intercept = attr(mt, "intercept") > 
            0L))

    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit$null.deviance <- eval(call(if (is.function(method)) "method" else method, 
            x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, 
            offset = offset, family = family, control = control, 
            intercept = TRUE))$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))

est=sum(cc * fit$coefficients)
cov.beta = summary.glm(fit)$cov.scaled
df = fit$df.residual
std.err = sqrt(sum(cc*cov.beta%*%cc))
t.stat = (est-c)/std.err
p.val= 2*pt(-abs(t.stat),df)
list(est=est, std.err=std.err, df=df, t.stat = t.stat, p.val=p.val)
}

#S3 method for class "glm"
test.lc.glm = function (f, cc, c, ...){

# tests hypothesis ccTb=c0
# cc: vector of coefficients specifying the linear comb in the hypothesis
# c0: hypothetical value of the linear comb
# f: the glm object storing the regression fit

# returns list consisting of the estimated linear comb, its stadnard error, 
# the residual df, the test statistic and a two-sided p-value

est=sum(cc*coef(f))
cov.beta = summary(f)$cov.scaled
std.err = sqrt(sum(cc*cov.beta%*%cc))
t.stat = (est-c)/std.err
p.val= 2*pnorm(-abs(t.stat))
list(est=est, std.err=std.err, t.stat = t.stat, p.val=p.val)
}

############################################################

WB.test <- function(f,n.rep=1000,...)
{
UseMethod("WB.test")
}


#S3 method formula
WB.test.formula<-function(f,n.rep=1000, data,  subset, weights, na.action, method = "qr", 
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
    contrasts = NULL, offset, ...) {

# Calculates the Weisberg-Bingham test
#  and its approximate p-value

#Code from lm
    			cl <- match.call()
    			mf <- match.call(expand.dots = FALSE)
   			m <- match(c("f", "data", "subset", "weights", "na.action", 
    			    "offset"), names(mf), 0L)
    			mf <- mf[c(1L, m)]
    			mf$drop.unused.levels <- TRUE
    			mf[[1L]] <- as.name("model.frame")
    			mf <- eval(mf, parent.frame())
    			if (method == "model.frame") 
        			return(mf)
    			else if (method != "qr") 
        			warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            			method), domain = NA)
    			mt <- attr(mf, "terms")
    			y <- model.response(mf, "numeric")
		if(any(y<=0)) stop("Responses must be positive")
    			w <- as.vector(model.weights(mf))
    			if (!is.null(w) && !is.numeric(w)) 
        			stop("'weights' must be a numeric vector")
    			offset <- as.vector(model.offset(mf))
    			if (!is.null(offset)) {
        			if (length(offset) != NROW(y)) 
            			stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                			length(offset), NROW(y)), domain = NA)
    			}

    			if (is.empty.model(mt)) stop("Model is empty") 
    
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }


fit=lm.fit(x,y,offset)
			

sorted.res<-sort(fit$residuals)   #residuals
n<-length(sorted.res)
normal.quantiles<-qnorm(((1:n)-0.375)/(n+0.25))
WB<-cor(normal.quantiles,sorted.res)
WB.vec<-numeric(n.rep)
for(i in 1:n.rep)WB.vec[i]<-cor(normal.quantiles,sort(rnorm(n)))
cat("WB test statistic = ",round(WB,3),"\n")
cat("p = ",round(mean(WB.vec<WB),2),"\n")
}


#S3 method lm
WB.test.lm<-function(f,n.rep=1000, ...){

# Calculates the Weisberg-Bingham test
#  and its approximate p-value
if (class(f)[1] != "lm") 
        stop("must be lm object")
sorted.res<-sort(residuals(f))
n<-length(sorted.res)
normal.quantiles<-qnorm(((1:n)-0.375)/(n+0.25))
WB<-cor(normal.quantiles,sorted.res)
WB.vec<-numeric(n.rep)
for(i in 1:n.rep)WB.vec[i]<-cor(normal.quantiles,sort(rnorm(n)))
cat("WB test statistic = ",round(WB,3),"\n")
cat("p = ",round(mean(WB.vec<WB),2),"\n")
}



#######################################################


influenceplots<-function(f, cex.lab = 0.7,...)
{
	UseMethod("influenceplots")
}

#S3 method for class "glm"
influenceplots.glm<-function(f, cex.lab = 0.7, ...)
{

oldpar= par(no.readonly = TRUE)

par(mfrow=c(2,2))

p<-f$rank
n = 1 + f$df.null
p.res<-residuals(f,type="pearson")
d.res<-residuals(f,type="deviance")

plot(d.res,type="h",xlab="Observation number",
     ylab="Deviance Residuals",
     main="Index plot of deviance residuals")
plt.text = abs(d.res)>2
if(any(plt.text))text((1:n)[plt.text], d.res[plt.text], (1:n)[plt.text], cex=cex.lab)
abline(h=2, lty=2)
abline(h=-2, lty=2)

hats<-hatvalues(f)
plot(hats,type="h",ylim=c(0,max(hats)),
     xlab="Observation Number" ,ylab="Leverage",
     main="Leverage plot")
plt.text = hats>3*p/n
if(any(plt.text))text((1:n)[plt.text], hats[plt.text], (1:n)[plt.text], cex=cex.lab)
abline(h=3*p/n, lty=2)



cooks<-(p.res^2*hats)/(p*(1-hats)^2)
plot(cooks,type="h",ylim=c(0,max(cooks)),
     xlab="Observation number",ylab="Cook's Distance", 
     main="Cook's Distance Plot")
plt.text = cooks>qchisq(0.1,p)/p
if(any(plt.text))text((1:n)[plt.text], cooks[plt.text],(1:n)[plt.text], cex=cex.lab)
abline(h=qchisq(0.1,p)/p, lty=2)

dev.changes<-d.res^2+p.res^2*(hats/(1-hats))
plot(dev.changes,type="h",ylim=c(0,max(dev.changes)),
     xlab="Observation number",ylab="Deviance changes",
     main="Deviance Changes Plot")
plt.text = dev.changes>4
if(any(plt.text))text((1:n)[plt.text], dev.changes[plt.text], (1:n)[plt.text], cex=cex.lab)
abline(h=4, lty=2)

par(oldpar)

invisible()
}





#S3 method for class "lm"
influenceplots.lm<-function(f,cex.lab=0.7, ...){
# draws index plots of influence measures
op <- par(ask = TRUE)
on.exit(par(op))
infl <- lm.influence(f)
p <- f$rank
e <- weighted.residuals(f)
s <- sqrt(sum(e^2, na.rm = TRUE)/df.residual(f))
xxi <- chol2inv(f$qr$qr, f$qr$rank)
si <- infl$sigma
h <- infl$hat
n<-length(h)
dfbetas <- abs(infl$coefficients/outer(infl$sigma, sqrt(diag(xxi))))

vn <- variable.names(f)
vn[vn == "(Intercept)"] <- "1_"
colnames(dfbetas) <- paste("dfb", abbreviate(vn), sep = ".")
dffits <- abs(e * sqrt(h)/(si * (1 - h)))
cov.ratio <- (si/s)^(2 * p)/(1 - h)
CR1<-abs(cov.ratio-1)
cooks.d <- ((e/(s * (1 - h)))^2 * h)/p
no.of.plots<-p+4

maintext<-c(colnames(dfbetas),"DFFITS"," ABS(COV RATIO-1)", "Cook's D", "Hats")
show.label<-matrix(FALSE,n,no.of.plots)
thresh = numeric(no.of.plots)
for(j in 1:p) show.label[,j]<-dfbetas[,j]>1
show.label[,p+1]<-dffits>(3*sqrt(p/(n-p))) # DFFITS
show.label[,p+2]<-CR1>3*p/n  # COV Ratio
show.label[,p+3]<-cooks.d>qf(0.1,p,n-p)  # Cooks D
show.label[,p+4]<-h>3*p/n   # Hats

for(j in 1:p) thresh[j]<-1
thresh[p+1]<-(3*sqrt(p/(n-p))) # DFFITS
thresh[p+2]<-3*p/n  # COV Ratio
thresh[p+3]<-qf(0.1,p,n-p)  # Cooks D
thresh[p+4]<-3*p/n   # Hats



plotmat<-cbind(dfbetas,dffits,CR1,cooks.d,h)

for(j in 1:(p+4)){
        show.i<-(1:n)[show.label[,j]]
        plot(c(0,n*1.05), range(plotmat[,j]), ylim=c(0,max(plotmat[,j])),
        main=maintext[j], xlab = "Obs. number", ylab = maintext[j], type="n")
        points(1:n, plotmat[,j],type = "h")
        abline(h=thresh[j], lty=2)
        if(length(show.i)>0)
              text(show.i+0.5,plotmat[show.i,j] , show.i, cex=cex.lab, adj=0)
}
}





#S3 method for class "formula"
influenceplots.formula<-function(f, cex.lab=0.7, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...)
{
  call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) 
        data <- environment(f)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("f", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) 
        return(mf)
    if (!is.character(method) && !is.function(method)) 
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit")) 
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- eval(call(if (is.function(method)) "method" else method, 
        x = X, y = Y, weights = weights, start = start, etastart = etastart, 
        mustart = mustart, offset = offset, family = family, 
        control = control, intercept = attr(mt, "intercept") > 
            0L))
    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit$null.deviance <- eval(call(if (is.function(method)) "method" else method, 
            x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, 
            offset = offset, family = family, control = control, 
            intercept = TRUE))$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = f, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))
    fit





if (family$family == "gaussian"){
	influenceplots.lm(fit)
} else

if (family$family == "binomial"){
	influenceplots.glm(fit)
} else

if (family$family == "poisson"){
	influenceplots.glm(fit)
}
else stop("family not supported")



}



###########################################################

# function to draw 3-d scatterplot and fitted regression plane
# assumes that first variable in data frame is the response
reg3d = function(data, wire=FALSE){
# wire = TRUE uses a wireframe surface
wire.val = if(wire) "lines" else "fill"
require(rgl)
X = data[,2]
X=(X-mean(X))/sd(X)
Y = data[,1]
Y=(Y-mean(Y))/sd(Y)
Z=data[,3]
Z=(Z-mean(Z))/sd(Z)
rgl.open()
rgl.bg(color="white")
rgl.material(lit = FALSE)

rgl.lines( c(0,5), c(0,0), c(0,0),lwd=2, col="black")
rgl.lines( c(0,0), c(0,5), c(0,0),lwd=2, col="black")
rgl.lines( c(0,0), c(0,0), c(0,5),lwd=2, col="black")
rgl.spheres(X, Y, Z, 
radius=0.15, col="red")

xg = seq(-5,5, length=20)
zg = seq(-5,5, length=20)
new.data = data.frame(X=rep(xg, 20), Z=rep(zg, rep(20,20)))
zp = matrix(predict(lm(Y~X+Z), newdata=new.data),20,20)
rgl.surface(xg,zg,zp, col="lightgrey", alpha=0.9, 
 front=wire.val,back=wire.val, lit=FALSE)
rgl.texts(c(5.5,0,0), text=colnames(data)[2], col="blue")
rgl.texts(c(0,5.5,0), text=colnames(data)[1], col="blue")
rgl.texts(c(0,0,5.5), text=colnames(data)[3], col="blue")
}