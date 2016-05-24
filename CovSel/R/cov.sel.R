`cov.sel` <-function(T, Y, X, type=c("dr","np"), alg = 3, scope=NULL, alpha = 0.1,thru=0.5,thro=0.25,thrc=100,...){
type<- match.arg(type)

 if(is.null(scope)){
    scope.dat <- NULL
    scope.0 <- NULL
    scope.1 <- NULL
 }else{
    scope.dat <- as.formula(paste("~", paste(paste("dat$",scope), collapse= "+")))  
    scope.0 <- as.formula(paste("~", paste(paste("data.0$",scope), collapse= "+")))
    scope.1 <- as.formula(paste("~", paste(paste("data.1$",scope), collapse= "+")))
 }
 
t <- length(T)
c <- nrow(X)
o <- length(Y)

  ## checking if of same length
 if(t!=c  || t!=o  ||  c!=o) stop("Data is not of same length")
  ## checking if enough df
 if((t - ncol(X)) <= 1 ) stop("Not Enough Observations or To Many Predictors")

dat <- as.data.frame(cbind(X,T,Y))
data.0 <- subset(dat,dat$T == 0 )
data.1 <- subset(dat,dat$T == 1 )

  ## checking missing data  and numericness
d <- sum(as.integer(is.na(cbind(T,Y)) == TRUE))
  if (d != 0) stop("Missing data in response and/or treatment")

namn <- colnames(dat)                          ##vector med namn
covar <- namn[1:(length(namn)-2)]              ##covariates
#Checking if any discrete covariates
if(any(sapply(X, is.factor)) && type=="dr") stop("'np' should be used when discrete covariates are present")
  
  
if(type=="np"){#source("cov.sel.np.R")
										l<-cov.sel.np(T, Y, X, alg, scope, thru, thro, thrc, dat, data.0, data.1, covar, ...)
										class(l) <- "cov.sel"
										invisible(return(l))

										
										}
if(type=="dr"){
				
if (alg == 1){
## Algoritm A
## Step 1
## creating dr object  and backward elimination
f1 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",covar), collapse= "+")))
step1 <- dr(f1, ...)
beA1 <- dr.step(step1, numdir=0, stop=alpha,scope=scope.dat,...)
X.T <- rownames(beA1[[9]])                        ##subset X.T, vector med namn
substring(X.T, 1) <- c("    ")
X.T <- sub('    ', '', X.T)
## Step 2
## creating dr object and BE    for control group
f2 <- as.formula(paste("data.0$Y ~ ", paste(paste("data.0$",X.T), collapse= "+")))
step2 <- dr(f2, ...)
beA2 <- dr.step(step2, numdir=0, stop=alpha,scope=scope.0,...)
Q.0 <- rownames(beA2[[9]])                     ##subset Q.0, vector med namn
substring(Q.0, 1) <- c("       ")
Q.0 <- sub('       ', '', Q.0)

## creating dr object and BE    for treated
f3 <- as.formula(paste("data.1$Y ~ ", paste(paste("data.1$",X.T), collapse= "+")))
step3 <- dr(f3, ...)
b <- dr.step(step3, numdir=0, stop=alpha,scope=scope.1,...)
Q.1 <- rownames(b[[9]])                     ##subset Q.1, vector med namn
substring(Q.1, 1) <- c("       ")
Q.1 <- sub('       ', '', Q.1)

##extract info
eq0 <- beA2[[9]]
rownames(eq0) <- Q.0
eq1 <- b[[9]]
rownames(eq1) <- Q.1
metoden <- beA2[[4]]

## Return Value
l <- list(X.T = X.T, Q.0 = Q.0, Q.1 = Q.1,
  evectorsQ.0 = eq0, evectorsQ.1 = eq1, method=metoden, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))

}else if (alg == 2){
## Algoritm B
## Step 1
  ## creating dr object and backward elimination     for control
f1 <- as.formula(paste("data.0$Y ~ ", paste(paste("data.0$",covar), collapse= "+")))
step1 <- dr(f1, ...)
beB1 <- dr.step(step1, numdir=0, stop=alpha, scope=scope.0,...)
X.0 <- rownames(beB1[[9]])                     ##subset X.0, name vector
substring(X.0, 1) <- c("       ")
X.0 <- sub('       ', '', X.0)

  ## creating dr object and BE    for treated
f2 <- as.formula(paste("data.1$Y ~ ", paste(paste("data.1$",covar), collapse= "+")))
step2 <- dr(f2, ...)
beB11 <- dr.step(step2, numdir=0, stop=alpha, scope=scope.1,...)
X.1 <- rownames(beB11[[9]])                     ##subset X.1, name vector
substring(X.1, 1) <- c("       ")
X.1 <- sub('       ', '', X.1)

## Step 2
  ## creating dr object    for X.0
f3 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",X.0), collapse= "+")))
step3 <- dr(f3, ...)
beB2 <- dr.step(step3, numdir=0, stop=alpha, scope=scope.dat,...)
Z.0 <- rownames(beB2[[9]])                        ##subset Z.0, name vector
substring(Z.0, 1) <- c("    ")
Z.0 <- sub('    ', '', Z.0)

  ## creating dr object   for X.1
f4 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",X.1), collapse= "+")))
step4 <- dr(f4, ...)
obj <- dr.step(step4, numdir=0, stop=alpha, scope=scope.dat,...)
Z.1 <- rownames(obj[[9]])                        ##subset Z.1, name vector
substring(Z.1, 1) <- c("    ")
Z.1 <- sub('    ', '', Z.1)

## extract info
ez0 <- beB2[[9]]
rownames(ez0) <- Z.0
ez1 <- obj[[9]]
rownames(ez1) <- Z.1
metoden <- beB2[[4]]

## Return Values
l <- list(X.0 = X.0, X.1 = X.1, Z.0 = Z.0, Z.1 = Z.1,
 evectorsZ.0 = ez0, evectorsZ.1 = ez1, method=metoden, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))

}else if (alg == 3){
## Algoritm A & B
## Step 1
## creating dr object
f1 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",covar), collapse= "+")))
step1 <- dr(f1, ...)
obj1 <- dr.step(step1, numdir=0, stop=alpha, scope=scope.dat,...)
X.T <- rownames(obj1[[9]])                        ##subset X.T, vector med namn
substring(X.T, 1) <- c("    ")
X.T <- sub('    ', '', X.T)

## Step 2
## creating dr object  for control group
f2 <- as.formula(paste("data.0$Y ~ ", paste(paste("data.0$",X.T), collapse= "+")))
step2 <- dr(f2, ...)
obj2 <- dr.step(step2, numdir=0, stop=alpha, scope=scope.0,...)
Q.0 <- rownames(obj2[[9]])                     ##subset Q.0, vector med namn
substring(Q.0, 1) <- c("       ")
Q.0 <- sub('       ', '', Q.0)

## creating dr object  for treated
f3 <- as.formula(paste("data.1$Y ~ ", paste(paste("data.1$",X.T), collapse= "+")))
step3 <- dr(f3,  ...)
obj3 <- dr.step(step3, numdir=0, stop=alpha, scope=scope.1,...)
Q.1 <- rownames(obj3[[9]])                     ##subset Q.1, vector med namn
substring(Q.1, 1) <- c("       ")
Q.1 <- sub('       ', '', Q.1)


## Step 1 in Algoritm B
  ## creating dr object for control
f4 <- as.formula(paste("data.0$Y ~ ", paste(paste("data.0$",covar), collapse= "+")))
step4 <- dr(f4, ...)
obj4 <- dr.step(step4, numdir=0, stop=alpha, scope=scope.0,...)
X.0 <- rownames(obj4[[9]])                     ##subset X.0, name vector
substring(X.0, 1) <- c("       ")
X.0 <- sub('       ', '', X.0)
                   
  ## creating dr object   for treated
f5 <- as.formula(paste("data.1$Y ~ ", paste(paste("data.1$",covar), collapse= "+")))
step5 <- dr(f5, ...)
obj5 <- dr.step(step5, numdir=0, stop=alpha, scope=scope.1,...)
X.1 <- rownames(obj5[[9]])                     ##subset X.1, name vector
substring(X.1, 1) <- c("       ")
X.1 <- sub('       ', '', X.1)

## Step 2
  ## creating dr object for 0
f6 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",X.0), collapse= "+")))
step6 <- dr(f6, ...)
obj6 <- dr.step(step6, numdir=0, stop=alpha, scope=scope.dat,...)
Z.0 <- rownames(obj6[[9]])                        ##subset Z.0, name vector
substring(Z.0, 1) <- c("    ")
Z.0 <- sub('    ', '', Z.0)

  ## creating dr object for 1
f7 <- as.formula(paste("dat$T ~ ", paste(paste("dat$",X.1), collapse= "+")))
step7 <- dr(f7, ...)
obj7 <- dr.step(step7, numdir=0, stop=alpha, scope=scope.dat,...)
Z.1 <- rownames(obj7[[9]])                        ##subset Z.1, name vector
substring(Z.1, 1) <- c("    ")
Z.1 <- sub('    ', '', Z.1)

##extract info
ez0 <- obj6[[9]]
rownames(ez0) <- Z.0
ez1 <- obj7[[9]]
rownames(ez1) <- Z.1
eq0 <- obj2[[9]]
rownames(eq0) <- Q.0
eq1 <- obj3[[9]]
rownames(eq1) <- Q.1
metoden <- obj3[[4]]

## Return Values
l <- list(X.T = X.T, Q.0 = Q.0, Q.1 = Q.1,
  X.0 = X.0, X.1 = X.1, Z.0 = Z.0, Z.1 = Z.1,
  evectorsQ.0 = eq0, evectorsQ.1 = eq1,
  evectorsZ.0 = ez0, evectorsZ.1 = ez1, method=metoden, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))

}else{
stop("Wrong Selected Algorithm")
}}
}

