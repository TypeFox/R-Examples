#' @useDynLib mpbart


mpbart_call <- function(formula, data, base = NULL,test.data = NULL, 
                        Prior = NULL, Mcmc = NULL, 
                        varying = NULL, sep = '.')
{

  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  formula <- callF$formula <- mFormula(formula)
  nframe <- length(sys.calls())
  
  response.name <- paste(deparse(attr(formula, "lhs")[[1L]]))
  m <- match(c("formula", "data"), names(callT), 0L)
  mf <- callT
  mf <- mf[c(1L, m)]
  
  mf$formula <- formula
  attr(mf$formula, "rhs") <- attr(formula,"rhs")[2L]
  
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf)
  Y <- as.factor(Y)
lev <- levels(Y)
p <- length(lev)
if (p < 3){
  stop("The number of alternatives should be at least 3.")
}

if (!is.null(base))
{  
  base <- base
} else {
  base <- lev[p]
}

if (base %in% lev) {
  Y <- relevel(Y, ref = base)
  lev <- levels(Y)
} else {
  stop(paste("Error: `base' does not exist in the response variable."))
}

base <- lev[1]
counts <- table(Y)
if (any(counts == 0)) {
  warning(paste("group(s)", paste(lev[counts == 0], collapse = " "), "are empty"))
  Y <- factor(Y, levels  = lev[counts > 0])
  lev <- lev[counts > 0]
}
Y <- as.numeric(unclass(Y)) - 1
Y <- ifelse(Y==0, p,Y)

cat("The base level is: '", lev[1], "'.\n\n", sep="") 

relvelved <- c(lev[2:length(lev)], lev[1])

Terms <- attr(mf, "terms")

X <- model.matrix.default(Terms, mf)


xcolnames <- colnames(X)
xcolnames <- xcolnames[-1]

if(length(xcolnames) == 1 ){
  X <- data.frame(X[,xcolnames])
  names(X) <- xcolnames[1]
  
} else {
  
  X <- X[,xcolnames]
}
if (!is.null(test.data)){
  
  if(length(xcolnames) == 1 ){
    Xtest <- data.frame(test.data[,xcolnames])
    names(Xtest) <- xcolnames[1]
  } else {
    Xtest <- test.data[,xcolnames]  
  }
  
  testXEx = NULL;
    for(i in 1:nrow(Xtest)){
      testXEx = rbind(testXEx, matrix(rep(Xtest[i,], p-1), byrow = TRUE, ncol = ncol(Xtest) ) )
    }
   
  
} else {
  testXEx = 0
}


XEx = NULL;
for(i in 1:nrow(X)){
  XEx = rbind(XEx, matrix(rep(X[i,], p-1), byrow = TRUE, ncol = ncol(X) ) )
}


Data = list(p=p,y=Y,X= XEx)
testData = list(p=p,X= testXEx)


cat("Table of y values",fill=TRUE)
print(table(model.response(mf) ))
            

n=length(Data$y)

pm1=p-1

if (!is.null(test.data)){
  testn <- nrow(testData$X)/(p-1)
} else {
  testn <- 0
}


#reading alternate specific variables
if (!is.null(varying)) {
  
  varying.names <- names(data[,varying])
  
  alt.names <- NULL
  for(vv in 1:(length(varying))){
    alt.names <- c(alt.names, unlist(strsplit(varying.names[vv], sep, fixed = TRUE))[1L])
  }
  
  

 alt.names <- unique(alt.names)


  if(length(alt.names) != length(varying)/p){
    stop("alternative variables names mismatch. Check names of alternative variables.")
  }
  
  reordered_names <- NULL
  for(nn in 1:(length(varying)/p)){
    reordered_names <- c(reordered_names,  paste0(alt.names[nn], sep, relvelved) )
  }
  
  XChSp <- createX(p,na=length(varying)/p ,
                  nd=NULL,Xa= data[,reordered_names],Xd=NULL, 
                  INT = FALSE,DIFF=TRUE,base=p)
  
  Data$X <- cbind(Data$X, XChSp) 
  
  if (!is.null(test.data)){  
  testXChSp = createX(p,na=length(varying)/p ,
                  nd=NULL,Xa= test.data[,reordered_names],Xd=NULL, 
                  INT = FALSE,DIFF=TRUE,base=p)  
  testData$X <- cbind(testData$X, testXChSp)
  }
  
}




if(missing(Prior)) 
{nu=pm1+3; V=nu*diag(pm1);
 ntrees=200; kfac=2.0;pbd=1.0;pb=0.5;beta = 2.0;alpha = 0.95; nc = 100; priorindep = 0; minobsnode = 10;
}
else 
{if(is.null(Prior$nu)) {nu=pm1+3} else {nu=Prior$nu}
 if(is.null(Prior$V)) {V=nu*diag(pm1)} else {V=Prior$V}
 if(is.null(Prior$ntrees)) {ntrees=200} else {ntrees=Prior$ntrees}
 if(is.null(Prior$kfac)) {kfac=2.0} else {kfac=Prior$kfac}
 if(is.null(Prior$pbd)) {pbd=1.0} else {pbd=Prior$pbd}
 if(is.null(Prior$pb)) {pb=0.5} else {pb=Prior$pb}
 if(is.null(Prior$beta)) {beta = 2.0} else {beta=Prior$beta}
 if(is.null(Prior$alpha)) {alpha = 0.95} else {alpha=Prior$alpha}
 if(is.null(Prior$nc)) {nc=100} else {nc=Prior$nc}
 if(is.null(Prior$priorindep)) {priorindep= FALSE} else {priorindep=Prior$priorindep}
 if(is.null(Prior$minobsnode)) {minobsnode= 10} else {minobsnode=Prior$minobsnode}
 
 
}

if(is.null(Mcmc$sigma0)) {sigma0=diag(pm1)} else {sigma0=Mcmc$sigma0}

if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
if(is.null(Mcmc$burn)) {burn=100} else {burn=Mcmc$burn}
if(is.null(Mcmc$ndraws)) {ndraws=1000} else {ndraws=Mcmc$ndraws}
if(is.null(Mcmc$keep_sigma_draws)) {keep_sigma_draws=FALSE} else {keep_sigma_draws=Mcmc$keep_sigma_draws}






C=chol(solve(sigma0))
#
#  C is upper triangular root of sigma^-1 (G) = C'C
#
sigmai=crossprod(C)

if( (priorindep ==TRUE) || (keep_sigma_draws==FALSE)){
  sigmasample = as.double(0);
  savesigma = 0;
} else {
  sigmasample = as.double(rep(sigma0, ndraws+burn));
  savesigma = 1;
}

cat("Number of trees: ", ntrees, "'.\n\n", sep="") 
cat("Number of draws: ", ndraws, "'.\n\n", sep="") 
cat("burn-in: ", burn, "'.\n\n", sep="") 



res =   .C("rmnpMDA",w=as.double(rep(0,nrow(Data$X))),
           trainx= as.double(t(Data$X)), 
           testx= as.double(t(testData$X)),
           mu = as.double(rep(0,nrow(Data$X))),
           sigmai = as.double(sigmai),
           V = as.double(V),
           n = as.integer(length(Data$y)),
           n_dim = as.integer(ncol(sigmai)),
           y = as.integer(Data$y), 
           n_cov = as.integer(ncol(Data$X)), 
           nu = as.integer(nu), 
           trainpred = as.double(rep(0,p*n)) , 
           testn = as.integer(testn), 
           testpred = as.double(rep(0,p*testn)), 
           ndraws = as.integer(ndraws), 
           burn = as.integer(burn),
           ntrees = as.integer(ntrees),
           kfac = as.double(kfac), 
           pbd = as.double(pbd), 
           pb = as.double(pb), 
           alpha = as.double(alpha),  
           beta =  as.double(beta),
           nc = as.integer(nc),
           savesigma = as.integer(savesigma),
           minobsnode = as.integer(minobsnode),
           sigmasample = sigmasample,
           PACKAGE="mpbart")      


class_prob_train <- matrix(res$trainpred,ncol = p, byrow = TRUE)
class_prob_train <- data.frame(class_prob_train)
names(class_prob_train) <- relvelved

predicted_class_train <- apply(class_prob_train,1,function(x) relvelved[which.max(x)] )
predicted_class_train <- as.factor(predicted_class_train)

if (!is.null(test.data)){
  class_prob_test <- matrix(res$testpred,ncol = p, byrow = TRUE)
  class_prob_test <- data.frame(class_prob_test)
  names(class_prob_test) <- relvelved
  predicted_class_test <- apply(class_prob_test,1,function(x) relvelved[which.max(x)] )
  predicted_class_test <- as.factor(predicted_class_test)
  
} else {
  class_prob_test <- NULL
  predicted_class_test <- NULL
}

ret = list(class_prob_train = class_prob_train, 
           predicted_class_train = predicted_class_train,
           class_prob_test = class_prob_test, 
           predicted_class_test = predicted_class_test, 
           sigmasample = res$sigmasample);



class(ret) = "mpbart"

return(ret)
}