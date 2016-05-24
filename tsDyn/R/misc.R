#.onAttach <- function(...) { 
#	cat("Version 0.7 of package tsDyn will have probably a few minor revisions rapidly.\n User should #update the package often.\n Minor revisions will be announced only on the tsdyn mailing list: #tsdyn@googlegroups.com\n")
#}

extend <- function(...)
	UseMethod("extend")

extend.list <- function(this, subclass, ..., listV) {
	class(this) <- c(subclass, class(this))
	if(missing(listV))
		return(structure(c(this, list(...)), class=class(this)))
	else
		return(structure(c(this, listV), class=class(this)))
}

MyEnv <- function(...) {
	this <- new.env()
	vars <- list(...)
	nms <- names(vars)
	if(length(nms)) for(i in 1:length(nms))
		assign(nms[i], vars[[i]], envir=this)
	structure(this, class=c("MyEnv",class(this)))
}

extend.MyEnv <- function(this, subclass, ...) {
	class(this) <- c(subclass, class(this))
	newvars <- list(...)
	if(length(names(newvars))) for(nm in names(newvars))
		assign(nm, newvars[[nm]], envir=this)
	return(this)
}



#'Available models
#'
#'Available built-in time series models
#'
#'Return the list of built-in available \sQuote{nlar} time series models
#'
#'@return A character vector containing built-in time series models. For help
#'on a specific model, type: \code{help(modelName)}.
#'@author Antonio, Fabio Di Narzo
#'@keywords ts
#'@export
#'@examples
#'
#'availableModels()
#'
availableModels <- function()
	fitters

latex <- function(obj, ...)
	UseMethod("latex")

formatSignedNum <- function(x, ...) {
  signChar <- c("-","+")[(x>=0)+1]
  nm <- names(x)
  res <- paste(signChar,format(abs(x), ...) )
  names(res) <- nm
  return(res)
}

build <- function(...)
	UseMethod("build")

add <- function(...)
	UseMethod("add")



#'sigmoid functions
#'
#'Some sigmoid functions. See R sources for their definition
#'
#'
#'@aliases sigmoid dsigmoid d2sigmoid
#'@param x numeric vector
#'@author J. L. Aznarte
#'@keywords ts
#'@export
sigmoid <- function(x) 1/(1 + exp(-x))

#' @export
#' @rdname sigmoid
dsigmoid <- function(x) x * (1 - x)

repmat <- function(a, b, c) kronecker(matrix(1,b,c), a)

norma <- function(v) sqrt(sum(v^2))

###Function to create the threshold in TVAR
TVAR_thresh<-function(mTh,thDelay,thVar=NULL,y, p){
T <- nrow(y) 
k <- ncol(y) 
if (!is.null(thVar)) {		
        if (length(thVar) > T) {
		z <- thVar[seq_len(T)]
		warning("The external threshold variable is not of same length as the original variable")
        }
        else
		z <- thVar
	z<-embed(z,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
} ###Combination (or single value indicating position) of contemporaneous variables
else {
	if (!length(mTh)%in%c(1,k))
		stop("length of 'mTh' should be equal to the number of variables, or just one")
	if(!all(mTh%in%seq_len(k)))
		stop("Unable to select the variable ",mTh[which(mTh%in%seq_len(k)==FALSE)], " for the threshold. Please see again mTh ")
	if(length(mTh)==1) {
		combin <- matrix(0,ncol=1, nrow=k)
		combin[mTh,]<-1
	}
	else 
		combin<-matrix(mTh,ncol=1, nrow=k)
	zcombin <- y %*% combin
	z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]		#if thDelay=2, ncol(z)=2
}
zcombin <- y %*% combin
z <- embed(zcombin,p+1)[,seq_len(max(thDelay))+1]	
trans<-as.matrix(z)
list(trans=trans, combin=combin)
}

#Function to obtain the number of digits of a value
getndp <- function(x, tol=2*.Machine$double.eps)
{
internal<-function(x){
  	ndp <- 0
  	while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1 
  	if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly underestimated.")
  	ndp 
}
if(!is.null(dim(x)))
	x<-c(x)
if(length(x)==1)
	return(internal(x))
else if(length(x) %in% c(2:20))
	return(max(apply(matrix(x),1,internal)))
else
	return(max(apply(matrix(sample(x,size=20)),1,internal)))
}

### Parameter matrix of tar with 1 threshold,  given the thresh and the delay
TAR1t_B<-function(thDelay,gamma,yy, xxl,xxh,z,m) {#
        isL <- ifelse(z[, thDelay + 1]<= gamma,1,0)	### isL: dummy 
	ndown<-mean(isL)
        xxthresh <- cbind(xxl * isL,xxh * (1 - isL))	### Lower matrix
	B<-round(matrix(solve(crossprod(xxthresh))%*%crossprod(xxthresh,yy), nrow=1),5)

	Bcolnames <- c("Trend", c(paste("t -", seq_len(m))))
	colnames(B)<-rep(Bcolnames,2)
	Bdown <- B[,seq_len(ncol(B)/2)]
	Bup <- B[,-seq_len(ncol(B)/2)]
	nobs <- c(ndown=ndown, nup=1-ndown)
	
	list(Bdown=Bdown, Bup=Bup, nobs=nobs)
}

### Parameter matrix of tar with 2 thresholds,  given the 2 thresh and the delay
TAR2t_B <- function(gam1,gam2,thDelay, yy, xx,z,m){
	##Threshold dummies
	dummydown <- ifelse(z[, thDelay + 1]<=gam1, 1, 0)
	regimedown <- dummydown*xx
	ndown <- mean(dummydown)
	dummyup <- ifelse(z[, thDelay + 1]>gam2, 1, 0)
	regimeup <- dummyup*xx
	nup <- mean(dummyup)
	##SSR from TAR(3)
	XX <- cbind(regimedown, (1-dummydown-dummyup)*xx, regimeup)		# dim k(p+1) x t	
	B <- round(matrix(solve(crossprod(XX))%*%crossprod(XX,yy),nrow=1),5)	#SSR
	Bcolnames <- c("Trend", c(paste("t -", seq_len(m))))
 	colnames(B)<-rep(Bcolnames,3)
	npar<-ncol(B)/3
	Bdown <- B[,c(1:npar)]
	Bmiddle <- B[,c(1:npar)+npar]
	Bup <- B[,c(1:npar)+2*npar]
	nobs <- c(ndown=ndown, nmiddle=1-ndown-nup,nup=nup)	
	list(Bdown=Bdown, Bmiddle=Bmiddle, Bup=Bup, nobs=nobs)
}


###Check if the AR coefficients in each regime lie inside the unit circle
is.InUnitCircle<-function(B,ninc,m, nthresh){

  a<-ninc
  B<-as.vector(B)
  
  para<-B[a+seq_len(m)]
  charPol<-c(1, -para)
  
  val1<-polyroot(charPol)
  root<-Mod(val1)
  
  for(i in 1:(nthresh+1))
  #extract coef
  #build poly
  #extract roots
  #name it
  
  if(nthresh==1|nthresh==2){
	  para2<-B[a+seq_len(m)+m+a]
	  charPol2<-c(1, -para2)
	  val2<-polyroot(charPol2)
	  root<-matrix(c(root,Mod(val2)), nrow=1)
	  colnames(root)<-c(paste("reg", rep(1,m)),paste("reg", rep(2,m)))
	  if(nthresh==2){
		  para3<-B[a+seq_len(m)+2*(m+a)]
		  charPol3<-c(1, -para3)
		  val3<-polyroot(charPol3)
		  root<-matrix(c(root, Mod(val3)), nrow=1)
		  colnames(root)<-c(paste("reg", rep(1,m)),paste("reg", rep(2,m)),paste("reg", rep(3,m)))
	  }
  }
  
  if(any(root<=1))
	  list(warn=TRUE, root=root)
  else
	  list(warn=FALSE, root=root)
}
  
isRoot<-function(coef, regime=c("L", "M", "H", "."), lags){
  regime<-match.arg(regime)
  coefName<-paste("phi", regime, sep="")
  phi<-coef[grep(coefName,names(coef))]
  vec<-rep(0, max(lags))
  vec[lags]<-phi
  charPol<-c(1, -vec)
  val<-polyroot(charPol)
  root<-Mod(val)
  if(any(root<=1)){
    regimeName<-switch(regime, "L"="low", "M"="medium", "H"="high", "."="")
    message<-paste("Possible unit root in the", regimeName, " regime. Roots are:", paste(round(root,4), collapse=" "))
    warn<-TRUE
    warning(message, call.=FALSE)
  }
   else{
    message<-NA
    warn<-FALSE
   }
   return(list(warn=FALSE, root=root, message=message))
}
  
  
percent<-function(x,digits=3,by100=FALSE){
	  a<-ifelse(by100,100,1)
	  paste(a*round(x,digits),"%",sep="")
}

myformat<-function(x,digits, toLatex=FALSE){
	r<-x
	littlex<-abs(x)<10^-(digits)
	r[!littlex]<-formatC(x[!littlex],digits=digits, format="f")
	r[littlex]<-format(x[littlex],digits=min(digits,2), scientific=TRUE)
	if(toLatex)
		r<-gsub("(e.*)","slashtext{\\1}",r)
	if(class(x)=="numeric")
		return(noquote(r))
	if(class(x)=="matrix")
		return(matrix(noquote(r), ncol=ncol(x), nrow=nrow(x)))
}

###Assign to class list is class matrix
asListIfMat<-function(x){
  if(class(x)=="matrix")
      return(list(x))
    else
      return(x)
}


myInsertCol<-function (m, c, v = NA) {
#m: matrix
#c: columns numbers
#v: value to add (scalar only)
    nr <- nrow(m)
    nc <- ncol(m)
    #first
    m2 <- if (1 %in% c) cbind(matrix(v, nrow = nr), m) else m
    #inter
    for(i in c[!c%in%c(1, length(c)+nc)])
        m2 <- cbind(m2[, 1:(i - 1)], matrix(v, nrow = nr), m2[,i:ncol(m2)])
    #last
    if (eval(ncol(m2) + 1) %in% c) 
      m2 <- cbind(m2, matrix(v, nrow = nr))

    return(m2)
}

if(FALSE){
X<-freeny.x
myInsertCol(X, 1, 2)
myInsertCol(X, c(1,2,4,5), 1)
myInsertCol(X, c(1,2,3,5,9), 1)
}
