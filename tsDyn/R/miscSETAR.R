###Build the xx matrix with 1 thresh and common=TRUE
buildXth1Common <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,ML]*isL,xx[,MH]*(1-isL))
}

buildXth1LagsIncCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const,xx[,1]*isL,xx[,1]*(1-isL), xx[,-1])
}

buildXth1LagsCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
  LH<-cbind(const*isL,xx[,1]*isL,xx[,1]*(1-isL),const*isL, xx[,-1])
}

###Build the xx matrix with 1 thresh and common=FALSE
buildXth1NoCommon <- function(gam1, thDelay, xx,trans, ML, MH,const, trim, temp=FALSE) {
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
  isL <- ifelse(trans[, thDelay + 1]<= gam1,1,0)	### isL: dummy variable
	xxL <- cbind(const,xx[,ML])*isL
	xxH <- cbind(const,xx[,MH])*(1-isL)
	xxLH<-cbind(xxL,xxH)
	nisL<-mean(isL)
	if(min(nisL, 1-nisL)<trim){
		cat("\n 1 T: Trim not respected: ", c(nisL, 1-nisL), "from th:", gam1)
	}
	return(xxLH)
}





###Build the xx matrix with 2 thresh and common=TRUE
buildXth2Common<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
	nup <- mean(dummyup)
  
	##Construction of the matrix
	xxLMH<-cbind(const,xx[,ML]*dummydown,xx[,MM]*(1-dummydown-dummyup),xx[,MH]*(dummyup))

  ##return result
	if(min(nup, ndown, 1-nup-ndown)>=trim){
		res <- xxLMH	#SSR
	}	else{
	  res <- NA 
	}

	return(res)
}

buildXth2LagsIncCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
      ML<-seq_len(ML)
      MH<-seq_len(MH)
    }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
  
	##Construction of the matrix
	xxLMH<-cbind(const,xx[,1]*dummydown,xx[,1]*(1-dummydown-dummyup),xx[,1]*(dummyup), xx[,-1])
	return(xxLMH)
}


buildXth2LagsCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
	trans<-as.matrix(trans)

	##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
	
  ##Construction of the matrix
	xxLMH<-cbind(const*dummydown,xx[,1]*dummydown,const*(1-dummydown-dummyup), xx[,1]*(1-dummydown-dummyup),const*dummyup, xx[,1]*(dummyup), xx[,-1])
	return(xxLMH)
}


###Build the xx matrix with 2 thresh and common=FALSE
buildXth2NoCommon<-function(gam1,gam2,thDelay,xx,trans, ML, MH,MM, const,trim, temp=FALSE){
  if(temp){
    ML<-seq_len(ML)
    MH<-seq_len(MH)
  }
	
  ##Threshold dummies
	dummydown <- ifelse(trans[, thDelay + 1]<=gam1, 1, 0)
	ndown <- mean(dummydown)
	dummyup <- ifelse(trans[, thDelay + 1]>gam2, 1, 0)
	nup <- mean(dummyup)
  
	##Construction of the matrix
	xxL <- cbind(const,xx[,ML])*dummydown
	xxM<-cbind(const, xx[,MM])*(1-dummydown-dummyup)
	xxH <- cbind(const,xx[,MH])*(dummyup)
	xxLMH<-cbind(xxL,xxM,xxH)

	##return result
	if(min(nup, ndown, 1-nup-ndown)>=trim-0.01){
		res <- xxLMH	#SSR
	}	else {
		cat("\nTrim not respected: ", c( ndown,1-nup-ndown, nup), "from", c(gam1, gam2))
		res <- xxLMH
	}
	return(res)
}



### Wrapper to just get  the SSR :

getSSR <- function(X,y){
  base::crossprod(lm.fit(X,y)$residuals) ## specify base if Matrix loaded
}

### Functions to assemble and return the SSR :
SSR_1thresh<- function(gam1,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH,const=const,trim, fun=buildXth1Common){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH,const, trim)
	if(any(is.na(XX))){
		res<-NA
	}	else {
		res <- getSSR(XX,yy)
	}
	return(res)
}

AIC.matrices<-function(X,y, T, k=2){
	SSR <- getSSR(X,y)
	res<-T*log(SSR/T)+k*(ncol(X)+2)
	return(res)
}




SSR_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim){
	XX<-buildXth2Common(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA
	}	else {
		res <- getSSR(XX,yy)
	}
	return(res)
}

# SSR_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=z, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2Common,trim=trim){
#   SSR_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,fun=buildXth2NoCommon,trim=trim)
# }

SSR_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim){
  	XX<-buildXth2NoCommon(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA
	}	else{
		res <- getSSR(XX, yy)
	}
	return(res)
}

AIC_1thresh<-function(gam1,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH,const=const,trim=trim,fun=buildXth1Common , k=2, T, temp=FALSE){
	XX<-fun(gam1,thDelay, xx,trans=trans, ML=ML, MH=MH, const, trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA
	}	else{
		SSR <- getSSR(XX,yy)
		res<-T*log(SSR/T)+k*(ncol(XX)+1)
	}
	return(res)
}

AIC_2threshCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2, T, temp=FALSE){
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA
	}	else {
		SSR <- getSSR(XX,yy)
		res<-T*log(SSR/T)+k*(ncol(XX)+2)
	}
	return(res)
}

AIC_2th <-function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2, T, temp=FALSE){
	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, temp=temp)
	if(any(is.na(XX))){
		res<-NA
	}	else {
		SSR <- getSSR(XX,yy)
		res<-T*log(SSR/T)+k*(ncol(XX)+2)
	}
	return(res)
}

SSR_2th<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon){
  	XX<-fun(gam1,gam2,thDelay,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim)
	if(any(is.na(XX))){
		res<-NA
	}	else {
		res <- getSSR(XX,yy)
	}
	return(res)
}



AIC_2threshNoCommon<- function(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2Common, k=2,T, temp=FALSE){
  AIC_2threshCommon(gam1,gam2,thDelay, yy=yy,xx=xx,trans=trans, ML=ML, MH=MH, MM=MM,const=const,trim=trim, fun=buildXth2NoCommon, k=k,T, temp=temp)
}


buildConstants<-function(include=c("const", "trend","none", "both"), n){
  include<-match.arg(include)
  if(include=="const"){
	  const <- rep(1,n)
	  incNames<-"const"
  } else if(include=="trend"){
	  const<-seq_len(n)
	  incNames<-"trend"
  }  else if(include=="both"){
	  const<-cbind(rep(1,n),seq_len(n))
	  incNames<-c("const","trend")
  } else {
	  const<-NULL
	  incNames<-NULL
  }
  ninc<-length(incNames)
  res<-list(const=const, incNames=incNames, ninc=ninc)
  return(res)
}





#'Specification of the threshold search
#'
#'This optional function allows the user to set different restrictions for the
#'threshold grid search in function \code{\link{selectSETAR}}.
#'
#'This function is just to check the inputs for the specification of the grid
#'search. If not provided, the search will be in the biggest interval
#'(\verb{ngrid ="All"}) between the minimum and maximum values. The user can
#'reduce it by giving setting "Half" (only every two points is taken) and so
#'on, or setting a number.
#'
#'The search can also be made around a point, or between two points. When
#'between a point, the argument \code{ngrid} is still used, whereas for around,
#'a value of 30 is taken as default value if \code{ngrid} is not specified by
#'user.
#'
#'@aliases MakeThSpec makeThSpec
#'@param exact The user give an exact threshold value
#'@param int The user gives an interval to search inside
#'@param around The user gives an point to search around
#'@param ngrid The number of values to search for
#'@param ... currently unused
#'@return The input values are given as output after checking for consistency
#'(only one of exact/int/around should be given).
#'@author Matthieu Stigler
#'@seealso \code{\link{selectSETAR}}
#'@export
#'@examples
#'
#'sun<-(sqrt(sunspot.year+1)-1)*2		
#'selectSETAR(sun, m=3, th=MakeThSpec(exact=10.40967),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=1)
#'#when pre-sepcified value does not correspond, function will search nearest value
#'selectSETAR(sun, m=3, th=MakeThSpec(exact=10.4),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=1)
#'#search around:
#'selectSETAR(sun, m=3, th=MakeThSpec(around=10.40967, ngrid=20),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=1)
#'#search in an interval
#'selectSETAR(sun, m=3, th=MakeThSpec(int=c(10, 11), ngrid=20),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=1)
#'#reduce size of the grid:
#'selectSETAR(sun, m=3, th=MakeThSpec(ngrid="Half"),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=1)
#'
#'
#'# 2 thresholds:
#'selectSETAR(sun, m=3, th=MakeThSpec(ngrid="Half"),criterion="SSR", d=1, thDelay=0:2,
#'            plot=FALSE, nthresh=2)
#'
#'
MakeThSpec<-function(ngrid=c("All", "Half", "Third", "Quarter"), exact=NULL, int=c("from","to"), around="val",...){
  if(is.character(ngrid))
    ngrid<-match.arg(ngrid)
  #check if only one solution exact/int/around is choosen
  exCheck<-ifelse(is.null(exact),0,1)
  inCheck<-ifelse(is.character(int),0,1)
  aroundCheck<-ifelse(around=="val",0,1)
  if(sum(exCheck, inCheck, aroundCheck)>1)
    stop("Only one of the makeThSpec args should be specified")
  
  #if around, user should give value for ngrid
  if(aroundCheck==1&&is.character(ngrid)){
    ngrid<-20
    cat("When setting arg MakeThSpec(around), a numeric ngrid should be given. Set to 20\n")
  }
  #if exact, ngrid is useless
  if(exCheck==1&&ngrid!="All"){
    ngrid<-20
    cat("When setting arg MakeThSpec(exact), it is useless to specify arg ngrid \n")
  }  
  return(list(exact=exact, int=int, around=around, ngrid=ngrid))
} 

makeTh<-function(allTh, trim, th=list(exact=NULL, int=c("from","to"), around="val",ngrid=c("All", "Half", "Third", "Quarter")), thSteps = 7, trace=FALSE, nthresh=1){
  ng <- length(allTh)
  down<-ceiling(trim*ng)
  up<-floor(ng*(1-trim))
  allin<-up-down
  ngrid<-th$ngrid

#threshold pre-specified
if(!is.null(th$exact)){
  if(length(th$exact)==1){
    if(nthresh==2)
      stop("Please either provide two pre-specified threshold values or set nthresh to 1")
    th<-unique(allTh[which.min(abs(allTh-th$exact))])
    if(length(th)>1){
      th<-th[1]
      cat("Many values correspond to the threshold you gave. The first one",th, "was taken")
    }
    ngrid<-1
  }
  else if(length(th$exact)==2){
    th1<-unique(allTh[which.min(abs(allTh-th$exact[1]))])
    th2<-unique(allTh[which.min(abs(allTh-th$exact[2]))])
    if(length(c(th1, th2))>2){
      th1<-th1[1]
      th2<-th2[2]
      cat("Many values correspond to the threshold you gave. The first ones",c(th1, th2), "were taken")
    }
  }
  else
    warning("Too many threshold values given")
}	  
#interval to search inside given by user
else if(is.numeric(th$int)){
	intDown<-which.min(abs(allTh-th$int[1]))
	intUp<-which.min(abs(allTh-th$int[2]))
	if(length(intDown)>1|length(intUp)>1){
	  intDown<-intDown[1]
	  intUp<-intUp[1]
	}
	#specify how many in the int
	nInt<-intUp-intDown
	if(is.character(ngrid))
	  lengthInt<-nInt*switch(ngrid, "All"=1, "Half"=1/2, "Third"=1/3, "Quarter"=1/4)
	else
	  lengthInt<-min(ngrid,nInt)
	if(trace)
		cat("Searching within",lengthInt, "values between",allTh[intDown], "and", allTh[intUp],"\n")
	th<-allTh[seq(from=intDown, to=intUp, length.out=lengthInt)]
	}
#value to search around	given by user
else if(is.numeric(th$around)){
	if(trace)
		cat("Searching within", ngrid, "values around", th$around,"\n")
	th<-aroundGrid(th$around,allvalues=allTh,ngrid=ngrid,trim=trim, trace=trace) #fun stored in TVECM.R
}

#Default method: grid from lower to higher point
else{
	if(is.character(ngrid))
	  ninGrid<-allin*switch(ngrid, "All"=1, "Half"=1/2, "Third"=1/3, "Quarter"=1/4)
	else 
	  ninGrid<-min(allin, ngrid)
	ths<-allTh[seq(from=down, to=up, length.out=ninGrid)]
	th <-unique(ths)
	if(trace)
		cat("Searching on",length(th), "possible threshold values within regimes with sufficient (",percent(trim*100,2),") number of observations\n")
}
# th<-round(th, getndp(x)) bad idea, rather use format in print and summary
return(th)
}

if(FALSE){

MakeThSpec(ngrid=20)
MakeThSpec(exact=30)
MakeThSpec(exact=30, ngrid=20)
MakeThSpec(int=c(30, 50), ngrid=40)
MakeThSpec(around=30)


environment(makeTh)<-environment(star)
x<-unique(embed(lynx, 2)[,2, drop=FALSE])
a<-makeTh(x, trim=0.15, th=list(ngrid="All"))
a
aHalf<-makeTh(x, trim=0.15, th=list(ngrid="Half"))
aNum<-makeTh(x, trim=0.15, th=list(ngrid=20))
length(a)
length(aHalf)
length(aNum)
length(unique(a))
length(lynx)
length(unique(lynx))
}


