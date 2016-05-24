#########################################################################
# ss.mar1 uses the model.object$bestfit to determine the zeros and fit a state-space MAR(1)
#########################################################################
ss.mar1=function(aggregated.data, MAR.obj=NULL, model=list(), control=list(),silent=FALSE){
if(!require(MARSS))
  stop("The MARSS R package is required to run state-space MAR(1) models.\nPlease install and retry.\n",call.=FALSE)

if(!missing(MAR.obj))
  if(class(MAR.obj) != "MAR")
    stop("MAR.obj must be a MAR object as output by run.mar().\n",call.=FALSE)

#Set up the default model
  if(is.null(model$B)){
    if(missing(MAR.obj)){ stop("The B matrix must be specified.\n",call.=FALSE)
    }else{ model$B=MAR.obj$bestfit$B }
  }
  if(is.null(model$Q)) model$Q="unconstrained"
  if(is.null(model$R)) model$R="diagonal and unequal"
  if(is.null(model$C)){
    if(missing(MAR.obj)){ stop("The C matrix must be specified.\n",call.=FALSE)
  }else{ model$C=MAR.obj$bestfit$C }
  }

#B and C must be matrices
if(!is.matrix(model$B)) stop("B must be passed in in the model list as a matrix.\n",call.=FALSE)
if(!is.numeric(model$B)) stop("B must be passed in as a numeric matrix. Locations of 0s should be 0.\n",call.=FALSE)
if(!is.matrix(model$C)) stop("C must be passed in in the model list as a matrix.\n",call.=FALSE)
if(!is.numeric(model$C)) stop("C must be passed in as a numeric matrix. Locations of 0s should be 0.\n",call.=FALSE)

#Determine the increment used for the data
inc=min(diff(aggregated.data$date))
if(inc < 32 & inc > 28){ by="month"; date.format="%Y-%m" }
if(inc > 0.5 & inc < 2){ by="day"; date.format="%Y-%d" } 
if(inc > 5 & inc < 8){ by="week";  date.format="%Y-%W" }
if(inc > 360 & inc < 370){ by="year";  date.format="%Y" }

#create a full list of dates separated by increment above
full.dates = seq(aggregated.data$date[1],max(aggregated.data$date)+1,by=by)
n=length(full.dates) #number of months
n.col=dim(aggregated.data)[2]
#create a data.frame to hold the full data
aggregated.full=data.frame(matrix(NA,n,n.col))
#set the column names of the full matrix
colnames(aggregated.full)=colnames(aggregated.data)
#set the data column to be the full dates
aggregated.full$date=full.dates
#find the locations of the rows in the full data.frame that match
#the months in the data.frame with missing months
match.rows=sapply(format(aggregated.data$date,date.format),
                  match,format(aggregated.full$date,date.format))
#fill in the months with values
aggregated.full[match.rows,]=aggregated.data

#specify which columns of variate data we need
if(is.null(rownames(model$B)))
  stop("The B matrix must have rownames corresponding the the variate names in the data.\n",call.=FALSE)
var.names=rownames(model$B)
n.var=length(var.names)
if(n.var==0) stop("There are no variates.\n",call.=FALSE)
#Subset based on the species and transpose so time goes across columns
marss.var=t(aggregated.full[,var.names])

#specify which columns of covariate data we need
cov.names=colnames(model$C)
n.cov=length(cov.names)
#Subset, transpose so time goes across columns, and stop R from
#dropping the matrix dimensions
marss.cov=t(aggregated.full[,cov.names,drop=FALSE])

#test that data are demeaned
if(!all(apply(marss.var,1,function(x){all.equal(mean(x,na.rm=TRUE),0)})))
   stop("Fittting a state-space MAR(1) model requires that the data be demeaned.\n",call.=FALSE)

B.names=as.list(paste("B(",rep(1:n.var,n.var),",",rep(1:n.var,each=n.var),")",sep=""))
B.names[model$B==0]=0
B=matrix(B.names,n.var,n.var)
C.names=as.list(paste("C(",rep(1:n.var,n.cov),",",rep(1:n.var,each=n.cov),")",sep=""))
C.names[model$C==0]=0
C=matrix(C.names,n.var,n.cov)

R.type=model$R
#set up the R matrix
if(identical(R.type,"diagonal and unequal")){
  R.names=as.list(paste("R(",1:n.var,",",1:n.var,")",sep=""))
  R=matrix(list(0),n.var,n.var)
  diag(R)=R.names
}
if(identical(R.type,"diagonal and equal")){
  R.names=as.list(rep("R",n.var))
  R=matrix(list(0),n.var,n.var)
  diag(R)=R.names
}
if(identical(R.type,"zero")) R=matrix(0,n.var,n.var)
if(identical(R.type,"identity")) R=diag(1,n.var)
if(!MARSS:::is.diagonal(R)) stop("R needs to be a diagonal matrix.\n",call.=FALSE)

Q.type=model$Q
#set up the Q matrix
if(identical(Q.type,"diagonal and unequal")){
  Q.names=as.list(paste("Q(",1:n.var,",",1:n.var,")",sep=""))
  Q=matrix(list(0),n.var,n.var)
  diag(Q)=Q.names
}
if(identical(Q.type,"diagonal and equal")){
  Q.names=as.list(rep("Q",n.var))
  Q=matrix(list(0),n.var,n.var)
  diag(Q)=Q.names
}
if(identical(Q.type,"identity")) Q=diag(1,n.var)
if(identical(Q.type,"unconstrained")){
  Q.names=as.list(paste("Q(",rep(1:n.var,n.var),",",rep(1:n.var,each=n.var),")",sep=""))
  Q=matrix(Q.names,n.var,n.var)
  Q[lower.tri(Q)]=t(Q)[lower.tri(Q)] #make symmetric 
}

#Check that R and Q are now matrices
if(!is.matrix(R)) stop("R must be either a matrix or the text \"identity\", \"diagonal and equal\", or \"diagonal and unequal\".", call.=FALSE)
if(!is.matrix(Q)) stop("Q must be either a matrix or the text \"identity\", \"diagonal and equal\", \"diagonal and unequal\", \"equalvarcov\" or \"unconstrained\".", call.=FALSE)

#Set A and Q
U="zero"

#the the parameters that are not in the model
Z="identity"
A="zero"
D="zero"

if(n.cov>0){ marss.model=list(B=B, U=U, C=C, Q=Q, Z=Z, A=A, D=D, R=R, tinitx=1,c=marss.cov)
}else{ marss.model=list(B=B, U=U, Q=Q, Z=Z, A=A, D=D, R=R, tinitx=1) } #C is zero
marss.dat = marss.var

if(any(is.na(marss.cov))){
  #there are missing values in the covariates so we need to model them
  #must be covs too
    B.cmd=matrix(list(0),n.var+n.cov,n.var+n.cov)
  B.cmd[1:n.var,]=cbind(B,C)
  B.cmd[(n.var+1):(n.var+n.cov),(n.var+1):(n.var+n.cov)]=diag(1,n.cov)

    Q.cmd=matrix(list(0),n.var+n.cov,n.var+n.cov)
    Q.cmd[1:n.var,1:n.var]=Q
Q.cmd[(n.var+1):(n.var+n.cov),(n.var+1):(n.var+n.cov)]=diag(1,n.cov)

R.cmd=matrix(list(0),n.var+n.cov,n.var+n.cov)
R.cmd[1:n.var,1:n.var]=R

  x0.cmd=matrix(list(),n.var+n.cov,1)
x0.names=paste("x0(",1:n.var,")",sep="")
x0.cmd[1:n.var,1]=x0.names
x0.cov=marss.cov[,1]
x0.cov[is.na(x0.cov)]=0 #set missing values to have x0 of 0
x0.cmd[(n.var+1):(n.var+n.cov),1]=x0.cov
  
  marss.model=list(B=B.cmd, U=U, C="zero", Q=Q.cmd, Z=Z, A=A, D=D, R=R.cmd, x0=x0.cmd, tinitx=1)
  marss.dat=rbind(marss.var, marss.cov)
  kem=MARSS(marss.dat, model=marss.model,control=control,silent=TRUE)
  B.est=print(kem,what="B",silent=TRUE)
  B.rtn=B.est[1:n.var,1:n.var,drop=FALSE]
  C.rtn=B.est[1:n.var,(n.var+1):(n.var+n.cov),drop=FALSE]
  Q.rtn=print(kem,what="Q",silent=TRUE)[1:n.var,1:n.var,drop=FALSE]
  R.rtn=print(kem,what="R",silent=TRUE)[1:n.var,1:n.var,drop=FALSE]
}else{ 
    kem=MARSS(marss.dat, model=marss.model,control=control,silent=silent)
    B.rtn=print(kem,what="B",silent=TRUE)
    C.rtn=print(kem,what="C",silent=TRUE)
    Q.rtn=print(kem,what="Q",silent=TRUE)
    R.rtn=print(kem,what="R",silent=TRUE)
    }
A.rtn=matrix(0,n.var,1)
rownames(B.rtn)=var.names; colnames(B.rtn)=var.names
rownames(Q.rtn)=var.names; colnames(Q.rtn)=var.names
rownames(C.rtn)=var.names; colnames(C.rtn)=cov.names
rownames(R.rtn)=var.names; colnames(R.rtn)=var.names
rownames(A.rtn)=var.names; colnames(A.rtn)="a"


rtn.list=list(ssfit=kem, A=A.rtn,B=B.rtn,C=C.rtn,process.errors=list(covariance=Q.rtn),observation.errors=list(covariance=R.rtn),AIC=kem$AIC,AICc=kem$AICc,log.likelihood=kem$logLik)

return(rtn.list)
}






