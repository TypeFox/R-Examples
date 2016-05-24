matrix.derv <-
function(Y,X0,order=1,rel.tol=.Machine$double.eps^0.5){
### Extract useful info
 if(length(Y)==1){Y<-c(Y)}                                                                  #create a vector for Y if Y has length of one
 nY<-length(Y);nX<-length(X0)                                                               #number of functions in Y & number of variables in X
 if(!(is.matrix(X0))){X0<-matrix(X0,ncol=length(X0))};ncol.X<-ncol(X0);nrow.X<-nrow(X0)     #matrixise X
 if(!(is.matrix(Y))){Y<-matrix(Y,nrow=length(Y))};ncol.Y<-ncol(Y);nrow.Y<-nrow(Y)           #matrixise Y
 vnames<-NULL;for(i in 1:length(Y)){vnames<-union(vnames,names(formals(Y[i][[1]])))}        #collect variable names
### Type classification
 if(nY==1 & nX==1){T<-0;Type<-"scalar-by-scalar"}
 if(nY==1 & nX>1){T<-1;Type<-"scalar-by-vector (gradient)"}
 if(nY>1 & nX==1){T<-2;Type<-"vector-by-scalar"}
 if(nrow.Y>1 & ncol.X>1){T<-3;Type<-"vector-by-vector (Jacobian)"}
 if(nY==1 & ncol.X>1 & nrow.X>1){T<-4;Type<-"scalar-by-matrix"}
 if(nX==1 & ncol.Y>1 & nrow.Y>1){T<-5;Type<-"matrix-by-scalar"}
### Input coherence check
 if(ncol.X>1 & nrow.X>1 & length(Y)>1){stop("If X0 is a matrix, Y can only be a multivariate function")}   #matrix-type multivariate evaluation points on one multivariate function
 if(ncol.Y>1 & nrow.Y>1 & length(X0)>1){stop("If Y is a matrix, X0 can only be a scalar")}                  #a scalar evaluation point on a matrix-valued function 
 if(length(Y)>1){for(i in 1:(length(Y)-1)){if(length(unique(names(formals(Y[i+1][[1]]))==names(formals(Y[i][[1]]))))>1){stop("Incoherent variables in some of the functions in Y")}}}                     #vaiables across all functions in Y
 if(length(vnames)!=length(X0)){stop("Number of variables in Y does not match number of evaluation points in X0")}  #Number of evaluation points in X and number of variables in Y
### Initialise computation
Y0<-Y
if(ncol.X>1){for(i in 1:(ncol.X-1)){Y0<-cbind(Y0,Y)}}
Y<-Y0
if(nrow.X>1){for(i in 1:(nrow.X-1)){Y<-rbind(Y,Y0)}}
vnames<-matrix(vnames,nrow=nrow.X,ncol=ncol.X)
nr<-nrow(Y);nc<-ncol(Y);D<-result<-matrix(NA,nr,nc)
### computation
for(i in 1:nr){
for(j in 1:nc){
D[i,j]<-body2string(Y[i,j][[1]])
  if(T==0|T==1|T==2|T==3){
    for(jj in 1:nc){
    if(j!=jj){D[i,j]<-object.substitute(D[i,j],vnames[1,jj],X0[1,jj])}
    }
    Df<-NULL
    eval(parse(text=paste("Df<-function(",vnames[1,j],"){",D[i,j],"}",collapse="",sep="")))
    result[i,j]<-derv(Df,X0[1,j],n=order,tol=rel.tol)}
  if(T==4){
    for(ii in 1:nr){
    for(jj in 1:nc){
    if(j!=jj | i!=ii){D[i,j]<-object.substitute(D[i,j],vnames[ii,jj],X0[ii,jj])}
    }}
    eval(parse(text=paste("Df<-function(",vnames[i,j],"){",D[i,j],"}",collapse="",sep="")))
    result[i,j]<-derv(Df,X0[i,j],n=order,tol=rel.tol)}
  if(T==5){
    eval(parse(text=paste("Df<-function(",vnames[1,1],"){",D[i,j],"}",collapse="",sep="")))
    result[i,j]<-derv(Df,X0[1,1],n=order,tol=rel.tol)}
}}
return(list(result=result,type=T))}
