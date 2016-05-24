StVal <-
function(fileDir,m=20,N.s=150,V=67749)
{
  #	added following lines to avoid package notes
	#PC=get("PC");
	PC=NULL;
	#files=get("files");
	#X=get("X");
	#N.s=get("N.s");
	#m=get("m");
	#n=get("n");
	
  #	library(fastICA)
  #	m=20
  #	N.s=150
  alpha=0.5
  #V=67749
  
  #	fileDir = "~aeloyan/DataMat/fICA"
  #	fileDir ="/home/bst/other/aeloyan/DataMat/fICA/"
  #	fileDir = "/home/bst/student/shachen/DataMat"
  files = dir(fileDir, pattern = "*", full.names = TRUE)
  X.full=matrix(0,V,m*N.s)
  
  for(i in 1:N.s){
    load(files[i])
    X.full[,(m*(i-1)+1):(m*i)]=PC
  }
  
  
  XtX=0
  ind=500
  mind=floor(V/ind)
  for(i in 1:mind){			
    XtX=XtX+t(X.full[((i-1)*ind+1):(i*ind),])%*%X.full[((i-1)*ind+1):(i*ind),]}
  if((mind*ind)<V){
    XtX=XtX+t(X.full[(mind*ind+1):V,])%*%X.full[(mind*ind+1):V,]}
  
  sv=svd(XtX)
  Sigma=diag(sqrt(sv$d))
  SigmaInv=diag(1/sqrt(sv$d))
  U=sv$u
  V=X.full%*%U%*%SigmaInv
  ## X.full==V%*%Sigma%*%t(U)
  ## a rank m approximation will be given by V[,1:m]%*%Sigma[1:m,1:m]%*%t(U)[1:m,]
  
  V.l=V[,1:m]
  Sigma.l=Sigma[1:m,1:m]
  U.l=t(U)[1:m,]
  X.app=V.l%*%Sigma.l%*%U.l
  
  ### The starting values of the mixing matrices
  
  W0=c()
  for(i in 1:N.s){
    W0=c(W0,list(t(solve(Sigma[1:m,1:m]%*%U.l[,((i-1)*m+1):(i*m)]))))
  }
  save(W0,file="W0.rda")
  
  for(i in 1:N.s){
    X=t(X.app[,((i-1)*m+1):(i*m)])
    inum=paste(paste(rep(0,nchar(as.character(N.s))-nchar(as.character(i))),collapse=""),i,sep="")
    save(X,file=paste(paste("app",inum,sep=""),".rda",sep=""))
  }
  
  
}
