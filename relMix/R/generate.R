generate <- function (R, K, x = 1){ 
   #if(length(K)==0) return(.gxMix(R,NULL,x))
   if(x==0) return(NULL)
   RI <- 1:length(sort(unique(R)))
   allgenos = allGenotypes(max(RI))
   allgenos_in_R = which(allgenos[, 1] %in% RI & allgenos[, 2] %in% RI)
   mixgrid = expand.grid(rep(list(allgenos_in_R), x))
   R_not_masked = setdiff(RI, K)
   if (!all(K %in% RI))
      stop("Input error. Known contributors with alleles outside mixture")
   if (length(R_not_masked) > 2*x)
      stop("Input error. Not enough unknown contributors")
   if (length(R_not_masked) > 0)
      mixgrid = mixgrid[apply(mixgrid, 1, function(r) all(R_not_masked %in% allgenos[r, ])),  ,drop=F]
   M <- matrix(apply(mixgrid, 1, function(allg_rows) allgenos[allg_rows, ]),nrow=x)
   dd <- dim(M)
   for (i in 1:dd[1])
     for (j in 1:dd[2])
	   M[i,j] <- R[M[i,j]]
   M
}
.gxMix=function(E,K,x=1){
if(x==1) {
set=.g1Mix(E,K)
dim1=nrow(set)
index=1:(2*dim1)
cc1=seq(1,2*dim1,by=2)
cc2=seq(2,2*dim1,by=2)
index[cc1]=1:dim1
index[cc2]=(dim1+1):(2*dim1)
set=c(set[,1],set[,2])[index]
set=matrix(set,nrow=1)
return(set)
}

U=setdiff(E,K)
lu=length(U)
z=max(lu-2*(x-1),0) #required alleles to occupy in U
if(z>2) return(NULL)
if(z==2) {
S=matrix(apply(expand.grid(U,U),1,sort),ncol=2,byrow=TRUE)
S=S[apply(S,1,function(x) x[1]!=x[2]),]
}
if(z==1) S=matrix(apply(expand.grid(E,U),1,sort),ncol=2,byrow=TRUE)
if(z==0) S=matrix(apply(expand.grid(E,E),1,sort),ncol=2,byrow=TRUE)
S=as.data.frame(unique(S))
ns=nrow(S)
set=vector("list",ns)
set=NULL
for (i in 1:ns){
aa=matrix(nrow=x-1,.gxMix(E,union(K,unlist(S[i,])),x-1))
line=rep(as.integer(S[i,]),dim(aa)[2]/2)
aa=rbind(aa,line)
set=cbind(set,aa)
}
dimnames(set)=NULL
set
}

.g1Mix=function(E,K){
if(length(setdiff(K,E))>0) return(NULL)
U=setdiff(E,K)
lu=length(U)
if(lu>2|lu<0) return(NULL)
if(lu==2) S=matrix(U,ncol=2,byrow=TRUE)
if(lu==1) S=expand.grid(U,E)
if(lu==0) S=matrix(apply(expand.grid(E,E),1,sort),ncol=2,byrow=TRUE)
S=unique(S)
S
}

  
.restAllel <- function(y,R) {
  #Finds the alleles in y also in the evidence R. 
  #Unused alleles are merged into a rest allele.
  #Only works for one marker!
  aa <- as.numeric(attr(y$markerdata[[1]], "alleles"))
  aint <- intersect(aa,R) 
  if(length(aint) < length(aa)){
    aanew <- as.integer(c(aint, max(R)+1))
    afreq <- attr(y$markerdata[[1]], "afreq")
    afreqnew <- c(afreq[aint], sum(afreq[-aint]))
    #names(afreqnew) <- aanew Thore Feb 7 2013
    for(i in 1:y$nInd)
      y <- modifyMarker(y, 1, ids=i, genotype=y$markerdata[[1]][i,], alleles=aanew, afreq=afreqnew)
  }  
  y
}  
 
.convert=function(x){
res=NULL
for (i in 1:(length(x)/2))
res=c(res,paste(sort(x[(2*i-1):(2*i)]),collapse="",sep=""))
res
}
   
   