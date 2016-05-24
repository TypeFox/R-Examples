reference <-
function(gen,ref=NULL){
  anyNA = function(x) any(is.na(x))
  if(is.null(ref)){
    maf=function(z){
      Z=z
      z[z==2]="A"; nA=length(which(z=="A"))
      z[z==0]="B"; nB=length(which(z=="B"))
      if(nA==nB){x=Z}
      if(nA>nB){z[z=="A"]=2;z[z=="B"]=0;x=z}
      if(nA<nB){z[z=="B"]=2;z[z=="A"]=0;x=z}
      return(x)}
    W=apply(gen,2,maf);W=(as.numeric(W));W=matrix(W,ncol=ncol(gen))
    colnames(W)=colnames(gen)
    
  }else{
    if(ncol(gen)!=length(ref)) stop("Reference parent and matrix of genotypes display non compatible dimensions")
    if(any(is.na(ref))||length(which(ref==5))>0) {"Reference parent must have no missing values"}
    gen[is.na(gen)]=5    
    CA = which(ref==0) # Changing Alleles
    W = gen-1
    W[,CA] = W[,CA]*-1
    W[W==(-4)]=4
    W = W+1
    W[W==(5)]=NA
  } 
  return(W)}
