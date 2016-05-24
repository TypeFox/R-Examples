#--< Function to check even/odd values >--
"is.odd"<-function(x) {
  as.logical(x%%2!=0)
}

#--< Function to rounds X.XXX5 values as X.XX1 >--
"roundup"<-function(x,digits=0){
  round(x+sign(x)*10^(-(digits+5)),digits)
}
#---< End of function >---

#--< Function returns the number of decimals corresponding to a given number of significant figures >--
"getdigits"<-function(x,sigf=6) {
  x[!is.na(x) & x>(-100) & x<100]<-signif(x[!is.na(x) & x>(-100) & x<100],sigf)   # for [-99...0...99] apply number of significant figures
  x[!is.na(x) & (x<=(-100) | x>=100)]<-roundup(x[!is.na(x) & (x<=(-100) | x>=100)]) # for [-Inf...-100,100...+Inf] around with no decimals
  om<-floor(log10(abs(x)))           # calculate order of magnitude for each vec[ ]
  om[which(om==-Inf | om==Inf)]<-0   # for vec[ ]=0, replace all results of log10(0)=-Inf with zero
  dp<-sigf-om-1                      # calculate number of decimal places
  dp[which(dp<0 | is.na(dp))]<-0     # use zero when dp is negative
  return(dp)
}
#---< End of function >---

#--< Function returns the difference between vector and its mean >--
#--  Clone of "subtractMeans" function in S-Plus
"subtractMeans"<-function(x){
  diff<-scale(x,center=T,scale=F)
  return(diff)
}
#---< End of function >---

#--< Function to set a character strin to the 'sentence' case >--
"tosent"<-function(x) {
  y<-x
  substring(y,1,1)<-toupper(substring(x,1,1))
  substring(y,2)<-tolower(substring(x,2))
  return(y)
}
#---< End of function >---

#--< Function to reverse a character >--
"revstr"<-function(x) {
  y<-sapply(strsplit(x,split=""), function(str) {paste(rev(str),collapse="")})
  return(y)
}
#---< End of function >---

#--< Function to convert list to matrix >--
"list2mat"<-function(x) {
  if(!is.null(x)){
    y<-t(sapply(x,'[',1:max(sapply(x,length))))
    if(length(x)>1 & nrow(y)==1) y<-t(y)
    return(y)
  }else{
    return(NULL)
  }
}
#---< End of function >---

#--< Function to convert symmetric matrix -> triangular matrics -> vector >--
"simmat2vect"<-function(x,pref="",...){
  if(!is.null(x)){
    if(is.vector(x)) x<-matrix(x)
    matnm<-colnames(x)
    ni<-nrow(x)
    vecnm<-NULL
    for(i in 1:ni){
      for(j in 1:ni){
        if(j>i) x[i,j]<-NA
      }
      vecnm<-c(vecnm,paste(pref,"(",matnm[i],",",matnm[1:i],")",sep=""))
    }
    x<-as.vector(t(x))
    x<-x[!is.na(x)]
    names(x)<-vecnm
    return(x)
  }else{
    return(NULL)
  }
}
#---< End of function >---

#--< Function to shift up/down a vector or matrix column >--
"vshift"<-function(x,k=1) {
  i<-is.vector(x)
  if(i) x<-as.matrix(x) else x<-as.matrix(x,nrow(x))
  if(nrow(x)>1){
    if(k>0) {
      x<-rbind(matrix(rep(NA,k*ncol(x)),ncol=ncol(x)),matrix(x[1:(nrow(x)-k),],ncol=ncol(x)))
    }
    else {
      x<-rbind(matrix(x[(-k+1):(nrow(x)),],ncol=ncol(x)),matrix(rep(NA,-k*ncol(x)),ncol=ncol(x)))
    }
  }else{
    if(k!=0) x<-matrix(rep(NA,length(x)),nrow=1) 
  }
  if(i) x<-x[1:length(x)]
  return(x)
}
#---< End of function >---
