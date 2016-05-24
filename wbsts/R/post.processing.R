post.processing <-
function(z,br,del=-1,epp=-1,C_i=NULL,scales=NULL){ 
  n = nz = dim(z)[1] 
  if (is.null(br)) {
    return(NA)
    stop
  }
  if(del<0){     del<-floor(n)^(2/3)
  }
  br=sort(br, decreasing=F)
  if((br[1]-1) < del) br=br[-1]
  if (is.null(br)) {
    return(NA)
    stop
  }
  if (sum(epp)<0) {
  epp<-c()
  for (j in 1:length(scales)) {
    ep = round(max(2*n/2^scales[j], ceiling(sqrt(n)/2)))
    epp = c(epp,ep)
  }
  epp=round(epp/2) 
  }
  TT=1
  sbr<-fbr<-br
  B<-L<-length(br)
  pp<-NULL
  l=dim(z)[2]
  pp<-matrix(fbr, B, 1)
  criterion<- log(n)  
  Cbr<-c(0, fbr, nz)
  temp<-rep(1, B)
  
  while(TT>0){
    cbr<-c(0, fbr, nz)
    
    for(i in 1:B) {
      b<-cbr[i+1] 
      if (cbr[i]==0) {
        s=1
      } else {
        ind=i
        while (temp[ind-1]==0) {
          ind=ind-1
          if (ind==1) break
        }
        
        s<-cbr[ind]+1
        
      }
      
      if (cbr[i+2]==nz) {
        e=nz
      } else {
        ind=i
        while(temp[ind+1]==0) {
          ind=ind+1
          if (ind==B) break
        }
        e<-cbr[ind+2]
      }
      cr.ip=0
      for (j in 1:l) {
        dis<-c(round((n-n/2^(scales[j])+1)):n)
        e.epp=e-epp[j]
        s.epp=s+epp[j]
        if (s.epp>b | b>e.epp) break
        v = abs(sqrt((e.epp-b)/(e.epp-s.epp+1)/(b-s.epp+1))*sum(z[s.epp:b,j])-sqrt((b-s.epp+1)/(e.epp-s.epp+1)/(e.epp-b))*sum(z[(b+1):e.epp,j]))
        v<-v/mean(z[s.epp:e.epp,j])
        cstar=max((b-s.epp)/(e.epp-s.epp),(e.epp-b)/(e.epp-s.epp))
        cr.ip=abs(cr.ip)+ifelse(v>criterion*C_i[j],v,0)
      }
      temp[sbr==b]<-cr.ip
    }
    pp<-cbind(pp, temp)
    TT=TT-1    
  }
  list(pp=pp)
  cn = ncol(pp)
  if(sum(pp[,cn])==0) {
    return(NA)
  } else return((pp[,1])[pp[,cn]>0])
}
