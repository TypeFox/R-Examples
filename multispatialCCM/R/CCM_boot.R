CCM_boot <-
function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), iterations=100) {
  #Note - input must have "NA" between any gaps in the library (e.g., if you are "stacking" plots)
  ## Initialize parameters
  out_tmp<-NULL
  plengtht=length(A[!is.na(A)]); pLibLength=length(A[!is.na(A)]); Aest=rep(0, length(A)); rho=Aest; varrho=Aest
  
  # Check to make sure library is not too long
  if(plengtht>pLibLength) {plengtht=pLibLength}
  
  # Make list of acceptable starting positions
  gapdist<-tau*(E-1)
  acceptablelib<-as.numeric(is.finite(A))
  lA<-length(A)
  for(i in 1:gapdist) {
    acceptablelib<-acceptablelib*as.numeric(is.finite(c(rep(NA, i),A[-c((lA-i+1):lA)])))
  }
  acceptablelib<-which(acceptablelib>0)-1 #Convert into positions in C array
  acceptablelib<-acceptablelib[acceptablelib<((plengtht-1)-(tau))]
  lengthacceptablelib<-length(acceptablelib)
  
  DesiredL<-DesiredL+E-2 #Convert to positions in C array
  
  for(i in 1:length(DesiredL)) { #Load nearby points from acceptablelib vector
    DesiredL[i]<-acceptablelib[which.min(abs(acceptablelib-DesiredL[i]))] 
  }
  DesiredL<-unique(DesiredL)
  # Update input to match actual indices
  A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"
  
  if(tau*(E+1)>lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E, "and tau =", tau))
    return(out=list(A=A, Aest=NA, B=B, rho=NA, varrho=NA, sdevrho=NA, Lobs=NA, E=out$E, tau=tau, FULLinfo=NA, rejectedL=NA))
  } else {
    
    lpos=NULL
    for(itlst in 1:iterations) {
    out<-.C("CCM_bootstrap", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
            plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)),
            acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
    out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
    
    out_tmp[[itlst]]<-list(Aest=out$Aest, rho=out$rho[out$rho!=0], Lobs=(1:length(A))[out$rho!=0]-E+1)
    lpos=sort(unique(c(lpos, out_tmp[[itlst]]$Lobs))) #Find all library lengths used
    }
    
    #Extract values from list
    Aest_mat<-matrix(nrow=length(A), ncol=iterations)
    rho_mat<-matrix(nrow=length(lpos), ncol=iterations)
    rhosq_mat<-matrix(nrow=length(lpos), ncol=iterations)
    for(itlst in 1:iterations) {
      Aest_mat[,itlst]<-out_tmp[[itlst]]$Aest
      rho_mat[,itlst]<-out_tmp[[itlst]]$rho[match(lpos, out_tmp[[itlst]]$Lobs)]
    }
    
    
    return(list(A=out$A, Aest=rowMeans(Aest_mat, na.rm=T), B=out$B, rho=rowMeans(rho_mat, na.rm=T), sdevrho=apply(rho_mat, 1, function(x) sd(x, na.rm=T)), Lobs=lpos, E=out$E, tau=out$tau, FULLinfo=rho_mat))
  }
}
