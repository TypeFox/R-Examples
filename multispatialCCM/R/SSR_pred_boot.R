SSR_pred_boot <-
function(A, B=A, E, tau=1, predstep=1, matchSugi=0) {
  repvec=as.numeric((sum(A[is.finite(A)]==B[is.finite(B)])==length(A[is.finite(A)]))&(length(A[is.finite(A)])==length(B[is.finite(B)])))
  
  #Predict elements of A using B
  #If A=B, uses cross-validation
  #matchSugi=1 removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
  #repvec=0 if A and B are not the same vector (i.e, not using "leave one out cross validation")
  # Make list of acceptable starting positions
  gapdist<-tau*(E-1)+predstep
  acceptablelib<-as.numeric(is.finite(A))
  lA<-length(A)
  for(i in 1:gapdist) {
    acceptablelib<-acceptablelib*as.numeric(is.finite(c(rep(NA, i),A[-c((lA-i+1):lA)])))
  }
  acceptablelib<-which(acceptablelib>0)-1 #Convert into positions in C array
  lengthacceptablelib<-length(acceptablelib)
  
  if(tau*(E+1)+predstep>=lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E, "tau =", tau, "and prestep =", predstep))
    return(out=list(A=A, Aest=NA, B=B, E=E, tau=tau, pBlength=length(B), pAlength=length(A), predstep=predstep,
                    prepvec=repvec, pmatchSugi=matchSugi, acceptablelib=acceptablelib, plengthacceptablelib=lengthacceptablelib, rho=NA))
  } else { # Don't attempt to run algorithm using more lags than datapoints
    A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"  
    
    out<-.C("SSR_predict_boot", A=as.double(A), Aest=as.double(rep(0, length(A))), B=as.double(B), E=as.integer(E),
            tau=as.integer(tau),pBlength=as.integer(length(B)), pAlength=as.integer(length(A)), predstep=as.integer(predstep),
            prepvec=as.integer(repvec), pmatchSugi=as.integer(matchSugi), acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
    #out$Aest is full of estimates for n predicted points ahead - fix this
    out$Aest[(1+predstep):length(out$Aest)]<-out$Aest[1:(length(out$Aest)-predstep)];
    out$Aest[out$Aest==0]<-NA
    
    out$rho<-cor(out$A[!is.na(out$Aest)], out$Aest[!is.na(out$Aest)])
    out$Aest[out$Aest==0]<-NA
    return(out)
  }
}
