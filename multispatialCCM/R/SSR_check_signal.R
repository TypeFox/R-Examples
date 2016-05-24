SSR_check_signal <-
function(A, E, tau=1,
         predsteplist=1:10, matchSugi=0) {
  repvec=1; B<-A; # Only a single vector: A = B
  
  #Predict future elements of A using historical observations of A
  #using cross-validation.
  #matchSugi=1 removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
  #repvec=0 if A and B are not the same vector (i.e, not using "leave one out cross validation")
  # Make list of acceptable starting positions
  predatout<-data.frame(predstep=predsteplist, rho=NA)
  npredstep<-1
  for(predstep in predsteplist) {
    pastdist<-max(c(tau*(E-1), E+1))
    futdist<-predstep
    acceptablelib<-as.numeric(is.finite(A))
    lA<-length(A)
    
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
      
      predatout[npredstep,2]<-out$rho
    }
    npredstep<-npredstep+1
  }
  cor.test(out$Aest, out$A, conf.level = 0.95)
  return(list(predatout=predatout,
    rho_pre_slope=summary(
      lm(predatout[,2]~predatout[,1]))$coefficients[2,c(1,4)],
    rho_predmaxCI=cor.test(
      out$Aest, out$A, conf.level = 0.95)$conf.int[1:2]))
}
