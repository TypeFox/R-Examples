##################################################
#This function computes the posterior summaries  #
#for an 'sir' class object                       #
##################################################

summary.sir <- function(object,...) {
   Beta <- object$Beta
   Q <- object$Q
   G <- object$G
   Tau_F <- object$tauF
   Tau_N <- object$tauN
   GSDF <- exp(sqrt(Tau_F))
   GSDN <- exp(sqrt(Tau_N))

   cred <- object$cred
   if(object$indep){
      suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_F, GSDN, GSDF),1,summary.out,cred=cred))
      PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F))
      rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "GSD(Tau_N)",  "GSD(Tau_F)")
   }
   else{
      Tau_NF <- object$tauNF
      suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_NF,Tau_F, GSDN, GSDF),1,summary.out, cred=cred))
      PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F, Tau_NF))
      rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_NF", "Tau_F", "GSD(Tau_N)", "GSD(Tau_F)")
   }

   colnames(suma)[c(4,5,2)] <- c("Mean", "SD", "Median")

   dic <- round(as.numeric(object$DIC),5)
   prop <- object$prop
   ESS <- object$ESS
   maxw <- object$maxw
   pD <- round(as.numeric(object$pD),5)
   Dbar <- round(as.numeric(object$Dbar),5)

   ans <- list(summary = suma,
            PostCovMat=PostCovMat, 
            DIC = dic,
            pD = pD,
            Dbar = Dbar,
            maxw = maxw, 
            ESS = ESS, 
            prop = prop,
            indep = object$indep)

   class(ans) <- 'summary.sir'
   return(ans)
}

