##################################################
#This function computes the posterior summaries  #
#for an 'mcmc' class object                      #
##################################################

summary.mcmc <- function(object,...) {
   index <- seq(object$burnin, object$NUpd, by = object$lag)
   Beta <- object$Beta[index]
   Q <- object$Q[index]
   G <- object$G[index]
   Tau_F <- object$tauF[index]
   Tau_N <- object$tauN[index]
   GSDF <- exp(sqrt(Tau_F))
   GSDN <- exp(sqrt(Tau_N))
   
   cred <- object$cred
   if(object$indep){
      suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_F,GSDN, GSDF),1,summary.out,cred=cred))
      PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F))
      rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_F", "GSD(Tau_N)",  "GSD(Tau_F)")
   }
   else{
      Tau_NF <- object$tauNF[index]
      suma <- t(apply(rbind(Beta,Q,G,Tau_N,Tau_NF,Tau_F, GSDN, GSDF),1,summary.out, cred=cred))
      PostCovMat <- cov(cbind(Beta,Q,G,Tau_N,Tau_F, Tau_NF))
      rownames(suma) <- c("Beta", "Q", "G", "Tau_N", "Tau_NF", "Tau_F", "GSD(Tau_N)", "GSD(Tau_F)")
   }

   colnames(suma)[c(4,5,2)] <- c("Mean", "SD", "Median")

   dic <- round(as.numeric(object$DIC),5)
   pD <- round(as.numeric(object$pD),5)
   Dbar <- round(as.numeric(object$Dbar),5)

   AR <- object$AR
   ESS <- matrix(object$ESS,1,3)
   rownames(ESS) <- ""
   colnames(ESS) <- c("Beta", "Q", "G")
   

   ans <- list(summary = suma,
            PostCovMat=PostCovMat, 
            DIC = dic, 
            pD = pD,
            Dbar = Dbar,
            ESS = ESS, 
            AcceptRate = AR,
            indep = object$indep)

   class(ans) <- 'summary.mcmc'
   
   return(ans)
}

