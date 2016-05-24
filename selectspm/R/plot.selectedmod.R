plot.selectedmod <-
function(x,...){
   selmod<- x
   aiccs <- sapply(selmod$aics, function(x) x$AICc)

    cualesHPP<- grep("HPP_", names(aiccs))
    cuantosHPP <- length(cualesHPP)
    cualHPP <- cualesHPP[which.min(aiccs[cualesHPP])]

     cualesHPC<- grep("HPC_", names(aiccs))
     cualHPC <- cualesHPC[which.min(aiccs[cualesHPC])]
     
     # for ploting homogeneous Poisson model
     rP<-selmod$HPPs[[cualHPP]]$r
     KP<-selmod$Kas[length(selmod$models)][[1]]
     
      par(mfrow=c(2,2))
      plot(selmod$models[[cualHPC]], main="HPC")
      plot(selmod$HPPs[[cualHPP]], sqrt(./pi)-r~r, xlab="r", ylab="L(r)", legend=F, main="HPP")
      plot(selmod$models[length(selmod$models)-1]$PC, main="PC")
      plot(rP, sqrt(KP/pi)-rP,main="P", xlab="r", ylab="L(r)", type="l")
            abline(h=0, lty=2,col=2)

}
