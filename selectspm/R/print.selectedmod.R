print.selectedmod <-
function(x,...){
      selmod<-x
     # AICcs de cada modelo
     aiccs <- sapply(selmod$aics, function(x) x$AICc)

    cualesHPP<- grep("HPP_", names(aiccs))
    cuantosHPP <- length(cualesHPP)
    cualHPP <- cualesHPP[which.min(aiccs[cualesHPP])]

     cualesHPC<- grep("HPC_", names(aiccs))
     cualHPC <- cualesHPC[which.min(aiccs[cualesHPC])]

     HPP.bw <- selmod$sigmas[cualHPP]
     names(HPP.bw)<- "HPP.bw"
 
    HPC.bw <- selmod$sigmas[cualHPC-cuantosHPP]
    names(HPC.bw)<- "HPC.bw"

    PCsigma2<-selmod$models["PC"][[1]]$sigma2
    names(PCsigma2)="PC.sigma2"
    PCrho<-selmod$models["PC"][[1]]$rho
    names(PCrho)="PC.rho"

    HPCsigma2<-selmod$models[cualHPC][[1]]$sigma2
    names(HPCsigma2)="HPC.sigma2"

   HPCrho<-selmod$models[cualHPC][[1]]$rho
   names(HPCrho)="HPC.rho"

    AICcs <-c(aiccs["P"],aiccs[cualHPP],aiccs["PC"],aiccs[cualHPC])
    names(AICcs)=c("P","HPP","PC","HPC")

  # BEST <- c("P","HPP","PC","HPC")[which.min(AICcs)]
    print(c(AICcs, PCsigma2,PCrho,HPCsigma2,HPCrho, HPP.bw, HPC.bw))
   return(c(AICcs, PCsigma2,PCrho,HPCsigma2,HPCrho, HPP.bw, HPC.bw))

}
