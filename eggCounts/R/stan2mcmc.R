stan2mcmc <- function(stanFit){
  modelName <- stanFit@model_name
  switch(modelName,
         "Zero-inflated Bayesian model for paired design"=,
         "zipaired"={
           meanEPG.untreated<-rowMeans(extract(stanFit,"mub")[[1]])*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-rowMeans(extract(stanFit,"mub")[[1]])*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,extract(stanFit,c("kappa","mu","phi","delta"),permuted=FALSE)[,1,])
         },
         "Zero-inflated Bayesian model for unpaired design"=,
         "ziunpaired"={
           meanEPG.untreated<-rowMeans(extract(stanFit,"mub")[[1]])*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-rowMeans(extract(stanFit,"mua")[[1]])*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,extract(stanFit,c("kappa","mu","phi","delta"),permuted=FALSE)[,1,])
         },
         "Bayesian model without zero-inflation for paired design"=,
         "paired"=={
           meanEPG.untreated<-rowMeans(extract(stanFit,"mub")[[1]])
           meanEPG.treated<-rowMeans(extract(stanFit,"mub")[[1]])*extract(stanFit,"delta")$delta
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,extract(stanFit,c("kappa","mu","delta"),permuted=FALSE)[,1,])
         },
         "Bayesian model without zero-inflation for unpaired design"=,
         "unpaired"={
           meanEPG.untreated<-rowMeans(extract(stanFit,"mub")[[1]])
           meanEPG.treated<-rowMeans(extract(stanFit,"mua")[[1]])*extract(stanFit,"delta")$delta
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,extract(stanFit,c("kappa","mu","delta"),permuted=FALSE)[,1,])
         },
         "Bayesian model without zero-inflation"=,
         "nb"={
           meanEPG<-rowMeans(extract(stanFit,"mui")[[1]])
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa)
         },
         "Zero-inflated Bayesian model"=,
         "zinb"={
           meanEPG<-rowMeans(extract(stanFit,"mui")[[1]])*extract(stanFit,"phi")[[1]]
           phi<-extract(stanFit,"phi")$phi
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa,phi=phi)
         }
        )
  return(invisible(mcmc(output,start=stanFit@sim$warmup+1,thin=stanFit@sim$thin)))
}