
mimwrap <-
function (dataset, estimator="pearson", disc="equalwidth", bins=sqrt(ncol(dataset))) {

   mi1=c("pearson", "spearman","kendall","gauss")
   mi2=c("emp","mm","shrink","sg")
   dataset=t(as.matrix(dataset))
   
   if(estimator %in% mi1){
      mim = cor(dataset, method=estimator, use="complete.obs")^2
      diag(mim)=0
      m=0.999999
      mim[mim>m]=m
      mim= -0.5 * log(1-mim)
      mim[is.na(mim)] = 0
   } else if(estimator %in% mi2){      
      dataset=discretize(dataset, disc, bins)
      mim=mutinformation(dataset, method=estimator)
      mim[is.na(mim)] = 0     
   } else if(estimator=="cor"){
      mim=cor(dataset)
   }
return(mim)
}

