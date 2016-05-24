cpg.everything.logical <-
function(x,perm=FALSE,levin=FALSE,...) {
  fdr<-x
      if(!perm) {
        info<-c("CPG.Labels","T.statistic","P.value","Holm.sig","FDR","F.statistic")
        if(!levin) {nametest<-info[1:5]}
        if(levin) { nametest<-info[c(1,6,3:5)]           
           }
 
      nametest } 
      else {
      permname<-c("Min.P.Value","Number.of.Holm.Significant",
                      "Number.of.FDR.Significant") 
     
      permname }
      }
