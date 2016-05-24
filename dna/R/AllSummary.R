setMethod("summary","modules",
function(object){
 cat("Class: modules\n")
 F=slot(object,"module")
 mF=length(levels(F))-as.numeric("0"%in%levels(F))
 if (mF>0){
  for (i in 1:mF)
   cat(sum(F==i),"genes in Module",i,"\n")
 }
 else{
  cat("No modules\n")
 }
})

setMethod("summary","resultsIndTest",
function(object){
 cat("Tests for differential connectivity of individual genes\n\n")
 p.val=slot(object,"p.values")
 cat(sum(p.val<.001),"genes are significant at level 0.001\n")
 cat(sum(p.val<.005),"genes are significant at level 0.005\n")
 cat(sum(p.val<.01),"genes are significant at level 0.01\n")
 cat(sum(p.val<.05),"genes are significant at level 0.05\n")
})

setMethod("summary","resultsModTest",
function(object){
 cat("Tests for differential modular structure in two networks of genes\n\n")
 p.val=slot(object,"p.value")
 N=slot(object,"N")
 mod1=slot(object,"modules1")
 mod2=slot(object,"modules2")
 cat("Network 1:\n")
 summary(mod1)
 cat("\n")
 cat("Network 2:\n")
 summary(mod2)
 cat("\n")
 cat("Test statistic: N=",N,"\n")
 cat("P-value=",p.val,"\n") 
})


