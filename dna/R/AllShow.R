setMethod("show","modules",
function(object){
 F=slot(object,"module")
 mF=length(levels(F))-(0%in%levels(F))
 if (is.null(names(F)))
  names(F)=paste("Gene ",1:length(F),sep="")
 else{
  NAnames=which(is.na(names(F)))
  names(F)[NAnames]=paste("Gene ",NAnames,sep="")
 }
 if (mF>0){
  for (i in 1:mF){
   cat(paste("Module ",i,":",sep=""),"\n")
   cat(names(F[F==i]),sep=",",fill=80)
   cat("\n")
  }
 }
 else{
  cat("No modules\n")
 }
})

setMethod("show","pairOfNetworks",
function(object){
 cat("Class: pairOfNetworks\n")
 N1=slot(object,"network1")
 N2=slot(object,"network2")
 cat("Network 1:",nrow(N1),"subjects and",ncol(N1),"genes.\n")
 cat("Network 2:",nrow(N2),"subjects and",ncol(N2),"genes.\n")
 if ((!is.null(colnames(N1)))&(!is.null(colnames(N2)))){
  n.common=length(intersect(colnames(N1),colnames(N2)))
  if ((n.common<ncol(N1))|(n.common<ncol(N2)))
   if (n.common!=1){
    cat("The networks have",length(intersect(colnames(N1),colnames(N2))),"genes in common.\n")
   }
   else{
    cat("The networks have 1 gene in common.\n")  
   }
 }
})

setMethod("show","resultsClassTest",
function(object){
 cat("Tests for differential connectivity of a class of genes\n\n")
 p.val=slot(object,"p.value")
 delta=slot(object,"delta")
 genelist=slot(object,"class.genes")
 cat("Class of genes:\n")
 cat(genelist,sep=",",fill=80)
 cat("\n")
 cat("Test statistic: delta=",delta,"\n")
 cat("P-value=",p.val,"\n")
})

setMethod("show","resultsIndTest",
function(object){
 cat("Tests for differential connectivity of individual genes\n")
 cat("(Up to 20 most significant genes)\n\n")
 p.val=slot(object,"p.values")
 test.stat=slot(object,"d")
 gene.names=names(p.val)
 p=length(p.val)
 out=data.frame(gene=gene.names,d=test.stat,p.value=p.val)
 print(out[order(p.val)[1:min(p,20)],],row.names=FALSE)
})

setMethod("show","resultsModTest",
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

