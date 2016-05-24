MAP.genes<-function(resx, value.dis, files=TRUE)
{
 hit <- resx
 probs <- list()
 probX <- list()
 Nprob <- nrow(value.dis) 
 n.entity <- ncol(value.dis) 
 probe.names<-rownames(value.dis)
 hitX <- patternmatrix(rownames(hit),n.entity)
 for ( i in 1:nrow(hitX))
 {
    comp.called <- which(hitX[i,] == 1)
    subX <- value.dis[,comp.called]
    probs[[i]] <- probe.names[which(apply(subX,1,sum)==
    length(comp.called))]
  }
if (files) {
 for (i in 1:nrow(hitX))
 write.table(probs[[i]], paste("probs_",i,".txt", sep=""),
 sep="\t")
}
return(probs)
}