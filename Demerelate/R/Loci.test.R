Loci.test <- function(tab.pop, ref.pop="NA", object=FALSE, value="rxy", bt=1000, file.output=FALSE)
  
{
  

  # Function tests on differences in mean relatedness based on number of loci used for calculation
  
  # Data are loaded for input
    if (object==FALSE) {
                          tab.pop <- input.txt(tab.pop, "pop")
                          if (ref.pop=="NA") {ref.pop <- tab.pop} else {ref.pop <- input.txt(ref.pop, "ref.pop")}
                        }
    
    if (object==TRUE) {
                          if (is(ref.pop)[1]!="data.frame") {ref.pop <- tab.pop}
                       }
  
  
  loci <- seq(3,length(tab.pop[1,]),2)
  boots <- vector("list",bt)
  loc <- vector("list", length(loci))
  
  # Calculate each loci in list
  lis.boot <- lapply(loci,function(x){allele.sharing(tab.pop, tab.pop, x, FALSE, value, ref.pop)})
  lis.boot <- lapply(lis.boot,as.dist)
  names(lis.boot) <- colnames(tab.pop)[loci]
  
    
  # Draw loci random for i loci from bt samples
  for (i in 1:length(loc))
  {
  # Draw random bt samples
  rsamp <- lapply(boots, function(x){sample(seq(1:length(as.vector(lis.boot[[1]]))),replace=T)})
  
  # Take loci information from bt samples
  list.bt.loc <- lapply(rsamp,function(rsamp){lapply(vector("list",i), function(x){lis.boot[[sample(1:length(loc),1)]][rsamp]})})
  
  # Mean over loci
  list.bt.loc <- lapply(list.bt.loc,function(x){do.call("cbind",lapply(x,array))})
  list.bt <- lapply(list.bt.loc,function(x){matrix(rowMeans(x,na.rm=TRUE))})
  loc[[i]] <- dist(lapply(list.bt,function(x){mean(x,na.rm=T)})) 
  }

if (file.output==TRUE)
{  
  pdf(paste("Loci.test",Sys.Date(),".pdf"))
  
  y <- unlist(lapply(loc,mean))
  errbar(x=c(1:length(loc)),y,yplus=y+unlist(lapply(loc,sd)),yminus=y-unlist(lapply(loc,sd)),xlab="Number of Random Loci", ylab=paste("Mean Difference in Relatedness [",value,"]"))
  lines(unlist(lapply(loc,mean))~c(1:length(loc)))
  
  dev.off()
}
  
  return(loc)
  
}
