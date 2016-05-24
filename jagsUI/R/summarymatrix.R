
summary.matrix <- function(output,samples,n.chains,codaOnly){
  
  hold <- unlist(output$mean[!names(output$mean)%in%codaOnly])
  toremove <- which(!is.na(hold))
  
  cleanup <- function(input,codaOnly){
    
    out.raw <- unlist(input[!names(input)%in%codaOnly])
    
    out <- out.raw[toremove]
    return(out)
  }
  
  y = data.frame(cleanup(output$mean,codaOnly),cleanup(output$sd,codaOnly),
                 cleanup(output$q2.5,codaOnly),cleanup(output$q25,codaOnly),
                 cleanup(output$q50,codaOnly),cleanup(output$q75,codaOnly),
                 cleanup(output$q97.5,codaOnly),
                 cleanup(output$Rhat,codaOnly),cleanup(output$n.eff,codaOnly),
                 cleanup(output$overlap0,codaOnly),cleanup(output$f,codaOnly))
 
  p <- colnames(samples[[1]])
  expand <- sapply(strsplit(p, "\\["), "[", 1)  
  row.names(y) = p[!expand%in%codaOnly]
  names(y) = c('mean','sd','2.5%','25%','50%','75%','97.5%','Rhat','n.eff','overlap0','f')
  if(n.chains==1){
    y = y[,-c(8,9)]
  }
  
  return(as.matrix(y))
}