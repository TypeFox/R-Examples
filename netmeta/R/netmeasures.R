netmeasures <- function(x){
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  if (x$comb.random==FALSE & length(x$tau.preset)==0)                           
    nmak <- nma.krahn(x)
  ##                                                                            
  if (length(x$tau.preset)==1)                                                  
    nmak <- nma.krahn(x, tau.preset=x$tau.preset)                                   
  ##                                                                            
  if (x$comb.random==TRUE & length(x$tau.preset)==0)                            
    nmak <- nma.krahn(x, tau.preset=x$tau)
  
  
  comparisons <- nmak$direct[,"comparison"]
  direct  <- nmak$direct
  network <- nmak$network
  ##
  H <- nmak$H
  H.studies <- nmak$H.studies
  
  
  ##
  ## Direct evidence proportion (Koenig 2013, subsection 3.4.1 and 3.4.2)
  ##
  proportion <- rep(0, nrow(H))
  names(proportion) <- rownames(H)
  proportion[comparisons] <- network[comparisons,"seTE"]^2 / direct$seTE^2
  
  
  ##
  ## Minimal parallelism (Koenig 2013, subsection 3.4.3)
  ## Mean path length (Koenig 2013, subsection 3.4.4)
  ##
  l <- lapply(split(H,
                    matrix(colnames(H), nrow=nrow(H),
                           ncol=length(colnames(H)), byrow=TRUE)),
              function(x) matrix(x, nrow=nrow(H))
              )
  ##
  H.tilde <- sapply(l, function(x){ apply(x, 1, function(x) 0.5*(abs(sum(x))+sum(abs(x))) )})
  ##
  rownames(H.tilde) <- rownames(H)
  ##  
  minpar  <- apply(H.tilde, 1, function(x) 1/max(abs(x)))
  meanpath <- apply(H.tilde, 1, sum)
  
  
  ##
  ## Minimal parallelism (study-level)
  ##
  l.studies <- lapply(split(H.studies,
                            matrix(colnames(H.studies), nrow=nrow(H.studies),
                                   ncol=length(colnames(H.studies)), byrow=TRUE)),
                      function(x) matrix(x,nrow=nrow(H.studies))
                      )
  ##
  H.tilde.studies <- sapply(l.studies,
                            function(x){ apply(x,1,function(x) 0.5*(abs(sum(x))+sum(abs(x))))}
                            )
  ##
  rownames(H.tilde.studies) <- rownames(H)
  ##  
  minpar.study <- apply(H.tilde.studies, 1, function(x) 1/max(abs(x)))
  
  
  res <- list(proportion=round(proportion, 4),
              meanpath=round(meanpath, 4),
              minpar=round(minpar, 4),
              minpar.study=round(minpar.study, 4),
              H.tilde=H.tilde)
  ##
  res
}
