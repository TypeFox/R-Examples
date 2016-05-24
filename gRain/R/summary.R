summary.grain <- function(object, type='std', ...){

  type <- match.arg(type,c("std","cliques","rip","configurations"))

  cat("Independence network: Compiled:", object$isCompiled, "Propagated:", object$isPropagated, "\n")

  if (length(object$evidence))
    getEvidence(object)

  cat(" Nodes :")
  utils::str(nodeNames(object)) ## $universe$nodes)

  if (object$isCompiled){
    rip <- object$rip
    cl  <- rip$clique
    se  <- rip$separators
    pa  <- rip$pa

    cl2 <- sapply(object$rip$cliques,length)

    cat(sprintf(" Number of cliques:              %4d \n",  length(cl2)))
    cat(sprintf(" Maximal clique size:            %4d \n",  max(cl2)))
    cat(sprintf(" Maximal state space in cliques: %4d \n",
                max(unlistPrim(lapply(object$equipot, length)))))

    if(length(e<-getEvidence(object))){
        print(e)
    }

    switch(type,
           "rip"={
               cat("\nRIP ordering:\n")
               print(rip)
           },
           "cliques"={
               cat("\nCliques:\n")
               .printList(object$rip$cliques)
           },
           "configurations"={
               cat("\nConfigurations:\n")
               for (i in 1:length(cl)){
                   cat(" ", i, ":", object$equipot[[i]]$ncells, "\n")
               }
           })
}
  return(invisible(object))
}

