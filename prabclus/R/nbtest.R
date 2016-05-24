"nbtest" <-
  function(nblist, n.regions=length(nblist)){
    nbok <- TRUE
    if (!is.list(nblist)){
      cat("Neighborhood is not a list.\n")
      nbok <- FALSE
    }
    lnb <- length(nblist)
    if (lnb!=n.regions){
      cat("Neighborhood length",lnb,"does not match regions number",
          n.regions,".\n")
      nbok <- FALSE
    }
    for (i in 1:lnb){
      if (length(nblist[[i]])>0){
        if (min(nblist[[i]]<1))
          stop("Neighborhood list contains elements smaller than 1.")
        if (max(nblist[[i]]>n.regions))
          stop("Neighborhood list contains elements larger than number of regions.")
        for (n in nblist[[i]]){
          if (n==i){
            cat(i,"is neighbor of itself.\n")
            nbok <- FALSE
          }
          if (all(nblist[[n]]!=i)){
            cat(i,"is neighbor of",n,"but",n,"is not neighbor of",
                i,".\n")
            nbok <- FALSE
          }
        }
      }
    }
    if (!nbok) stop("Improper neighborhood list.")
    invisible(nbok)
  }


