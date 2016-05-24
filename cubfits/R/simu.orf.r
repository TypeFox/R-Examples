### This file contains functions generating sequences for simulation studies.

simu.orf <- function(n, b.Init, phi.Obs = NULL, AA.prob = NULL,
    orf.length = NULL, orf.names = NULL, model = .CF.CT$model){
  ### Check n.
  if(n <= 0){
    stop("n should be a positive integer.")
  }

  ### Check b.Init.
  if(!is.list(b.Init)){
    stop("b.Init should be a list.")
  }
  if(!all(names(b.Init) %in% .CF.GV$amino.acid.split)){
    stop("names of b.Init is incorrect.")
  }
  tmp <- do.call("c", lapply(b.Init, function(i.b) is.matrix(i.b$coef.mat)))
  if(!all(tmp)){
    stop("b.Init[[i]]$coef.mat should be in matrix format for all i.")
  }

  ### Check phi.Obs
  if(is.null(phi.Obs)){
    phi.Obs <- rlnorm(n)
  }
  if(!all(phi.Obs > 0)){
    stop("all phi.Obs should be positive.")
  }
  if(length(phi.Obs) != n){
    stop("length(phi.Obs) != n.")
  }

  ### Check AA.prob.
  if(is.null(AA.prob)){
    AA.prob <- rep(1, length(b.Init))
  }
  if(length(AA.prob) != length(b.Init)){
    stop("AA.prob should be corresponding to b.Init.")
  }
  if(all(AA.prob > 0)){
    AA.prob <- AA.prob / sum(AA.prob)
  } else{
    stop("AA.prob should be all positilve.")
  }

  ### Check orf.length.
  if(is.null(orf.length)){
    orf.length <- sample(10:20, n, replace = TRUE)
  }
  if(length(orf.length) != n){
    stop("orf.length is incorrect")
  }

  ### Check orf.names.
  if(length(orf.names) != n){
    orf.names <- paste("ORF", 1:n, sep = "")
  }

  ### Check model.
  if(model[1] == "rocnsef"){
    simu.orf.model <- simu.orf.rocnsef
  } else if(model[1] == "roc"){
    simu.orf.model <- simu.orf.roc
  } else if(model[1] == "nsef"){
    simu.orf.model <- simu.orf.nsef
  } else{
    stop("model is not found.")
  }

  ### Run the function.
  ret <- lapply(1:n, function(i.n){
                       simu.orf.model(orf.length[i.n], b.Init, phi.Obs[i.n],
                                      AA.prob)
                     })
  names(ret) <- orf.names
  ret
} # End of simu.orf().


simu.orf.rocnsef <- function(orf.length, b.Init, phi.Obs, AA.prob){
  aa <- names(b.Init)
  if(length(aa) == 0){
    stop("AA names are not found.")
  }
  orf <- sample(aa, orf.length, replace = TRUE, prob = AA.prob)

  if("Z" %in% aa){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- rep(NA, orf.length) 
  for(i.aa in 1:orf.length){
    scodon <- synonymous.codon[[orf[i.aa]]]
    if(length(scodon) == 1){
      ret[i.aa] <- scodon
    } else{
      x <- cbind(1, phi.Obs, phi.Obs * i.aa)
      exponent <- x %*% b.Init[[orf[i.aa]]]$coef.mat
      scodon.prob <- my.inverse.mlogit(exponent)
      ret[i.aa] <- sample(scodon, 1, prob = scodon.prob)
    }
  }

  ret <- paste(c("ATG", ret, "TAA"), collapse = "")
  ret
} # End of simu.orf.rocnsef().

simu.orf.roc <- function(orf.length, b.Init, phi.Obs, AA.prob){
  aa <- names(b.Init)
  if(length(aa) == 0){
    stop("AA names are not found.")
  }
  orf <- sample(aa, orf.length, replace = TRUE, prob = AA.prob)

  if("Z" %in% names(b.Init)){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- rep(NA, orf.length)
  for(i.aa in 1:length(aa)){
    scodon <- synonymous.codon[[aa[i.aa]]]
    id <- orf == aa[i.aa]
    if(sum(id) > 0){
      if(length(scodon) == 1){
        ret[id] <- scodon
      } else{
        x <- cbind(1, phi.Obs)
        exponent <- x %*% b.Init[[aa[i.aa]]]$coef.mat
        scodon.prob <- my.inverse.mlogit(exponent) 
        ret[id] <- sample(scodon, sum(id), replace = TRUE, prob = scodon.prob)
      }
    }
  }

  ret <- paste(c("ATG", ret, "TAA"), collapse = "")
  ret
} # End of simu.orf.roc().

simu.orf.nsef <- function(orf.length, b.Init, phi.Obs, AA.prob){
  aa <- names(b.Init)
  if(length(aa) == 0){
    stop("AA names are not found.")
  }
  orf <- sample(aa, orf.length, replace = TRUE, prob = AA.prob)

  if("Z" %in% aa){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- rep(NA, orf.length) 
  for(i.aa in 1:orf.length){
    scodon <- synonymous.codon[[orf[i.aa]]]
    if(length(scodon) == 1){
      ret[i.aa] <- scodon
    } else{
      x <- cbind(1, phi.Obs * i.aa)
      exponent <- x %*% b.Init[[orf[i.aa]]]$coef.mat
      scodon.prob <- my.inverse.mlogit(exponent)
      ret[i.aa] <- sample(scodon, 1, prob = scodon.prob)
    }
  }

  ret <- paste(c("ATG", ret, "TAA"), collapse = "")
  ret
} # End of simu.orf.nsef().

