print.qtlnet <- function(x, cutoff = 0.01, digits = 3, ...)
{
  cat("\nModel averaged probabilities for edge direction (row -> col):\n")
  mav <- get.model.average(x)
  print(round(mav, digits))

  cat("\nPosterior probabilities by causal model:\n")
  pp <- get.posterior.prob(x)
  cutoff <- min(cutoff, max(pp$prob))
  wh <- which(pp$post.prob >= cutoff)
  print(pp[wh,, drop = FALSE])
  
  invisible(list(mav = mav, pp = pp))
}  
######################################################################
threshold.net <- function(qtlnet.object,
                          mav = get.model.average(qtlnet.object),
                          min.prob = 0.9, ...)
{
  ## This finds the directed edge (or none) that has highest post.prob.
  ## There may be two problems:
  ## 1. the highest post.prob may be modest (.33 to .67, say).
  ## 2. the mav may not be a DAG. Need to check.
  n <- nrow(mav)

  ## Make sure at least one edge is detected.
  max.mav <- max(mav)
  if(min.prob > max.mav)
    min.prob <- max.mav
  
  new <- mav
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      m1 <- mav[i,j]
      m2 <- mav[j,i]
      m3 <- 1 - m1 - m2
      aux <- which.max(c(m1,m2,m3))
      if(aux == 1) {
        new[j,i] <- 0
        if(new[i,j] < min.prob)
          new[i,j] <- 0
      }
      if(aux == 2) {
        new[i,j] <- 0
        if(new[j,i] < min.prob)
          new[j,i] <- 0
      }
      if(aux == 3) { 
        new[j,i] <- 0
        new[i,j] <- 0
      }
    }
  }
  attr(new, "min.prob") <- min.prob
  new
}
######################################################################
get.averaged.net <- function(qtlnet.object,
                             pheno.names = attr(qtlnet.object, "pheno.names"),
                             min.prob = 0.5, ...)
{
  new <- threshold.net(qtlnet.object, min.prob = min.prob, ...)
  min.prob <- attr(new, "min.prob")
  out <- M.2.lista(new, pheno.names)
  attr(out, "min.prob") <- min.prob

  out
}
######################################################################
summary.qtlnet <- function(object, parent.patterns = FALSE, ...)
{
  mav <- get.model.average(object)
  pheno.names <- attr(object, "pheno.names")

  freq.accept <- object$freq.accept
  if(length(freq.accept) > 1)
    freq.accept <- weighted.mean(freq.accept, attr(object, "nSamples"))
  
  out <- list(freq.accept = freq.accept,
              averaged.net = get.averaged.net(object, pheno.names, ...),
              posterior.table = averaged.posterior.table(mav, pheno.names))

  if(parent.patterns)
    out$parent.patterns <- parent.qtlnet(object)
  
  class(out) <- c("summary.qtlnet", "list")
  out
}
######################################################################
print.summary.qtlnet <- function(x, ...)
{
  min.prob <- attr(x$averaged.net, "min.prob")
  cat("\nModel-averaged network:")
  if(min.prob > 1/3)
    cat(paste(" (min.prob = ", round(min(min.prob, 1), 3), ")", sep = ""))
  cat("\n")
  print(x$averaged.net)

  cat("\nPosterior probabilities by direction:\n")
  print(x$posterior.table)

  cat("\nAcceptance frequency for MCMC:", x$freq.accept, "\n")

  if(!is.null(x$parent.patterns))
    print(x$parent.patterns, ...)
  
  invisible(x)
}
######################################################################
parent.qtlnet <- function(qtlnet.object)
{
  ## Split model into trait|depends), with "|" and "depends" optional.
  ## This is slow. Too many apply of apply operations.

  ## Model patterns.
  model <- qtlnet.object$post.model
  model <- unlist(strsplit(model, "(", fixed = TRUE))
  model <- model[model != ""]

  ## Parent patterns.
  parents <- substring(sub("^[0-9]*", "", model), 2)

  ## Size of parent patterns.
  tmp <- parents[parents != ""]
  n.parents <- rep(0, length(parents))
  n.parents[parents != ""] <- sapply(strsplit(tmp, ",", fixed = TRUE), length)

  ## How many different traits share the same parent set?
  tmp <- sapply(tapply(model, parents, table), length)
  tmp2 <- sapply(strsplit(names(tmp), ",", fixed = TRUE), length)
  scans <- sapply(tapply(tmp, tmp2,
                       function(x) c(solo = sum(x == 1), more = sum(x > 1))),
                function(x) x)
  scans <- rbind(scans, total = apply(scans, 2, sum))
  scans <- cbind(scans, total = apply(scans, 1, sum))
  
  out <- list(model = model, parents = parents, n.parents = n.parents,
              scans = scans)
  class(out) <- c("parent.qtlnet", "list")
  out
}
######################################################################
print.parent.qtlnet <- function(x, freq.max = 10, ...)
{
  mytable <- function(parents, freq.max = 10) {
    tmp <- table(parents)
    tmp2 <- "rest" ## paste(">=", freq.max, sep = "")
    tmp[tmp >= freq.max] <- tmp2
    tmp <- ordered(tmp, c(seq(freq.max - 1), tmp2))
    tmp <- c(table(tmp))
    tmp["total"] <- sum(tmp)
    tmp
  }
  
  cat("\nHow frequently are model patterns sampled?\n")
  print(mytable(x$model, freq.max = freq.max))
  
  cat("\nHow frequently are parent patterns sampled?\n")
  print(mytable(x$parents, freq.max = freq.max))
  
  cat("\nWhat are sizes of parent pattern sets?\n")
  tmp <- table(x$n.parents)
  tmp["total"] <- sum(tmp)
  print(tmp)
  
  cat("\nHow frequently are parent sets sampled (by set size)?\n")
  parents.size <- tapply(x$parents[x$n.parents > 0],
                         x$n.parents[x$n.parents > 0],
                         mytable, freq.max = freq.max)
  print(parents.size)

  cat("\nWhat are frequencies of each single-parent pattern?\n")
  tmp <- table(x$parents[x$n.parents == 1])
  names(tmp) <- substring(names(tmp), 1, nchar(names(tmp)) - 1)
  tmp <- tmp[order(as.numeric(names(tmp)))]
  print(tmp)
  
  cat("\nHow many scans are solo or more for each size parent pattern?\n")
  print(x$scans)
  
  ## Evidence suggest it is not worth doing multiple traits
  ## for n > le.pheno / 2.
}
######################################################################
######################################################################
averaged.posterior.table <- function(maM,nms)
{
  nn <- nrow(maM)
  np <- choose(nn,2)
  out <- data.frame(matrix(NA,np,5))
  k <- 1
  for(i in 1:(nn-1)){
    for(j in (i+1):nn){
      out[k,1] <- i
      out[k,2] <- j
      out[k,3] <- round(maM[i,j],3)
      out[k,4] <- round(maM[j,i],3)
      out[k,5] <- round(1 - out[k,3] - out[k,4], 3)
      k <- k + 1
    }
  }
  new.out <- out
  for(k in 1:np){
    new.out[k,1] <- nms[out[k,1]]
    new.out[k,2] <- nms[out[k,2]]
  }
  names(new.out) <- c("node1","node2","-->","<--","no")
  new.out
}
################################################################
M.2.lista <- function(M,pheno.names)
{
  le <- nrow(M)
  k <- 1
  le2 <- le*(le-1)/2
  out <- data.frame(cause = factor(rep(NA, le2), pheno.names),
                    effect = factor(rep(NA, le2), pheno.names),
                    prob = rep(NA, le2))

  for(i in 1:(le-1)){
    for(j in (i+1):le){
      if(M[i,j] != 0){
        out[k,1] <- pheno.names[i]
        out[k,2] <- pheno.names[j]
        out[k,3] <- M[i,j]
        k <- k + 1
      }
      if(M[j,i] != 0){
        out[k,2] <- pheno.names[i]
        out[k,1] <- pheno.names[j]
        out[k,3] <- M[j,i]
        k <- k + 1
      }
    }
  }
  out[1:(k-1),]
}
