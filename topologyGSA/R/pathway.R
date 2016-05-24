pathway.var.test <- function(y1, y2, dag, alpha,
                             variance=FALSE, s1=NULL, s2=NULL) {
  vs <- c(substitute(y1), substitute(y2), substitute(dag))

  l  <- .procParams(y1, y2, dag)
  y1 <- l$y1
  y2 <- l$y2

  res           <- .runPathwayVarTest(y1, y2, l$graph, alpha, variance, s1, s2)
  ns            <- nodes(res$graph)
  res$cli.moral <- lapply(res$cli.moral, function(ixs) ns[ixs])

  res$var.names <- vs
  attr(res, "class") <- "pathway.var.test"
  return(res)
}

pathway.mean.test <- function(y1, y2, dag, alpha,
                              perm.num=10000, variance=TRUE, paired=FALSE) {
  vs <- c(substitute(y1), substitute(y2), substitute(dag))

  l  <- .procParams(y1, y2, dag)
  y1 <- l$y1
  y2 <- l$y2

  y      <- rbind(y1, y2)
  y.num  <- nrow(y1) + nrow(y2)

  l <- .runPathwayVarTest(y1, y2, l$graph, alpha, variance)
  cli.moral <- l$cli.moral
  var.equal <- l$var.equal

  if (paired) {

    t.obs <- .hotePaired(y1, y2, cli.moral)

    stat.perm <- vector("numeric", perm.num)
    for (i in seq_len(perm.num))
      stat.perm[i] <- .hotePaired(y1, y2, cli.moral, perm=TRUE)

    p.value <- sum(stat.perm >= t.obs) / perm.num

    l <- list(p.value=p.value,
              cli.moral=cli.moral,
              graph=l$graph,
              t.value=t.obs)

  } else {

    s <- .hote(y1, y2, FALSE, cli.moral)

    y1.num <- nrow(y1)
    stat.perm <- vector("numeric", perm.num)
    for (i in seq_len(perm.num)) {
      ind          <- sample(y.num)
      y1.perm      <- y[ind[1:y1.num],]
      y2.perm      <- y[ind[(y1.num+1):y.num],]
      stat.perm[i] <- .hote(y1.perm, y2.perm, FALSE, cli.moral)$t.obs
    }

    p.value <- sum(stat.perm >= s$t.obs) / perm.num

    l <- list(t.value=s$t.obs,
              df.mean=s$df,
              p.value=p.value,
              lambda.value=l$lambda.value,
              df.var=l$df,
              p.value.var=l$p.value,
              qchisq.value=l$qchisq.value,
              var.equal=var.equal,
              cli.moral=cli.moral,
              graph=l$graph)
    }

  l$var.names <- vs
  attr(l, "class") <- "pathway.mean.test"
  return(l)
}

.runPathwayVarTest <- function(y1, y2, graph, alpha, variance,
                               s1=NULL, s2=NULL) {
  cliques <- graph$cli.moral
  maxCliqueSize <- max(sapply(cliques, length))
  if (nrow(y1) <= maxCliqueSize)
    stop("y1 should have more than ", maxCliqueSize, " rows (samples)")
  else if (nrow(y2) <= maxCliqueSize)
    stop("y2 should have more than ", maxCliqueSize, " rows (samples)")

  if (is.null(s1) != is.null(s2)) {
    stop("You must provide both s1 and s2 or neither.")

  } else if (is.null(s1) && is.null(s2)) {
    cov <- .estimateCov(y1, y2)

    s1.hat <- qpIPF(cov$s1, cliques)
    s2.hat <- qpIPF(cov$s2, cliques)
    s.hat  <- qpIPF(cov$s,  cliques)

  } else {

    s1.hat <- s1
    s2.hat <- s2
    s.hat  <- .estimateCovPool(nrow(y1), nrow(y2), s1, s2)
  }

  s1.det <- det(s1.hat)
  s2.det <- det(s2.hat)
  s.det  <- det(s.hat)

  lambda.value <- nrow(y1)*log(s.det/s1.det) + nrow(y2)*log(s.det/s2.det)
  df           <- (sum(graph$adj.moral)/2) + ncol(y1)
  qchisq.value <- qchisq(1 - alpha, df)
  p.value      <- 1 - pchisq(lambda.value, df)
  var.equal    <- p.value <= alpha

  res <- list(lambda.value=lambda.value,
              p.value=p.value,
              df=df,
              qchisq.value=qchisq.value,
              cli.moral=cliques,
              var.equal=var.equal)

  if (variance) {
    res$s1 <- s1.hat
    res$s2 <- s2.hat
  }

  res$graph <- graph$moral
  return(res)
}
