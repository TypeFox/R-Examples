clique.var.test <- function(y1, y2, dag, alpha) {
  vs <- c(substitute(y1), substitute(y2), substitute(dag))

  p <- .procParams(y1, y2, dag)
  l <- .runCliqueVarTest(p$y1, p$y2, p$graph, alpha)

  l$var.names <- vs
  attr(l, "class") <- "clique.var.test"
  return(l)
}

clique.mean.test <- function(y1, y2, dag, alpha, perm.num=1000, paired=FALSE) {
  vs <- c(substitute(y1), substitute(y2), substitute(dag))

  l <- .procParams(y1, y2, dag)
  y1 <- l$y1
  y2 <- l$y2

  cli.test   <- .runCliqueVarTest(y1, y2, l$graph, alpha)
  check      <- cli.test$var.equal
  cliques    <- cli.test$cliques
  clique.num <- length(cliques)

  alpha.obs  <- vector("numeric", clique.num)
  t.obs      <- vector("numeric", clique.num)
  for (i in seq_len(clique.num)) {
    cli    <- unlist(cliques[i])
    y1.cli <- y1[,cli, drop=FALSE]
    y2.cli <- y2[,cli, drop=FALSE]

    if (length(cli) != 1) {
      if (paired) {

        y.diff <- y1 - y2
        y.bar <- colMeans(y.diff)
        y.centr <- y.diff - y.bar
        y.num <- nrow(y1)
        y.s <- (t(y.centr) %*% y.centr) / y.num
        t2 <- y.num * (t(y.bar) %*% solve(y.s) %*% y.bar)

        p <- ncol(y1)
        np <- y.num - p
        t.value <- t2 * np / (p * (y.num-1))

        r <- list(alpha.obs=1-pf(t.value, p, np),
                  t.obs=t.value)

      } else if (check[i]) {
        r <- .mult.test(y1.cli, y2.cli, perm.num)

      } else {
        r <- .hote(y1.cli, y2.cli, TRUE)
      }

      alpha.obs[i] <- r$alpha.obs
      t.obs[i]     <- r$t.obs

    } else {
      r            <- t.test(y1.cli, y2.cli, paired=paired)
      alpha.obs[i] <- r$p.value
      t.obs[i]     <- r$statistic
    }
  }

  l <- list(t.value=t.obs,
            p.value=as.numeric(alpha.obs))

  if (!paired) {
    l$lambda.value <- cli.test$lambda.value
    l$p.value.var  <- cli.test$p.value
    l$var.equal    <- check
  }

  l$cliques <- cliques
  l$graph   <- cli.test$graph

  l$var.names <- vs
  attr(l, "class") <- "clique.mean.test"
  return(l)
}

.runCliqueVarTest <- function(y1, y2, graph, alpha) {
  cliques <- graph$cli.tg
  maxCliqueSize <- max(sapply(cliques, length))
  if (nrow(y1) <= maxCliqueSize)
    stop("y1 should have more than ", maxCliqueSize, " rows (samples)")
  else if (nrow(y2) <= maxCliqueSize)
    stop("y2 should have more than ", maxCliqueSize, " rows (samples)")

  cov <- .estimateCov(y1, y2)

  clique.num <- length(cliques)

  alpha.obs  <- rep(0,     clique.num)
  lambda.obs <- rep(0,     clique.num)
  check      <- rep(FALSE, clique.num)

  for (i in seq_along(cliques)) {
    cli <- unlist(cliques[i])
    p   <- length(cli)

    s1.hat <- cov$s1[cli, cli, drop=FALSE]
    s2.hat <- cov$s2[cli, cli, drop=FALSE]
    s.hat  <- cov$s[cli, cli, drop=FALSE]

    s1.det <- det(s1.hat)
    s2.det <- det(s2.hat)
    s.det  <- det(s.hat)

    lambda.obs[i] <- nrow(y1)*log(s.det/s1.det) + nrow(y2)*log(s.det/s2.det)
    alpha.obs[i]  <- 1 - pchisq(lambda.obs[i], p*(p+1)/2)

    if (alpha.obs[i] <= alpha)
      check[i] <- TRUE
  }

  list(lambda.value=lambda.obs,
       p.value=alpha.obs,
       cliques=cliques,
       var.equal=check,
       graph=graph$tg)
}
