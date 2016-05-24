# $Id: cscores.R,v 1.14.4.1 2003/11/03 15:16:50 hothorn Exp $

cscores <- function(y, ...) UseMethod("cscores")

cscores.default <- function(y, type=c("Data", "Wilcoxon", 
                            "NormalQuantile", "AnsariBradley", "Median",
                            "Savage", "ConSal"), int = FALSE, 
                            maxs=length(y), ... ) {
  type <- match.arg(type)
  if (!(all.equal(floor(maxs),maxs))  || maxs < 1) 
    stop("maxs is not an positiv integer")
  N <- length(y)
  RET <- switch(type, "Data" = y,
           "Wilcoxon" = rank(y),
           "NormalQuantile" = qnorm(rank(y)/(N+1)),
           "AnsariBradley" = { 
             r <- rank(y)
             pmin(r, N - r + 1) },
           "Median" = {
             r <- rank(y)
             r[r <= (N+1)/2] <- 0
             r[r > 0] <- 1
             r},
           "Savage" = {
             cscores.Surv(cbind(y, 1)) },
           "ConSal" = {
             (rank(y)/(N+1))^4
           }
  )
  attr(RET, "scores") <- type
  if (int) {
    fscore <- RET - floor(RET)
    if (all(fscore[fscore != 0] == 0.5)) 
      RET <- 2*RET
    else RET <- round(RET*maxs/max(RET))
  }
  RET
}

cscores.factor <- function(y, ...) {
  if (nlevels(y) > 2) stop("cannot compute scores for more than 2 levels")
  RET <- as.integer(y)
  attr(RET, "scores") <- "Median"
  RET
}

irank <- function(x, ox=NULL) {
  if (is.null(ox))
    .Call("irank", as.double(x), as.integer(order(x)-1), 
          PACKAGE="exactRankTests")
  else 
    .Call("irank", as.double(x), as.integer(ox-1), 
          PACKAGE="exactRankTests")
}

cscores.Surv <- function(y, type="LogRank", int=FALSE, maxs=nrow(y), ...) {
  type <- match.arg(type)
  if (!(all.equal(floor(maxs),maxs))  || maxs < 1)
    stop("maxs is not an positiv integer")
  time <- y[,1]
  event <- y[,2]
  N <- length(time)
  ot <- order(time)
  rt <- irank(time, ot)
  fact <- event/(N - rt + 1)
  scores <- event - cumsum(fact[ot])[rt]
  if (int)
    RET <- round(scores*maxs/max(scores))
  else
    RET <- scores
  attr(RET, "scores") <- "LogRank"
  RET
}

