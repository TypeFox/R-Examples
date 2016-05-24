getroc <- function(yhat, y) {
  TPR = apply(yhat, 2, f <- function(x) {
    sum((x == y) & (y == 1))/sum(y == 1)
  })
  FPR = apply(yhat, 2, f <- function(x) {
    sum((x != y) & (y == 0))/sum(y == 0)
  })
  cbind(FPR, TPR)
}

get.vio <- function(roc.th, roc.em, nproc.em, alphalist) {
  alphaall = sort(c(roc.th[, 1], roc.em[, 1]))
  t1 = approx(roc.th[, 1], roc.th[, 2], alphaall, rule = 2)$y
  t2 = approx(roc.em[, 1], roc.em[, 2], alphaall, rule = 2)$y
  loc = which(t2 > t1)
  vio1 = sum(diff(alphaall)[loc[loc < length(alphaall)]])
  alphaall2 = sort(c(roc.th[, 1], alphalist))
  t3 = approx(roc.th[, 1], roc.th[, 2], alphaall2, rule = 2)$y
  t4 = approx(alphalist, nproc.em[, 2], alphaall2, rule = 2)$y
  loc = which(t4 > t3)
  vio2 = sum(diff(alphaall)[loc[loc < length(alphaall2)]])
  return(c(vio1, vio2))
}

getalpha <- function(roc1, roc2) {
  alphaall = sort(c(roc1[, 1], roc2[, 1]))
  t1 = approx(roc1[, 1], roc1[, 2], alphaall, rule = 2)$y
  t2 = approx(roc2[, 1], roc2[, 2], alphaall, rule = 2)$y
  alphaall[which.min(diff(sign(t1 - t2) != 0)) + 1]
}









## try replace bootstrap with subsampling interpolation
bs <- function(prob, method = c("bs", "ss"), alpha = c(0.01, 0.05), delta = 0.05,
               B = 1000, n.cores = 1) {
  ### prob: probability of assigning to 1 method: bs: boostrap, ss: subsampling n/2
  ### of n alpha: a vector of type-I error delta: violation rate control B: number of
  ### bootstraps set.seed(0)
  n = length(prob)
  cutoff.list = as.vector(sort(prob))
  if (method == "bs") {
    cutoff.bs = matrix(cutoff.list[sample(x = 1:n, size = n * B, replace = TRUE)],
                       ncol = B)  # store the bootstrapped lr.Hat into a matrix. Every column is a bootstrap sample
  } else if (method == "ss") {
    cutoff.bs = matrix(0, round(n/2), B)
    for (i in 1:B) {
      cutoff.bs[, i] = cutoff.list[sample(1:n, size = round(n/2), replace = FALSE)]
    }
  }
  if (n.cores > 1) {
    percents <- mcmapply(cutoff.list, FUN = function(cutoff) {
      # for every possible cutoff, find the percentage of empirical type I errors >
      # alpha
      temp1 = cutoff.bs - cutoff
      temp2 = matrix(as.numeric(temp1 > 0), ncol = B)
      temp3 = colSums(temp2)/n
      apply(outer(temp3, alpha, ">"), 2, sum)/B
    }, mc.cores = n.cores)
  } else {
    percents <- sapply(cutoff.list, FUN = function(cutoff) {
      # for every possible cutoff, find the percentage of empirical type I errors >
      # alpha
      temp1 = cutoff.bs - cutoff
      temp2 = matrix(as.numeric(temp1 > 0), ncol = B)
      temp3 = colSums(temp2)/n
      apply(outer(temp3, alpha, ">"), 2, sum)/B
    })
  }

  if (is.matrix(percents) == FALSE) {
    percents = matrix(percents, 1, length(percents))
  }
  ind1vec = apply(percents, 1, f <- function(x) {
    min(which(x <= delta))
  })
  #ind2vec = apply(percents, 1, f <- function(x) {
  #  max(which(x > delta))
  #})
  w1 = abs(percents[ind1vec] - delta)
  #w2 = abs(percents[ind2vec] - delta)
  # cutoff = (w1* cutoff.list[ind2] + w2* cutoff.list[ind1])/(w1+w2)

  cutoff = cutoff.list[ind1vec]

  return(cutoff)

}
