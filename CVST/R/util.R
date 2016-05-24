# data is a list
# x: either a list or a matrix containing the data rowwise
# y: vector of labels/values
constructData = function(x, y) {
  stopifnot(is.list(x) || is.vector(x) || is.matrix(x))
  stopifnot(is.list(y) || is.vector(y) || is.factor(y))
  data = list(x=x, y=y)
  class(data) = "CVST.data"
  return(data)
}

getN = function(data) {
  stopifnot(class(data) == "CVST.data")
  if (is.list(data$x) || is.vector(data$x)) {
    N = length(data$x)
  }
  else {
    N = nrow(data$x)
  }
  return(N)
}

shuffleData = function(data) {
  stopifnot(class(data) == "CVST.data")
  shuffle = sample.int(getN(data))
  return(getSubset(data, shuffle))
}

getSubset = function(data, subset) {
  stopifnot(class(data) == "CVST.data")
  x = getX(data, subset)
  y = data$y[subset]
  ret = constructData(x=x, y=y)
  return(ret)
}

getX = function(data, subset=NULL) {
  stopifnot(class(data) == "CVST.data")
  if (is.null(subset)) {
    ret = data$x
  }
  else {
    if (is.list(data$x) || is.vector(data$x)) {
      ret = data$x[subset]
    }
    else {
      ret = data$x[subset, ,drop=FALSE]
    }
  }
  return(ret)
}

isClassification = function(data) {
  stopifnot(class(data) == "CVST.data")
  return(is.factor(data$y))
}

isRegression = function(data) {
  stopifnot(class(data) == "CVST.data")
  return(!isClassification(data))
}

constructLearner = function(learn, predict) {
  stopifnot(is.function(learn) && is.function(predict))
  learner = list(learn=learn, predict=predict)
  class(learner) = "CVST.learner"
  return(learner)
}

constructCVSTModel = function(steps=10, beta=.1, alpha=.01, similaritySignificance=.05, earlyStoppingSignificance=.05, earlyStoppingWindow=3, regressionSimilarityViaOutliers=FALSE) {
  ret = list(steps=steps,
    beta=beta,
    alpha=alpha,
    similaritySignificance=similaritySignificance,
    earlyStoppingSignificance=earlyStoppingSignificance,
    earlyStoppingWindow=earlyStoppingWindow,
    regressionSimilarityViaOutliers=regressionSimilarityViaOutliers)
  class(ret) = "CVST.setup"
  return(ret)
}

constructParams = function(...) {
  pn = names(substitute(c(...)))[-1]
  ret = expand.grid(..., stringsAsFactors=FALSE, KEEP.OUT.ATTRS = FALSE)
  params = lapply(1:nrow(ret), function(ind) as.list(ret[ind, ]))
  paramNames = lapply(1:nrow(ret), function(ind) paste(pn, ret[ind, ], sep="=", collapse=" "))
  names(params) = paramNames
  class(params) = "CVST.params"
  return(params)
}


.getResult = function(train, test, learner, param, squared=TRUE) {
  stopifnot(class(learner) == "CVST.learner" && class(train) == "CVST.data" && class(test) == "CVST.data")
  model = try(learner$learn(train, param))
  if (class(model) == "try-error") {
    pred = rep(NA, length(test$y))
  }
  else {
    pred = try(learner$predict(model, test))
    if (class(pred) == "try-error") {
      pred = rep(NA, length(test$y))
    }
  }
  if (isClassification(test)) {
    res = (test$y != pred)
  }
  else {
    if (squared) {
      res = (pred - test$y)^2
    }
    else {
      res = (pred - test$y)
    }
  }
  return(res)
}

cochranq.test = function(mat) {
  cochransQtest = list(statistic = 0, parameter = 0, p.value = 1,
    method = "Cochran's Q Test",
    data.name = deparse(substitute(mat)))
  class(cochransQtest) = "htest"

  if (is.vector(mat) || any(dim(mat) <= 1)) {
    return(cochransQtest)
  }

  # we expect the individuals in the rows, repetitions/treatments in the columns
  m = ncol(mat)
  df = m - 1
  L = apply(mat, 1, sum)
  index = (L > 0 & L < m)
  if (sum(index) <= 1) {
    # all rows are either one or zero... no effect!
    return(cochransQtest)
  }
  
  if (sum(index) * m <= 24) {
    return(.perm.cochranq.test(mat[index, ]))
  }

  L = L[index]
  T = apply(mat[index, ], 2, sum)
  Q = ((m-1) * (m * sum(T^2) - sum(T)^2)) / (m * sum(L) - sum(L^2))
  names(df) = "df"
  names(Q) = "Cochran's Q"

  if (is.nan(Q)) {
    p.val = 1.0
  }
  else {
    p.val = pchisq(Q, df, lower.tail=FALSE)
  }
  cochransQtest$statistic = Q
  cochransQtest$parameter = df
  cochransQtest$p.value = p.val
  return(cochransQtest)
} 

.perm.cochranq.test = function(mat, nperm=1000) {
  if (is.vector(mat) || any(dim(mat) <= 1)) {
    cochransQtest = list(statistic = 0, parameter = 0, p.value = 1,
      method = "Cochran's Q Test",
      data.name = deparse(substitute(mat)))
    class(cochransQtest) = "htest"
    return(cochransQtest)
  }
  # we expect no straight zero or one-rows in mat
  m = ncol(mat)
  df = m - 1
  L = apply(mat, 1, sum)
  T = apply(mat, 2, sum)
  quot = (m * sum(L) - sum(L^2))
  Q = ((m-1) * (m * sum(T^2) - sum(T)^2)) / quot
  names(df) = "df"
  names(Q) = "Cochran's Q"
  
  permFun = function() {
    newPerm = mat
    for (i in 1:nrow(mat)) {
        newPerm[i, ] = mat[i, sample(m)]
      }
    T = apply(newPerm, 2, sum)
    Q = ((m-1) * (m * sum(T^2) - sum(T)^2)) / quot
    return(Q)
  }
  
  QS = replicate(nperm, permFun())
  p.value = mean(QS >= Q)
  cochransQtest = list(statistic = Q, parameter = df, p.value = p.value,
    method = "Cochran's Q Test (monte-carlo)",
    data.name = deparse(substitute(mat)))
  class(cochransQtest) = "htest"
  return(cochransQtest)
}

constructSequentialTest = function(piH0=.5, piH1=.9, beta, alpha) {
  a1 = log((1 - beta) / alpha) / (log(piH1 / piH0) + log((1 - piH0) / (1 - piH1)))
  a0 = -log(beta / (1 - alpha)) / (log(piH1 / piH0) + log((1 - piH0) / (1 - piH1)))
  b = log((1 - piH0) / (1 - piH1)) / (log(piH1 / piH0) + log((1 - piH0) / (1 - piH1)))
  ret = list(a1=a1, a0=a0, b=b, piH0=piH0, piH1=piH1, alpha=alpha, beta=beta)
  class(ret) = "CVST.sequentialTest"
  return(ret)
}

plotSequence = function(st, s) {
  y = cumsum(s)
  if (!is.null(st$steps)) {
    plot(y, xlim=c(1, st$steps), ylim=c(1, st$steps))
  }
  else {
    plot(y)
  }
  abline(a=st$a1, b=st$b, col="red")
  abline(a=-st$a0, b=st$b, col="red", lty=2)

  abline(h=0)
  abline(a=0, b=1)
  title(sprintf("one-sided H0:%0.2f; H1:%0.2f", st$piH0, st$piH1))
}


testSequence = function(st, s) {
  stopifnot(class(st) == "CVST.sequentialTest")  
  n = length(s)
  y = cumsum(s)
  ret = 0
  if (y[n] >= st$b * n + st$a1) {
    ret = 1
  }
  else if (y[n] <= st$b * n - st$a0) {
    ret = -1
  }
  return(ret)
}

noisySinc = function(n, dim=2, sigma=0.1) {
  if (length(n) > 1) {
    x = n
  }
  else {
    x = runif(n, -pi, pi)
  }
  sinc = function(d) sin(d) / (d)
  y = sinc(4 * x) + 0.2 * sin(15 * x * dim) + sigma*rnorm(n)
  y[is.nan(y)] = 1
  return(constructData(x=as.matrix(x), y=y))
}

noisySine = function(n, dim=5, sigma=.25) {
  x = runif(n, 0, 2 * pi * dim)
  y = sin(x)
  if (!is.null(sigma) && sigma > 0) {
    y = y + rnorm(n, sd=sigma)
  }
  label = factor(y == abs(y))  
  return(constructData(x=as.matrix(x), y=label))
}

noisyDonoho = function(n, fun=doppler, sigma=1) {
  x = matrix(runif(n, 0, 1), n, 1)
  y = as.vector(fun(x)) + rnorm(n, sd=sigma)
  return(constructData(x=x, y=y))
}

blocks = function(x, scale=3.656993) {
  t = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  ret = t(sapply(x, function(xx) (1 + sign(xx - t)) / 2)) %*% h
  
  ret = ret * scale
  return(ret)
}

bumps = function(x, scale=10.52884) {
  t = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  h = c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  w = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  ret = t(sapply(x, function(xx) (1 + abs((xx - t) / w))^-4 )) %*% h
  ret = ret * scale
  return(ret)
}

heavisine = function(x, scale=2.356934) {
  ret = 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - x)
  ret = ret * scale
  return(ret)
}

doppler = function(x, scale=24.22172) {
  ret = sqrt(x * (1 - x)) * sin((2.1 * pi) / (x + 0.05)) 
  ret = ret * scale
  return(ret)
}
