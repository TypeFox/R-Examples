library(doMC)

registerDoMC()

x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000

ptime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})[3]

cat(sprintf('Parallel time using doMC on %d workers: %f\n',
            getDoParWorkers(), ptime))

stime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %do% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})[3]

cat(sprintf('Sequential time: %f\n', stime))
cat(sprintf('Speed up for %d workers: %f\n',
            getDoParWorkers(), round(stime / ptime, digits=2)))
