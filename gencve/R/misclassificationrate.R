`misclassificationrate` <-
  function(y, yp){
    tb<-table(y, yp)
    1-sum(diag(tb))/sum(tb)
  }

logloss <- function(y, yp) {
  tb<-table(y, yp)
  adjprob <- function(a,b) {
    eps <- 10^(-15)
    prob <- a/b
    prob <- ifelse(prob<eps, eps, prob)
    prob <- ifelse(prob>1-eps, 1-eps, prob)
    prob
  }
  probs <- sweep(tb, 1, rowSums(tb), FUN=adjprob)
  -sum(tb*log(probs))
}


#logloss <- function(CT) {
#  adjprob <- function(a,b) {
#    eps <- 10^(-15)
#    prob <- a/b
#    prob <- ifelse(prob<eps, eps, prob)
#    prob <- ifelse(prob>1-eps, 1-eps, prob)
#    prob
#  }
#  probs <- sweep(CT, 1, rowSums(CT), FUN=adjprob)
#  -sum(CT*log(probs))
#}
#logloss <- function(y, pr) {
#    eps <- 1e-15; #pr=confusion matrix
#    nr <- nrow(pr)
#    pr <-  matrix(sapply(pr, function(x) max(eps,x)), nrow = nr)
#    pr <- matrix(sapply(pr, function(x) min(1-eps,x)), nrow = nr)
#    -sum(y*log(pr) + (1-y)*log(1-pr))/nrow(y)
#  }


