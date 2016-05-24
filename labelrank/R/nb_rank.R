#' Predicting label rankings based on the naive Bayes ranking model
#' 
#' This function predicts the rankings given prior and conditional probabilities obtained from \code{model_nbr}
#' @param x is  \code{n x p} matrix of \code{n} observations and \code{p} training attributes and can have continuous or nominal values.
#' @param y is \code{n x j} matrix of label rankings
#' @param n is a parameter of 'memory'; that is, how fast past gets forgotten. (see details of \link{time_weights}).
#' @param new.x is a vector of new attributes
#' @return a numeric vector of ranking
#' @details
#' This function predicts a ranking for \code{test.x} attributes. It initially builds a model for naive Bayes algorithm that calculates priors and conditional label ranking probabilities and then use them to predict rankings. The attributes can be nominal or continuous data.
#' @examples
#' train.x <- lr.nom[1:16,]
#' test.x <- lr.nom[17,]
#' predrank <- nb_rank(train.x,y,test.x,n=1)
#' @importFrom stats dnorm
#' @export
nb_rank <- function(x, y, new.x, n=1) {
  model <- model_nbr(x, y, n)
  if (all(is.numeric(new.x))) {
    prod <- -log(model$priors) + sapply(1:nrow(y), function(r) {
      cond <- sapply(seq_along(new.x), function(i) {
        dnorm(new.x[i], mean = model$cond$mean[r, i], sd = model$cond$sdev[r, i])
      })
      sum(-log(cond), na.rm = T)
    })
    y[which.min(prod), ]
  } else {
    prod <- -log(model$priors) + apply(-log(sapply(new.x, function(i) {
      cond <- model$cond[rownames(model$cond) == i, ]
      if (is.matrix(cond)) {
        apply(cond, 2, sum, na.rm = T)
      } else {
        cond
      }
    }
    )), 1, sum)
    y[which.min(prod), ]
  }
}
