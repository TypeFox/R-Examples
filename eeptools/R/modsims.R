##' Generate prediction intervals for model functions
##'
##' Generate prediction intervals from R models following Gelman and Hill
##'
##' @param mod Name of a model object such as \code{\link{lm}}, \code{\link{glm}}, or \code{merMod}
##' @param newdata Sets of new data to generate predictions for
##' @param n.sims Number of simulations per case
##' @param na.omit Logical indicating whether to remove NAs from \code{newdata}
##' @return A dataframe with newdata and prediction intervals
##' @references Modified from Gelman and Hill 2006. Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge University Press.
##' @details Currently gelmansim does not work for \code{\link{lm}} objects because of the way \code{\link{sim}} in the 
##' \code{arm} package handles variable names for these objects. It is recommended users use \code{\link{glm}} in these cases.
##' @export
##' @import arm
##' @import memisc
##' @examples
##'  #Examples of "sim" 
##' set.seed (1)
##' J <- 15
##' n <- J*(J+1)/2
##' group <- rep (1:J, 1:J)
##' mu.a <- 5
##' sigma.a <- 2
##' a <- rnorm (J, mu.a, sigma.a)
##' b <- -3
##' x <- rnorm (n, 2, 1)
##' sigma.y <- 6
##' y <- rnorm (n, a[group] + b*x, sigma.y)
##' u <- runif (J, 0, 3)
##' y123.dat <- cbind (y, x, group)
##' # Linear regression 
##' x1 <- y123.dat[,2]
##' y1 <- y123.dat[,1]
##' M1 <- glm (y1 ~ x1)
##'
##' cases <- data.frame(x1 = seq(-2, 2, by=0.1))
##' sim.results <- gelmansim(M1, newdata=cases, n.sims=200, na.omit=TRUE)
##' \dontrun{
##' 
##' dat <- as.data.frame(y123.dat)
##' M2 <- glm (y1 ~ x1 + group, data=dat)
##' 
##' cases <- expand.grid(x1 = seq(-2, 2, by=0.1), 
##'                      group=seq(1, 14, by=2))
##' 
##' sim.results <- gelmansim(M2, newdata=cases, n.sims=200, na.omit=TRUE)
##'  
##' }
gelmansim <- function(mod, newdata, n.sims, na.omit=TRUE){
    sims.tmp <- arm::sim(mod, n.sims=n.sims)
  if("lm" %in% class(mod)[1]){
    form.tmp <- as.formula(paste("~",as.character(formula(mod)[3])))
  } else{
    form.tmp <- mod$formula[-2]
  }
 
  
  # generate factors matrix from model
  not.numeric <- function(x) ifelse(is.numeric(x), FALSE, TRUE)
  
  # Get factors
  facs <- mod$model[,sapply(mod$model, not.numeric)]
  if(length(facs) != 0){
  if(class(facs) == "data.frame"){
    contr.tmp <- lapply(facs, unique)
    for(i in names(contr.tmp)){
      newdata[, i] <- factor(newdata[, i], levels=unlist(contr.tmp[i]))
      X.tilde <- model.matrix(form.tmp, newdata, xlev=contr.tmp, 
                              contrasts.arg = lapply(mod$model[,sapply(mod$model, 
                                                                       is.factor)],
                                                     contrasts, contrasts=FALSE))
    }
  } else {
    contr.tmp <- levels(facs)
    names.tmp <- sapply(mod$model, not.numeric)
    facname <- names(names.tmp[names.tmp==TRUE])
    newdata[, facname] <- factor(newdata[, facname], 
                                 levels=contr.tmp)
    contr.list <- list(contrasts(mod$model[, facname], 
                                 contrasts=FALSE))
    names(contr.list) <- facname
    X.tilde <- model.matrix(form.tmp, newdata, xlev=contr.tmp, 
                            contrasts.arg = contr.list)
  }
  } else {
    X.tilde <- model.matrix(form.tmp, newdata)
  }
  
  if(na.omit==TRUE){
    newdata <- na.omit(newdata)
  }
  
  pred.tmp <- unlist(dimnames(sims.tmp@coef)[2])
  X.tilde <- X.tilde[, unlist(colnames(X.tilde)) %in% pred.tmp]
  X.tilde <- reorder(X.tilde, dim=2, names=pred.tmp)
  n.tilde <- nrow(X.tilde)
  y.tilde <- array(NA, c(n.sims, n.tilde))
  
  for(s in 1:n.sims){
    y.tilde[s,] <- rnorm(n.tilde, X.tilde %*% sims.tmp@coef[s,], 
                         sims.tmp@sigma[s])
  }
  
  yhats <- colMeans(y.tilde)
  yhatMin <- apply(y.tilde, 2, quantile, probs=0.2)
  yhatMax <- apply(y.tilde, 2,  quantile, probs=0.8)
  
  sim.results <- cbind(newdata, yhats)
  sim.results <- cbind(sim.results, cbind(yhatMin, yhatMax))
  sim.results <- as.data.frame(sim.results)
  return(sim.results)
}