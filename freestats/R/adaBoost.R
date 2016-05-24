

#' @title Adaboost algorithm
#' @description Do classification using adaboost algorithm with decisionStump as weak learner
#' @export adaBoost
#' @details
#' Train function can be any weak learner algorithm. For now, train function must has form train(X,w,y,...). see more in \code{\link{decisionStump}}
#' 
#' If you have any good weak learner but can't use it in this function, feel free to let me know.
#' 
#' @return \item{alpha}{The weight for different weak learners}
#' \item{allPars}{A list of parameters for different weak learners}
#' @author Xiaoyao Yang
#' @param train Function of weak learner that would be used in adaboost, must have form train(dat.train,w,y.train)
#' @param dat.train Training data set
#' @param y.train Label for training data set
#' @param B Number of weak learners that will used
#' @param \dots Other parameters that need to passed in train function
#' @examples
#' set.seed(1024)
#' z <- runif(n=5)
#' mydata <- fakedata(w=z,n=100)
#' X<- mydata$S[,1:4]
#' y <- mydata$y
#' res <- adaBoost(dat.train=X,y.train=y,B=3)





adaBoost <- function(train=decisionStump,dat.train,y.train,B=10,...)
{
    #implement boosting
    if (!is.data.frame(y.train)) {y.train<-data.frame(y.train=y.train)}
    allPars <- matrix(list())
    n<-dim(y.train)[1]
    if (n!=dim(dat.train)[1]){stop('data and label must have same length')}
    
    w=rep(1/n,n)
    alpha <- vector()
    for (b in 1:B){
#         pars <- do.call(train,c(list(dat.train,w,y.train),list(...)))
        pars <- train(dat.train, w, y.train, ...)
        yhat <- classify(pars,dat.train)
        error <- t(w)%*%(yhat!=y.train)/sum(w)
        vote <- as.numeric(log((1-error)/error))
        w <- w*exp(vote*(yhat!=y.train))
        alpha[b] <- vote
        allPars[[b]] <- pars
    }
    res <- list(alpha=alpha,allPars=allPars)
    class(res) <- 'ab'
    return(res)
}
