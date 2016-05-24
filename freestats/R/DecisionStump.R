
#' @title Decision Stump Algorithm
#' @description Do classification with tree method in one step
#' @export decisionStump
#' @return \item{j}{The best dimention to cut the tree}
#' \item{theta}{Value that seperate tree in the best dimention}
#' \item{m}{the routine label value (for now only 1)}
#' @author Xiaoyao Yang
#' @param X Data matrix / Data frame
#' @param w Weight that given to each observation. Used in calculate cost function. 
#' @param y Class label for data points in X, must be -1 or 1
#' @examples
#' 
#'set.seed(1024)
#'z <- runif(n=5)
#'mydata <- fakedata(w=z,n=100)
#'X<- mydata$S[,1:4]
#'y <- mydata$y
#'w <- rep(1/100,100)
#'pars <- decisionStump(X=X,w=w,y=y)
#' 

#consider factor!!!


# given predictor matrix X, with label y, and weight of data. return weak learner pars(j,theta,m)
decisionStump <- function(X,w,y)
{
    if(class(w)!='matrix'){
        w <- as.matrix(w)
    }
    #     DecisionStump <- function(X,w,y) 
    #     decision stump
    #     w must be 1 column, y=-1/1
    #     find best theta for each dim
    theta <- vector()
    best.cost <- vector()
    Compare <- 10
    for(d in 1:dim(X)[2]){
        if (!is.factor(X[,d]))
            {
            for(n in 1:dim(X)[1])
                {
                yhat <- 2*(X[,d]>X[n,d])-1
                cost.temp <- t(w)%*%(y!=yhat)/sum(w)
                if (cost.temp<Compare) 
                    {
                    Compare <- cost.temp
                    best.theta <- X[n,d] 
                    best.j<-d
                    best.n<-n    
                    }
                 }
             }
        else stop('only for numeric input now!')
    }
    res <-list(j=best.j,theta=best.theta,m=1)
    class(res) <- 'ds'
    return(res)
}
