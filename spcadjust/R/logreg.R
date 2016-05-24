#######################################
## Risk adjusted logistic regression ##
#######################################

#' @include main.R model.R
NULL

#' Data Model for Binary Responses using a logarithmic model and likelihood ratio updates.
#'
#' @param formula The formula of the model.
#' @param Delta This value will be added to the log odds ratio in the
#' out-of-control model in the likelihood ratio between the
#' out-of-control and the in-control model constituting the updates. 
#'
#' @examples
#' n <- 1000
#' Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
#' xbeta <- -1+Xlogreg$x1*100+Xlogreg$x2+Xlogreg$x3
#' Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))
#' chartlogreg <- new("SPCCUSUM",
#'                    model=SPCModellogregLikRatio(Delta= 1, formula="y~x1+x2+x3"))
#' SPCproperty(data=Xlogreg,nrep=10,property="calARL",
#'             chart=chartlogreg,params=list(target=100))
#' #increase nrep for real applications.
#' @export
SPCModellogregLikRatio <- function(formula,Delta=1){
    SPCModelNonpar(
        updates=function(xi,data){
              xbeta <- predict.glm(xi,newdata=data)
              response <-  model.response( model.frame( formula,data=data))
              Delta*response + log(1+exp(xbeta)) - log(1+exp(Delta+xbeta))
        },
        xiofP=function(P)
            glm(formula,data=P,family=binomial("logit"))
    )
}


#' Data Model for Binary Responses using a Logarithmic Model and observed minus expected updates.
#'
#' @param formula The formula of the model.
#' @param Delta Half of this value will be subtracted for every update.
#' @examples
#' n <- 1000
#' Xlogreg <- data.frame(x1=rbinom(n,1,0.4), x2=runif(n,0,1), x3=rnorm(n))
#' xbeta <- -1+Xlogreg$x1*100+Xlogreg$x2+Xlogreg$x3
#' Xlogreg$y <- rbinom(n,1,exp(xbeta)/(1+exp(xbeta)))
#' chartlogreg <- new("SPCEWMA",
#'                    model=SPCModellogregOE(Delta= 0, formula="y~x1+x2+x3"), lambda=0.8)
#' SPCproperty(data=Xlogreg,nrep=10,property="calARL",
#'             chart=chartlogreg,params=list(target=100))
#' #increase nrep for real applications.
#' @export
SPCModellogregOE <- function(formula,Delta=0){
    SPCModelNonpar(
        updates=function(xi,data){
          response <-  model.response( model.frame( formula,data=data))
          response - predict.glm(xi,newdata=data,type="response")-Delta/2
         },
        xiofP=function(P)
            glm(formula,data=P,family=binomial("logit"))
    )
}


