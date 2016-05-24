################################
## Risk Adjusted Linear Model ##
################################

#' @include main.R model.R
NULL

#' Data model based on a linear model
#'
#' The parameters needed for running the chart are the fitted linear
#' model. Resampled data sets are created by resampling cases with
#' replacement (i.e. keeping observations together).
#'
#'
#' @param formula the linear model specified as a formula.
#' @param Delta Object of class \code{"numeric"}.
#'
#' @examples
#' n <- 1000
#' Xlinreg <- data.frame(x1= rbinom(n,1,0.4), x2= runif(n,0,1), x3= rnorm(n))
#' Xlinreg$y <- 2 + Xlinreg$x1 + Xlinreg$x2 + Xlinreg$x3 + rnorm(n)
#'\dontrun{
#' chartlinregCUSUM <- new("SPCCUSUM", model=SPCModellm(Delta=1,formula="y~x1+x2+x3"))
#' SPCproperty(data=Xlinreg,nrep=10,property="calARL",
#'             chart=chartlinregCUSUM,params=list(target=100))
#' #increase nrep in real applications.
#'#' chartlinregCUSUM2 <- new("SPCCUSUM",model=SPCModellm(Delta=1,formula="y~x1"))
#' SPCproperty(data=Xlinreg,nrep=10,property="calARL",
#'             chart=chartlinregCUSUM2,params=list(target=100))
#' #increase nrep in real applications.
#'
#' chartlinregEWMA <- new("SPCEWMA", model=SPCModellm(Delta=0,formula="y~x1+x2+x3"),lambda=0.8)
#' SPCproperty(data=Xlinreg,nrep=10,property="calARL",
#'             chart=chartlinregEWMA,params=list(target=100))
#' #increase nrep in real applications.
#'
#' chartlinregEWMA2 <- new("SPCEWMA",model=SPCModellm(Delta=0,formula="y~x1"),lambda=0.8)
#' SPCproperty(data=Xlinreg,nrep=10,property="calARL",
#'             chart=chartlinregEWMA2,params=list(target=100))
#' }
#' #increase nrep in real applications.
#' @export
SPCModellm <- function(formula,Delta=0){
    SPCModelNonpar(
        updates=function(xi,data){
            response <-  model.response( model.frame( formula,data=data))
            response - predict(xi,newdata=data)-Delta/2
        },
        xiofP=function(P)
            lm(formula,data=P)
    )
}
