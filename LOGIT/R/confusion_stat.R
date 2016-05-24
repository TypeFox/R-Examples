# Hilbe, J.M., Practical Guide to Logistic Regression 2015
# Rafael de Souza, Eotvos Lorand Univ.
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#' @title  Display confusion or classification matrix following logistic regression using glm with binomial family.
#' @description Provides a confusion matrix of classification statistics following logistic regression.
#' @aliases confusion_stat
#' @usage confusion_stat(pred=pred,obs=obs)
#' @importFrom stats ftable addmargins
#' @format \describe{
#' \item{x}{
#' The function has two arguments: predicted values, response values}
#' }
#' @param pred Predicted values
#' @param obs  Observed values
#' @return confusion matrix
#' @import caret e1071
#' @note confusion_stat() must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user.
#'@examples
#' library(MASS)
#' library(LOGIT)
#' data(R84)
#' R84$cage <- R84$age - mean(R84$age)
#' R84$cdoc <- R84$docvis - mean(R84$docvis)
#' mylogit <- glm(outwork ~ cdoc + female + kids + cage + factor(edlevel),
#' family=binomial, data=R84)
#' mu <- predict(mylogit, type="response")
#' cutpoint<-ROCtest(mylogit, fold=10, type="Sensitivity")$cut
#' mu[mu>=cutpoint]<-1
#' mu[mu<cutpoint]<-0
#' confusion_stat(mu, R84$outwork)
#'
#' @seealso \code{\link{glm}}
#' @author Rafael de Souza, ELTE  University,  and Joseph M. Hilbe, Arizona State University
#'
#' @references Hilbe, Joseph M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#'
#' Hilbe, Joseph M. (2009), Logistic Regression Models, Chapman & Hall/CRC.
#' @keywords models
#' @export
#'
confusion_stat<-function(pred=pred,obs=obs)
{
  CM<-table(pred,obs)
  table=ftable(addmargins(CM))
  stat<-confusionMatrix(pred,obs)
  return(list(matrix=table,statistics=c(stat$overall[1],stat$byClass[1],stat$byClass[2] )))
}
