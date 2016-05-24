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
#' @title  Display Hosmer-Lemeshow statistic and table of probabilities following logistic
#' regression using glm with binomial family.
#' @description Provides a Hosmer-Lemeshow statistic and table following logistic regression.
#' @aliases hlGOF.test
#' @usage hlGOF.test(observed, predicted, breaks = 15)
#' @importFrom stats pchisq
#' @format \describe{
#' \item{x}{
#' The function has three arguments: observed term, predicted values, # groups}
#' }
#' @details hlGOF.test is a post-estimation function for logistic regression, following the use
#' of glm(). Usage displays a table of observed vs predicted groups and an overall
#' H-L goodness-of-fit statistic. The test is originally from Hilbe (2009).
#' @param observed response variable
#' @param predicted predicted statistic
#' @param breaks breaks or groups
#' @return numeric
#' @note hlGOF.test must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user. My thanks to Prof. Robert LaBudde for the initial
#' version of this function.
#'@examples
#'library(MASS)
#'library(LOGIT)
#'data(medpar)
#'mylogit <- glm( died ~  los + white + hmo, family=binomial, data=medpar)
#'summary(mylogit)
#'medpar2 <- na.omit(medpar)
#'hlGOF.test(medpar2$died, predict(mylogit,medpar2, type="response"), breaks=12)
#'
#' @seealso \code{\link{glm}}
#' @author Joseph M. Hilbe, Arizona State University, Robert LaBudde,
#' Institute for Statisical Education (Statistics.com), provided initial code
#' for this function for Hilbe, Logistic Regression Models, text.
#'
#' @references Hilbe, J. M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#'
#' Hilbe, J. M. (2009), Logistic Regression Models, Chapman & Hall/CRC.
#' @keywords models
#' @export
#'
hlGOF.test <- function (observed, predicted, breaks=15) {
  #H-L GOF test for logistic regression
  #observed and predicted should not have missing values and match by index
  cat('\n', 'Hosmer-Lemeshow GOF test', '\n')
  ndata <- length(predicted)
  cuts <- c(round(.75*breaks), breaks, round(1.25*breaks))
  pvals <- rep(1,cuts[3]) #p-values
  for (nCuts in cuts) {
    ip <- order(predicted)
    iq <- iQuantile(predicted, nCuts) #indices for cuts
    iqInd<- iq$index
    cat('\n','For # Cuts =',nCuts,'  # Data =',ndata,'\n')
    cat('Cut  # Total #Patterns # Resp.    # Pred.  Mean Resp. Mean Pred.','\n')

    x2 <- 0
    ntot <- 0
    for (i in 1:nCuts) {
      if (i==1) {
        isubs <- ip[1:iqInd[2]]
      } else {
        isubs <- ip[(iqInd[i]+1):iqInd[i+1]]
      }
      nsubs <- length(isubs)
      ntot <- ntot + nsubs
      aobs <- mean(observed[isubs])
      mobs <- sum(observed[isubs])
      ncvp <- length(unique(predicted[isubs]))
      apred <- mean(predicted[isubs])
      mpred <- apred*nsubs
      x2 <- x2 + (mobs-mpred)^2/mpred + ((nsubs-mobs) - (nsubs-mpred))^2/(nsubs-mpred)
      cat(sprintf('%3d',i), sprintf('%8d', nsubs), sprintf('%8d', ncvp), sprintf('%8d', mobs),
        sprintf('%10.2f',mpred), sprintf('%8.5f', aobs), sprintf('%8.5f',apred), '\n')
    }
    cat('Total # Data:',ndata,' Total over cuts:',ntot,'\n')
    pvals[nCuts] <- pchisq(x2,nCuts-2,lower.tail=FALSE)
    cat('Chisq:', x2, '  d.f.:', sprintf('%d',nCuts-2), ' P-value:',
      sprintf('%8.5f', pvals[nCuts]),'\n')
  }
  cat('\n','Minimum P-value: ',sprintf('%8.5f',min(pvals)),'\n')
}
