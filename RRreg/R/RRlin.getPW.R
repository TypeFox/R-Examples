# PWarray -> missclassification array for all RR variables and groups

getPWarray <- function(M, G, models, p.list, group, par2, Kukrep=1){
  # first: no groups, only one PW matrix
  # for first RR model separately
  PW <- getPW(models[1],p.list[[1]],group=1,par2[1], Kukrep=Kukrep)
  #  if more RR models: cronecker product for all RR models
  if (M>1){
    for (m in 2:M){
      PWtemp <- getPW(models[m],p.list[[m]],group=1,par2[m], Kukrep=Kukrep)
      PW <- kronecker(PW,PWtemp,make.dimnames = T)
    }
  }
  J <- nrow(PW)    # number of possible response patterns
  K <- ncol(PW)    # number of possible response patterns
  gcombs <- unique(group)
  PWarray <- array(dim=c(J,K,G),
                   dimnames=list("response"=rownames(PW),
                                 "true"=colnames(PW),
                                 "group"= apply(gcombs,1,paste,collapse="")))
  PWarray[,,1] <- PW
  
  if (G>1){  # if there are groups: fill array
    # loop for group combinations g
    for (g in 1:G){  
      PW <- getPW(models[1],p.list[[1]],group=gcombs[g,1],par2[1], Kukrep=Kukrep)
      # loop for models m 
      if (M>1){
        for (m in 2:M){
          PWtemp <- getPW(models[m],p.list[[m]],group=gcombs[g,m],par2[m], Kukrep=Kukrep)
          PW <- kronecker(PW,PWtemp,make.dimnames = T)
        }
      }
      PWarray[,,g] <- PW
    }
  }
  PWarray
}


#' Get Misclassification Matrices for RR Models
#' 
#' Given some randomization probabilities \code{p}, each RR design corresponds to a misclassification matrix PW.  This square matrix has entries defined as: \code{PW[i,j] = P(respond i | true state j)}. 
#' @param model one of the available models in the package \code{RRreg}: 
#' @param p randomization probability
#' @param group group index (1 or 2) for two-group designs such as  \code{"UQTunknown"} or \code{"SLD"}
#' @param par2 the second, estimated parameter in two-group designs (e.g., the unknown prevalence of the irrelevant question in \code{"UQTunknown"}, the t-parameter for truth in the \code{"SLD"})
#' @param Kukrep number of replications in Kuk's RR design (how many cards are drawn)
#' @details The method is used internally for estimation. Moreover, the method might be useful to check the exact definition of the RR designs.
#' 
#' Note that for two-group designs, the matrix dependends on a second parameter that is estimated from the data (e.g., the unknown prevalence of the unknown question in the unrelated question technique). Hence, the matrix itself is not constant, but an estimate of a random variable itself. 
#' @examples 
#' getPW(model = "Warner", p = 2/12)
#' getPW(model = "UQTknown", p = c(2/12, .3))
#' getPW(model = "UQTunknown", p = c(2/12, .10/12), group=2, par2=.4)
#' @references van den Hout, A., & Kooiman, P. (2006). Estimating the Linear Regression Model with Categorical Covariates Subject to Randomized Response. Computational Statistics & Data Analysis, 50(11), 3311–3323.

#' @export
getPW <- function (model, p, group = 1, par2=NULL, Kukrep=1){
  gr <- max(group,1)
  switch(model,
         "Warner" = PW <- matrix(c(p,       # true 0 -> 0 response
                                   1-p,          # 0 -> 1
                                   1-p,          # 1 -> 0
                                   p),           # 1 -> 1
                                 nrow=2,dimnames= list(0:1,0:1)) ,  
         "UQTknown" = PW <- matrix( c( 1-(1-p[1])*p[2],     #true 0 -> 0 response
                                       (1-p[1])*p[2],           # 0 -> 1
                                       (1-p[1])*(1-p[2]),       # 1 -> 0
                                       p[1] + (1-p[1])*p[2]),   # 1 -> 1
                                    nrow=2,dimnames= list(0:1,0:1)),   
         "UQTunknown" = PW <- matrix( c( 1-(1-p[gr])*par2,     #true 0 -> 0 response
                                         (1-p[gr])*par2,           # 0 -> 1
                                         (1-p[gr])*(1-par2),       # 1 -> 0
                                         p[gr] + (1-p[gr])*par2),      # 1 -> 1
                                  nrow=2,dimnames= list(0:1,0:1)), 
         "Mangat" = PW <- matrix(c(p,       # true 0 -> 0 response
                                   1-p,          # 0 -> 1
                                   0,            # 1 -> 0
                                   1),           # 1 -> 1
                                 nrow=2,dimnames= list(0:1,0:1)) , 
         "Kuk" = {
#             PW <- matrix(c(1-p[2],     # true 0 -> 0 response
#                                 p[2],            # 0 -> 1
#                                 1-p[1],          # 1 -> 0
#                                 p[1]),           # 1 -> 1
#                               nrow=2,dimnames= list(0:1,0:1)) 
#          
            ## Kukrep >1 : 
            PW <- cbind(dbinom(0:Kukrep,Kukrep,p[2]), dbinom(0:Kukrep,Kukrep,p[1]))
            dimnames(PW) <- list(0:Kukrep,0:1)
         },
         "FR" = {
              numcat <- length(p)
              PW <- matrix(rep(p,numcat),
                           nrow=numcat, 
                           dimnames= list(0:(numcat-1),0:(numcat-1)))
              for (i in 1:numcat){
                PW[i,i]<- 1-sum(p[-i])
              }
         },
         "Crosswise" = PW <- matrix(c(p,       # true 0 -> 0 response
                                      1-p,          # 0 -> 1
                                      1-p,          # 1 -> 0
                                      p),nrow=2,dimnames= list(0:1,0:1)),   # 1 -> 1
         "CDM" = PW <-  matrix(c(0,1,
                                 1-p[gr], p[gr], # non-user
                                 1, 0),# cheater
                               nrow=2, dimnames = list(0:1, c(1, 0, "cheater"))),
          # CDM possible??: 3 true underlying states but only 2 observed!
         "CDMsym" = PW <-  matrix(c(p[2*gr], 1-p[2*gr],
                                    1-p[2*gr-1], p[2*gr-1], # non-user
                                    1, 0),# cheater
                                  nrow=2, dimnames = list(0:1, c(1, 0, "cheater"))),
         "SLD" = PW <- matrix(c(p[gr],   # true 0 -> 0 response
                                1-p[gr],      # 0 -> 1
                                1-par2,          # 1 -> 0
                                par2),           # 1 -> 1
                              nrow=2,dimnames= list(0:1,0:1)),
         "custom" = {
              if (class(p) != "matrix" || nrow(p) != ncol(p) || any(colSums(p)!=1))
                stop("If the RR method 'custom' is used, a missclassification matrix 'p' 
                     must be provided, where p[i,j] gives the probability of responding
                     'i' (i-th row), while being in the true state 'j' (j-th column). 
                     Within a column, probabilities must sum up to one.")
              PW <- p  # own specification of missclassification matrix PW
              dimnames(PW) <- list(response=1:ncol(p)-1,true=1:ncol(p)-1)}
           )
  
  # for direct questioning: identity matrix, no missclassification
  # use loops to keep size and names of matrix
if (group==0){
  PW1 <- diag(nrow=nrow(PW),ncol=ncol(PW))
  dimnames(PW1) = dimnames(PW)
  PW <- PW1
}
  PW
}