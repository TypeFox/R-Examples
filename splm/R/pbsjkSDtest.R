`pbsjkSDtest` <-
function(formula, data, w, index=NULL, ...) {

  ## performs Baltagi, Song, Jung and Koh C.1 test
  ## for spatial dependence conditional on RE and serial corr.
  ## Giovanni Millo, Trieste, this version 2: 19/03/2009

  ## new interface for operation with pbsjktest2.R--> and
  ## -->spreml.R

  ## NB! not the same numbers as in the nlme-based version!!
  ## ...but these look like fitting better with a Wald test
  ## on semsrre

  ## depends: spreml()>ssrREmod(), fdHess{nlme} for numerical hessians

  #require(nlme) # not needed any more

  ## reorder data if needed
  if(!is.null(index)) {
    #require(plm)
    data <- plm.data(data, index)
    }

  gindex <- data[,1]
  tindex <- data[,2]

  ## for our purpose data has to be (re)ordered
  ## by time, then group
  data <- data[order(tindex, gindex),]

  ## est. MLE AR-RE model
  mymod <- spreml(formula=formula, data=data, w=w,
                  index=index, lag=FALSE, errors="srre", ...)

  ## def. 'trace' function
  tr<-function(x) sum(diag(x))

  nt. <- dim(data)[[1]]
  n. <- length(unique(gindex))
  t. <- length(unique(tindex))
  Jt<-matrix(1,ncol=t.,nrow=t.)

  ## make W matrix from listw object, if needed
  if("listw" %in% class(w)) w<-listw2mat(w) #if(require(spdep)) w<-listw2mat(w)

  ## retrieve restricted model's residuals ### substitute by a direct extraction
  X<-model.matrix(formula, data)
  y<-model.response(model.frame(formula,data))
  beta0<-mymod$coefficients
  u.hat<-as.numeric(y-X%*%beta0)

  ## retrieve AR(1) coefficient
  rho<-mymod$errcomp["psi"] # notice change in parm names in spreml

  ## henceforth notation as in Baltagi, Song, Jung, Koh (JE 2007)

  b<-tr(w%*%w+crossprod(w))

  d2<-(1+rho)/(1-rho)+t.-1

  ## retrieve variance components sigma.e and sigma.mu from lme object
  sigma2tot<-sum(u.hat^2)/length(u.hat)
  phi<- mymod$errcomp["phi"]  ## sigmatot^2=sigma.e^2*(phi+1)
  sigma2.e<-sigma2tot/(phi+1)
  sigma2.u<-sigma2tot-sigma2.e

  sigma.e<-sqrt(sigma2.e)
  sigma.u<-sqrt(sigma2.u)

  c. <- (sigma.e^2 * sigma.u^2) / (d2*(1-rho)^2*sigma.u^2+sigma.e^2)

  g. <- (1-rho)/sigma.e^2 * ( 2 + (t.-2)*(1-rho) )

  V1<-matrix(ncol=t.,nrow=t.)
  for(i in 1:t.) V1[i,]<-rho^abs(1:t.-i)

  V <- sigma.e^2 * (1/(1-rho^2)) * V1
  iV<-solve(V)
  VJt <- solve(V,Jt)

  bluestar<-(iV - 2*c. * VJt %*% iV +
             c.^2 * VJt%*%VJt %*% iV)

  bluespade<-kronecker(bluestar,(t(w)+w))

  Dhat <- 1/2 * crossprod(u.hat, bluespade) %*% u.hat

  LMl.rm <- (Dhat^2) / (b*(t. - 2*c.*g. + c.^2*g.^2))

  df.<-1
  pval <- pchisq(LMl.rm,df=df.,lower.tail=FALSE)

  names(LMl.rm)="LM"
  names(df.)<-"df"

  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = LMl.rm, parameter = df.,
               method = "Baltagi, Song, Jung and Koh C.1 conditional test",
               alternative = "spatial dependence in error terms, sub RE and serial corr.",
               p.value = pval,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}

