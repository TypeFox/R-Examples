`pbsjkARtest` <-
function(formula, data, w, index=NULL, ...) {

  ## performs Baltagi, Song, Jung and Koh C.2 test
  ## for serial correlation conditional on RE and spatial corr.
  ## Giovanni Millo, Trieste, this version compatible with spreml()
  ## 19/03/2009

  ## depends: spreml()>semREmod(), fdHess{nlme} for numerical hessians

  #require(nlme) #not needed any more


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

  ## est. MLE SEM-RE model
  mymod <- spreml(formula=formula, data=data, w=w,
                  index=index, lag=FALSE, errors="semre", ...)

  nt. <- dim(data)[[1]]
  n. <- length(unique(gindex))
  t. <- length(unique(tindex))

  ## def. 'trace' function
  tr<-function(x) sum(diag(x))

  ## def. 'matrix square' function
  msq<-function(x) x%*%x

  ## make W matrix from listw object, if needed
  if("listw" %in% class(w)) w<-listw2mat(w) #if(require(spdep)) w<-listw2mat(w)

  ## retrieve restricted model's residuals ### substitute by a direct extraction
  X<-model.matrix(formula, data)
  y<-model.response(model.frame(formula,data))
  beta0<-mymod$coefficients
  u.hat<-as.numeric(y-X%*%beta0)

  ## retrieve SEM coefficient from model coef
  lambda <- mymod$errcomp["rho"]

  ## retrieve variance components sigma.e and sigma.mu from lme object
  eta <- mymod$errcomp["phi"]  # pay attention to this renaming
  sigma2tot <- as.numeric(crossprod(u.hat)/nt.)
  sigma2e <- as.numeric(sigma2tot/(eta+1))
  sigma2mu <- as.numeric(eta*sigma2e)

  ## henceforth notation as in Baltagi, Song, Jung, Koh (JE 2007)

  JaT<-matrix(1,nrow=t.,ncol=t.)/t.
  It<-diag(1,t.)
  Et<-It-JaT


  B<-diag(1,n.)-lambda*w
  BB<-crossprod(B)
  BB.1 <- solve(BB)

  wBBw<-crossprod(w,B)+crossprod(B,w)

  Z0 <- solve( t. * sigma2mu * diag(1,n.) + sigma2e * BB.1 )

  G<-matrix(0,ncol=t.,nrow=t.)
  for(i in 2:t.) {
    G[i-1,i]<-1
    G[i,i-1]<-1
    }

  EGE <- Et%*%G%*%Et
  JGE <- JaT%*%G%*%Et
  EGJ <- Et%*%G%*%JaT
  JGJ <- JaT%*%G%*%JaT


  redspade <- 1/sigma2e^2*kronecker(EGE,BB) + 1/sigma2e*kronecker(JGE,Z0) +
              1/sigma2e*kronecker(EGJ,Z0) + kronecker(JGJ,Z0%*%BB.1%*%Z0)

  Dhat <- -(t.-1)/t. * (sigma2e * tr(Z0%*%BB.1) -n.) +
          1/2 * sigma2e * crossprod(u.hat, redspade) %*% u.hat

  ## information matrix:
  d1<-tr( msq(Z0%*%BB.1) )
  d2<-tr(Z0%*%BB.1%*%Z0)
  d3<-tr( wBBw%*%BB.1 )
  d4<-tr( Z0 %*% BB.1 %*% wBBw %*% BB.1 %*% Z0 %*% BB.1 )
  d5<-tr( Z0 %*% BB.1 %*% wBBw %*% BB.1 %*% Z0 )
  d6<-tr( msq( wBBw %*% BB.1 ) )
  d7<-tr( msq( Z0 %*% BB.1 %*% wBBw %*% BB.1 ) )

  j11<-(n.*(t.-1)/sigma2e^2 + d1)/2
  j12<-t./2*d2
  j13<-(t.-1)/t.*(sigma2e*d1-n./sigma2e)
  j14<-((t.-1)/sigma2e*d3 + sigma2e*d4)/2
  j22<-t.^2/2*tr(msq(Z0))
  j23<-(t.-1)*sigma2e*d2
  j24<-t./2*sigma2e*d5
  j33<-n./t.^2 * (t.^3-3*t.^2+2*t.+2) + (2*(t.-1)^2*sigma2e^2)/t.^2*d1
  j34<-(t.-1)/t. * (sigma2e^2*d4 - d3)
  j44<-((t.-1)*d6 + sigma2e^2*d7)/2

  Jtheta<-matrix(ncol=4,nrow=4)
  Jtheta[1,]<-c(j11,j12,j13,j14)
  Jtheta[2,]<-c(j12,j22,j23,j24)
  Jtheta[3,]<-c(j13,j23,j33,j34)
  Jtheta[4,]<-c(j14,j24,j34,j44)

  J33.1<-solve(Jtheta)[3,3]

  LMr.lm <- (Dhat^2) * J33.1

  df.<-1
  pval <- pchisq(LMr.lm,df=df.,lower.tail=FALSE)

  names(LMr.lm)="LM"
  names(df.)<-"df"

  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = LMr.lm, parameter = df.,
               method = "Baltagi, Song, Jung and Koh C.2 conditional test",
               alternative = "serial corr. in error terms, sub RE and spatial dependence",
               p.value = pval,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}

