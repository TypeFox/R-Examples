`pbsjkJtest` <-
function(formula, data, w, index=NULL, ...) {

  ## performs Baltagi, Song, Jung and Koh J(oint) test
  ## for RE, serial correlation and spatial corr.
  ## Giovanni Millo, Trieste, this version: 19/03/2009

  ## for our purpose data has to be (re)ordered
  ## by time, then group (but this is cared for just below)

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

  ## def. 'trace' function
  tr<-function(x) sum(diag(x))

  ## def. 'matrix square' function
  msq<-function(x) x%*%x

  ## make W matrix from listw object, if needed
  if("listw" %in% class(w)) w <- listw2mat(w) #if(require(spdep)) w<-listw2mat(w)

  ## retrieve restricted model's (OLS) residuals (ordered!)
  X<-model.matrix(formula, data)
  y<-model.response(model.frame(formula,data))
  beta0<-lm(y~X-1)$coef
  u.hat<-y-X%*%beta0

  ## calc. data numerosities (do it better)
  nt.<- length(y)
  n.<- dim(w)[[1]]
  t.<-nt./n.

  ## henceforth notation as in Baltagi, Song, Jung, Koh (JE 2007)
  Jt<-matrix(1,ncol=t.,nrow=t.)
  In<-diag(1,n.)
  It<-diag(1,t.)
  G<-matrix(0,ncol=t.,nrow=t.)
  for(i in 2:t.) {
    G[i-1,i]<-1
    G[i,i-1]<-1
    }

  ## NB do all this without Kronecker prods.!
  A <- (crossprod(u.hat, kronecker(Jt, In)) %*% u.hat)/crossprod(u.hat)-1
  F <- 1/2 * (crossprod(u.hat, kronecker(G, In)) %*% u.hat)/crossprod(u.hat)
  H <- 1/2 * (crossprod(u.hat, kronecker(It, (t(w)+w))) %*% u.hat)/crossprod(u.hat)
  b <- tr(msq(w+t(w)))/2

  LMj <- n.*t.^2 / (2*(t.-1)*(t.-2)) * (A^2 - 4*A*F + 2*t.*F^2) + (n.^2*t.)/b*H^2

  df.<-3
  pval <- pchisq(LMj,df=df.,lower.tail=FALSE)

  names(LMj)="LM"
  names(df.)<-"df"

  ##(insert usual htest features)
  dname <- deparse(formula)
  RVAL <- list(statistic = LMj, parameter = df.,
               method = "Baltagi, Song, Jung and Koh joint test (J)",
               alternative = "random effects or serial corr. or spatial dependence in error terms",
               p.value = pval,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}

