"mlmmmemex"<-function(test=1:20,print.it=F)
{
  y<-matrix(c(1.03,NA,
	     1.54,477,
	     1.82,444,
	     NA,370,
	     1.31,403,
	     2.16,NA,
	     2.13,450,
	     NA,393,
	     1.59,394,
	     2.53,499,
	     2.33,482,
	     1.80,317,
	     2.09,499,
	     NA,411,
	     2.21,391,
	     2.82,396,
	     1.66,371,
	     2.30,418,
	     2.65,486,
	     2.18,NA,
	     1.42,395,
	     NA,325,
	     1.58,316,
	     1.49,311,
	     1.41,414,
	     1.65,313,
	     NA,309,
	     1.34,323,
	     0.18,315,
	     0.64,376,
	     0.76,308,
	     0.70,439),ncol=2,nrow=32,byrow=T)
  subj<-c(1,1,1,1,
	 2,2,2,2,
	 3,3,3,3,
	 4,4,4,4,
	 5,5,5,5,
	 6,6,6,6,
	 7,7,7,7,
	 8,8,8,8)
  occ<-c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
# assign("occ",occ,frame=1)
  assign("occ",occ,env=.GlobalEnv)
  pred<-cbind(int=rep(1,32),
	     dummy1=rep(c(1,0,0,0),8),
	     dummy2=rep(c(0,1,0,0),8),
	     dummy3=rep(c(0,0,1,0),8))
  xcol<-1:4
  zcol<-1
########################################################################
  res2<-mlmmm.em(y,subj,pred,xcol,zcol,maxits=200,eps=0.0001)
  print(res2)
  res2.bd<-mlmmmbd.em(y,subj,pred,xcol,zcol,maxits=200,eps=0.0001)
  print(res2.bd)
########################################################################
  result <- rep(T, 20)
  ans <- F
  compare <- function(old, new, name)
    {
      if(is.character(all.equal(old, new)[1])) {
	browser()
	cat(name, " F\n")
	return(F)
      }
      else {
	return(T)
      }
    }
  if(any(test == 01)) {
    if(print.it) cat("== Test 01, compare `beta' for `em.pan' ====\n")
    old <- res$beta
    new <- res2$beta
    result[01] <- compare(old, new, 01)
  }
  if(any(test == 02)) {
    if(print.it) cat("== Test 02, compare `sigma' for `em.pan' ====\n")
    old <- res$sigma
    new <- res2$sigma
    result[02] <- compare(old, new, 02)
  }
  if(any(test == 03)) {
    if(print.it) cat("== Test 03, compare `psi' for `em.pan' ====\n")
    old <- res$psi
    new <- res2$psi
    result[03] <- compare(old, new, 03)
  }
  if(any(test == 04)) {
    if(print.it) cat("== Test 04, compare `converged' for `em.pan' ====\n")
    old <- res$converged
    new <- res2$converged
    result[04] <- compare(old, new, 04)
  }
  if(any(test == 05)) {
    if(print.it) cat("== Test 05, compare `iter' for `em.pan' ====\n")
    old <- res$iter
    new <- res2$iter
    result[05] <- compare(old, new, 05)
  }
  if(any(test == 06)) {
    if(print.it) cat("== Test 06, compare `logll' for `em.pan' ====\n")
    old <- res$logll
    new <- res2$logll
    result[06] <- compare(old, new, 06)
  }
  if(any(test == 07)) {
    if(print.it) cat("== Test 07, compare `logoll' for `em.pan' ====\n")
    old <- res$logoll
    new <- res2$logoll
    result[07] <- compare(old, new, 07)
  }
#
  if(any(test == 11)) {
    if(print.it) cat("== Test 11, compare `beta' for `embd.pan' ====\n")
    old <- res.bd$beta
    new <- res2.bd$beta
    result[11] <- compare(old, new, 11)
  }
  if(any(test == 12)) {
    if(print.it) cat("== Test 12, compare `sigma' for `embd.pan' ====\n")
    old <- res.bd$sigma
    new <- res2.bd$sigma
    result[12] <- compare(old, new, 12)
  }
  if(any(test == 13)) {
    if(print.it) cat("== Test 13, compare `psi' for `embd.pan' ====\n")
    old <- res.bd$psi
    new <- res2.bd$psi
    result[13] <- compare(old, new, 13)
  }
  if(any(test == 14)) {
    if(print.it) cat("== Test 14, compare `converged' for `embd.pan' ====\n")
    old <- res.bd$converged
    new <- res2.bd$converged
    result[14] <- compare(old, new, 14)
  }
  if(any(test == 15)) {
    if(print.it) cat("== Test 15, compare `iter' for `embd.pan' ====\n")
    old <- res.bd$iter
    new <- res2.bd$iter
    result[15] <- compare(old, new, 15)
  }
  if(any(test == 16)) {
    if(print.it) cat("== Test 16, compare `logll' for `embd.pan' ====\n")
    old <- res.bd$logll
    new <- res2.bd$logll
    result[16] <- compare(old, new, 16)
  }
  if(any(test == 17)) {
    if(print.it) cat("== Test 17, compare `logoll' for `embd.pan' ====\n")
    old <- res.bd$logoll
    new <- res2.bd$logoll
    result[17] <- compare(old, new, 17)
  }
#
  ans <- all(result)
  if(ans) ans
  else result
}

library(mlmmm)
mlmmmemex(20,TRUE)
