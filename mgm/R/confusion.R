
confusion <- function(tg, eg) {
  

  #check on input
  check <- apply(tg, 1:2, function(x) { x %in% c(0,1)})  
  if(sum(check)!=length(tg)){
    stop("Only zeros and ones allowed in the adjacency matrix!")
  }
  
  
  diag(tg) <- diag(eg) <- NA
  m_eval <- cbind(tg[which(!is.na(tg))], eg[which(!is.na(eg))])
  
  ev <- t(apply(m_eval, 1, function(x) {
    acc <- x[1]==x[2]
    if(x[1]==1 & x[2]==1) se<-1 else {se<-0}
    if(x[1]==0 & x[2]==0) sp<-1 else {sp<-0}
    return(cbind(acc, se, sp))
  }))
  
  m_eval2 <- cbind(m_eval, ev)
  
  #precision
  m_eval_top <- m_eval[which(m_eval[,2]==1),] #m_eval - test outcome positive (top)
  pre <- mean(m_eval_top[,1]==m_eval_top[,2]) #precision
  
  a <- mean(m_eval2[,3]) #accuracy
  se <- sum(m_eval2[,4])/sum(m_eval[,1]) #sensitivity
  sp <- sum(m_eval2[,5])/(nrow(m_eval)-sum(m_eval[,1])) #specifity  
  out.list <- list("accuracy"=a, "sensitivity"=se, "specificity"=sp, "precision"=pre)
  
  return(out.list)
}




