spmHrf <-function(RT,p=c(6, 16, 1, 1, 6, 0, 32)){
  fMRI_T<-16
  dt<- RT/fMRI_T
  u <-c(0:(p[7]/dt)) - p[6]/dt
  hrf <- dgamma(u,p[1]/p[3],dt/p[3]) - dgamma(u,p[2]/p[4],dt/p[4])/p[1]
  hrf <-hrf[(0:(p[7]/RT))*fMRI_T + 1];
  hrf <- hrf/sum(hrf)
  return(list(hrf=hrf,p=p))
}
