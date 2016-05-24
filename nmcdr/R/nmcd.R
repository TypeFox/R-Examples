#_______________________________ Standardized C-M test statistic
tcs <- function(nk, n ,y)              
 {
  cvm <- CramerVonMisesTwoSamples(y[1:nk], y[(nk+1):n])
  cvm^0.5
 }

#________________________________ Initial screening
isp <- function(x, n = length(x), an = as.integer(log(n)^1.5/2))
 {
  cn <- 0
  test_score <- numeric(n)                    
  cand_cp_set <- numeric(n)                
  #### Calculate the Cramer-Von Mise statistic
  no_need = sapply(an:(n-an), function(j) test_score[j] <<- 
             tcs(an, 2*an, x[(j-an+1):(j+an)]))
  l<-0             
  for(j in an:(n-an))           ####  Find the local maximum
   {
     i = max(an, j-an+1)
     flag = 0
     while(i <= (j+an) & flag == 0)
      {
        if(test_score[i] > test_score[j]) flag <- 1
        i <- i + 1
      }
     if(flag == 0)            # find one possible candiate!
      {
        l <- l + 1
        cand_cp_set[l] <- j        
      }
   }
  tt<-list(ncp = l, cpp = cand_cp_set[1:l]+1, data = x)
  class(tt) <- c('isp','nmcd')
  return(tt)
 }

#____nonparamatric multiple change-point detection procedure
nmcd <- function(x, kmax, cpp, ncp = length(cpp), n){
  if(missing(n)) n = length(x)
  else{
  	if(n != length(x)) x = x[1:n]
  } 
  x = as.numeric(x)
  n = as.integer(n) 
  if(missing(cpp)) cpp = isp(x)$cpp
  cpp = as.integer(cpp)
  if(missing(kmax)) kmax = ncp
  else{
  	if(kmax > ncp) kmax = ncp
  }
  kmax = as.integer(kmax)
  ncp = as.integer(ncp)
  ke = 0L
  cpe = integer(kmax)
  gc()
  w = matrix(0,n,n)
  bic=numeric(1)
  r <- .Fortran('nmcdp', n = n, y = x, steps = kmax, ncp = ncp, cpp = cpp,
                 ke = ke, cpe = cpe, w = w, tmin = bic)
  gc()
  result = list(ncp = r$ke, cpp = r$cpe[2:(r$ke+1)], data = x, bic=r$tmin)
  class(result) <- 'nmcd'
  return(result)
}

print.nmcd <- function(x,...){
  a = x$ncp
  b = x$cpp
  t = list(ncp = a, cpp = b)
  print(t)
} 

plot.nmcd <- function(x,xlab = "position",ylab = "data",main = 'NMCD',...){
  data = x$data
  n = length(data)
  position = 1:n
  nt = x$ncp + 1
  yy = integer(nt)
  yy[1] = 1
  yy[2:nt] = x$cpp
  yy[nt+1] = n+1
  zz <- 0
  xx <- 0
  for(i in 1:nt)
   {
     zz[i] <- median(data[yy[i]:(yy[i+1]-1)])
     xx[yy[i]:(yy[i+1]-1)]<-zz[i]
   }  
  plot(position, data, lty=1, ylab=ylab,xlab=xlab, pch = 15, cex=0.3,
    cex.lab=1, cex.axis=1, main = main)
  lines(position,xx,type="l",col="red",lty=1,lwd=2.0)
  }

plot.isp <- function(x,xlab = "position",ylab = "data",main = 'ISP',...){
  data = x$data
  n = length(data)
  position = 1:n
  nt = x$ncp + 1
  yy = integer(nt)
  yy[1] = 1
  yy[2:nt] = x$cpp
  yy[nt+1] = n+1
  zz <- 0
  xx <- 0
  for(i in 1:nt)
   {
     zz[i] <- median(data[yy[i]:(yy[i+1]-1)])
     xx[yy[i]:(yy[i+1]-1)]<-zz[i]
   }  
  plot(position, data, lty=1, ylab=ylab,xlab=xlab, pch = 15, cex=0.3,
    cex.lab=1, cex.axis=1, main = main)
  lines(position,xx,type="l",col="green",lty=1,lwd=2.0)
  }
  
summary.nmcd<-function(object,...){
  data=object$data
  ncp=object$ncp
  cpp=object$cpp
  bic=object$bic
  if(class(object)[1]=='isp') bic=Inf
  cat('number of observations:',length(data),'\n')
  cat('estimated number of change-points:',ncp,'\n')
  cat('estimated change-points:',cpp,'\n')
  cat('minimal BIC:',bic,'\n')
  }