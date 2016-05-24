IBMD <-
function(x,conf.levels = 0.95){
 d = dim(x)
if(is.null(d)) stop("argument is not a matrix or a dataframe")
if(conf.levels>1 || conf.levels<0){stop("magnitude of the returned confidence interval must be a single number between 0 and 1")}
if (is.na(table(sapply(x,is.numeric))["FALSE"])==FALSE){ stop("argument is not numeric")}

nsm = 1000
w = matrix(0,1,nsm)
ii = 1:d[1]

for (sm in 1:nsm){
 s = 0
 co = 0
 for (i in ii){
 for (j in 1:(d[2]-1)){
 if(is.na(x[i,j]) == FALSE){
 for (k in (j+1):d[2]){
 if (is.na(x[i,k]) == FALSE){
 diff <- x[i,j]-x[i,k]
 v <- ifelse(diff == 0, 0,log2(1+abs(diff)/max(x[i,k],x[i,j])))
 s <- s + sum(v, na.rm = FALSE)
 co = co+1
  }
  }
  }
  }
 }
w[1,sm] = s/co
ii = round(runif(d[1],1,d[1]))
 }
 q1 = (1-conf.levels)/2
 q2 = 1-q1
 q = quantile(w,c(q1,q2))
 dfo = data.frame(w[1,1],q[1],q[2],row.names=NULL)
 colnames(dfo) = c('value',paste('(',q1*100,"%"),paste('-',q2*100,"%",')'))
 return(list(subjects=d[1],observers=d[2],IBMD = dfo))
}
