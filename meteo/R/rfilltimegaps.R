rfilltimegaps <- function(stfdf,
                          tunits="day",
                          attrname=1,
                        ...){
    
  nt=dim(stfdf)[2]
  time = as.POSIXlt(index(stfdf@time))
  times= seq(from=time[1], to=time[nt],  by = tunits)
  zcol = attrname
  
mat=matrix(ncol=length(times),nrow=dim(stfdf)[1] )

idx <- sapply(times, function(i){
  any(as.character(i) == as.character(time)) })

tt=1
for(ii in which(idx) ) {
  mat[,ii] = stfdf[,tt]@data[,zcol] 
  tt=tt+1
}

pb<-txtProgressBar(style = 3, max=length(mat[,1]))
xxx= lapply(1:length(mat[,1]), function(i) {
  x=spline(1:length(times) ,mat[i,] , xout=1:length(times) , ...)$y
  setTxtProgressBar(pb, i )
  x } ) # 
close(pb)

mat <- do.call(rbind,  xxx)

df <- reshape(as.data.frame(mat), varying=list(1:length(times)), v.names="spl", 
              direction="long", 
              times=times, ids=1:dim(stfdf)[1])

row.names(stfdf@sp) = 1:dim(stfdf)[1]

ret <- STFDF(stfdf@sp, 
             time= times, 
             data=data.frame(spl=df[,'spl']))

  
  return(ret)
}


