`polylist2list` <-
function(data)
{
if(class(data)[1]!="polylist")
stop("Must be a polylist object")

n<-length(data[])
contours<-rbind(NA,NA,NA,data[[1]][,])

 for(i in 1:(n-1))
  {
   contours=rbind(contours,NA,NA,NA,data[[i+1]][,])   
  }

contours=rbind(contours,NA,NA,NA)
return(contours)

}

