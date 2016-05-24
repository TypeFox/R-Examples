denmean<-function(volume,obslkm,n)
{
#Calculates the value of the estimate at rectangle

#volum is >0
#obslkm is lkm-vector, pointers to the rows of x
#n is the sample size

#returns non-negative real number

#Value of the estimate is the number of observations
#in rec divided by the total number of obs, 
#and divided by the volume of the rectangle
#value=n_rec/(n*volume(rec))

ans<-obslkm/(n*volume)
return(ans)
}
