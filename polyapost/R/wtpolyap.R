#Given a vector of sampled  y values, ysamp, of size n
#from a population of size N this function uses the
#weighted Polya posterior to  generate one completed
#copy of the entire population. Each member of the sample
#is assigned a nonnegative weight. These weights are given
#in the vector wts which is of length n. With k=N-n it returns
#a vector of length N with ysamp being the first n elements.

wtpolyap<-function(ysamp,wts,k)
  {
    out<-cwpolya(ysamp,wts,k)
    return(out)
  }

