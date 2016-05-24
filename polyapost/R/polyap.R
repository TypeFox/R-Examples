#Given a sample of y values, say ysamp,  of size n
#from a population of size N this function uses the
#Polya posterior to  generate one completed copy of
#the entire population. With k=N-n it returns a vector
#of length N with ysamp being the first n elements.

polyap<-function(ysamp,k)
  {
    out<-cwpolya(ysamp,rep(1,length(ysamp)),k)
    return(out)
  }
