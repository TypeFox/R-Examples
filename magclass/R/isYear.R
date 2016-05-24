isYear<-function(x, with_y=TRUE)
{
  if(!is.vector(x)){stop("Year Object is no Vector")}
  return_vector<-rep(TRUE,length(x))
  for (i in 1:length(x))
    {
      if (with_y==FALSE)
        {
          if (nchar(x[i])!=4)                 {return_vector[i]=FALSE}
          if (grepl("[0-9]{4}",x[i])==FALSE)  {return_vector[i]=FALSE}
        } else
        { 
          if (nchar(x[i])!=5)                 {return_vector[i]=FALSE}
          if (grepl("y[0-9]{4}",x[i])==FALSE) {return_vector[i]=FALSE}
        }      
    }
  return(return_vector)
}
