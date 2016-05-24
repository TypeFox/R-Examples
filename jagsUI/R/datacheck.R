
data.check <- function(x,name,verbose=TRUE){
  
  test = FALSE
  
  if(is.data.frame(x)){
    if(!is.null(dim(x))){
      if(verbose){cat('\nConverting data frame \'',name,'\' to matrix.\n',sep="")}
      x = as.matrix(x)
    } else {
      if(verbose){cat('\nConverting data frame',name,'to vector.\n')}
      x = as.vector(x)}
  }
  
  
  if (is.numeric(x)&&is.matrix(x)){
    if(1%in%dim(x)){
      if(verbose){cat('\nConverting 1-column matrix \'',name,'\' to vector\n',sep="")}
      x = as.vector(x)
    }
    test = TRUE
  }
  
  if(is.numeric(x)&&is.array(x)&&!test){
    test = TRUE
  }
  
  if (is.numeric(x)&&is.vector(x)&&!test){
    test = TRUE
  }
  
  if(test){
    return(x)
  } else{return('error')}
  
}