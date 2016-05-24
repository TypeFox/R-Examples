
place_x_in_y<-function(x,y,expand=T){
  if(!all(getYears(x)%in%getYears(y))) {
    if (expand){
      print("x has years that dont exist in y. Expand y.")
      y<-add_columns(x = y,dim = 2.1,addnm = setdiff(getYears(x),getYears(y)))
      y<-magpiesort(y)
    }else{
      x<-x[,getYears(x)[getYears(x)%in%getYears(y)],]
    }
  }
  if(!all(getRegions(x)%in%getRegions(y))) {
    if (expand){
      print("x has regions that dont exist in y. Expand y.")
      y<-add_columns(x = y,dim = 2.1,addnm = setdiff(getYears(x),getYears(y)))
      y<-magpiesort(y)
    } else {
      x<-x[getRegions(x)[getRegions(x)%in%getRegions(y)],,]
    }
  }
  if(!all(getNames(x)%in%getNames(y))) {
    if (expand){
      stop("x has names that dont exist in y. Cannot handle this yet. Please improve me!")}
  } else {
    x<-x[,,getNames(x)[getNames(x)%in%getNames(y)]]
  }
  y[getRegions(x),getYears(x),getNames(x)]<-x
  return(y)
}
