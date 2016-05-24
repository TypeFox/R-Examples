eqcat <-
function(x){
#check for variable names: long,lat,z,time,magn1
	  a<-c("time",  "lat",   "long",  "z",     "magn1")
	  n<-5-length(intersect(a,names(x)))
if (n>0) {
	  print("WRONG EARTHQUAKE CATALOG DEFINITION")
	  print(c(n," wrong variable names"))
	  print("an object of class eqcat (earthquake catalog) must contains the five names")
	  print (a)
	  return(list(cat=x,ok=FALSE));
	  } 
else
	  {
	  class(x)<-c("eqcat","data.frame")
	  return(list(cat=x,ok=TRUE));
	  }
}
