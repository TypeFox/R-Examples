get.weight.function <-
function(type){ #To get the weight function from the chracter string specified in type. 
  if(is.character(type)) { 
  if(!type%in%c("Trio","Gaussian", "Uniform", "Triweight", "Triangle", "Epanechnikov", "Quartic"))
      stop("type  not in known types.\n")
  f<-switch(EXPR=type, 
        "Trio"=  function(x)        (1 - abs(x)^3)^3 * (abs(x) < 1),       
        "Gaussian" =function(x)     exp(-x^2/2)/sqrt(2 * pi), 
        "Uniform"=function(x)       (abs(x) <= 1)/2,
        "Triweight"=  function(x)   (1 - x^2)^3 * 35/32 * (abs(x) <= 1), 
        "Triangle"= function(x)     (1 - abs(x))*(abs(x) <= 1),   
        "Epanechnikov"= function(x) (1 - x^2) * 3/4 * (abs(x) <= 1),
        "Quartic"=function(x)       (1 - x^2)^2 * 15/16 * (abs(x) <= 1)) 
  attr(f,"name")=type 
  return(f)
  }else if(is.function(type)) return(type)  else stop("Wrong type specified.")
}
