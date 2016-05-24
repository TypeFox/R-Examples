soft.threshold <- structure(function # Soft-thresholding
### Apply the soft-threshold function to a vector.
(x,
### A vector of numeric data.
 lambda=1
### The largest absolute value that will be mapped to zero.
 ){
  stopifnot(lambda>=0)
  ifelse(abs(x)<lambda,0,x-sign(x)*lambda)
### The vector of observations after applying 
### the soft-thresholding function.
},ex=function(){
  x <- seq(-5,5,l=50)
  y <- soft.threshold(x)
  plot(x,y)
})

.result <- 
 list(soft.threshold = list(definition = "soft.threshold <- structure(function # Soft-thresholding\n### Apply the soft-threshold function to a vector.\n(x,\n### A vector of numeric data.\n lambda=1\n### The largest absolute value that will be mapped to zero.\n ){\n  stopifnot(lambda>=0)\n  ifelse(abs(x)<lambda,0,x-sign(x)*lambda)\n### The vector of observations after applying \n### the soft-thresholding function.\n},ex=function(){\n  x <- seq(-5,5,l=50)\n  y <- soft.threshold(x)\n  plot(x,y)\n})",
     title = "Soft-thresholding",
        format="",
     description = "Apply the soft-threshold function to a vector.",  
     `item{x}` = "A vector of numeric data.",
     `item{lambda}`= "The largest absolute value that will be mapped to zero.",
     examples = "\nx <- seq(-5,5,l=50)\ny <- soft.threshold(x)\nplot(x,y)\n",  
     value = "The vector of observations after applying \nthe soft-thresholding function.")) 
