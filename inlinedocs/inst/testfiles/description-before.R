### A multi-line comment before a data item goes into the description
### section of its Rd file. The items of this list are all norms,
### which are functions that take a vector of N elements and return a
### non-negative real number.
norms <-
  structure(list(l0=function(x)sum(x!=0),
                 l1=function(x)sum(abs(x)),
                 l2=function(x)sqrt(sum(x^2)),
                 linfty=function(x)max(abs(x))),
            ex=function(){
              m <- replicate(2,rnorm(10))
              cbind(m,sapply(norms,function(norm)apply(m,1,norm)))
            })
.result <- 
 list(norms = list(definition = "norms <-\n  structure(list(l0=function(x)sum(x!=0),\n                 l1=function(x)sum(abs(x)),\n                 l2=function(x)sqrt(sum(x^2)),\n                 linfty=function(x)max(abs(x))),\n            ex=function(){\n              m <- replicate(2,rnorm(10))\n              cbind(m,sapply(norms,function(norm)apply(m,1,norm)))\n            })",  
     description = "A multi-line comment before a data item goes into the description\nsection of its Rd file. The items of this list are all norms,\nwhich are functions that take a vector of N elements and return a\nnon-negative real number.",  
     format = "", title = "norms", examples = "\nm <- replicate(2,rnorm(10))\ncbind(m,sapply(norms,function(norm)apply(m,1,norm)))\n")) 
