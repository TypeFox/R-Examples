
# coerce to "DistrList"
setAs(from = "MultivariateDistribution", to = "MVDistrList",
    def = function(from){ do.call("MultivarDistrList", args = list(from)) })

# generating function for class 'MultivarDistrList'
MultivarDistrList <- function(..., Dlist){
    ldots <- list(...)
    if(!missing(Dlist)){
        Dlist.L <- as(Dlist, "list")
        if(!is(try(do.call(MultivarDistrList,args=Dlist.L),silent=TRUE),"try-error"))
            ldots <- c(ldots, Dlist.L)
       }
    l <- length(ldots)
    if(l==0) return(NULL) else{
   #    print(dim(ldots[[1]]))
       if(dim(ldots[[1]])==1) return(do.call(UnivarDistrList, args = ldots))
       else do.call(new, args=list("MVDistrList", c(ldots)))
    }
}


# coerce to "MultivarDistrList"
setAs(from = "MultivariateDistribution", to = "MultivarDistrList",
    def = function(from){ do.call("MultivarDistrList", args = list(from)) })

setMethod("dim", "MultivarDistrList", function(x) x[[1]]@img@dimension)
setMethod("dimension", "MultivarDistrList", function(object) object[[1]]@img@dimension)

