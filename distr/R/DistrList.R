# generating function for class 'DistrList'


DistrList <- function(..., Dlist){
    ldots <- list(...)
    if(!missing(Dlist)){
        Dlist.L <- as(Dlist, "list")
        if(!is(try(do.call(DistrList,args=Dlist.L),silent=TRUE),"try-error"))
            ldots <- c(ldots, Dlist.L)
       }
    l <- length(ldots)
    if(l==0) return(NULL) else{
       do.call(new, args=list("DistrList", c(ldots)))
    }
}


# coerce to "DistrList"
setAs(from = "Distribution", to = "DistrList", 
    def = function(from){ new("DistrList", list(from)) })


# generating function for class 'UnivarDistrList'
UnivarDistrList <- function(..., Dlist){
    ldots <- list(...)
    if(!missing(Dlist)){
        Dlist.L <- as(Dlist, "list")
        if(!is(try(do.call(UnivarDistrList,args=Dlist.L),silent=TRUE),"try-error"))
            ldots <- c(ldots, Dlist.L)
       }
    l <- length(ldots)
    if(l==0) return(NULL) else{
       do.call(new, args=list("UnivarDistrList", c(ldots)))
    }
}


# coerce to "UnivarDistrList"
setAs(from = "UnivariateDistribution", to = "UnivarDistrList", 
    def = function(from){ new("UnivarDistrList", list(from)) })

