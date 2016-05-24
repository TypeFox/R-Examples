
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Jan 08 10:59:20 EST 2014 -0500 (Week 01)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


##' A class inherited directly from envRefClass
##'
##' 
##' @title A class inherited directly from envRefClass
##' @field data list, a base::list
##' @author Xiaobei Zhao
##' 
List <- setRefClass(
  'List',
  list(
    data='list'
    ),
  methods=list(
    initialize=function(...){
      data <<- list(...)
    },
    get_names=function(){
      base::names(data)
    },
    get_data=function(){
      data
    },
    as.list=function(){
      data
    },
    has_key=function(x){
      .names <- names(data)
      if (!length(.names)){
        ret <- FALSE
      } else {
        ret <- (x %in% .names)
      }
      ret
    },
    getone=function(name,index,type,default,msg=TRUE){
      if (missing(name) & missing(index)){
        stop("List | getone | missing: name and index")
      }
      if (!missing(name)) {
        .value <- ..(getone)(data,name,msg=msg)
      } else {
        .value <- NULL
      }
      ## logme(.value,"List | getone.1",log='DEBUG') 
      
      if (is.null(.value)) {
        if (!missing(index)) {
          .value <- ..(getone)(data,index,msg=msg)
        }
      }
      ## logme(.value,"List | getone.2",log='DEBUG')
      
      if (is.null(.value) & !missing(default)){
        .value <- default
      }
      ## logme(.value,"List | getone.3",log='DEBUG')
      
      if (missing(type)){
        type <- 'ANY'
      }
      ret <- ValueParser$new(value=.value,type=type)$get_value()   
      ## logme(ret,"List | getone",log='DEBUG')

      return(ret)
    },
    popone=function(x, warn=TRUE, error=TRUE){
      'Pop the one at the given index/position (or name) in the list, and return it. If no index is specified, obj$popone() removes and returns the last one in the list.'
      .envir <- as.environment(.self)
      ..(popone)(data,x,warn=warn,error=error,envir=.envir)
    },
    popmany=function(x){
      'Pop many by indexes.'
      .envir <- as.environment(.self)
      ..(popmany)(data,x,envir=.envir)
    },
    removeone=function(x){
      'Remove the first matched element whose value is x. Display an error if it does not exist.'
      .envir <- as.environment(.self)
      ..(removeone)(data,x,envir=.envir)
    }
    )
  )




