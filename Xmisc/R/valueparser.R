
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Feb 24 10:14:08 EST 2014 -0500 (Week 08)
## 
## 
## ************************************************************************


##' Parser for values
##'
##' 
##' @title Parser for values
##' @author Xiaobei Zhao
##' @field value ANY
##' @field type character
##' 
##' @examples
##' ValueParser$new(value="",type='character')$get_value()
##' ValueParser$new(type='character')$get_value()
##' 
ValueParser <- 
  setRefClass(
    'ValueParser',
    list(
      value='ANY',
      type='character'
      ),
    methods=list(
      initialize=function(value,type){
        if(missing(type)){
          .type <- 'ANY'
        } else {
          if(!length(type)){
            .type <- 'ANY'
          } else {
            .type <- type
          }
        }
        type <<- .type

        if (missing(value)) {
          .value <- NULL
        } else {
          .value <- value
        }
        ## logme(.value,'ValueParser | initialize',log='DEBUG')
        .self$value <- R5.value.parse(.value,type)
        ## logme(.self$value,'ValueParser | initialize',log='DEBUG')
      },
      get_value=function(){
        ## logme(.self$value,'ValueParser | get_value',log='DEBUG')
        return(.self$value)
      }
      )
    )


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

##' R5.value.default
##'
##' 
##' @title R5.value.default
##' @param type the type of an R object
##' @param default the default value
##' @return ANY
##' @author Xiaobei Zhao
##' @examples
##' require(Xmisc)
##' R5.value.default('character')
##' try(R5.value.default(NULL))
##' R5.value.default('environment')
##' 
##' R5.value.default('hclust')
##' R5.value.default('dendrogram')
##' R5.value.default('formula')
##' R5.value.default('lm')
##' 
R5.value.default <- function(type,default='list'){
  ## logme('START','R5.value.default','DEBUG') ##
  if(!is.character(type)){
    stop("R5.value.default | type must be character of length one.")
  }
  .envir <- new.env()
  .s <- sprintf("%.12f",Sys.Epoch())
  ret <- suppressWarnings(
    switch(
      type,
      'data.table'={
        if (check.packages(data.table)){
          tmp=as.data.table(
            setRefClass(sprintf('.TMP%s',.s),list(.value='matrix'),where=.envir)$new()$.value)
          return(tmp)          
        } else {
          stop('Install `data.table` for full functionality.')
        }
      },
      'expression'=setRefClass(sprintf('.TMP%s',.s),list(.value='character'),
        where=.envir)$new()$.value,
      'RData'=setRefClass(sprintf('.TMP%s',.s),list(.value='ANY'),
        where=.envir)$new()$.value,
      'environment'=emptyenv(), # empty environment
      ## ## 'ANY'=NULL,
      {
        tmp=tryCatch(
          {setRefClass(
            sprintf('.TMP%s',.s),list(.value=type),where=.envir)$new()$.value},
          error=function(e){
            ##cat(as.character(e));
            return(NULL)}
          )
        if (is.null(tmp)){
          msg <- sprintf('R5.value.default | Undefined type (`%s`).',type)
          if (length(default)){
            if (default=='list'){
              x <- list()
              class(x)=c(type,class(x))
              return(x)
            }
          }
          stop(msg)
        } else {
          return(tmp)
        }
      }
      )
    )
  rm(.envir)
  return(ret)
}


##' R5.value.parse
##'
##' 
##' @title R5.value.parse
##' @param value the value passed to the parameter 
##' @param type the type of an R object
##' @param default the default value
##' @return ANY
##' @author Xiaobei Zhao
##' @examples
##' R5.value.parse(NULL,'logical')
##' R5.value.parse(1,'logical')
##' 
##' R5.value.parse(NULL,'hclust')
##' R5.value.parse(NULL,'dendrogram')
##' R5.value.parse(NULL,'formula')
##' R5.value.parse(NULL,'lm')
##' 
##' R5.value.parse(NULL,'character')
##' R5.value.parse("",'character')
##' ## [1] ""
##' 
R5.value.parse <- function(value,type,default='list'){
  ## logme('START','R5.value.parse','DEBUG') ## 
  if (missing(type)) {
    type <- 'ANY'
  }
  if(!is.character(type)){
    stop("R5.value.parse | type must be character of length one.")
  }
  if(missing(value)){
    return(R5.value.default(type))
  }
  ## logme(value,'R5.value.parse','DEBUG') ## 
  
  if (is.null(value)){
    return(R5.value.default(type))
  }
  ## logme(value,'R5.value.parse','DEBUG') ## 
  if(is.function(value) & type!='function'){
    value <- eval(parse(text=sprintf("value()")))
  }
  ret <- suppressWarnings(
    switch(
      type,
      'data.frame'={
        if (is.data.frame(value)) value else as.data.frame(value)
      },
      'data.table'={
        if (check.packages(data.table)){
          if (is.data.table(value)) {
            value
          } else {
            as.data.table(value)
          }
        } else {
          stop('Install `data.table` for full functionality.')
        }
      },
      'expression'={
        eval(base::parse(text=value))
      },
      'RData'={
        .envir <- new.env()
        base::load(file=value,envir=.envir) #.variable will not be loaded
        .names <- ls(.envir)
        if (!length(.names)) { # empty
          .value <- NULL 
        } else if (length(.names)==1){
          .name <- .names[1]
          .value <- base::get(.name,envir=.envir)
          rm(.envir)
          .value
        } else {
          stop("ValueParser | must be sinlge-object RData !")
        }
      },
      logical={
        .value <- value
        if (!is.logical(.value)){
          if (any(c(0,1,"0","1") %in% .value)) {
            .value <- as.logical(as.numeric(.value))
          } else if (any( c("TRUE","FALSE") %in% toupper(.value) ) ) {
            .value <- as.logical(toupper(.value))
          } else {
            stop("ValueParser | invalid logical value")
          }
        }
        .value
      },
      character={
        if (is.character(value)) value else as.character(value)
      },
      factor={
        if (is.factor(value)) value else as.factor(value)
      },
      numeric={
        if (is.numeric(value)) value else as.numeric(value)
      },
      integer={
        if (is.integer(value)) value else as.integer(value)
      },
      matrix={
        if (is.matrix(value)) value else as.matrix(value)
      },
      list={
        if (is.list(value)) value else as.list(value)
      },
      "function"={
        if (is.function(value)) value else as.function(value)
      },
      environment={
        if(is.null(value)){
          .value <- emptyenv()
        } else {
          if (is.environment(value)) {
            .value <- value
          } else {
            .value <- as.environment(value)
          }
        }
        .value
      },
      'ANY'=value,
      {
        msg <- sprintf('R5.value.parse | Undefined type (`%s`).',type)
        if (length(default)){
          if (default=='list'){
            if (length(value)) {
              return(value)
            } else {
              return(R5.value.default(type))
            }
          }
          stop(msg)
        } else {
          stop(sprintf('ValueParser | .parse | Invalid type: %s',type))
        }
      }
      )
    )
  return(ret)
}
