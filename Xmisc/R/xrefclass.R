
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Fri Jan 10 13:36:51 EST 2014 -0500 (Week 01)
## 
## 
## Reference: 
## 
## ************************************************************************



##' Extended Reference Class
##'
##' 
##' @title Extended Reference Class
##' @description
##' The Extended Reference Class (xRefClass) inherits directly from \code{envRefClass}. Listed are some of its key features: 
##' \itemize{
##' \item Method \code{initialize} passes arguments by position or name. 
##' \item Method \code{copy2}, a modified version of \code{copy}, 
##'   is tolerant to \code{activeBindingFunction} as fields. 
##' \item Method \code{update} updates a class instance's methods
##'   according to any update of the class.
##' }
##' @author Xiaobei Zhao
##' @field .index named_numeric indexes of args
##' @field .default named_list default values of args
##' @field .meta named_list additional args (meta information)
##' @field .envir environment. Default: as.environment(.self)
##' @field .tmp.list list for temporary storage
##' @field .out.list list for outputting
##' @examples
##' \dontrun{
##' MyClass <- 
##'   setRefClass(
##'     "MyClass",
##'     list(
##'       x="numeric",
##'       y="numeric",
##'       z=function(){x+y}
##'       ),
##'     contains="xRefClass",
##'     methods=list(
##'       initialize=function(...){
##'         .idx <- c(x=1,y=2)
##'         callSuper(...,.index=.idx)
##'       },
##'       printme=function(){
##'         cat('Hello World!','\n')
##'       }
##'       )
##'     )
##'
##' ## Method initialize - pass by position
##' obj <- MyClass$new(1,2)
##' obj$x
##' obj$y
##'
##' ## Method initialize - pass by name
##' obj <- MyClass$new(y=2)
##' obj$x
##' obj$y
##'
##' ## Method copy
##' ## obj <- MyClass$new(1,2)
##' ## obk <- obj$copy()    # Fail!
##' ## ## Error in (function ()  : unused argument (quote("myclass"))
##'
##' ## Method copy2
##' obj <- MyClass$new(1,2) # No such error!
##' obk <- obj$copy2()
##' obk$z
##'
##' ## Method update
##' obj <- MyClass$new()
##' obj$printme()
##' MyClass <- # To modify one of the original functions
##'   setRefClass(
##'     "MyClass",
##'     list(
##'       x="numeric",
##'       y="numeric",
##'       z=function(){x+y}
##'       ),
##'     contains="xRefClass",
##'     methods=list(
##'       initialize=function(...){
##'         .idx <- c(x=1,y=2)
##'         callSuper(...,.index=.idx)
##'       },
##'       printme=function(){ # This function is modified
##'         cat('Hello R!','\n')
##'       }
##'       )
##'     )
##' obj$printme() # The function is yet not modified
##' ## Hello World!
##' obj$update("printme") # update the function
##' obj$printme() # The function is modified
##' ## Hello R!
##' }
##' @seealso \code{methods::ReferenceClasses}
##' @exportClass xRefClass
xRefClass <- 
  setRefClass(
    'xRefClass',
    list(
      .index='numeric',
      .default='list',
      .meta='list',
      .envir='environment',
      .lib=function(){
        tolower(as.character(.self$get_def()@className))
      },
      .class=function(){
        as.character(class(.self))[1]
        ## tolower(as.character(class(.self)))[1]
      },
      .package=function(){
        utils::packageName()
      },
      .tmp.list='list',
      .out.list='list'
      ),
    methods=list(
      initialize=function(...){
        ''
        ## logme('START','xRefClass | initialize','DEBUG') ##
        
        .args <- List$new(...)
        ## logme(.args,'xRefClass | initialize','DEBUG') ##
        
        .field.classes <- get_fieldclasses()
        .names <- names(.field.classes)
        if (".index" %in% .names){
          .field.classes <- c(
            list(.index=.field.classes[[".index"]]),
            .field.classes[-which(.names %in% ".index")]
            )
        }
        .fields <- names(.field.classes)
        ## logme(.field.classes) ##
        ## logme(.field) ##
        .idx.to.pop <- c() 
        for (i in seq_along(.fields)) {
          .field <- .fields[i]
          .fclass <- .field.classes[[.field]]
          if(!is.activeBindingFunction(.fclass)){
            if (.field %in% names(.args)){
              .idx <- which(names(.args)==.field)[1]
            } else {
              if (.field %in% names(.index)){
                .idx <- .index[.field]
              } else {
                .idx <- numeric()
              }
            }
            if (.field %in% names(.default)){
              .dflt <- .default[[.field]]
            } else {
              .dflt <- NULL
            }
            ## if (.field=='default'){ ##
            ##   ## logme(.index) ##
            ##   logme(.field,"xRefClass | initialize",'DEBUG') ##
            ##   logme(.idx,"xRefClass | initialize",'DEBUG') ##
            ##   logme(.fclass,"xRefClass | initialize",'DEBUG') ##
            ##   logme(.dflt,"xRefClass | initialize",'DEBUG') ##
            ## }
            .res <- .args$getone(.field,.idx,.fclass,.dflt,msg=FALSE)
            base::assign(.field, .res, envir=as.environment(.self))
            
            ## if(.field=='.envir'){ ##
            ##   logme(.envir,'xRefClass | initialize')
            ## }
            .idx0 <- valid.arg.index(.args$as.list(),.field,safe=TRUE)
            if(!length(.idx0)){
              .idx0 <- valid.arg.index(.args$as.list(),.idx,safe=TRUE)
            }
            if(length(.idx0)){
              .idx.to.pop <- c(.idx.to.pop,.idx0)
            }
          }
        }
        ## logme(.idx.to.pop) ##
        if (length(.idx.to.pop)){
          .args$popmany(.idx.to.pop)
        }
        ## save nonrequired args
        .meta <<- .args$as.list()

        ## 
        if (identical(.envir, emptyenv())){
          .envir <<- as.environment(.self) # for local evaluation
          ## NOT globalenv()
        }
        
        ## logme('END','xRefClass | initialize','DEBUG') ##
      },
      setme=function(){        
      },
      get_envir=function(){
        .envir
      },
      get_def=function(){
        ## logme('START','xRefClass | get_def','DEBUG') ##
        .def <- .refClassDef
        ## logme('END','xRefClass | get_def','DEBUG') ##
        return(.def)
      },
      get_fieldclasses=function(){
        ## logme('START','xRefClass | get_fieldclasses','DEBUG') ##
        .def <- get_def()
        .field.classes <- .def@fieldClasses
        ## logme('END','xRefClass | get_fieldclasses','DEBUG') ##
        return(.field.classes)
      },
      copy2=function(shallow=FALSE){
        'Modified version of `copy\' to allow `activeBindingFunction\' as fields.
'
        .def <- get_def()
        .field.classes <- get_fieldclasses()
        .fields <- names(.field.classes)
        .new <- new(.def)
        .new.env <- as.environment(.new)
        .self.env <- as.environment(.self)
        for (i in seq_along(.fields)) {
          .field <- .fields[i]
          .fclass <- .field.classes[[.field]]
          if (.fclass != "activeBindingFunction"){
            ## logme(.field)
            ## logme(.fclass)
            if (shallow){
              base::assign(.field,get(.field,envir=.self.env),envir=.new.env)
            } else {
              current <- get(.field,envir=.self.env)
              if (is(current,"envRefClass")){
                current <- current$copy2(FALSE)
              }
              ## logme(.field)
              ## logme(current)
              ## logme(.new.env)
              base::assign(.field, current, envir=.new.env) # copy won't work due to override of `assign`
            }
          }
        }
        .new
      },
      saveme=function(outFpre=""){
        .name <- paste(tolower(.class),".self",sep='')
        assign(.name,.self)        
        base::save(file=sprintf("%s__%s.RData",outFpre,.name),list=.name)   
      },
      ##' @param x character, methods to be updated.
      update=function(x){
        'Modify method definition without re-create the class instance. x: character, methods to be updated.
'
        ## See discussions: http://stackoverflow.com/questions/22510100/manual-modifications-of-the-class-definition-of-a-reference-class-instance/22845207
        selfEnv <- as.environment(.self)
        ## .fields <- as.vector(utils::lsf.str(selfEnv))
        for (e in x){
          v <- get(e,eval(parse(text=sprintf("%s@generator$def@refMethods",.class))))
          assign(e,v,envir=selfEnv)
        }
      }
      )
    )


