## ------------------------------------------------------------------------
library(settings)
par('lwd')
par(lwd=4)
par('lwd')
reset_par()
par('lwd')

## -----------------------------------------------------------------------------
reset_options()

## -----------------------------------------------------------------------------
my_options <- options_manager(foo = 1, bar = 2, baz = 'hello')
my_options()

## -----------------------------------------------------------------------------
my_options('foo')
my_options('foo','baz')

## -----------------------------------------------------------------------------
my_options(foo=7)
my_options()
# or multiple options at once
my_options(foo=7,bar=0)
my_options()

## -----------------------------------------------------------------------------
reset(my_options)
my_options()

## -----------------------------------------------------------------------------
opt <- options_manager(foo="up", bar=2
  , .allowed = list(
      foo = inlist("up","down")
    , bar = inrange(min=0, max=3)
  )
)

## ---- eval=FALSE--------------------------------------------------------------
#  > opt(foo="middle")
#  Error: Option value out of range. Allowed values are up, down
#  > opt(bar=7)
#  Error: Option value out of range. Allowed values are in [0, 3]

## -----------------------------------------------------------------------------
my_options <- options_manager(a=2,b=3)

## -----------------------------------------------------------------------------
f <- function(x,...){
  # create local copy of options, merged with the global options.
  local_opts <- clone_and_merge(my_options,...)
  # local options can be used
  local_opts('a') + local_opts('b') * x 
}

## -----------------------------------------------------------------------------
# a and b are taken from global option set.
f(1)         # 2 + 3 * 1
# specify 'a'
f(1,a=10)    # 10 + 3 * 1
#specify 'a' and 'b'
f(1,a=10,b=100) # 10 + 100 * 1

# global options are unaltered, as expected.
my_options()

## ----eval=FALSE---------------------------------------------------------------
#  # Variable, global to package's namespace.
#  # This function is not exported to user space and does not need to be documented.
#  MYPKGOPTIONS <- options_manager(a=1, b=2)
#  
#  # User function that gets exported:
#  
#  #' Set or get options for my package
#  #'
#  #' @param ... Option names to retrieve option values or \code{[key]=[value]} pairs to set options.
#  #'
#  #' @section Supported options:
#  #' The following options are supported
#  #' \itemize{
#  #'  \item{\code{a}}{(\code{numeric};1) The value of a }
#  #'  \item{\code{b}}{(\code{numeric};2) The value of b }
#  #' }
#  #'
#  #' @export
#  pkg_options <- function(...){
#    # protect against the use of reserved words.
#    stop_if_reserved(...)
#    MYPKGOPTIONS(...)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  #' Reset global options for pkg
#  #'
#  #' @export
#  pkg_reset() reset(MYPKGOPTIONS)

## -----------------------------------------------------------------------------
# general options manager, will be invisible to user.
opt <- options_manager(foo=1,bar=2)

## -----------------------------------------------------------------------------
# class definition containing default options in prototype.
TestClass <- setClass("TestClass"
  , slots=list(options='function',value='numeric')
  , prototype = list(
     options = opt
     , value = 0
    )
)

## -----------------------------------------------------------------------------
setGeneric("test_options",function(where=NULL,...) standardGeneric("test_options"))

# method for accessing global options
setMethod("test_options","ANY",function(where=NULL,...){
  do.call(opt,c(where,list(...)))
})

# method for getting/setting functions in a slot.
setMethod("test_options","TestClass", function(where=NULL,...){
  if (is_setting(...)){
    where@options <- clone_and_merge(where@options,...)
    where
  } else {
    where@options(...)
  }
})

## -----------------------------------------------------------------------------
# instantiate a class; with global options as currently set.
test <- TestClass()

# get global options
test_options()

# set a global option
test_options(foo=2)
test_options('foo')
# check that 'test' uses global option
test_options(test)

# set local option
test <- test_options(test,bar=3)
test_options(test)
# check global option
test_options()

## -----------------------------------------------------------------------------
opt <- options_manager(foo=1,bar=2)

## -----------------------------------------------------------------------------
RefTest <- setRefClass("RefTest"
  , fields =  list(.options='function',value='numeric')
  , methods = list(
    initialize = function(){
      .self$.options <- opt
      .self$value <- 0
    }
    , options = function(...){
        if(is_setting(...)){
          .self$.options <- clone_and_merge(.self$.options,...)
        } else {
          .self$.options(...)
        }
      }
    , reset = function(){
        # explicitly reference the 'settings' package here to avoid recursion.
        settings::reset(.self$.options) 
    }
    )
)

## -----------------------------------------------------------------------------
reftest <- RefTest()

reftest$options()

# set global options
opt(foo=10)
reftest$options()

# set local options
reftest$options(bar=3)
reftest$options()
opt()

# reset local options
reftest$reset()
reftest$options()

