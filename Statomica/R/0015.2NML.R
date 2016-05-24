setClass( "DensityFunNML"
        , representation( "function"
                        , family     = "Family"
                        , COMP       = "Complexity"
                        )
        )
setValidity("DensityFunNML", function(object)
{
       len.ok <- length(object@COMP) == 1
       ok <- len.ok
       if(!ok)
       { printInvalid(object); browser()}
       ok
})
new.DensityFunNML <- function( object
                             , family
                             , COMP
                             ) #, weight = blank.Weight())
{
       new( "DensityFunNML"
          , object
          , family = family
          , COMP   = COMP
          ) # , weight = weight)
}
DensityFunNML <- function( object
                         , COMP               = NULL
                         , truncation.tol     = numeric(0)
                         , incidental.x       = numeric(0)
                         , weight             = numeric(0)
                         , allow.0.Complexity = FALSE
                         , verbose
                         , allow.unused.args  = FALSE
                         , ...
                         )
{
       assert.is(object, c("extendedFamily", "finiteFamily"))
       if("x" %in% names(list(...)))
               stop("x was specified for DensityFunNML")
       if("param0" %in% names(list(...)))
               stop("param0 was specified for DensityFunNML")
       if(is(incidental.x, "xprnSetObject") && is.nothing(weight)) #  && !missing(i) && length(i) == 1)
       {
               weight <- Weight(x = incidental.x, incidental = TRUE) #as(incidental.x, "Weight")
               es <- incidental.x # [-i, ]
               #XXX: from statistic to nstatistic
               incidental.x <- statistic(x = es, object = object)
               stopifnot(length(incidental.x) == length(weight@incidental.weight))
       }
#       else
#               stopifnot(missing(i))
       assert.is(incidental.x, "numeric")
       assert.is(weight, c("numeric", "Weight"))
       COMP <- if(is.nothing(COMP))
       {
          if(missing(verbose)) verbose <- default(TRUE, "verbose")
          Complexity( object             = object
                    , truncation.tol     = truncation.tol
                    , incidental.x       = incidental.x
                    , weight             = weight
                    , allow.0.Complexity = allow.0.Complexity
                    , ...
                    )
       }
       else if(length(list(...)) == 0 || allow.unused.args)
       {
          if(missing(verbose)) verbose <- FALSE
          if(is(COMP, "Complexity"))
            COMP
          else if(is(COMP, "numeric"))
            stop("numeric non-Complexity COMP is not accepted as of 100916") #Scalar(COMP)
          else
            stop("non-numeric COMP")
       }
       else
          stop("DensityFunNML cannot use the COMP provided as the parametric complexity of object")
       if(verbose)
          message("COMP: ", COMP)
       assert.is(COMP, "Complexity", "outside DensityFun")
       DensityFun <- function(x, ...)
       {
          assert.is(COMP, "Complexity", "inside DensityFun")
          NML( object             = object
             , x                  = x
             , COMP               = COMP
             , truncation.tol     = truncation.tol
             , incidental.x       = incidental.x
             , weight             = weight
             , allow.0.Complexity = allow.0.Complexity
             , ...
             )
       }
       new.DensityFunNML(  object = DensityFun
                         , family = object
                         , COMP   = COMP
                         )
} # end DensityFunNML