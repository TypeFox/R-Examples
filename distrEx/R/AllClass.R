.onLoad <- function(lib, pkg){
#    require("methods", character = TRUE, quietly = TRUE)
}
.ULC.cast <- function(x){
         if( is(x,"AbscontDistribution"))
             x <- as(as(x,"AbscontDistribution"), "UnivarLebDecDistribution")
         if(is(x,"DiscreteDistribution"))
             x <- as(as(x,"DiscreteDistribution"), "UnivarLebDecDistribution")
         if(!is(x,"UnivarLebDecDistribution"))
            x <- as(x,"UnivarLebDecDistribution")
         return(x)
}


.onAttach <- function(library, pkg)
{
unlockBinding(".distrExOptions", asNamespace("distrEx"))
msga <- gettext("Note: Packages \"e1071\", \"moments\", \"fBasics\" should be attached ")
msgb <- gettext("/before/ package \"distrEx\". See distrExMASK().")
msgc <- gettext("Note: Extreme value distribution functionality has been moved to
      package \"RobExtremes\". See distrExMOVED().")
buildStartupMessage(pkg = "distrEx", msga, msgb, msgc, library = library, packageHelp = TRUE,
#                    MANUAL="http://www.uni-bayreuth.de/departments/math/org/mathe7/DISTR/distr.pdf",
VIGNETTE = gettext("Package \"distrDoc\" provides a vignette to this package as well as to several related packages; try vignette(\"distr\")."))
  invisible()
}

distrExMASK <- function(library = NULL) 
{
    infoShow(pkg = "distrEx", filename = "MASKING", library = library)
}

distrExMOVED <- function(library = NULL)
{
    infoShow(pkg = "distrEx", filename = "MOVED", library = library)
}

.onUnload <- function(libpath)
{
    library.dynam.unload("distrEx", libpath)
}


# multivariate distribution
setClass("MultivariateDistribution", 
            prototype = prototype(r = function(n){ matrix(rep(c(0,0), n), ncol=2) }, 
                                  d = NULL, p = NULL, q = NULL, param = NULL,
                                  img = new("EuclideanSpace", dimension = 2),
                                  support = matrix(c(0,0), ncol = 2)),
            contains = "Distribution")

# discrete mulitvariate distribution
setClass("DiscreteMVDistribution", representation(support = "matrix"),
            prototype(r = function(n){ matrix(rep(c(0,0), n), ncol=2) }, 
                      d = NULL, p = NULL, q = NULL, param = NULL,
                      img = new("EuclideanSpace", dimension = 2),
                      support = matrix(c(0,0), ncol = 2)),
            contains = "MultivariateDistribution")

# condition
setClass("Condition", representation(name = "character"),
            prototype(name = "a condition"))

# conditioning by an Euclidean space
setClass("EuclCondition", 
            representation(Range = "EuclideanSpace"), 
            prototype(name = gettext("conditioning by an Euclidean space"),
                      Range = new("EuclideanSpace")),
            contains = "Condition")

## conditional univariate distribution
setClass("UnivariateCondDistribution", 
            representation(cond = "Condition"), 
            prototype(r = function(n, cond){ rnorm(n, mean = 0, sd = 1) },
                      d = NULL, p = NULL, q = NULL, img = new("Reals"), 
                      param = NULL, cond = new("Condition")),
            contains = "Distribution")

# absolutely continuous conditional distribution
setClass("AbscontCondDistribution", 
            representation(cond = "Condition"), 
            prototype(r = function(n, cond){ rnorm(n, mean = 0, sd = 1) },
                      d = NULL, p = NULL, q = NULL, img = new("Reals"), 
                      param = NULL, cond = new("Condition")),
            contains = "UnivariateCondDistribution")

# discrete conditional distribution
setClass("DiscreteCondDistribution", 
            representation(support = "function", 
                           cond = "Condition"), 
            prototype(r = function(n, cond){ rep(0, n) },
                      d = NULL, p = NULL, q = NULL, img = new("Reals"), 
                      param = NULL, support = function(cond){0},
                      cond = new("Condition")),
            contains = "UnivariateCondDistribution")


# Parameter of a linear regression model (with intercept and scale)
setClass("LMParameter", 
            representation(theta = "numeric",
                           intercept = "numeric",
                           scale = "numeric"), 
            prototype(name = gettext("parameter of a linear regression model"),
                      theta = 0, intercept = 0, scale = 1),
            contains = "Parameter",
            validity = function(object){
                if(any(!is.finite(object@theta)))
                    stop("inifinite or missing values in 'theta'")
                if(length(object@intercept) != 1)
                    stop("'intercept' has to be of length 1")
                if(!is.finite(object@intercept))
                    stop("inifinite or missing value in 'intercept'")
                if(length(object@scale) != 1)
                    stop("'scale' has to be of length 1")
                if(!is.finite(object@scale))
                    stop("inifinite or missing value in 'scale'")
                return(TRUE)
            })


