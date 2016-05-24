#' Generate a virtual species distributions with responses to environmental variables
#' 
#' This function generates a virtual species distribution from a RasterStack of environmental
#' variables and a defined set of responses to each environmental parameter.
#' 
#' @param raster.stack a RasterStack object, in which each layer represent an environmental 
#' variable.
#' @param parameters a list containing the functions of response of the species to 
#' environmental variables with their parameters. See details.
#' @param rescale \code{TRUE} or \code{FALSE}. If \code{TRUE}, the final probability of presence is rescaled between 0 and 1.
#' @param formula a character string or \code{NULL}. The formula used to combine partial responses into the final
#' environmental suitability value (e.g., \code{"layername1 + 2 * layername2 +
#' layername3 * layername4 etc."}). If \code{NULL} then partial responses will be added or multiplied according to
#' \code{species.type}
#' @param species.type \code{"additive"} or \code{"multiplicative"}. Only used if \code{formula = NULL}. 
#' Defines how the final environmental suitability is calculated: if \code{"additive"}, responses to each
#' variable are summed; if \code{"multiplicative"}, responses are multiplicated.
#' @param rescale.each.response \code{TRUE} or \code{FALSE}. If \code{TRUE}, the individual responses to
#' each environmental variable are rescaled between 0 and 1 (see details).
#' @param plot \code{TRUE} or \code{FALSE}. If \code{TRUE}, the generated virtual species will be plotted.
#' @return a \code{list} with 3 elements:
#' \itemize{
#' \item{\code{approach}: the approach used to generate the species, \emph{i.e.}, \code{"response"}}
#' \item{\code{details}: the details and parameters used to generate the species}
#' \item{\code{suitab.raster}: the raster containing the environmental suitability of the virtual species}
#' }
#' The structure of the virtualspecies object can be seen using str()
#' @seealso \code{\link{generateSpFromPCA}} to generate a virtual species with a PCA approach
#' @details
#' This functions proceeds into several steps:
#' \enumerate{
#' \item{The response to each environmental variable is calculated with the functions provided
#' in \code{parameters}. This results in a probability of presence for each variable.}
#' \item{If \code{rescale.each.response} is \code{TRUE}, each probability of presence is rescaled between 0 and 1.}
#' \item{The final probability of presence is calculated according to the chosen \code{species.type}.}
#' \item{If \code{rescale} is \code{TRUE}, the final probability of presence is rescaled between 0 and 1,
#' with the formula (val - min) / (max - min).}
#' }
#' The RasterStack containing environmental variables must have consistent names, 
#' because they will be checked with the \code{parameters}. For example, they can be named
#' var1, var2, etc. Names can be checked and set with \code{names(my.stack)}.
#' 
#' The \code{parameters} have to be carefully created, otherwise the function will not work:
#' \itemize{
#' \item{Either see \code{\link{formatFunctions}} to easily create your list of parameters}
#' \item{Or create a \code{list} defined according to the following template: \cr
#' \code{list(
#'            var1 = list(fun = 'fun1', args = list(arg1 = ..., arg2 = ..., etc.)),
#'            var2 = list(fun = 'fun2', args = list(arg1 = ..., arg2 = ..., etc.)))}\cr
#' It is important to keep the same names in the parameters as in the stack of environmental
#' variables. Similarly, argument names must be identical to argument names in the associated 
#' function (e.g., if you use \code{fun = 'dnorm'}, then args should look like \code{list(mean = 0, sd = 1)}).
#' 
#' See the example section below for more examples.}}
#'            
#' 
#' Any response function that can be applied to the environmental variables can
#' be chosen here. Several functions are proposed in this package:
#' \code{\link{linearFun}}, \code{\link{logisticFun}} and \code{\link{quadraticFun}}.
#' Another classical example is the normal distribution: \code{\link[stats]{dnorm}}.
#' Ther users can also create and use their own functions.
#' 
#'   
#' If \code{rescale.each.response = TRUE}, then the probability response to each
#' variable will be normalised between 0 and 1 according to the following formula:
#' P.rescaled = (P - min(P)) / (max(P) - min (P))
#' This rescaling has a strong impact on response functions, so users may prefer to
#' use \code{rescale.each.response = FALSE} and apply their own rescaling within
#' their response functions.
#' 
#' 
#' @export
#' @import raster
#' @author
#' Boris Leroy \email{leroy.boris@@gmail.com}
#' 
#' with help from C. N. Meynard, C. Bellard & F. Courchamp
#' @examples
#' # Create an example stack with two environmental variables
#' a <- matrix(rep(dnorm(1:100, 50, sd = 25)), 
#'             nrow = 100, ncol = 100, byrow = TRUE)
#' env <- stack(raster(a * dnorm(1:100, 50, sd = 25)),
#'              raster(a * 1:100))
#' names(env) <- c("variable1", "variable2")
#' plot(env) # Illustration of the variables
#' 
#' # Easy creation of the parameter list:
#' # see in real time the shape of the response functions
#' parameters <- formatFunctions(variable1 = c(fun = 'dnorm', mean = 1e-04, 
#'                                              sd = 1e-04),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#'                               
#' # If you provide env, then you can see the shape of response functions:
#' parameters <- formatFunctions(x = env,
#'                               variable1 = c(fun = 'dnorm', mean = 1e-04, 
#'                                              sd = 1e-04),
#'                               variable2 = c(fun = 'linearFun', a = 1, b = 0))
#' 
#' # Generation of the virtual species
#' sp1 <- generateSpFromFun(env, parameters)
#' sp1
#' par(mfrow = c(1, 1))
#' plot(sp1)
#' 
#' 
#' # Manual creation of the parameter list
#' # Note that the variable names are the same as above
#' parameters <- list(variable1 = list(fun = 'dnorm',
#'                                     args = list(mean = 0.00012,
#'                                                 sd = 0.0001)),
#'                    variable2 = list(fun = 'linearFun',
#'                                     args = list(a = 1, b = 0)))
#' # Generation of the virtual species
#' sp1 <- generateSpFromFun(env, parameters, plot = TRUE)
#' sp1
#' plot(sp1)


generateSpFromFun <- function(raster.stack, parameters, 
                              rescale = TRUE, formula = NULL, 
                              species.type = "multiplicative", rescale.each.response = TRUE,
                              plot = FALSE)
{
  approach <- "response"
  if(!(is(raster.stack, "Raster")))
  {
    stop("raster.stack must be a raster stack object")
  }
  if(any(is.na(maxValue(raster.stack))))
  {
    raster.stack <- setMinMax(raster.stack)
  }
  n.l <- nlayers(raster.stack)
  if(n.l != length(parameters)) 
  {stop("Provide as many layers in raster.stack as functions on parameters")}
  if(any(!(names(parameters) %in% names(raster.stack)) |
           !(names(raster.stack) %in% names(parameters))))
     {stop("Layer names and names of parameters must be identical")}
  # Checking the structure and consistency of parameters
  for (i in 1:length(parameters))
  {
    if(any(!(c("fun", "args") %in% names(parameters[[i]]))))
    {stop("The structure of parameters does not seem correct. 
          Please provide function and arguments for variable '",
          names(parameters)[i], "'. See help(generateSpFromFun) for more details.",
          sep = "")}
    test <- tryCatch(match.fun(parameters[[i]]$fun), error = function(c) "error")
    if(class(test) != "function")
    {
      stop(paste("The function ", parameters[[i]]$fun, " does not exist, please verify spelling.", sep = ""))
    }
    if(any(!(names(parameters[[i]]$args) %in% names(formals(fun = test)))))
    {
      stop(paste("Arguments of variable '", names(parameters)[i], "' (", 
                 paste(names(parameters[[i]]$args), collapse = ", "), 
                 ") do not match arguments of the associated function\n
                 List of possible arguments for this function: ",
                 paste(names(formals(fun = test)), collapse = ", "), sep = ""))
    }
    rm(test)
  }
  suitab.raster <- stack(sapply(names(raster.stack), FUN = function(y)
  {
    calc(raster.stack[[y]], fun = function(x)
    {
      do.call(match.fun(parameters[[y]]$fun), args = c(list(x), parameters[[y]]$args))
    }
    )
  }))
  
  for (var in names(raster.stack))
  {
    parameters[[var]]$min <- raster.stack[[var]]@data@min
    parameters[[var]]$max <- raster.stack[[var]]@data@max
  }
  
  if(rescale.each.response)
  {
    suitab.raster <- stack(sapply(names(suitab.raster), function(y)
      {
        (suitab.raster[[y]] - suitab.raster[[y]]@data@min) / (suitab.raster[[y]]@data@max - suitab.raster[[y]]@data@min)
      }))
  }

  
  if(is.null(formula))
  {
    if(species.type == "multiplicative")
    {
      formula <- paste(names(suitab.raster), collapse = " * ")
      suitab.raster <- raster::overlay(suitab.raster, fun = prod)
    } else if (species.type == "additive")
    {
      formula <- paste(names(suitab.raster), collapse = " + ")
      suitab.raster <- raster::overlay(suitab.raster, fun = sum)
    } else stop("If you do not provide a formula, please choose either species.type = 'additive' or 'multiplicative'")
  } else
  {
    if(any(!(all.vars(reformulate(formula)) %in% names(suitab.raster))))
    {
      stop("Please verify that the variable names in your formula are correctly spelled") 
    } else if(any(!(names(suitab.raster) %in% all.vars(reformulate(formula)))))
    {
      stop("Please verify that your formula contains all the variables of your input raster stack")
    } else
    {
      custom.fun <- NULL # To remove the note in rcheck
      eval(parse(text = paste("custom.fun <- function(",
                              paste(names(suitab.raster), collapse = ", "),
                              ") {",
                              formula,
                              "}"
      )))
      suitab.raster <- raster::overlay(suitab.raster, fun = custom.fun)
      print(formula)
    }
  }

  if(rescale)
  {
    suitab.raster <- (suitab.raster - suitab.raster@data@min) / (suitab.raster@data@max - suitab.raster@data@min)
  }
    

  results <- list(approach = approach,
                  
                  details = list(variables = names(parameters),
                                 formula = formula,
                                 rescale.each.response = rescale.each.response,
                                 rescale = rescale,
                                 parameters = parameters),
                  suitab.raster = suitab.raster
  )

  if(plot)
  {
    plot(results$suitab.raster, main = "Environmental suitability of the virtual species")
  }
  
  class(results) <- append(class(results), "virtualspecies")

  return(results)
}
