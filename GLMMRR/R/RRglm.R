
#' Fitting Generalized Linear Models with binary Randomized Response data
#'
#' Fit a generalized linear model (GLM) with binary Randomized Response data.
#' Implemented as a wrapper for \code{\link{glm}}. Reference: Fox, J-P, Klotzke, K. and Veen, D. (2016).
#' \emph{Generalized Linear Mixed Models for Randomized Responses.} Manuscript submitted for publication.
#'
#' @param formula
#' a two-sided linear formula object describing the model to be fitted,
#' with the response on the left of a ~ operator and the terms, separated by + operators, on the right.
#' @param link
#' a glm link function for binary outcomes. Must be a function name.
#' Available options: "RRlink.logit", "RRlink.probit", "RRlink.cloglog" and "RRlink.cauchit"
#' @param RRmodel
#' the Randomized Response model. Can be a single value or a vector of models.
#' Available options: "DQ", "Warner", "Forced", "UQM", "Crosswise", "Triangular" and "Kuk"
#' @param p1
#' the Randomized Response parameter p1. Must be 0 <= p1 <= 1.
#' @param p2
#' the Randomized Response parameter p2. Must be 0 <= p2 <= 1.
#' @param data
#' a data frame containing the variables named in \code{\link{formula}} as well as the Randomized Response model and parameters.
#' If the required information cannot be found in the data frame, or if no data frame is given, then the variables are taken
#' from the environment from which RRglm is called.
#' @param na.action
#' a function that indicates what should happen when the data contain NAs.
#' The default action (\code{\link{na.omit}}, as given by \code{getOption("na.action"))})
#' strips any observations with any missing values in any variables.
#' @param ...
#' other potential arguments to be passed to \code{\link{glm}}.
#'
#' @return
#' An object of class RRglm. Extends the class \code{glm} with Randomize Response data.
#' @export
#' @seealso \code{\link{glm}}
#'
#' @examples
#' # Fit the model with fixed effects for gender, RR, pp and age using the logit link function.
#' # The Randomized Response parameters p1, p2 and model
#' # are specified for each observation in the dataset.
#' out <- RRglm(response ~ Gender + RR + pp + age, link="RRlink.logit", RRmodel=RRmodel,
#'          p1=RRp1, p2=RRp2, data=Plagiarism, etastart=rep(0.01, nrow(Plagiarism)))
#' summary(out)
RRglm <- function (formula, link, RRmodel, p1, p2, data, na.action = "na.omit", ...) {

  # Try to find the RR parameters in the data set first
  tryCatch(RRmodel <- eval(substitute(RRmodel), data), error=function(e) NULL)
  tryCatch(p1 <- eval(substitute(p1), data), error=function(e) NULL)

  if(!missing(p2))
    tryCatch(p2 <- eval(substitute(p2), data), error=function(e) NULL)
  else
    p2 <- rep(1, length(p1))

  if(!(length(RRmodel) == length(p1) && length(RRmodel) == length(p2)))
    stop("RR parameter vectors are not of same length.")

  if (!(length(RRmodel) == 1 || length(RRmodel) == nrow(data)))
    stop("Length of RR parameters does not match 1 or number of elements in the data.")

  if (any(is.na(RRmodel) | is.na(p1) | is.na(p2)))
    stop("RR parameters must be specified for each case.")

  # Translate p1, p2 to c, d for the chosen RR model
  RRparameters <- getRRparameters(RRmodel, p1, p2)

  # Create a dataset containing the variables of the formula and the the RR parameters
  varnames <- all.vars(formula)
  RRdata <- data.frame(data[,varnames], "RRmodel" = RRmodel, "c" = RRparameters$c, "d" = RRparameters$d, "p1" = p1, "p2" = p2)

  # Debug
  #View(RRdata)

  # If NA's shall be deleted by glm(), we must do the same for the RR paramters
  if (na.action == "na.omit" || na.action == "na.exclude")
    RRdata <- na.omit(RRdata)

  # Must be a data frame
  df <- as.data.frame(data)

  # Get a reference to the correct link function by name
  glmlink <- match.fun(link)

  # Create call object with given arguments
  cl <- call("glm", formula = formula, family = quote(binomial(link = glmlink(RRdata$c, RRdata$d))), data = quote(df),
             na.action = na.action)

  # Make sure that additional arguments are passed
  m <- match.call(expand.dots = FALSE)
  dots <- m$...
  for (ii in 1:(length(dots))) {
    cl[names(dots[ii])] <- dots[ii]
  }

  # Evaluate the call and fit the model
  output <- eval(cl)

  # Add RR parameters to the output object
  output$RRmodel <- RRdata$RRmodel
  output$RRp1 <- RRdata$p1
  output$RRp2 <- RRdata$p2
  output$RRc <- RRdata$c
  output$RRd <- RRdata$d

  # Return an object of class RRglm
  class(output) <- c("RRglm", "glm", "lm")
  return(output)

}
