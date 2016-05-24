#' attribute.profile.class
#' @slot results A data frame of probability of mastering each attribute profile
#' @slot attribute.profile.matrix A matrix of all possible attribute profiles
#' @exportClass attribute.profile.class
attribute.profile.class <- setClass("attribute.profile.class",
                          slots = c(results = "data.frame",
                                    attribute.profile.matrix = "matrix"),
                          validity = function(object) return(TRUE)
)

#' attribute.class
#' @slot results A data frame of probability of mastering each attribute
#' @exportClass attribute.class
attribute.class <- setClass("attribute.class", 
                              slots = c(results = "data.frame"),
                              validity = function(object) return(TRUE)
)

#' parameter.class
#' @slot results A data frame of parameter estimates 
#' @exportClass parameter.class
parameter.class <- setClass("parameter.class",
                              slots = c(results = "data.frame"),
                              validity = function(object) return(TRUE)
)

#' all.results.class
#' @slot attribute.profile.result An object of \code{\link{attribute.profile.class}}
#' @slot attribute.result An object of \code{\link{attribute.class}}
#' @slot parameter.result An object of \code{\link{parameter.class}}
#' @exportClass all.results.class
all.results.class <- setClass("all.results.class",
                              slots = c(attribute.profile.result = "attribute.profile.class",
                                        attribute.result = "attribute.class",
                                        parameter.result = "parameter.class")
)

#' dcm.scorer.class
#' @slot inputs A list of inputs provided by user
#' @slot mcmc.inputs A list of inputs for MCMC estimation
#' @slot results An object of \code{\link{all.results.class}}
#' @slot mcmc.output A list of MCMC runs
#' @exportClass dcm.scorer.class
dcm.scorer.class <- setClass("dcm.scorer.class",
                             slots = c(inputs = 'list',
                                       mcmc.inputs = 'list',
                                       results = 'all.results.class',
                                       mcmc.output = 'list')
)


