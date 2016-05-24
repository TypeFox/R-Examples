##' Extract variable names from latent variable model
##'
##' Extract exogenous variables (predictors), endogenous variables (outcomes),
##' latent variables (random effects), manifest (observed) variables from a
##' \code{lvm} object.
##'
##' \code{vars} returns all variables of the \code{lvm}-object including
##' manifest and latent variables. Similarily \code{manifest} and \code{latent}
##' returns the observered resp. latent variables of the model.
##' \code{exogenous} returns all manifest variables without parents, e.g.
##' covariates in the model, however the argument \code{latent=TRUE} can be used
##' to also include latent variables without parents in the result. Pr. default
##' \code{lava} will not include the parameters of the exogenous variables in
##' the optimisation routine during estimation (likelihood of the remaining
##' observered variables conditional on the covariates), however this behaviour
##' can be altered via the assignment function \code{exogenous<-} telling
##' \code{lava} which subset of (valid) variables to condition on.  Finally
##' \code{latent} returns a vector with the names of the latent variables in
##' \code{x}. The assigment function \code{latent<-} can be used to change the
##' latent status of variables in the model.
##'
##' @aliases vars vars.lvm vars.lvmfit latent latent<- latent.lvm latent<-.lvm
##' latent.lvmfit latent.multigroup manifest manifest.lvm manifest.lvmfit
##' manifest.multigroup exogenous exogenous<- exogenous.lvm exogenous<-.lvm
##' exogenous.lvmfit exogenous.multigroup endogenous endogenous.lvm
##' endogenous.lvmfit endogenous.multigroup
##' @usage
##'
##' vars(x,...)
##'
##' endogenous(x,...)
##'
##' exogenous(x,...)
##'
##' manifest(x,...)
##'
##' latent(x,...)
##'
##' \method{exogenous}{lvm}(x,silent = FALSE, xfree = TRUE,...) <- value
##'
##' \method{exogenous}{lvm}(x,latent=FALSE,index=TRUE,...)
##'
##' \method{latent}{lvm}(x,clear=FALSE,...) <- value
##'
##' @param x \code{lvm}-object
##' @param latent Logical defining whether latent variables without parents
##' should be included in the result
##' @param index For internal use only
##' @param clear Logical indicating whether to add or remove latent variable
##' status
##' @param silent Suppress messages
##' @param xfree For internal use only
##' @param value Formula or character vector of variable names.
##' @param \dots Additional arguments to be passed to the low level functions
##' @return Vector of variable names.
##' @author Klaus K. Holst
##' @seealso \code{\link{endogenous}}, \code{\link{manifest}},
##' \code{\link{latent}}, \code{\link{exogenous}}, \code{\link{vars}}
##' @keywords models regression
##' @examples
##'
##' g <- lvm(eta1 ~ x1+x2)
##' regression(g) <- c(y1,y2,y3) ~ eta1
##' latent(g) <- ~eta1
##' endogenous(g)
##' exogenous(g)
##' identical(latent(g), setdiff(vars(g),manifest(g)))
##'
##' @export
`vars` <-
function(x,...) UseMethod("vars")

##' @export
`vars.graph` <-
  function(x,...) {
    graph::nodes(x)
  }

##' @export
`vars.lvm` <-
  function(x,...) {
    colnames(x$M)
  }

##' @export
`vars.lvmfit` <-
  function(x,...) {
    vars(Model(x),...)
  }

##' @export
vars.list <- function(x,...) {
  varlist <- c()
  for (i in seq_along(x)) {
    varlist <- c(varlist, vars(x[[i]]))
  }
  varlist <- unique(varlist)
  return(varlist)
}
