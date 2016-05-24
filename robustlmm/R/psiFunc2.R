##' Create psi_func_cached object using cached numerical integration for
##' E... slots.
##'
##' @title psiFuncCached constructor
##' @param rho rho-function
##' @param psi psi-function
##' @param wgt wgt-function
##' @param Dwgt derivative of weight function
##' @param Dpsi derivative of psi
##' @param name descriptor of this function family
##' @param ... default values for tuning constants
##' @return psi_func_cached-class object
##' @section Warning: the E... slots will not be fully functional: they just
##' return the value for the current defaults and ignore their
##' arguments.
##' @seealso \code{\link{psi_func_cached-class}}
##' @keywords utilities
##' @examples
##' ## re-define cPsi as psiFuncCached.
##' F0 <- function(x=1, .) rep.int(0, length(x))
##' F1 <- function(x=1, .) rep.int(1, length(x))
##' cPsi2 <- psiFuncCached(rho = function(x, .) x^2 / 2,
##'                        psi = function(x, .) x,
##'                        wgt = F1, Dwgt = F0, Dpsi = F1, 
##'                        name = "classic (x^2/2)",
##'                        . = Inf ## dummy, need at least one parameter
##'                        )
##' stopifnot(all.equal(cPsi@@Erho(), cPsi2@@Erho()),
##'           all.equal(cPsi@@Epsi2(), cPsi2@@Epsi2()),
##'           all.equal(cPsi@@EDpsi(), cPsi2@@EDpsi()))
##' @export
psiFuncCached <- function(rho,psi,wgt,Dwgt,Dpsi,name=NULL, ...) {
    lent <- length(dotsargs <- list(...))
    ## '...'  must contain all tuning parameters and their defaults:
    stopifnot(length(nt <- names(dotsargs)) == lent,
              all(nchar(nt)) >= 1)
    if(lent >= 1) {
        ## rho, psi,... checking: must have argument names
        argn <- c("x", nt)
        for(fnam in list("rho", "psi", "wgt", "Dwgt", "Dpsi")) {
            f <- get(fnam, inherits = FALSE)
            ef <- environment(f)
            nf <- names(ff <- formals(f)) # "x" and "k" for Huber's
            if(!identical(nf, argn))
                stop("arguments of function '",fnam,"' are (",
                     paste(nf,  collapse=","),") but should be (",
                     paste(argn,collapse=","),").")
            
            formals(f)[-1] <- dotsargs
            environment(f) <- ef
            assign(fnam, f, inherits = FALSE)
        }
    }

    Erho.val <- integrate(function(x) rho(x)*dnorm(x),-Inf, Inf,
                          rel.tol = .Machine$double.eps^0.5)$value
    Epsi2.val <- integrate(function(x) psi(x)^2*dnorm(x),-Inf, Inf,
                           rel.tol = .Machine$double.eps^0.5)$value
    EDpsi.val <- integrate(function(x) Dpsi(x)*dnorm(x),-Inf, Inf,
                           rel.tol = .Machine$double.eps^0.5)$value
    
    new("psi_func_cached",
        rho = new("functionX", rho),
        psi = new("functionX", psi),
        wgt = new("functionX", wgt),
        Dpsi= new("functionX", Dpsi),
        Dwgt= new("functionX", Dwgt),
        ## tNams = if(lent) nt else character(0),
        tDefs = if(lent) unlist(dotsargs) else numeric(0),
        Erho= Erho <- new("functionXal", function(arg=1) rep(Erho.val, length(arg))),
        Epsi2= Epsi2 <- new("functionXal", function(arg=1) rep(Epsi2.val, length(arg))),
        EDpsi= EDpsi <- new("functionXal", function(arg=1) rep(EDpsi.val, length(arg))),
        name= name
        )
}

##' Change the default arguments for a psi_func_cached object
##'
##' @title Change default arguments
##' @param ... arguments to change
##' @keywords utilities
##' @examples
##' hPsi <- chgDefaults(huberPsi, k=2)
##' curve(huberPsi@@psi(x), 0, 3)
##' curve(hPsi@@psi(x), 0, 3, color="blue", add=TRUE)
##' @export
setMethod("chgDefaults", signature("psi_func_cached"),
          function(object, ...) {
              ##cat("~~~~ chgDefaults of psi_func_cached ~~~~~\n")
              lent <- length(dotsargs <- list(...))
              ## '...'  must contain all tuning parameters and their defaults:
              stopifnot(length(nt <- names(dotsargs)) == lent,
                        all(nchar(nt)) >= 1)
              if(lent >= 1) {
                  ## rho "..." must conform to rho, etc:
                  nf <- names(ff <- formals(object@rho))
                  if(!identical(nf[-1], nt))
                     stop("invalid tuning parameter names: ",
                          paste(nt,    collapse=",")," instead of ",
                          paste(nf[-1],collapse=","),".")

                  for(fnam in list("rho", "psi", "wgt", "Dwgt", "Dpsi")) {
                      f <- slot(object, fnam)
                      ef <- environment(f)
                      formals(f)[-1] <- dotsargs
                      environment(f) <- ef
                      ## lowlevel {faster than}: slot(..) <- new("functionX", f)
                      slot(object, fnam)@.Data <- f
                  }
                  object@tDefs <- unlist(dotsargs)
              }

              Erho.val <- integrate(function(x) object@rho(x)*dnorm(x),-Inf, Inf,
                                    rel.tol = .Machine$double.eps^0.5)$value
              Epsi2.val <- integrate(function(x) object@psi(x)^2*dnorm(x),-Inf, Inf,
                                     rel.tol = .Machine$double.eps^0.5)$value
              EDpsi.val <- integrate(function(x) object@Dpsi(x)*dnorm(x),-Inf, Inf,
                                     rel.tol = .Machine$double.eps^0.5)$value
              object@Erho <- new("functionXal", function(arg=1) rep(Erho.val, length(arg)))
              object@Epsi2 <- new("functionXal", function(arg=1) rep(Epsi2.val, length(arg)))
              object@EDpsi <- new("functionXal", function(arg=1) rep(EDpsi.val, length(arg)))
              
              object
          })

.sprintPsiFunc <- function(x, short=FALSE) {
    v <- x@tDefs
    n <- names(v)
    ## do not print a single dummy parameter "."
    if (length(n) == 1 && n == ".") {
        v <- numeric(0)
        n <- character(0)
    }
    name <- x@name
    if (short) name <- gsub('\\s?(psi|function|\\(.*\\))', '', name)
    if (length(v) >= 1) {
        paste(name, " (",
              paste(n, round(v, 3), sep = " = ", collapse = ", "), ")",
              sep="")
    } else name
}

## from example(psiFunc)
F0 <- function(x=1, .) rep.int(0, length(x))
F1 <- function(x=1, .) rep.int(1, length(x))
FF1 <- function(.) rep.int(1, length(.))
FF1.2 <- function(.) rep.int(1/2, length(.))

##' \eqn{\psi}{Psi}-functions are used by \code{\link{rlmer}}
##' in the estimating equations and to compute robustness
##' weights. Change tuning parameters using \code{\link{chgDefaults}}
##' and convert to squared robustness weights using the
##' \code{\link{psi2propII}} function.
##' 
##' The \bold{\dQuote{classical} \eqn{\psi}{psi}-function \code{cPsi}}
##' can be used to get a non-robust, i.e., classical, fit.
##' The \code{psi} slot equals the identity function, and
##' the \code{rho} slot equals quadratic function. Accordingly,
##' the robustness weights will always be 1 when using \code{cPsi}.
##'
##' The \bold{Huber \eqn{\psi}{psi}-function \code{huberPsi}} is identical to
##' the one in the package \code{robustbase}. The \code{psi} slot equals
##' the identity function within \eqn{\pm k}{+-k} (where \eqn{k}{k} is
##' the tuning parameter). Outside this interval it is equal to
##' \eqn{\pm k}{+-k}. The \code{rho} slot equals the quadratic
##' function within \eqn{\pm k}{+-k} and a linear function outside.
##'
##' The \bold{smoothed Huber \eqn{\psi}{psi}-function} is very similar to
##' the regular Huber \eqn{\psi}{psi}-function.
##' Instead of a sharp bend like the Huber function,
##' the smoothe Huber function bends smoothly. The first tuning
##' contant, k, can be compared to the tuning constant
##' of the original Huber function. The second tuning
##' constant, s, determines the smoothness of the bend.
##'
##' @title Classical, Huber and smoother Huber psi- or rho-functions
##' @name psi-functions
##' @rdname psi-functions
##' @aliases cPsi huberPsi smoothPsi
##' @usage ## see examples
##' @seealso \code{\link{chgDefaults}} and \code{\link{psi2propII}}
##' for changing tuning parameters;
##' \code{\link{psi_func_cached-class}} and
##' \code{\link{psi_func-class}} for a more detailed description of the
##' slots; \code{\link{psiFuncCached}} for a constructor function to
##' create custom \eqn{\psi}{psi}-functions.
##' @examples
##' plot(cPsi)
##' plot(huberPsi)
##' plot(smoothPsi)
##' curve(cPsi@@psi(x), -3, 3)
##' curve(smoothPsi@@psi(x, 1.345, 10), -3, 3, add=TRUE, col="red")
##' curve(huberPsi@@psi(x, 1.345), -3, 3, add=TRUE, col="blue")
##' @export cPsi
cPsi <- psiFunc(rho = function(x, .) x^2 / 2, psi = function(x, .) x,
                 wgt = F1, Dwgt = F0, Dpsi = F1, Erho = FF1.2,
                 Epsi2 = FF1, EDpsi = FF1,
                 name = "classic (x^2/2)", . = Inf)

##' @exportMethod plot
##' @export huberPsi
##' @export
smoothPsi <- psiFuncCached(rho = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, x^2/2, c^2/2 + k*(ax-c) -
                                       ((ax-d)^(1-s) - a^(1-s))/(1-s))
                            },
                            psi = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, x, sign(x)*(k - (ax-d)^(-s)))
                            },
                            Dpsi = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 1, s*(ax-d)^(-s-1))
                            },
                            wgt = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 1, (k - (ax-d)^(-s))/ax)
                            },
                            Dwgt = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 0,
                                       (ax - d)^(-s-1)*s/x -
                                       (k - (ax-d)^(-s))/(x*ax))
                            },
                            k = 1.345, s = 10,
                            name = "smoothed Huber")


.psi2propII <- function(object, ...) {
    ## do not do anything for cPsi
    if (identical(object, cPsi)) return(object)
    
    ## Convert a regular psi-function into a proposal II psi function
    ## (with squared weights)
    f <- formals(object@psi)
    nf <- names(f)
    args <- paste(nf, collapse=",")
    x <- nf[1]

    ## wgt
    fun <- paste("function(",args,") object@wgt(", args, ")^2")
    wgt <- eval(parse(text=fun))
    formals(wgt) <- f
    ## Dwgt
    fun <- paste("function(",args,") 2*object@wgt(", args, ")*object@Dwgt(",args,")")
    Dwgt <- eval(parse(text=fun))
    formals(Dwgt) <- f
    ## psi
    fun <- paste("function(",args,") object@wgt(", args, ")*object@psi(",args,")")
    psi <- eval(parse(text=fun))
    formals(psi) <- f
    ## Dpsi
    fun <- paste("function(",args,") object@wgt(", args, ")*object@Dpsi(",args,
                 ") + object@Dwgt(", args, ")*object@psi(",args,")")
    Dpsi <- eval(parse(text=fun))
    formals(Dpsi) <- f
    ## rho
    intRho <- function(psi, x, ...) {
        ret <- x
        for (i in seq_along(x)) {
            if (is.infinite(x[i])) next
            ret[i] <- integrate(psi, 0, x[i], ..., rel.tol = .Machine$double.eps^0.5)$value
        }
        ret
    }
    fun <- paste("function(",args,") intRho(psi,",args,")")
    rho <- eval(parse(text=fun))
    formals(rho) <- f

    ret <- do.call(psiFuncCached, c(list(wgt=wgt, Dwgt=Dwgt, psi=psi, Dpsi=Dpsi, rho=rho),
                                    f[-1], name=paste(object@name, ", Proposal II", sep="")))
    ## if ... is given: pass it to chgDefaults
    chgArgs <- list(...)
    if (length(chgArgs) > 0) {
        if (is.null(names(chgArgs))) stop("Extra arguments in ... need to be named")
        ## extend list, all arguments need to be passed to chgDefaults
        for (name in names(ret@tDefs))
            if (is.null(chgArgs[[name]])) chgArgs[[name]] <- ret@tDefs[[name]]
        ret <- do.call("chgDefaults", c(list(ret), chgArgs))
    }
    return(ret)
}

## hP2 <- psi2propII(huberPsi)
## x <- -3:10
## all.equal(hP2@wgt(x), huberPsi@wgt(x)^2)
## all.equal(hP2@Dwgt(x), 2*huberPsi@wgt(x)*huberPsi@Dwgt(x))
## all.equal(hP2@psi(x), huberPsi@wgt(x)^2*x)
## all.equal(hP2@Dpsi(x), 2*huberPsi@wgt(x)*huberPsi@Dwgt(x)*x + huberPsi@wgt(x)^2)
## all.equal(hP2@Dpsi(x), huberPsi@Dwgt(x)*huberPsi@psi(x) + huberPsi@wgt(x)*huberPsi@Dpsi(x))

##' Converts the psi_func object into a function that corresponds
##' to Proposal II, i.e., a function of the squared weights.
##' The other elements of the psi_func object are adapted accordingly.
##'
##' @title Convert to Propsal II weight function
##' @param object psi_func object to convert
##' @param ... optional, new default arguments passed to chgDefaults.
##' @aliases psi2propII,psi_func-method
##' @keywords utilities
##' @examples
##' par(mfrow=c(2,1))
##' plot(smoothPsi)
##' plot(psi2propII(smoothPsi))
##' @export
setGeneric("psi2propII", function(object, ...) standardGeneric("psi2propII"))
##' @exportMethod psi2propII
setMethod("psi2propII", signature("psi_func_cached"), .psi2propII)
##' @exportMethod psi2propII
setMethod("psi2propII", signature("psi_func"), .psi2propII)
