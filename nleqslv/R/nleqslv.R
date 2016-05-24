#
# Interface to Fortran library for solving a system of nonlinear equations
# with either a Broyden or a full Newton method
# There a six global search methods:
#   cubic, quadratic and geometric linesearch
#   double dogleg trust region a la Dennis Schnabel
#   powell single dogleg a la Minpack
#   so-called hook step (Levenberg-Marquardt)
#

nleqslv <- function(x, fn, jac = NULL, ...,
                    method = c("Broyden", "Newton"),
                    global = c("dbldog", "pwldog", "cline", "qline", "gline", "hook", "none"),
                    xscalm = c("fixed","auto"),
                    jacobian=FALSE,
                    control = list())
{
    fn1  <- function(par) fn(par, ...)

    if( is.null(jac ) ) {
        jacfunc <- NULL
    } else {
        if(!is.function(jac)) stop("argument 'jac' is not a function!")
        jacfunc <- function(par) jac(par, ...)
    }

    method <- match.arg(method)
    global <- match.arg(global)
    xscalm <- match.arg(xscalm)

    ## Defaults
    con <- list(ftol=1e-8, xtol=1e-8,
                btol=1e-3,
                stepmax=-1.0, delta="newton", sigma=0.5,
                scalex = rep(1, length(x)),
                maxit=150,
                trace=0,
                chkjac=FALSE,
                cndtol=1e-12,
                allowSingular=FALSE,
                dsub=-1L,
                dsuper=-1L
               )

    # limit maximum number of iterations for pure local strategy
    if( global == "none" ) con$maxit=20

    # strict validity test of control
    # based on test of control argument in nlminb
    if( length(control) ) {
        namc <- names(control)
        if( !is.list(control) || is.null(namc) )
            stop("'control' argument must be a named list")
        # check names of control argument
        if( !all(namc %in% names(con)) )
            stop("unknown names in control: ", paste(sQuote(namc[!(namc %in% names(con))]), collapse=", "))
        con[namc] <- control
    }

    tmp <- con[["delta"]]
    if( is.character(tmp) ) {
        k <- match(tolower(tmp), c("cauchy","newton"))
        con[["delta"]] <- as.numeric(-k)
    }
    else if( !is.numeric(tmp) ) stop('control["delta"] should be either character or numeric')

    # to reset flag for checking recursive calls (not allowed for now)
    on.exit(.C("deactivatenleq",PACKAGE="nleqslv"))
    out <- .Call("nleqslv", x, fn1, jacfunc, method, global, xscalm, jacobian, con, new.env(), PACKAGE = "nleqslv")

    out
}
