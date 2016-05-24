##' Calculates a "residual-centered" interaction regression.
##'
##' Given a fitted \code{lm}, this function scans for coefficients
##' estimated from "interaction terms" by checking for colon
##' symbols. The function then calculates the "residual centered"
##' estimate of the interaction term and replaces the interaction term
##' with that residual centered estimate. It works for any order of
##' interaction, unlike other implementations of the same
##' approach. The function \code{lmres} in the now-archived package
##' pequod was a similar function.
##' @param model A fitted lm object
##' @export 
##' @rdname residualCenter
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @references  Little, T. D., Bovaird, J. A.,
##' & Widaman, K. F. (2006). On the Merits of Orthogonalizing
##' Powered and Product Terms: Implications for Modeling
##' Interactions Among Latent Variables.
##' Structural Equation Modeling, 13(4), 497-519.
residualCenter <- function(model)
{
    UseMethod("residualCenter")
}

##' @rdname residualCenter
##' @export
##' @return a regression model of the type as the input model, with
##' the exception that the residualCentered predictor is used in place
##' of the original interaction. The return model includes new
##' variable centeringRegressions: a list including each of the
##' intermediate regressions that was calculated in order to create
##' the residual centered interaction terms. These latter objects may
##' be necessary for diagnostics and to calculate predicted values for
##' hypothetical values of the inputs. If there are no interactive
##' terms, then NULL is returned.
##' @method residualCenter default
##' @export
##' @example inst/examples/residualCenter-ex.R
residualCenter.default <- function (model)
{
    createRCinteraction <- function(x, data) {
        nterms <- length(strsplit( x, ":")[[1]])
        if (nterms == 1) return(NULL)
        dv <- paste("I(", gsub(":", "*", x), ")", sep = "")
        iv <- paste("(", gsub(":"," + ",x),")")
        if(nterms >= 3) iv <- paste(iv, "^", nterms-1, sep = "")
        myformula <- paste(dv, "~", iv)
        myformula <- as.formula(myformula)
        rcterms <-  attr(terms(myformula), "term.labels")
        rcterms <- gsub(":", ".X.", rcterms)
        finalFormula <- paste(dv, "~", paste(rcterms, collapse = " + "))
        aReg <- lm(finalFormula, as.data.frame(data))
        list(residuals = as.numeric(resid(aReg)), reg = aReg)
    }

    dat <- model.matrix(model)

    if (colnames(dat)[1] == "(Intercept)")  dat <- dat[ , -1] # remove intercept
    datn <- colnames(dat)

    frmla <- paste(datn, collapse = " + ")

    rcRegressions <- list()
    for (i in seq_along(datn)) {
        newname <- gsub(":", ".X.", datn[i])
        res <- createRCinteraction( datn[i] , data = dat)
        if (!is.null(res)){
            dat[ , datn[i]] <- res[["residuals"]]
            rcRegressions[[newname]] <- res[["reg"]]
            dimnames(dat)[[2]][i] <- newname
        }
    }

    modelt <- terms(model)
    frmla <- paste(colnames(dat), collapse = " + ")
    frmlanew <- as.formula(paste( modelt[[2L]], "~" , frmla))

    mfnew <- model.frame(frmlanew, data=cbind(dat, model$model))

    mc <- model$call
    mc$formula <- frmlanew
    mc$data <- quote(mfnew)
    res <- eval(mc)

    class(res) <- c("rcreg", class(model))
    res$rcRegressions <- rcRegressions
    res
}
NULL

##' predict method for rcreg objects
##'
##' Calculates predicted values of
##' residual centered interaction regressions estimated in
##' any type of regression framework (lm, glm, etc).
##' @method predict rcreg
##' @export 
##' @rdname residualCenter
##' @example inst/examples/predict.rcreg-ex.R
##' @param object Fitted residual-centered regression from residualCenter
##' @param ... Other named arguments. May include newdata, a dataframe of 
##' predictors. That should include values for individual predictor, need 
##' not include interactions that are constructed by residualCenter. 
##' These parameters that will be passed to the predict method of the model.
predict.rcreg <-
    function (object, ...)
{
    if ( ! c("rcreg") %in% class(object) ) stop("predict.rcreg is intended for rcreg objects, which are created by residualCenter in the rockchalk package")

    dots <- list(...)
    newdata <- NULL
    if (!is.null(dots[["newdata"]])) {
        newdata <- dots[["newdata"]]
        dots[["newdata"]] <- NULL
    } else {
        newdata <- model.frame(object)
    }

    objectTerms <- terms(object)
    dvname <- names(attr(objectTerms, "dataClasses"))[1]
    rcRegs <- object$rcRegressions
    nams <- names(rcRegs)
    ## remove user's residual centered variables, if any
    ## for(i in nams)  newdata[[i]] <- NULL
    for(i in seq_along(rcRegs)){
        aReg <- rcRegs[[i]]
        prodVars <-  unlist(strsplit(nams[i], ".X."))
        predvals <-  predict.lm(aReg, newdata = newdata )
        actualProduct <- apply(newdata[ ,prodVars], 1, prod)
        myresids <-  actualProduct - predvals
        newdata[[nams[i]]] <- myresids
    }
    NextMethod(object, newdata = newdata, dots)
}
