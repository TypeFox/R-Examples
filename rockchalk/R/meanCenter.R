##' Estimate standardized regression coefficients for all variables
##'
##' This is brain-dead standardization of all variables in the design
##' matrix.  It mimics the silly output of SPSS, which standardizes
##' all regressors, even if they represent categorical variables.
##'
##' @param model a fitted lm object
##' @return an lm fitted with the standardized variables
##' @export
##' @author Paul Johnson <pauljohn@@ku.edu>
##' @seealso \code{\link[rockchalk]{meanCenter}} which will center or
##' re-scale only numberic variables
standardize <-
    function(model)
{
    UseMethod("standardize")
}

##' @return a standardized regression object
##' @rdname standardize
##' @export
##' @method standardize lm
##' @example inst/examples/standardize-ex.R
standardize.lm <-
    function(model)
{
    formulaReplace <- function(fmla, xname, newname){
        do.call("substitute", list(fmla, setNames(list(as.name(newname)), xname)))
    }

    mt <- terms(model)
    mdata <- model.frame(model)
    ys  <- drop(scale(mdata[, 1]))
    ##dm = design matrix, columns of predictors as numerically coded
    dm <- model.matrix(model)[ , -1, drop=FALSE] #no intercept
    dmnames <- paste(colnames(dm),"s", sep = "")
    dmnamesticked <- paste("`",dmnames,"`", sep = "")
    dmnamesticked <- gsub("``","`", dmnamesticked)
    dvname <- paste(colnames(mdata)[1],"s", sep = "")
    dvnameticked <-  paste("`", dvname,"`", sep = "")
    dvnameticked <- gsub("``","`", dvnameticked)
    std <- function(x) if(is.numeric(x)) scale(x) else x
    stddat <- apply(dm, 2, std)  ##standardize numeric vars
    colnames(stddat) <- paste(colnames(stddat), "s", sep = "")
    stddat <- cbind(ys, stddat)
    stddat <- as.data.frame(stddat)
    colnames(stddat) <- c(dvname, dmnames)
    colnames(stddat) <- gsub("`","", colnames(stddat))

    mc <- model$call
    mc$data <- quote(stddat)
    fmla <- paste(dvnameticked, " ~ ", "-1 + ", paste(dmnamesticked, collapse= " + "))
    mc$formula <- formula(fmla)
    res <- eval(mc)
    class(res) <- c("stdreg", class(model))
    res
}


##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @method summary stdreg
##' @export
summary.stdreg <-
    function(object, ...)
{
    dm <- model.matrix(object)
    dm <- dm[ , which(attr(dm, "assign") != 0)] #remove intercept, if any
    dm <- cbind( model.frame(object)[ , deparse(terms(object)[[2]])], dm)
    colnames(dm)[1] <- deparse(terms(object)[[2]])
    dmmeans <- apply(dm, 2, mean)
    dmstds <- apply(dm, 2, sd)
    summstat <- zapsmall(data.frame("mean" = dmmeans, "std.dev." = dmstds))
    summ <- NextMethod(generic = "summary", object = object, ...)
    summ$summstat <- summstat
    class(summ) <- paste("summary.", class(object), sep = "")
    summ
}
NULL


##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @method print stdreg
##' @export 
print.stdreg <- function(x, ...){
    cat("The standardized variables are suffixed with the letter \"s\" \n")
    NextMethod(generic = "print", object = x, ...)
}
NULL

##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @method print summary.stdreg
##' @export 
print.summary.stdreg <-
    function (x, ...)
{
    cat("All variables in the model matrix and the dependent variable
were centered. The centered variables have the letter \"s\" appended to their
non-centered counterparts, even constructed
variables like `x1:x2` and poly(x1,2). We agree, that's probably
ill-advised, but you asked for it by running standardize().\n
The rockchalk function meanCenter is a smarter option, probably. \n
The summary statistics of the variables in the design matrix. \n")
    print(x$summstat)
    NextMethod()
}
NULL





##' meanCenter selectively centers or standarizes variables in a regression model.
##'
##' Works with "lm" class objects, objects estimated by \code{glm()}. This
##' centers some or all of the the predictors and then re-fits the
##' original model with the new variables. This is a convenience to
##' researchers who are often urged to center their predictors.  This
##' is sometimes suggested as a way to ameliorate multi-collinearity
##' in models that include interaction terms (Aiken and West, 1991;
##' Cohen, et al 2002). Mean-centering may enhance interpretation of
##' the regression intercept, but it actually does not help with
##' multicollinearity.  (Echambadi and Hess, 2007). This function
##' facilitates comparison of mean-centered models with others by
##' calculating centered variables.  The defaults will cause a
##' regression's numeric interactive variables to be mean
##' centered. Variations on the arguments are discussed in details.
##'
##' Suppose the user's formula that fits the original model is
##' \code{m1 <- lm(y ~ x1*x2 + x3 + x4, data = dat)}. The fitted model
##' will include estimates for predictors \code{x1}, \code{x2},
##' \code{x1:x2}, \code{x3} and \code{x4}. By default,
##' \code{meanCenter(m1)} scans the output to see if there are
##' interaction terms of the form \code{x1:x2}. If so, then x1 and x2
##' are replaced by centered versions (m1-mean(m1)) and
##' (m2-mean(m2)). The model is re-estimated with those new variables.
##' model (the main effect and the interaction). The resulting thing
##' is "just another regression model", which can be analyzed or
##' plotted like any R regression object.
##'
##' The user can claim control over which variables are centered in
##' several ways. Most directly, by specifying a vector of variable
##' names, the user can claim direct control. For example, the
##' argument \code{terms=c("x1","x2","x3")} would cause 3 predictors
##' to be centered. If one wants all predictors to be centered, the
##' argument \code{centerOnlyInteractors} should be set to
##' FALSE. Please note, this WILL NOT center factor variables. But it
##' will find all numeric predictors and center them.
##'
##' The dependent variable will not be centered, unless the user
##' explicitly requests it by setting centerDV = TRUE.
##'
##' As an additional convenience to the user, the argument
##' \code{standardize = TRUE} can be used.  This will divide each
##' centered variable by its observed standard deviation. For people
##' who like standardized regression, I suggest this is a better
##' approach than the \code{standardize} function (which is brain-dead
##' in the style of SPSS). meanCenter with \code{standardize = TRUE}
##' will only try to standardize the numeric predictors.
##'
##' To be completely clear, I believe mean-centering is not helpful
##' with the multicollinearity problem. It doesn't help, it doesn't
##' hurt.  Only a misunderstanding leads its proponents to claim
##' otherwise. This is emphasized in the vignette "rockchalk" that is
##' distributed with this package.
##'
##' @title meanCenter
##' @param model a fitted regression model (presumably from lm)
##' @param centerOnlyInteractors Default TRUE. If FALSE, all numeric
##'     predictors in the regression data frame are centered before
##'     the regression is conducted.
##' @param centerDV Default FALSE. Should the dependent variable be
##'     centered? Do not set this option to TRUE unless the dependent
##'     variable is a numeric variable. Otherwise, it is an error.
##' @param standardize Default FALSE. Instead of simply mean-centering
##'     the variables, should they also be "standardized" by first
##'     mean-centering and then dividing by the estimated standard
##'     deviation.
##' @param terms Optional. A vector of variable names to be
##'     centered. Supplying this argument will stop meanCenter from
##'     searching for interaction terms that might need to be
##'     centered.
##' @export meanCenter
##' @rdname meanCenter
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @seealso
##'     \code{\link[rockchalk]{standardize}}
##'     \code{\link[rockchalk]{residualCenter}}
##' @references Aiken, L. S. and West, S.G. (1991). Multiple
##'     Regression: Testing and Interpreting Interactions. Newbury
##'     Park, Calif: Sage Publications.
##'
##' Cohen, J., Cohen, P., West, S. G., and Aiken, L. S. (2002). Applied
##' Multiple Regression/Correlation Analysis for the Behavioral Sciences
##' (Third.). Routledge Academic.
##'
##' Echambadi, R., and Hess, J. D. (2007). Mean-Centering Does Not Alleviate
##' Collinearity Problems in Moderated Multiple Regression Models.
##' Marketing Science, 26(3), 438-445.
##' @example inst/examples/meanCenter-ex.R
meanCenter <-
    function(model, centerOnlyInteractors = TRUE, centerDV = FALSE,
             standardize=FALSE, terms = NULL)
{
    UseMethod("meanCenter")
}

##' @return A regression model of the same type as the input model,
##' with attributes representing the names of the centered variables.
##' @rdname meanCenter
##' @export
##' @method meanCenter default
##' @export 
meanCenter.default <-
    function(model, centerOnlyInteractors = TRUE, centerDV = FALSE,
             standardize = FALSE, terms = NULL)
{
    std <- function(x) {
        if(!is.numeric(x)) stop("can't center a factor variable. No Can Do!")
        xmean <- mean(x, na.rm = TRUE)
        if (standardize) {
            xsd <- sd(x, na.rm = TRUE)
        } else {
            xsd <- 1
        }
        x <- (x-xmean)/xsd
        list(x = x, xmean = xmean, xsd = xsd)
    }

    ## rdf <- get_all_vars(formula(model), model$model) #raw data frame
    rdf <- model.data(model)
    t <- terms(model)
    ## TODO 20140417: look at using na.action attribute of model to be more delicate here.
    tl <- attr(t, "term.labels")
    tmdc <- attr(t, "dataClasses") ##term model data classes

    isNumeric <- names(tmdc)[ which(tmdc %in% c("numeric"))]
    isFac <-  names(tmdc)[ which(tmdc %in% c("factor"))]

    if (centerDV & tmdc[1] != "numeric")
        stop("Sorry, the DV is not a numeric column, it does not make sense to center it.")

    ##Build "nc", a vector of variable names that "need centering"
    ##
    if (!centerDV) {
        if (!is.null(terms)){
            nc <- as.vector(terms)
            nc <- unique(nc)
        } else if (centerOnlyInteractors == FALSE){
            nc <- isNumeric[-1] #-1 excludes response
            nc <- unique(nc)
        } else {
            interactTerms <- tl[grep(":", tl)]
            nc <- unique(unlist(strsplit( interactTerms, ":")))
            nc <-  nc[which(nc %in% isNumeric)]
        }
    } else {
        if (!is.null(terms)){
            nc <- as.vector(terms)
            nc <- c(names(tmdc)[1] , nc)
        } else if (centerOnlyInteractors == FALSE){
            nc <- isNumeric
        } else {
            interactTerms <- tl[grep(":", tl)]
            nc <- unique(unlist(strsplit( interactTerms, ":")))
            nc <- nc[which(nc %in% isNumeric)]
            nc <- c(names(tmdc)[1] , nc)
        }
    }

    mc <- model$call
    ## run same model call, replacing non centered data with centered data.
    ##
    stddat <- rdf
    centeredVars <- matrix(NA, nrow=2, ncol=length(nc))
    colnames(centeredVars) <- nc
    rownames(centeredVars) <- c("mean","scale")

    formulaReplace <- function(fmla, xname, newname){
        do.call("substitute", list(fmla, setNames(list(as.name(newname)), xname)))
    }

    newFmla <- mc$formula
    for (i in seq_along(nc)){
        icenter <- std(stddat[, nc[i]])
        centeredVars[1, nc[i]] <- icenter$xmean
        centeredVars[2, nc[i]] <- icenter$xsd
        newname <- paste(as.character(nc[i]), "c", sep = "")
        if (isTRUE(standardize)) newname <- paste(newname, "s", sep = "")
        stddat[ ,newname] <- icenter$x
        newFmla <- formulaReplace(newFmla,  as.character(nc[i]), newname)
        nc[i] <- newname
    }
    colnames(centeredVars) <- nc
    mc$formula <- newFmla
    mc$data <- quote(stddat)
    res <- eval(mc)
    class(res) <- c("mcreg", class(model))
    attr(res, "centeredVars") <- centeredVars
    attr(res, "centerCall") <-  match.call()
    res
}

##' @author <pauljohn@@ku.edu>
##' @export 
##' @method summary mcreg
summary.mcreg <-
    function(object, ...)
{
    centeredVars <- attr(object, "centeredVars")
    dm <- model.matrix(object)
    dm <- dm[ , which(attr(dm, "assign") != 0), drop=FALSE] #remove intercept, if any
    dm <- cbind( model.frame(object)[ , deparse(terms(object)[[2]])], dm)
    colnames(dm)[1] <- deparse(terms(object)[[2]])
    dmmeans <- apply(dm, 2, mean)
    dmstds <- apply(dm, 2, sd)
    summstat <- zapsmall(data.frame("mean" = dmmeans, "std.dev." = dmstds))
    summ <- NextMethod(generic = "summary", object = object, ...)
    summ$summstat <- summstat
    summ$centeredVars <- centeredVars
    class(summ) <- paste("summary.", class(object), sep="")
    summ$mc <- attr(object, "centerCall")
    summ
}
NULL

##' @author <pauljohn@@ku.edu>
##' @method print mcreg
##' @export 
print.mcreg <- function(x, ...){
    centeredVars <- attr(x, "centeredVars")
    cat("The centered variables are: \n")
    print(centeredVars)
    mc <- attr(x, "centerCall")
    cat("The call that requested centering was: \n")
    print(mc)
    NextMethod(generic = "print", object = x, ...)
}
NULL


##' @author <pauljohn@@ku.edu>
##' @method print summary.mcreg
##' @export 
print.summary.mcreg <-
    function (x, ...)
{
    cat("These variables were mean-centered before any transformations were made on the design matrix.\n")
    print(colnames(x$centeredVars))
    cat("The centers and scale factors were \n")
    print(x$centeredVars)

    cat("The summary statistics of the variables in the design matrix (after centering). \n")
    print(x$summstat)
    cat("\nThe following results were produced from: \n")
    print(x$mc)
    ##NextMethod(generic = "print", x = x, ...)
    NextMethod()
}
NULL


##' @author <pauljohn@@ku.edu>
##' @method predict mcreg
##' @export
predict.mcreg <-
    function (object, ...)
{
    originalCall <- object$call
    dots <- list(...)
    newdata <- NULL
    if (! is.null(dots[["newdata"]])) {
        newdata <- dots[["newdata"]]
        dots[["newdata"]] <- NULL
    }
    centeredVars <- attr(object, "centeredVars")
    nc <- colnames(centeredVars) #need centering
    dvname <- parse(text = formula(originalCall)[[2]])
    nc <- setdiff(nc, dvname) #remove dv name if present
    if (is.null(newdata)) {
        newdata <- model.frame(object)##should be centered already
        tmeans <- sapply(newdata[ , nc, drop = FALSE], mean, na.rm = TRUE)
        if (! isTRUE(all.equal(abs(tmeans), rep(0, length(nc)), check.attributes = FALSE)))
            stop(paste("In predict.mcreg, the fitted regression claims to have centered variables,",
                       paste(nc, collapse=" "), "but not all of those centered values have",
                       "observed means very close to 0. Something's wrong"))
    }
    NextMethod(object, newdata = newdata, dots)
}
NULL


##' Find numeric columns, center them, re-name them, and join them with the original data.
##'
##' The meanCentered regression function requires centered-inputs when
##' calculations are predicted. For comparison with ordinary
##' regression, it is convenient to have both centered and the
##' original data side-by-side.  This function handles that.  If the
##' input data has columns c("x1","x2","x3"), then the centered result
##' will have columns c("x1","x2","x3","x1c","x2c","x3c"), where "c"
##' indicates "mean-centered". If standardize=TRUE, then the result
##' will have columns c("x1","x2","x3","x1cs","x2cs","x3cs"), where "cs"
##' indicate "centered and scaled".
##' @param data Required. data frame or matrix.
##' @param center Optional. If nc is NOT supplied, then all numeric columns
##' in data will be centered (possiblly scaled).  Can be specified in 2 formats. 1) Vector of column names that are to be centered, 2) Vector named elements giving values of means to be used in centering.  Values must be named, as in c("x1" = 17, "x2" = 44).
##' (possibly scaled).
##' @param standardize  Default FALSE. If TRUE, the variables are
##' first mean-centered, and then divided by their standard deviations
##' (scaled). User can supply a named vector of scale values by which
##' to divide each variable (otherwise sd is used). Vector must have same
##' names and length as center argument. Variables can be entered in any order (will be resorted inside function).
##' @return A data frame with 1) All original columns 2) additional
##' columns with centered/scaled data, variables renamed "c" or "cs"
##' to indicate the data is centered or centered and
##' scaled. Attributes "centers" and "scales" are created for "record
##' keeping" on centering and scaling values.
##' @author <pauljohn@@ku.edu>
##' @export
##' @examples
##' set.seed(12345)
##' dat <- data.frame(x1=rnorm(100, m = 50), x2 = rnorm(100, m = 50),
##'     x3 = rnorm(100, m = 50), y = rnorm(100),
##'     x4 = gl(2, 50, labels = c("Male","Female")))
##' datc1 <- centerNumerics(dat)
##' head(datc1)
##' summarize(datc1)
##' datc2 <- centerNumerics(dat, center=c("x1", "x2"))
##' head(datc2)
##' summarize(datc2)
##' attributes(datc2)
##' datc3 <- centerNumerics(dat, center = c("x1" = 30, "x2" = 40))
##' head(datc3)
##' summarize(datc3)
##' attributes(datc3)
##' datc4 <- centerNumerics(dat, center=c("x1", "x2"), standardize = TRUE)
##' head(datc3)
##' summarize(datc4)
##' attributes(datc4)
##' datc5 <- centerNumerics(dat, center=c("x1"=30, "x2"=40),
##' standardize = c("x2" = 5, "x1" = 7))
##' head(datc5)
##' summarize(datc5)
##' attributes(datc5)
centerNumerics <- function(data, center, standardize = FALSE){
    if (!is.data.frame(data))
        data <- as.data.frame(data)
    isN <- sapply(data, is.numeric)
    if (sum(isN) == 0)  return(data)

    if (missing(center)) {
        center <- TRUE
        nc <- colnames(data)[isN] ##all need centering
    } else if (is.character(center)){
        if(sum (!center %in% colnames(data)[isN]) != 0)
            stop(paste("centerNumerics failed. Argument center includes",
                       "column names that are not numeric variables in the"
                       , "data frame", deparse(substitute(data))))
        nc <- center
        center <- TRUE
    } else if (is.numeric(center)) {
        if(is.null(names(center)))
            stop(paste("centerNumerics failed.",
                       "The center vector must be a named vector so that",
                       "centerNumerics can decide which columns need centering"))
        if(sum (!names(center) %in% colnames(data)[isN]) != 0)
            stop(paste("centerNumerics failed. Argument center includes column names that are not numeric variables in the data frame", deparse(substitute(data))))
        nc <- names(center) ##names
    }

    if (is.numeric(standardize)) {
        if (!setequal(names(standardize), nc))
            stop("centerNumerics failed. Names of standardize argument must be identical to names of center argument")
        standardize <- standardize[nc] ##sorts
    }

    datas <- scale(data[ , nc], center = center, scale = standardize)
    if(!is.null(attr(datas, "scaled:center"))) centers <- attr(datas, "scaled:center")
    colnames(datas) <- paste(colnames(datas), "c", sep="")
    if (!is.null(attr(datas, "scaled:scale"))) {
        scales <- attr(datas, "scaled:scale")
        colnames(datas) <- paste(colnames(datas), "s", sep = "")
    }
    data <- as.data.frame(cbind(data, datas))
    attr(data, "centers") <- centers
    if (!is.null(attr(datas, "scaled:scale")))attr(data, "scales") <- scales

    data
}
NULL

