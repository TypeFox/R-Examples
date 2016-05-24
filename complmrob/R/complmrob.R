#' MM-type estimators for linear regression on compositional data
#'
#' Uses the \code{\link[robustbase]{lmrob}} method for robust linear regression models to fit
#' a linear regression models to compositional data.
#'
#' The variables on the right-hand-side of the formula will be transformed with the isometric log-ratio
#' transformation (\code{\link{isomLR}}) and then the robust linear regression model is applied to
#' those transformed variables. The orthonormal basis can be constructed in \code{p} different ways,
#' where \code{p} is the number of variables on the RHS of the formula.
#'
#' To get an interpretable estimate of the regression coefficient for each part of the composition,
#' the data has to be transformed according to each of these orthonormal basis and a regression model
#' has to be fit to every transformed data set.
#'
#' @param formula The formula for the regression model
#' @param data The data.frame to use
#' @return A list of type \code{complmrob} with fields
#'      \describe{
#'          \item{coefficients}{the estimated coefficients}
#'          \item{models}{the single regression models (one for each orthonormal basis)}
#'          \item{npred}{the number of predictor variables}
#'          \item{predictors}{the names of the predictor variables}
#'          \item{coefind}{the index of the relevent coefficient in the single regression models}
#'          \item{call}{how the function was called}
#'          \item{intercept}{if an intercept is included}
#'      }
#' @import robustbase
#' @importFrom boot boot
#' @import parallel
#' @importFrom stats .lm.fit as.formula coef confint density fitted formula lm.wfit model.frame
#' @importFrom stats model.matrix model.response na.omit predict printCoefmat pt
#' @importFrom stats residuals sd terms update
#'
#' @references K. Hron, P. Filzmoser & K. Thompson (2012): Linear regression with compositional explanatory
#'      variables, Journal of Applied Statistics, DOI:10.1080/02664763.2011.644268
#' @export
#' @examples
#' data <- data.frame(lifeExp = state.x77[, "Life Exp"], USArrests[ , -3])
#' mUSArr <- complmrob(lifeExp ~ ., data = data)
#' summary(mUSArr)
#'
complmrob <- function(formula, data) {
    #
    # Initialize auxiliary variables
    #
    mf <- match.call(expand.dots = FALSE);

    m <- match(c("formula", "data"), names(mf), 0);
    mf <- mf[c(1, m)];
    mf$drop.unused.levels <- TRUE;
    mf[[1]] <- as.name("model.frame");
    mf <- eval(mf, parent.frame());
    mt <- attr(mf, "terms");
    y <- model.response(mf, "numeric");

    compPred <- attr(mt, "term.labels");
    npred <- length(compPred);

    int <- (attr(mt, "intercept") == 1);
    coefind <- 1L + int;

    #
    # Initialize return object
    #
    ret <- list(
        coefficients = numeric(npred),
        models = vector("list", npred),
        npred = npred,
        predictors = compPred,
        coefind = coefind,
        call = match.call(),
        intercept = int
    );

    class(ret) <- "complmrob";

    names(ret$coefficients) <- compPred;
    names(ret$models) <- compPred;

    for(predInd in seq_along(compPred)) {
        tmpPred <- isomLR(as.matrix(mf[ , compPred]), predInd);
        tmp <- data.frame(y = y, tmpPred);

        tmpFormula <- as.formula(sprintf(" y ~ %s", paste(colnames(tmpPred), collapse = " + ")))

        if(int == FALSE) {
            tmpFormula <- update(tmpFormula, . ~ . - 1)
        }

        m <- robustbase::lmrob(tmpFormula, data = tmp);

        ret$models[[predInd]] <- m;
        ret$coefficients[predInd] <- coef(m)[coefind];
    }

    # If the intercept is included in the model, evaluate it as well, but only once!
    # The variable m holds the last fit model!
    if(int == TRUE) {
        ret$coefficients <- c(coef(m)["(Intercept)"], ret$coefficients);
    }

    return(ret);
}