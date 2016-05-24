orlm.lm <- function (model, ui, ci = NULL, index = 2:length(coef(model)), 
    meq = 0, orig.out = FALSE, boot = FALSE, B = 1000, fixed = FALSE, 
    tol = sqrt(.Machine$double.eps), ...) 
{
    ## check model
    if (!("lm" %in% class(model))) 
        stop("ERROR: model must be of class lm.")
    if (any(c("glm", "mlm", "rlm") %in% class(model))) 
        stop("orlm does not work on classes glm, mlm or rlm.")
    ## preliminary calculations
    so <- summary(model)
    b <- coef(model)
    V <- so$cov.unscaled  ## works also, if sigma^2=0
    g <- length(b)
    ## check inputs
    if (!(is.vector(index))) 
        stop("index must be a vector.")
    if (1 %in% index & boot & !fixed) 
        stop("no restrictions on intercept possible when linear model is bootstrapped and fixed=FALSE.")
    if (is.vector(ui)) 
        ui <- matrix(ui, 1, length(ui))
    if (!is.matrix(ui)) 
        stop("ui must be a matrix.")
    if (!length(index) == ncol(ui)) 
        stop("mismatch between number of columns for ui and length of index")
    if (is.null(ci)) 
        ci <- rep(0, nrow(ui))
    if (!is.vector(ci)) 
        stop("ci must be a vector.")
    if (!nrow(ui) == length(ci)) 
        stop("mismatch between number of rows in ui and elements in ci")
    hilf <- RREF(t(ui))
    if (hilf$rank < nrow(ui)) 
        stop(paste("Matrix ui must have full row-rank (choose e.g. rows", 
            paste(hilf$pivot, collapse = " "), ")."))
    ## expand ui by 0 columns, if necessary 
    uiw <- ui
    if (length(index) > g | max(index) > g) 
        stop(paste("index must be vector of index positions, at most of length ", 
            g))
    uiw <- matrix(0, nrow(ui), g)
    uiw[, index] <- ui
    ## inequality restrictions only, all fulfilled
    if (all(uiw %*% b - ci >= 0 * ci) & meq == 0) 
        aus <- list(b.unrestr = b, b.restr = b, R2 = so$r.squared, 
            residuals = model$residuals, fitted.values = model$fitted.values, 
            weights = weights(model), orig.R2 = so$r.squared, 
            df.error = model$df.residual, s2 = so$sigma^2, Sigma = vcov(model), 
            origmodel = NULL, ui = ui, ci = ci, iact = NULL, 
            restr.index = index, meq = meq, bootout = NULL)
    else {
        ## equality restrictions involved or some inequality restrictions violated
        ## calculate restricted estimate
        aus <- solve.QP(Dmat = solve(V), dvec = solve(V, b), 
            Amat = t(uiw), bvec = ci, meq = meq)
        y <- model$model[, attr(model$terms, "response")]
        names(aus$solution) <- names(b)
        aus$solution[abs(aus$solution) < tol] <- 0
        ## initialize output list
        aus <- list(b.restr = aus$solution, b.unrestr = b, 
            R2 = NULL, residuals = NULL, fitted.values = NULL, 
            weights = weights(model), orig.R2 = so$r.squared, 
            df.error = model$df.residual, s2 = so$sigma^2, Sigma = vcov(model), 
            origmodel = NULL, ui = ui, ci = ci, iact = aus$iact, 
            restr.index = index, meq = meq, bootout = NULL)
        aus$fitted.values <- model.matrix(model) %*% aus$b.restr
        aus$residuals <- y - aus$fitted.values
        ### R2 for models with and without weights and with and without intercept
        aus$R2 <- 1 - sum(aus$residuals^2)/sum((y - mean(y))^2)
        if (is.null(weights(model)) & !attr(model$terms, "intercept")) 
            aus$R2 <- 1 - sum(aus$residuals^2)/sum(y^2)
        if (attr(model$terms, "intercept") & !is.null(weights(model))) 
            aus$R2 <- 1 - sum(weights(model) * aus$residuals^2)/sum(weights(model) * 
                (y - weighted.mean(y, w = aus$weights))^2)
        if (!(attr(model$terms, "intercept") | is.null(weights(model)))) 
            aus$R2 <- 1 - sum(weights(model) * aus$residuals^2)/sum(weights(model) * 
                y^2)
    }
    if (orig.out) 
        aus$origmodel <- model
    if (boot) 
        aus$bootout <- boot.orlm(model, B = B, fixed = fixed, 
            ui = ui, ci = ci, index = index, meq = meq)
    class(aus) <- c("orlm", "orest")
    aus
}
