neuralnet <-
function (formula, data, hidden = 1, threshold = 0.01, stepmax = 1e+05, 
    rep = 1, startweights = NULL, learningrate.limit = NULL, 
    learningrate.factor = list(minus = 0.5, plus = 1.2), learningrate = NULL, 
    lifesign = "none", lifesign.step = 1000, algorithm = "rprop+", 
    err.fct = "sse", act.fct = "logistic", linear.output = TRUE, 
    exclude = NULL, constant.weights = NULL, likelihood = FALSE) 
{
    call <- match.call()
    options(scipen = 100, digits = 10)
    result <- varify.variables(data, formula, startweights, learningrate.limit, 
        learningrate.factor, learningrate, lifesign, algorithm, 
        threshold, lifesign.step, hidden, rep, stepmax, err.fct, 
        act.fct)
    data <- result$data
    formula <- result$formula
    startweights <- result$startweights
    learningrate.limit <- result$learningrate.limit
    learningrate.factor <- result$learningrate.factor
    learningrate.bp <- result$learningrate.bp
    lifesign <- result$lifesign
    algorithm <- result$algorithm
    threshold <- result$threshold
    lifesign.step <- result$lifesign.step
    hidden <- result$hidden
    rep <- result$rep
    stepmax <- result$stepmax
    model.list <- result$model.list
    matrix <- NULL
    list.result <- NULL
    result <- generate.initial.variables(data, model.list, hidden, 
        act.fct, err.fct, algorithm, linear.output, formula)
    covariate <- result$covariate
    response <- result$response
    err.fct <- result$err.fct
    err.deriv.fct <- result$err.deriv.fct
    act.fct <- result$act.fct
    act.deriv.fct <- result$act.deriv.fct
    for (i in 1:rep) {
        if (lifesign != "none") {
            lifesign <- display(hidden, threshold, rep, i, lifesign)
        }
        flush.console()
        result <- calculate.neuralnet(learningrate.limit = learningrate.limit, 
            learningrate.factor = learningrate.factor, covariate = covariate, 
            response = response, data = data, model.list = model.list, 
            threshold = threshold, lifesign.step = lifesign.step, 
            stepmax = stepmax, hidden = hidden, lifesign = lifesign, 
            startweights = startweights, algorithm = algorithm, 
            err.fct = err.fct, err.deriv.fct = err.deriv.fct, 
            act.fct = act.fct, act.deriv.fct = act.deriv.fct, 
            rep = i, linear.output = linear.output, exclude = exclude, 
            constant.weights = constant.weights, likelihood = likelihood, 
            learningrate.bp = learningrate.bp)
        if (!is.null(result$output.vector)) {
            list.result <- c(list.result, list(result))
            matrix <- cbind(matrix, result$output.vector)
        }
    }
    flush.console()
    if (!is.null(matrix)) {
        weight.count <- length(unlist(list.result[[1]]$weights)) - 
            length(exclude) + length(constant.weights) - sum(constant.weights == 
            0)
        if (!is.null(startweights) && length(startweights) < 
            (rep * weight.count)) {
            warning("some weights were randomly generated, because 'startweights' did not contain enough values", 
                call. = F)
        }
        ncol.matrix <- ncol(matrix)
    }
    else ncol.matrix <- 0
    if (ncol.matrix < rep) 
        warning(sprintf("algorithm did not converge in %s of %s repetition(s) within the stepmax", 
            (rep - ncol.matrix), rep), call. = FALSE)
    nn <- generate.output(covariate, call, rep, threshold, matrix, 
        startweights, model.list, response, err.fct, act.fct, 
        data, list.result, linear.output, exclude)
    return(nn)
}
varify.variables <-
function (data, formula, startweights, learningrate.limit, learningrate.factor, 
    learningrate.bp, lifesign, algorithm, threshold, lifesign.step, 
    hidden, rep, stepmax, err.fct, act.fct) 
{
    if (is.null(data)) 
        stop("'data' is missing", call. = FALSE)
    if (is.null(formula)) 
        stop("'formula' is missing", call. = FALSE)
    if (!is.null(startweights)) {
        startweights <- as.vector(unlist(startweights))
        if (any(is.na(startweights))) 
            startweights <- startweights[!is.na(startweights)]
    }
    data <- as.data.frame(data)
    formula <- as.formula(formula)
    model.vars <- attr(terms(formula), "term.labels")
    formula.reverse <- formula
    formula.reverse[[3]] <- formula[[2]]
    model.resp <- attr(terms(formula.reverse), "term.labels")
    model.list <- list(response = model.resp, variables = model.vars)
    if (!is.null(learningrate.limit)) {
        if (length(learningrate.limit) != 2) 
            stop("'learningrate.factor' must consist of two components", 
                call. = FALSE)
        learningrate.limit <- as.list(learningrate.limit)
        names(learningrate.limit) <- c("min", "max")
        learningrate.limit$min <- as.vector(as.numeric(learningrate.limit$min))
        learningrate.limit$max <- as.vector(as.numeric(learningrate.limit$max))
        if (is.na(learningrate.limit$min) || is.na(learningrate.limit$max)) 
            stop("'learningrate.limit' must be a numeric vector", 
                call. = FALSE)
    }
    if (!is.null(learningrate.factor)) {
        if (length(learningrate.factor) != 2) 
            stop("'learningrate.factor' must consist of two components", 
                call. = FALSE)
        learningrate.factor <- as.list(learningrate.factor)
        names(learningrate.factor) <- c("minus", "plus")
        learningrate.factor$minus <- as.vector(as.numeric(learningrate.factor$minus))
        learningrate.factor$plus <- as.vector(as.numeric(learningrate.factor$plus))
        if (is.na(learningrate.factor$minus) || is.na(learningrate.factor$plus)) 
            stop("'learningrate.factor' must be a numeric vector", 
                call. = FALSE)
    }
    else learningrate.factor <- list(minus = c(0.5), plus = c(1.2))
    if (is.null(lifesign)) 
        lifesign <- "none"
    lifesign <- as.character(lifesign)
    if (!((lifesign == "none") || (lifesign == "minimal") || 
        (lifesign == "full"))) 
        lifesign <- "minimal"
    if (is.na(lifesign)) 
        stop("'lifesign' must be a character", call. = FALSE)
    if (is.null(algorithm)) 
        algorithm <- "rprop+"
    algorithm <- as.character(algorithm)
    if (!((algorithm == "rprop+") || (algorithm == "rprop-") || 
        (algorithm == "slr") || (algorithm == "sag") || (algorithm == 
        "backprop"))) 
        stop("'algorithm' is not known", call. = FALSE)
    if (is.null(threshold)) 
        threshold <- 0.01
    threshold <- as.numeric(threshold)
    if (is.na(threshold)) 
        stop("'threshold' must be a numeric value", call. = FALSE)
    if (algorithm == "backprop") 
        if (is.null(learningrate.bp) || !is.numeric(learningrate.bp)) 
            stop("'learningrate' must be a numeric value, if the backpropagation algorithm is used", 
                call. = FALSE)
    if (is.null(lifesign.step)) 
        lifesign.step <- 1000
    lifesign.step <- as.integer(lifesign.step)
    if (is.na(lifesign.step)) 
        stop("'lifesign.step' must be an integer", call. = FALSE)
    if (lifesign.step < 1) 
        lifesign.step <- as.integer(100)
    if (is.null(hidden)) 
        hidden <- 0
    hidden <- as.vector(as.integer(hidden))
    if (prod(!is.na(hidden)) == 0) 
        stop("'hidden' must be an integer vector or a single integer", 
            call. = FALSE)
    if (length(hidden) > 1 && prod(hidden) == 0) 
        stop("'hidden' contains at least one 0", call. = FALSE)
    if (is.null(rep)) 
        rep <- 1
    rep <- as.integer(rep)
    if (is.na(rep)) 
        stop("'rep' must be an integer", call. = FALSE)
    if (is.null(stepmax)) 
        stepmax <- 10000
    stepmax <- as.integer(stepmax)
    if (is.na(stepmax)) 
        stop("'stepmax' must be an integer", call. = FALSE)
    if (stepmax < 1) 
        stepmax <- as.integer(1000)
    if (is.null(hidden)) {
        if (is.null(learningrate.limit)) 
            learningrate.limit <- list(min = c(1e-08), max = c(50))
    }
    else {
        if (is.null(learningrate.limit)) 
            learningrate.limit <- list(min = c(1e-10), max = c(0.1))
    }
    if (!is.function(act.fct) && act.fct != "logistic" && act.fct != 
        "tanh") 
        stop("''act.fct' is not known", call. = FALSE)
    if (!is.function(err.fct) && err.fct != "sse" && err.fct != 
        "ce") 
        stop("'err.fct' is not known", call. = FALSE)
    return(list(data = data, formula = formula, startweights = startweights, 
        learningrate.limit = learningrate.limit, learningrate.factor = learningrate.factor, 
        learningrate.bp = learningrate.bp, lifesign = lifesign, 
        algorithm = algorithm, threshold = threshold, lifesign.step = lifesign.step, 
        hidden = hidden, rep = rep, stepmax = stepmax, model.list = model.list))
}
generate.initial.variables <-
function (data, model.list, hidden, act.fct, err.fct, algorithm, 
    linear.output, formula) 
{
    formula.reverse <- formula
    formula.reverse[[2]] <- as.formula(paste(model.list$response[[1]], 
        "~", model.list$variables[[1]], sep = ""))[[2]]
    formula.reverse[[3]] <- formula[[2]]
    response <- as.matrix(model.frame(formula.reverse, data))
    formula.reverse[[3]] <- formula[[3]]
    covariate <- as.matrix(model.frame(formula.reverse, data))
    covariate[, 1] <- 1
    colnames(covariate)[1] <- "intercept"
    if (is.function(act.fct)) {
        act.deriv.fct <- differentiate(act.fct)
        attr(act.fct, "type") <- "function"
    }
    else {
        if (act.fct == "tanh") {
            act.fct <- function(x) {
                tanh(x)
            }
            attr(act.fct, "type") <- "tanh"
            act.deriv.fct <- function(x) {
                1 - x^2
            }
        }
        else if (act.fct == "logistic") {
            act.fct <- function(x) {
                1/(1 + exp(-x))
            }
            attr(act.fct, "type") <- "logistic"
            act.deriv.fct <- function(x) {
                x * (1 - x)
            }
        }
    }
    if (is.function(err.fct)) {
        err.deriv.fct <- differentiate(err.fct)
        attr(err.fct, "type") <- "function"
    }
    else {
        if (err.fct == "ce") {
            if (all(response == 0 | response == 1)) {
                err.fct <- function(x, y) {
                  -(y * log(x) + (1 - y) * log(1 - x))
                }
                attr(err.fct, "type") <- "ce"
                err.deriv.fct <- function(x, y) {
                  (1 - y)/(1 - x) - y/x
                }
            }
            else {
                err.fct <- function(x, y) {
                  1/2 * (y - x)^2
                }
                attr(err.fct, "type") <- "sse"
                err.deriv.fct <- function(x, y) {
                  x - y
                }
                warning("'err.fct' was automatically set to sum of squared error (sse), because the response is not binary", 
                  call. = F)
            }
        }
        else if (err.fct == "sse") {
            err.fct <- function(x, y) {
                1/2 * (y - x)^2
            }
            attr(err.fct, "type") <- "sse"
            err.deriv.fct <- function(x, y) {
                x - y
            }
        }
    }
    return(list(covariate = covariate, response = response, err.fct = err.fct, 
        err.deriv.fct = err.deriv.fct, act.fct = act.fct, act.deriv.fct = act.deriv.fct))
}
differentiate <-
function (orig.fct, hessian = FALSE) 
{
    body.fct <- deparse(body(orig.fct))
    if (body.fct[1] == "{") 
        body.fct <- body.fct[2]
    text <- paste("y~", body.fct, sep = "")
    text2 <- paste(deparse(orig.fct)[1], "{}")
    temp <- deriv(eval(parse(text = text)), "x", func = eval(parse(text = text2)), 
        hessian = hessian)
    temp <- deparse(temp)
    derivative <- NULL
    if (!hessian) 
        for (i in 1:length(temp)) {
            if (!any(grep("value", temp[i]))) 
                derivative <- c(derivative, temp[i])
        }
    else for (i in 1:length(temp)) {
        if (!any(grep("value", temp[i]), grep("grad", temp[i]), 
            grep(", c", temp[i]))) 
            derivative <- c(derivative, temp[i])
    }
    number <- NULL
    for (i in 1:length(derivative)) {
        if (any(grep("<-", derivative[i]))) 
            number <- i
    }
    if (is.null(number)) {
        return(function(x) {
            matrix(0, nrow(x), ncol(x))
        })
    }
    else {
        derivative[number] <- unlist(strsplit(derivative[number], 
            "<-"))[2]
        derivative <- eval(parse(text = derivative))
    }
    if (length(formals(derivative)) == 1 && length(derivative(c(1, 
        1))) == 1) 
        derivative <- eval(parse(text = paste("function(x){matrix(", 
            derivative(1), ", nrow(x), ncol(x))}")))
    if (length(formals(derivative)) == 2 && length(derivative(c(1, 
        1), c(1, 1))) == 1) 
        derivative <- eval(parse(text = paste("function(x, y){matrix(", 
            derivative(1, 1), ", nrow(x), ncol(x))}")))
    return(derivative)
}
display <-
function (hidden, threshold, rep, i.rep, lifesign) 
{
    text <- paste("    rep: %", nchar(rep) - nchar(i.rep), "s", 
        sep = "")
    cat("hidden: ", paste(hidden, collapse = ", "), "    thresh: ", 
        threshold, sprintf(eval(expression(text)), ""), i.rep, 
        "/", rep, "    steps: ", sep = "")
    if (lifesign == "full") 
        lifesign <- sum(nchar(hidden)) + 2 * length(hidden) - 
            2 + max(nchar(threshold)) + 2 * nchar(rep) + 41
    return(lifesign)
}
calculate.neuralnet <-
function (data, model.list, hidden, stepmax, rep, threshold, 
    learningrate.limit, learningrate.factor, lifesign, covariate, 
    response, lifesign.step, startweights, algorithm, act.fct, 
    act.deriv.fct, err.fct, err.deriv.fct, linear.output, likelihood, 
    exclude, constant.weights, learningrate.bp) 
{
    time.start.local <- Sys.time()
    result <- generate.startweights(model.list, hidden, startweights, 
        rep, exclude, constant.weights)
    weights <- result$weights
    exclude <- result$exclude
    nrow.weights <- sapply(weights, nrow)
    ncol.weights <- sapply(weights, ncol)
    result <- rprop(weights = weights, threshold = threshold, 
        response = response, covariate = covariate, learningrate.limit = learningrate.limit, 
        learningrate.factor = learningrate.factor, stepmax = stepmax, 
        lifesign = lifesign, lifesign.step = lifesign.step, act.fct = act.fct, 
        act.deriv.fct = act.deriv.fct, err.fct = err.fct, err.deriv.fct = err.deriv.fct, 
        algorithm = algorithm, linear.output = linear.output, 
        exclude = exclude, learningrate.bp = learningrate.bp)
    startweights <- weights
    weights <- result$weights
    step <- result$step
    reached.threshold <- result$reached.threshold
    net.result <- result$net.result
    error <- sum(err.fct(net.result, response))
    if (is.na(error) & type(err.fct) == "ce") 
        if (all(net.result <= 1, net.result >= 0)) 
            error <- sum(err.fct(net.result, response), na.rm = T)
    if (!is.null(constant.weights) && any(constant.weights != 
        0)) 
        exclude <- exclude[-which(constant.weights != 0)]
    if (length(exclude) == 0) 
        exclude <- NULL
    aic <- NULL
    bic <- NULL
    if (likelihood) {
        synapse.count <- length(unlist(weights)) - length(exclude)
        aic <- 2 * error + (2 * synapse.count)
        bic <- 2 * error + log(nrow(response)) * synapse.count
    }
    if (is.na(error)) 
        warning("'err.fct' does not fit 'data' or 'act.fct'", 
            call. = F)
    if (lifesign != "none") {
        if (reached.threshold <= threshold) {
            cat(rep(" ", (max(nchar(stepmax), nchar("stepmax")) - 
                nchar(step))), step, sep = "")
            cat("\terror: ", round(error, 5), rep(" ", 6 - (nchar(round(error, 
                5)) - nchar(round(error, 0)))), sep = "")
            if (!is.null(aic)) {
                cat("\taic: ", round(aic, 5), rep(" ", 6 - (nchar(round(aic, 
                  5)) - nchar(round(aic, 0)))), sep = "")
            }
            if (!is.null(bic)) {
                cat("\tbic: ", round(bic, 5), rep(" ", 6 - (nchar(round(bic, 
                  5)) - nchar(round(bic, 0)))), sep = "")
            }
            time <- difftime(Sys.time(), time.start.local)
            cat("\ttime: ", round(time, 2), " ", attr(time, "units"), 
                sep = "")
            cat("\n")
        }
    }
    if (reached.threshold > threshold) 
        return(result = list(output.vector = NULL, weights = NULL))
    output.vector <- c(error = error, reached.threshold = reached.threshold, 
        steps = step)
    if (!is.null(aic)) {
        output.vector <- c(output.vector, aic = aic)
    }
    if (!is.null(bic)) {
        output.vector <- c(output.vector, bic = bic)
    }
    for (w in 1:length(weights)) output.vector <- c(output.vector, 
        as.vector(weights[[w]]))
    generalized.weights <- calculate.generalized.weights(weights, 
        neuron.deriv = result$neuron.deriv, net.result = net.result)
    startweights <- unlist(startweights)
    weights <- unlist(weights)
    if (!is.null(exclude)) {
        startweights[exclude] <- NA
        weights[exclude] <- NA
    }
    startweights <- relist(startweights, nrow.weights, ncol.weights)
    weights <- relist(weights, nrow.weights, ncol.weights)
    return(list(generalized.weights = generalized.weights, weights = weights, 
        startweights = startweights, net.result = result$net.result, 
        output.vector = output.vector))
}
generate.startweights <-
function (model.list, hidden, startweights, rep, exclude, constant.weights) 
{
    input.count <- length(model.list$variables)
    output.count <- length(model.list$response)
    if (!(length(hidden) == 1 && hidden == 0)) {
        length.weights <- length(hidden) + 1
        nrow.weights <- array(0, dim = c(length.weights))
        ncol.weights <- array(0, dim = c(length.weights))
        nrow.weights[1] <- (input.count + 1)
        ncol.weights[1] <- hidden[1]
        if (length(hidden) > 1) 
            for (i in 2:length(hidden)) {
                nrow.weights[i] <- hidden[i - 1] + 1
                ncol.weights[i] <- hidden[i]
            }
        nrow.weights[length.weights] <- hidden[length.weights - 
            1] + 1
        ncol.weights[length.weights] <- output.count
    }
    else {
        length.weights <- 1
        nrow.weights <- array((input.count + 1), dim = c(1))
        ncol.weights <- array(output.count, dim = c(1))
    }
    length <- sum(ncol.weights * nrow.weights)
    vector <- rep(0, length)
    if (!is.null(exclude)) {
        if (is.matrix(exclude)) {
            exclude <- matrix(as.integer(exclude), ncol = ncol(exclude), 
                nrow = nrow(exclude))
            if (nrow(exclude) >= length || ncol(exclude) != 3) 
                stop("'exclude' has wrong dimensions", call. = FALSE)
            if (any(exclude < 1)) 
                stop("'exclude' contains at least one invalid weight", 
                  call. = FALSE)
            temp <- relist(vector, nrow.weights, ncol.weights)
            for (i in 1:nrow(exclude)) {
                if (exclude[i, 1] > length.weights || exclude[i, 
                  2] > nrow.weights[exclude[i, 1]] || exclude[i, 
                  3] > ncol.weights[exclude[i, 1]]) 
                  stop("'exclude' contains at least one invalid weight", 
                    call. = FALSE)
                temp[[exclude[i, 1]]][exclude[i, 2], exclude[i, 
                  3]] <- 1
            }
            exclude <- which(unlist(temp) == 1)
        }
        else if (is.vector(exclude)) {
            exclude <- as.integer(exclude)
            if (max(exclude) > length || min(exclude) < 1) {
                stop("'exclude' contains at least one invalid weight", 
                  call. = FALSE)
            }
        }
        else {
            stop("'exclude' must be a vector or matrix", call. = FALSE)
        }
        if (length(exclude) >= length) 
            stop("all weights are exluded", call. = FALSE)
    }
    length <- length - length(exclude)
    if (!is.null(exclude)) {
        if (is.null(startweights) || length(startweights) < (length * 
            rep)) 
            vector[-exclude] <- rnorm(length)
        else vector[-exclude] <- startweights[((rep - 1) * length + 
            1):(length * rep)]
    }
    else {
        if (is.null(startweights) || length(startweights) < (length * 
            rep)) 
            vector <- rnorm(length)
        else vector <- startweights[((rep - 1) * length + 1):(length * 
            rep)]
    }
    if (!is.null(exclude) && !is.null(constant.weights)) {
        if (length(exclude) < length(constant.weights)) 
            stop("constant.weights contains more weights than exclude", 
                call. = FALSE)
        else vector[exclude[1:length(constant.weights)]] <- constant.weights
    }
    weights <- relist(vector, nrow.weights, ncol.weights)
    return(list(weights = weights, exclude = exclude))
}
rprop <-
function (weights, response, covariate, threshold, learningrate.limit, 
    learningrate.factor, stepmax, lifesign, lifesign.step, act.fct, 
    act.deriv.fct, err.fct, err.deriv.fct, algorithm, linear.output, 
    exclude, learningrate.bp) 
{
    step <- 1
    nchar.stepmax <- max(nchar(stepmax), 7)
    length.weights <- length(weights)
    nrow.weights <- sapply(weights, nrow)
    ncol.weights <- sapply(weights, ncol)
    length.unlist <- length(unlist(weights)) - length(exclude)
    learningrate <- as.vector(matrix(0.1, nrow = 1, ncol = length.unlist))
    gradients.old <- as.vector(matrix(0, nrow = 1, ncol = length.unlist))
    if (is.null(exclude)) 
        exclude <- length(unlist(weights)) + 1
    if (type(act.fct) == "tanh" || type(act.fct) == "logistic") 
        special <- TRUE
    else special <- FALSE
    if (linear.output) {
        output.act.fct <- function(x) {
            x
        }
        output.act.deriv.fct <- function(x) {
            matrix(1, nrow(x), ncol(x))
        }
    }
    else {
        if (type(err.fct) == "ce" && type(act.fct) == "logistic") {
            err.deriv.fct <- function(x, y) {
                x * (1 - y) - y * (1 - x)
            }
            linear.output <- TRUE
        }
        output.act.fct <- act.fct
        output.act.deriv.fct <- act.deriv.fct
    }
    result <- compute.net(weights, length.weights, covariate = covariate, 
        act.fct = act.fct, act.deriv.fct = act.deriv.fct, output.act.fct = output.act.fct, 
        output.act.deriv.fct = output.act.deriv.fct, special)
    err.deriv <- err.deriv.fct(result$net.result, response)
    gradients <- calculate.gradients(weights = weights, length.weights = length.weights, 
        neurons = result$neurons, neuron.deriv = result$neuron.deriv, 
        err.deriv = err.deriv, exclude = exclude, linear.output = linear.output)
    reached.threshold <- max(abs(gradients))
    min.reached.threshold <- reached.threshold
    while (step < stepmax && reached.threshold > threshold) {
        if (!is.character(lifesign) && step%%lifesign.step == 
            0) {
            text <- paste("%", nchar.stepmax, "s", sep = "")
            cat(sprintf(eval(expression(text)), step), "\tmin thresh: ", 
                min.reached.threshold, "\n", rep(" ", lifesign), 
                sep = "")
            flush.console()
        }
        if (algorithm == "rprop+") 
            result <- plus(gradients, gradients.old, weights, 
                nrow.weights, ncol.weights, learningrate, learningrate.factor, 
                learningrate.limit, exclude)
        else if (algorithm == "backprop") 
            result <- backprop(gradients, weights, length.weights, 
                nrow.weights, ncol.weights, learningrate.bp, 
                exclude)
        else result <- minus(gradients, gradients.old, weights, 
            length.weights, nrow.weights, ncol.weights, learningrate, 
            learningrate.factor, learningrate.limit, algorithm, 
            exclude)
        gradients.old <- result$gradients.old
        weights <- result$weights
        learningrate <- result$learningrate
        result <- compute.net(weights, length.weights, covariate = covariate, 
            act.fct = act.fct, act.deriv.fct = act.deriv.fct, 
            output.act.fct = output.act.fct, output.act.deriv.fct = output.act.deriv.fct, 
            special)
        err.deriv <- err.deriv.fct(result$net.result, response)
        gradients <- calculate.gradients(weights = weights, length.weights = length.weights, 
            neurons = result$neurons, neuron.deriv = result$neuron.deriv, 
            err.deriv = err.deriv, exclude = exclude, linear.output = linear.output)
        reached.threshold <- max(abs(gradients))
        if (reached.threshold < min.reached.threshold) {
            min.reached.threshold <- reached.threshold
        }
        step <- step + 1
    }
    if (lifesign != "none" && reached.threshold > threshold) {
        cat("stepmax\tmin thresh: ", min.reached.threshold, "\n", 
            sep = "")
    }
    return(list(weights = weights, step = as.integer(step), reached.threshold = reached.threshold, 
        net.result = result$net.result, neuron.deriv = result$neuron.deriv))
}
compute.net <-
function (weights, length.weights, covariate, act.fct, act.deriv.fct, 
    output.act.fct, output.act.deriv.fct, special) 
{
    neuron.deriv <- NULL
    neurons <- list(covariate)
    if (length.weights > 1) 
        for (i in 1:(length.weights - 1)) {
            temp <- neurons[[i]] %*% weights[[i]]
            act.temp <- act.fct(temp)
            if (special) 
                neuron.deriv[[i]] <- act.deriv.fct(act.temp)
            else neuron.deriv[[i]] <- act.deriv.fct(temp)
            neurons[[i + 1]] <- cbind(1, act.temp)
        }
    if (!is.list(neuron.deriv)) 
        neuron.deriv <- list(neuron.deriv)
    temp <- neurons[[length.weights]] %*% weights[[length.weights]]
    net.result <- output.act.fct(temp)
    if (special) 
        neuron.deriv[[length.weights]] <- output.act.deriv.fct(net.result)
    else neuron.deriv[[length.weights]] <- output.act.deriv.fct(temp)
    if (any(is.na(neuron.deriv))) 
        stop("neuron derivatives contain a NA; varify that the derivative function does not divide by 0", 
            call. = FALSE)
    list(neurons = neurons, neuron.deriv = neuron.deriv, net.result = net.result)
}
calculate.gradients <-
function (weights, length.weights, neurons, neuron.deriv, err.deriv, 
    exclude, linear.output) 
{
    if (any(is.na(err.deriv))) 
        stop("the error derivative contains a NA; varify that the derivative function does not divide by 0 (e.g. cross entropy)", 
            call. = FALSE)
    if (!linear.output) 
        delta <- neuron.deriv[[length.weights]] * err.deriv
    else delta <- err.deriv
    gradients <- crossprod(neurons[[length.weights]], delta)
    if (length.weights > 1) 
        for (w in (length.weights - 1):1) {
            delta <- neuron.deriv[[w]] * tcrossprod(delta, remove.intercept(weights[[w + 
                1]]))
            gradients <- c(crossprod(neurons[[w]], delta), gradients)
        }
    gradients[-exclude]
}
plus <-
function (gradients, gradients.old, weights, nrow.weights, ncol.weights, 
    learningrate, learningrate.factor, learningrate.limit, exclude) 
{
    weights <- unlist(weights)
    sign.gradient <- sign(gradients)
    temp <- gradients.old * sign.gradient
    positive <- temp > 0
    negative <- temp < 0
    not.negative <- !negative
    if (any(positive)) {
        learningrate[positive] <- pmin.int(learningrate[positive] * 
            learningrate.factor$plus, learningrate.limit$max)
    }
    if (any(negative)) {
        weights[-exclude][negative] <- weights[-exclude][negative] + 
            gradients.old[negative] * learningrate[negative]
        learningrate[negative] <- pmax.int(learningrate[negative] * 
            learningrate.factor$minus, learningrate.limit$min)
        gradients.old[negative] <- 0
        if (any(not.negative)) {
            weights[-exclude][not.negative] <- weights[-exclude][not.negative] - 
                sign.gradient[not.negative] * learningrate[not.negative]
            gradients.old[not.negative] <- sign.gradient[not.negative]
        }
    }
    else {
        weights[-exclude] <- weights[-exclude] - sign.gradient * 
            learningrate
        gradients.old <- sign.gradient
    }
    list(gradients.old = gradients.old, weights = relist(weights, 
        nrow.weights, ncol.weights), learningrate = learningrate)
}
backprop <-
function (gradients, weights, length.weights, nrow.weights, ncol.weights, 
    learningrate.bp, exclude) 
{
    weights <- unlist(weights)
    if (!is.null(exclude)) 
        weights[-exclude] <- weights[-exclude] - gradients * 
            learningrate.bp
    else weights <- weights - gradients * learningrate.bp
    list(gradients.old = gradients, weights = relist(weights, 
        nrow.weights, ncol.weights), learningrate = learningrate.bp)
}
minus <-
function (gradients, gradients.old, weights, length.weights, 
    nrow.weights, ncol.weights, learningrate, learningrate.factor, 
    learningrate.limit, algorithm, exclude) 
{
    weights <- unlist(weights)
    temp <- gradients.old * gradients
    positive <- temp > 0
    negative <- temp < 0
    if (any(positive)) 
        learningrate[positive] <- pmin.int(learningrate[positive] * 
            learningrate.factor$plus, learningrate.limit$max)
    if (any(negative)) 
        learningrate[negative] <- pmax.int(learningrate[negative] * 
            learningrate.factor$minus, learningrate.limit$min)
    if (algorithm != "rprop-") {
        delta <- 10^-6
        notzero <- gradients != 0
        gradients.notzero <- gradients[notzero]
        if (algorithm == "slr") {
            min <- which.min(learningrate[notzero])
        }
        else if (algorithm == "sag") {
            min <- which.min(abs(gradients.notzero))
        }
        if (length(min) != 0) {
            temp <- learningrate[notzero] * gradients.notzero
            sum <- sum(temp[-min]) + delta
            learningrate[notzero][min] <- min(max(-sum/gradients.notzero[min], 
                learningrate.limit$min), learningrate.limit$max)
        }
    }
    weights[-exclude] <- weights[-exclude] - sign(gradients) * 
        learningrate
    list(gradients.old = gradients, weights = relist(weights, 
        nrow.weights, ncol.weights), learningrate = learningrate)
}
calculate.generalized.weights <-
function (weights, neuron.deriv, net.result) 
{
    for (w in 1:length(weights)) {
        weights[[w]] <- remove.intercept(weights[[w]])
    }
    generalized.weights <- NULL
    for (k in 1:ncol(net.result)) {
        for (w in length(weights):1) {
            if (w == length(weights)) {
                temp <- neuron.deriv[[length(weights)]][, k] * 
                  1/(net.result[, k] * (1 - (net.result[, k])))
                delta <- tcrossprod(temp, weights[[w]][, k])
            }
            else {
                delta <- tcrossprod(delta * neuron.deriv[[w]], 
                  weights[[w]])
            }
        }
        generalized.weights <- cbind(generalized.weights, delta)
    }
    return(generalized.weights)
}
generate.output <-
function (covariate, call, rep, threshold, matrix, startweights, 
    model.list, response, err.fct, act.fct, data, list.result, 
    linear.output, exclude) 
{
    covariate <- t(remove.intercept(t(covariate)))
    nn <- list(call = call)
    class(nn) <- c("nn")
    nn$response <- response
    nn$covariate <- covariate
    nn$model.list <- model.list
    nn$err.fct <- err.fct
    nn$act.fct <- act.fct
    nn$linear.output <- linear.output
    nn$data <- data
    nn$exclude <- exclude
    if (!is.null(matrix)) {
        nn$net.result <- NULL
        nn$weights <- NULL
        nn$generalized.weights <- NULL
        nn$startweights <- NULL
        for (i in 1:length(list.result)) {
            nn$net.result <- c(nn$net.result, list(list.result[[i]]$net.result))
            nn$weights <- c(nn$weights, list(list.result[[i]]$weights))
            nn$startweights <- c(nn$startweights, list(list.result[[i]]$startweights))
            nn$generalized.weights <- c(nn$generalized.weights, 
                list(list.result[[i]]$generalized.weights))
        }
        nn$result.matrix <- generate.rownames(matrix, nn$weights[[1]], 
            model.list)
    }
    return(nn)
}
generate.rownames <-
function (matrix, weights, model.list) 
{
    rownames <- rownames(matrix)[rownames(matrix) != ""]
    for (w in 1:length(weights)) {
        for (j in 1:ncol(weights[[w]])) {
            for (i in 1:nrow(weights[[w]])) {
                if (i == 1) {
                  if (w == length(weights)) {
                    rownames <- c(rownames, paste("Intercept.to.", 
                      model.list$response[j], sep = ""))
                  }
                  else {
                    rownames <- c(rownames, paste("Intercept.to.", 
                      w, "layhid", j, sep = ""))
                  }
                }
                else {
                  if (w == 1) {
                    if (w == length(weights)) {
                      rownames <- c(rownames, paste(model.list$variables[i - 
                        1], ".to.", model.list$response[j], sep = ""))
                    }
                    else {
                      rownames <- c(rownames, paste(model.list$variables[i - 
                        1], ".to.1layhid", j, sep = ""))
                    }
                  }
                  else {
                    if (w == length(weights)) {
                      rownames <- c(rownames, paste(w - 1, "layhid.", 
                        i - 1, ".to.", model.list$response[j], 
                        sep = ""))
                    }
                    else {
                      rownames <- c(rownames, paste(w - 1, "layhid.", 
                        i - 1, ".to.", w, "layhid", j, sep = ""))
                    }
                  }
                }
            }
        }
    }
    rownames(matrix) <- rownames
    colnames(matrix) <- 1:(ncol(matrix))
    return(matrix)
}
relist <-
function (x, nrow, ncol) 
{
    list.x <- NULL
    for (w in 1:length(nrow)) {
        length <- nrow[w] * ncol[w]
        list.x[[w]] <- matrix(x[1:length], nrow = nrow[w], ncol = ncol[w])
        x <- x[-(1:length)]
    }
    list.x
}
remove.intercept <-
function (matrix) 
{
    matrix(matrix[-1, ], ncol = ncol(matrix))
}
type <-
function (fct) 
{
    attr(fct, "type")
}
print.nn <-
function (x, ...) 
{
    matrix <- x$result.matrix
    cat("Call: ", deparse(x$call), "\n\n", sep = "")
    if (!is.null(matrix)) {
        if (ncol(matrix) > 1) {
            cat(ncol(matrix), " repetitions were calculated.\n\n", 
                sep = "")
            sorted.matrix <- matrix[, order(matrix["error", ])]
            if (any(rownames(sorted.matrix) == "aic")) {
                print(t(rbind(Error = sorted.matrix["error", 
                  ], AIC = sorted.matrix["aic", ], BIC = sorted.matrix["bic", 
                  ], `Reached Threshold` = sorted.matrix["reached.threshold", 
                  ], Steps = sorted.matrix["steps", ])))
            }
            else {
                print(t(rbind(Error = sorted.matrix["error", 
                  ], `Reached Threshold` = sorted.matrix["reached.threshold", 
                  ], Steps = sorted.matrix["steps", ])))
            }
        }
        else {
            cat(ncol(matrix), " repetition was calculated.\n\n", 
                sep = "")
            if (any(rownames(matrix) == "aic")) {
                print(t(matrix(c(matrix["error", ], matrix["aic", 
                  ], matrix["bic", ], matrix["reached.threshold", 
                  ], matrix["steps", ]), dimnames = list(c("Error", 
                  "AIC", "BIC", "Reached Threshold", "Steps"), 
                  c(1)))))
            }
            else {
                print(t(matrix(c(matrix["error", ], matrix["reached.threshold", 
                  ], matrix["steps", ]), dimnames = list(c("Error", 
                  "Reached Threshold", "Steps"), c(1)))))
            }
        }
    }
    cat("\n")
}
