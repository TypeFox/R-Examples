gevcdn.bag <-
function (x, y, iter.max = 1000, iter.step = 10, n.bootstrap = 30,
          n.hidden = 3, Th = gevcdn.logistic, fixed = NULL,
          init.range = c(-0.25, 0.25), scale.min = .Machine$double.eps,
          beta.p = 3.3, beta.q = 2, sd.norm = Inf,
          method = c("BFGS", "Nelder-Mead"),
          max.fails = 100, silent = TRUE, ...)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    method <- match.arg(method)
    if (identical(Th, gevcdn.identity)) n.hidden <- 3
    x.min <- apply(x, 2, min)
    x.max <- apply(x, 2, max)
    x <- sweep(sweep(x, 2, x.min, '-'), 2, x.max - x.min, '/')
    y.min <- min(y)
    y.max <- max(y)
    y <- (y - y.min)/(y.max - y.min)
    w.bootstrap <- list()
    for (i in seq_len(n.bootstrap)){
        cat("*** Bootstrap sample", i, "\n")
        cases.in <- sample(nrow(x), replace=TRUE)
        cases.out <- setdiff(1:nrow(x), cases.in)
        x.train <- x[cases.in,,drop=FALSE]
        x.valid <- x[cases.out,,drop=FALSE]
        y.train <- y[cases.in,,drop=FALSE]
        y.valid <- y[cases.out,,drop=FALSE]
        if (identical(names(init.range), c("W1", "W2"))){
            weights <- unlist(init.range) + gevcdn.initialize(x,
                           n.hidden, c(-0.25, 0.25))
        } else{
            weights <- gevcdn.initialize(x, n.hidden, init.range)
        }
        cost.best <- cost.train <- cost.valid <- Inf
        iter.best <- iter <- 0
        while(iter < iter.max){
            if (!silent) cat("Iter", iter, ";",
                             sprintf("%.6g", cost.train),
                             ";", sprintf("%.6g", cost.valid), '\n')
            exception <- TRUE
            exception.count <- 0
            while (exception){
                fit <- try(suppressWarnings(optim(weights, gevcdn.cost,
                               method = method,
                               control = list(maxit = iter.step, ...),
                               x = x.train, y = y.train,
                               n.hidden = n.hidden, Th = Th,
                               fixed = fixed, scale.min = scale.min,
                               beta.p = beta.p, beta.q = beta.q,
                               sd.norm = sd.norm)),
                       silent = TRUE)
                if (!class(fit) == "try-error"){
                    exception <- fit$value > 1e+308
                }
                if (exception){
                    exception.count <- exception.count + 1
                    weights <- gevcdn.initialize(x, n.hidden,
                                                 init.range)
                    w.best <- gevcdn.reshape(x, weights, n.hidden)
                    cost.best <- cost.train <- cost.valid <- Inf
                    iter.best <- iter <- 0
                }
                if (exception.count == max.fails){
                    stop("exception... check arguments")
                }
            }
            weights <- fit$par
            cost.prev <- cost.train
            cost.train <- gevcdn.cost(weights, x.train, y.train,
                                      n.hidden, Th, fixed, scale.min,
                                      beta.p, beta.q, sd.norm)
            cost.valid <- gevcdn.cost(weights, x.valid, y.valid,
                                      n.hidden, Th, fixed, scale.min,
                                      beta.p, beta.q, sd.norm)
            iter <- iter + iter.step
            if (cost.valid <= cost.best){
                w.best <- gevcdn.reshape(x, weights, n.hidden)
                cost.best <- cost.valid
                iter.best <- iter
            }
            if (abs(cost.train - cost.prev) < .Machine$double.eps){
                cat("local minimum\n")
                break()
            }
        }
        cat("* Best weights at iter", iter.best, ";",
            sprintf("%.6g", cost.best), "\n")
        attr(w.best, "x.min") <- x.min
        attr(w.best, "x.max") <- x.max
        attr(w.best, "y.min") <- y.min
        attr(w.best, "y.max") <- y.max
        attr(w.best, "Th") <- Th
        attr(w.best, "fixed") <- fixed
        attr(w.best, "scale.min") <- scale.min
        attr(w.best, "stopped.training") <- TRUE
        attr(w.best, "cost.valid") <- as.numeric(cost.best)/
                                                 length(cases.out)
        w.bootstrap[[i]] <- w.best
    }
    if (n.bootstrap==1) w.bootstrap <- w.bootstrap[[i]]
    w.bootstrap
}

