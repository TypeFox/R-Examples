exploreOutcomes <- function(outcomes, probs, name, ...) {
        maxvalue <- if (is.character(name) && name == "poisson") 50 else 1000
        minvalue <- if (is.character(name) == "poisson") -50 else -1000
    
        outcomes[outcomes == Inf] <- maxvalue
        outcomes[outcomes == -Inf] <- minvalue
        
        myprobs <- probs(outcomes[1]:outcomes[2], ...)
        myouts <- (outcomes[1]:outcomes[2])[!is.nan(myprobs)]
        myprobs <- myprobs[!is.nan(myprobs)]
        
        return(myouts[myprobs > .Machine$double.eps^0.5])
}

#' Make a random variable consisting of possible outcome values and their probabilities or odds
#' 
#' @name RV
#' @docType package
#' @param outcomes Vector of possible outcomes
#' @param probs Vector of probabilities or function defining probabilities
#' @param odds Vector of odds
#' @param fractions If TRUE, return the probabilities as fractions when printing
#' @param range If TRUE, outcomes specify a range of values in the form c(lower, upper)
#' @param verifyprobs If TRUE, verify that the probs sum to one
#' @param id Set the id of the random variable
#' @param ... Additional parameters passed to the function defining outcome probabilities
#' @return random variable as RV object.
#' @importFrom stats rnorm
#' @importFrom utils type.convert
#' @importFrom plyr ddply
#' @importFrom plyr .
#' @importFrom plyr summarise
#' @export
#' @examples
#' # Make a 50:50 Bernoulli random variable:
#' X.Bern <- RV(c(1,0), c(.5,.5))   
#' 
#' # An equivalent method
#' X.Bern <- RV("bernoulli")
#'   
#' # Make a fair coin flip game with payoffs +$1 and -$1:
#' X.fair.coin <- RV(c(1,-1), c(.5,.5))
#' 
#' # Make a biased coin flip game with odds 1:2 and with fair payoffs +$2 and -$1
#' X.biased.coin <- RV(c(2,-1), odds = c(1,2))
#' 
#' # Make a fair die
#' X.fair.die <- RV(1:6, 1/6)
#' 
#' # Make a loaded die, specifying odds 1:1:1:1:2:4 rather than probabilities:
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' 
#' # Make a Poisson random variable
#' pois.func <- function(x, lambda) { lambda^x * exp(-lambda) / factorial(x) }
#' X.pois <- RV(c(0, Inf), pois.func, lambda = 5)
#' 
#' # An equivalent method
#' X.pois <- RV("poisson")
RV <- function(outcomes, probs = NULL, odds = NULL, fractions = (class(probs) != "function"), range = any(is.infinite(outcomes)), verifyprobs = TRUE, id = rnorm(1), ...) {
    ##
    ## Special dists
    ##
    outtry <- outcomes
    if (length(outcomes) == 1 && outcomes == "bernoulli") {
        probs <- function(x, ...) { p <- ifelse(length(list(...)) == 0, .5, list(...)[[1]]); p^x * (1 - p)^(1 - x) }
        outtry <- 0:1
    } else if (length(outcomes) == 1 && outcomes == "poisson") {
        probs <- function(x, ...) { lambda <- ifelse(length(list(...)) == 0, 5, list(...)[[1]]); lambda^x * exp(-lambda) / factorial(x) }
        range <- TRUE
        fractions <- FALSE
        outtry <- c(0, Inf)
    } else if (length(outcomes) == 1 && outcomes == "geometric") {
        probs <- function(x, ...) { p <- ifelse(length(list(...)) == 0, .5, list(...)[[1]]); (1 - p)^(x - 1) * p }
        range <- TRUE
        fractions <- FALSE
        outtry <- c(1, Inf)
    }
    
    test <- fractions # TODO: Fix
    old.out <- outtry
    
    if (range) outcomes <- suppressWarnings(exploreOutcomes(outtry, probs, outcomes, ...))     else outcomes <- outtry
    
    if (class(probs) == "function") probs <- suppressWarnings(probs(outcomes, ...))
    
    pr <- probs
    if (is.null(pr)) pr <- odds
    
    probsSum <- sum(pr)
    
    if ((probsSum > (1 + .Machine$double.eps^0.5)) && is.null(odds) && verifyprobs) stop("Probabilities sum to over 1")
    if (any(pr < 0)) stop("Probabilities cannot be negative")
    
    isOdds <- !is.null(odds)
    
    if (length(outcomes) < length(pr)) {
        stop("More probabilities/odds than outcomes provided")
    } else if (length(outcomes) > length(pr)) {
        pr <- c(pr, rep(ifelse(isOdds, 1, (1 - probsSum) / (length(outcomes) - length(pr))), length(outcomes) - length(pr)))
    }
    
    ## Convert to probs
    if (verifyprobs) probs <- pr / sum(pr)
    names(probs) <- outcomes
    
    ## Remove zero prob events
    ind <- (probs > .Machine$double.eps^0.5)
    outcomes <- outcomes[ind]
    probs <- probs[ind]
    
    if (any(duplicated(outcomes))) {
        out <- NULL
        my.df <- data.frame(out = as.vector(outcomes), pr = as.numeric(probs))
        my.sum <- ddply(my.df, .(out), summarise, pr = sum(pr))
        outcomes <- type.convert(as.character(my.sum$out))
        probs <- as.numeric(my.sum$pr)
        names(probs) <- outcomes
    }
    
    if (length(grep(",", outcomes)) > 0) {
        ord1 <- as.numeric(unlist(lapply(strsplit(outcomes, ","), '[[', 1)))
        ord2 <- as.numeric(unlist(lapply(strsplit(outcomes, ","), '[[', 2)))
        ord3 <- ord1 * length(unique(ord2)) + ord2
        
        outcomes <- outcomes[order(ord3)]
        probs <- probs[order(ord3)]
        names(probs) <- outcomes
    }
    
    class(outcomes) <- "RV"
    
    attr(outcomes, "probs") <- probs
    attr(outcomes, "odds") <- isOdds
    attr(outcomes, "fractions") <- fractions
    attr(outcomes, "range") <- range
    attr(outcomes, "outcomes") <- old.out
    attr(outcomes, "id") <- id

    return(outcomes)
}

#' Make a joint random variable consisting 
#' 
#' @name jointRV
#' @param outcomes The possible outcomes of the joint random variable, as a list
#' @param probs The probabilities of each event, in the order (x1, y1, x1, y2, ..., x2, y1, x2, y2, ..., xn, yn)
#' @param ... Further arguments to be passed to the RV function
#' @return An RV object
#' @export
jointRV <- function(outcomes, probs = NULL, ...) {
    if (!is.list(outcomes)) stop("outcomes must be presented as a list of outcomes for each variable")
    
    myouts <- apply(expand.grid(outcomes)[with(expand.grid(outcomes), order(Var1)),], 1, paste, collapse = ",")
    
    RV(myouts, probs, ...)
}

ugh <- function(X, Y, op) {
    op <- get(op)
    jointdist <- attr(X, "joint")
    if (is.null(jointdist)) jointdist <- X * Y
    
    if (class(Y) == "RV") {
        result <- unlist(lapply(strsplit(jointdist, ","), function(hi){op(as.numeric(hi[attr(X, "num")]), as.numeric(hi[attr(Y, "num")]))}))
    } else {
        result <- unlist(lapply(strsplit(jointdist, ","), function(hi){op(hi[attr(X, "num")], Y)}))
    }
    
    class(result) <- "RVresult"
    
    attr(result, "outcomes") <- outcomes(jointdist)
    attr(result, "rv") <- "joint"
    attr(result, "probs") <- probs(jointdist)
    
    return(result)
}

unopset <- function(X, Xchar, cond, x) {
    if (is.character(x)) x <- paste("\"", x, "\"", sep = "")
    
    X.notrv <- X
    class(X.notrv) <- NULL
    
    result <- eval(parse(text = paste("X.notrv", cond, "c(", paste(x, collapse = ","), ")")))
    class(result) <- "RVresult"
    
    attr(result, "outcomes") <- as.vector(X)
    attr(result, "rv") <- Xchar
    attr(result, "probs") <- probs(X)
    
    return(result)
}

binopset <- function(X, Xchar, cond, Y) {    
    result <- eval(parse(text = paste("as.logical(X)", cond, "as.logical(Y)")))
    class(result) <- "RVresult"
    
    attr(result, "outcomes") <- attr(X, "outcomes")
    attr(result, "rv") <- Xchar
    attr(result, "probs") <- probs(X)
        
    return(result)
}

#' @export
"<.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, "<"))
    else return(unopset(X, deparse(substitute(X)), "<", x)) 
}
#' @export
"<=.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, "<="))
    else return(unopset(X, deparse(substitute(X)), "<=", x)) 
}
#' @export
"==.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, "=="))
    else return(unopset(X, deparse(substitute(X)), "==", x)) 
}
#' @export
"!=.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, "!="))
    else return(unopset(X, deparse(substitute(X)), "!=", x)) 
}
#' @export
">=.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, ">="))
    else return(unopset(X, deparse(substitute(X)), ">=", x)) 
}
#' @export
">.RV" <- function(X, x) { 
    if (class(x) == "RV" || !is.null(attr(X, "joint"))) return(ugh(X, x, ">"))
    else return(unopset(X, deparse(substitute(X)), ">", x)) 
}

#' @export
"+.RV" <- function(X, Y) {
    if (class(X) == "RV" && class(Y) == "RV" && attr(X, "id") != attr(Y, "id")) {
        return(SofI(X, Y, fractions = (attr(X, "fractions") && attr(Y, "fractions"))))
    } else if (class(X) == "RV" && class(Y) != "RV") {
        return(RV(as.numeric(outcomes(X)) + Y, probs(X), fractions = attr(X, "fractions"), id = attr(X, "id")))
    } else if (class(X) != "RV" && class(Y) == "RV") {
        return(RV(as.numeric(outcomes(Y)) + X, probs(Y), fractions = attr(Y, "fractions"), id = attr(Y, "id")))
    } else {
        return(RV(as.numeric(outcomes(X)) + as.numeric(outcomes(Y)), probs(X), fractions = attr(X, "fractions"), id = attr(X, "id")))
    }
}
#' @export
"-.RV" <- function(X, Y) {
    X + (-1 * Y)
}
#' @export
"*.RV" <- function(X, Y) { 
    if (class(X) == "RV" && class(Y) == "RV" && attr(X, "id") != attr(Y, "id")) {
        return(joint(X, Y))
    } else if (class(X) == "RV" && class(Y) != "RV") {
        return(RV(as.numeric(outcomes(X)) * Y, probs(X), fractions = attr(X, "fractions"), id = attr(X, "id")))
    } else if (class(X) != "RV" && class(Y) == "RV") {
        return(RV(as.numeric(outcomes(Y)) * X, probs(Y), fractions = attr(Y, "fractions"), id = attr(Y, "id")))
    } else {
        return(RV(as.numeric(outcomes(X)) * as.numeric(outcomes(Y)), probs(X), fractions = attr(X, "fractions"), id = attr(X, "id")))
    }
}
#' @export
"^.RV" <- function(X, Y) { return(RV(as.numeric(outcomes(X))^Y, probs(X), fractions = attr(X, "fractions"), id = attr(X, "id"))) }

#' Generic method for in operator function
#' 
#' @name %in%
#' @param e1 First vector
#' @param e2 Second vector
#' @return A logical vector indicating which elements of e1 are in e2
#' @export
"%in%" <- function(e1, e2) { UseMethod("%in%") } 

#' @export
"%in%.default" <- function(e1, e2) { base::`%in%`(e1, e2) }
#' @export
"%in%.RV" <- function(e1, e2) { 
    if (!is.null(attr(e1, "joint"))) return(ugh(e1, e2, "%in%"))
    else return(unopset(e1, deparse(substitute(e1)), "%in%", e2))
}

#' Compute the logical OR of two events
#' 
#' @name %OR%
#' @param X RVcond object
#' @param Y RVcond object
#' @return An RVresult object which is two events ORed together
#' @export
#' @examples
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' P((X.fair.die == 4) %OR% (X.fair.die == 3))
"%OR%" <- function(X, Y) { return(binopset(X, deparse(substitute(X)), "|", Y)) }

#' Compute the logical AND of two events
#' 
#' @name %AND%
#' @param X RVcond object
#' @param Y RVcond object
#' @return An RVresult object which is two events ANDed together
#' @export
#' @examples
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' P((X.fair.die == 4) %AND% (X.fair.die == 3))
"%AND%" <- function(X, Y) { return(binopset(X, deparse(substitute(X)), "&", Y)) }

#' @export
"|.RVresult" <- function(vec1, vec2) {
    if (class(vec1) == "RV") {
        jointdist <- attr(vec1, "joint")
        vec <- jointdist[vec2]
        mynum <- attr(vec1, "num")
        
        mysplit <- strsplit(jointdist, ",")
        myout <- as.numeric(unlist(lapply(mysplit, function(test){test[mynum]})))
        
        final.outcomes <- myout[vec2]
        final.probs <- as.numeric(probs(jointdist)[vec2] / sum(probs(jointdist)[vec2]))
        
        return(RV(final.outcomes, final.probs))
    }
    result <- P(vec1 %AND% vec2) / P(vec2)
    class(result) <- "RVcond"
    attr(result, "probs") <- attr(vec1, "probs")
    
    return(result)
}

#' @export
"print.RVcond" <- function(x, ...) {
    return(print.default(as.numeric(x)))
}

#' @export
"print.RVresult" <- function(x, ...) {
    vec <- (as.logical(x))
    names(vec) <- attr(x, "outcomes")
    
    return(print.default(vec, ...))
}

#' Outcomes of random variable X 
#'
#' Obtain the list of outcomes from a random variable
#'
#' @param X random variable
#' @return vector of outcomes of X
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' outcomes(X.Bern)
#' 
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' outcomes(X.fair.die)
#' 
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' outcomes(X.loaded.die)
outcomes <- function(X) {
    return(as.vector(X))
}

#' Probability mass function of random variable X 
#'
#' Obtain the list of probabilities from a random variable: p(x)
#'
#' @param X random variable
#' @return named vector of probablities for each element of the random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' probs(X.Bern)
#' 
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' probs(X.fair.die)
#' 
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' probs(X.loaded.die)
probs <- function(X) { 
    return(attr(X, "probs"))
}

#' Joint probability mass function of random variables X and Y
#'
#' @author Heike Hofmann \email{hofmann@@iastate.edu}
#' @param X random variable
#' @param Y random variable
#' @param sep separator between items from marginal distributions, by default set to ","
#' @param fractions If TRUE, return the probabilities as fractions
#' @export
#' @examples
#' d <- RV(c("A","B","C"), odds = c(3,5,11))
#' d2 <- joint(d,d)
#' probs(d2)
joint <- function(X, Y, sep=",", fractions = (attr(X, "fractions") & attr(Y, "fractions"))) {
    S <- X
    tmp <- tapply(outer(probs(S), probs(Y), FUN="*"),
                  outer(S, Y, FUN="paste", sep=sep), paste, sep=sep)
    S <- names(tmp)

    return(RV(outcomes = as.character(S), probs = as.numeric(tmp), fractions = fractions))
}

#' Probability mass function of  X^n
#'
#' @author Heike Hofmann \email{hofmann@@iastate.edu}
#' @param X random variable
#' @param n power
#' @param sep separator between items from marginal distributions, by default set to ","
#' @param fractions If TRUE, return the probabilities as fractions
#' @export
#' @examples
#' d <- RV(c("A","B","C"), odds = c(3,5,11))
#' d2 <- iid(d)
#' probs(d2)
iid <- function(X, n=2, sep=",", fractions=attr(X, "fractions")) {
    S <- X;  i <- 2
    while(i<=n) {
        tmp <- tapply(outer(probs(S), probs(X), FUN="*"),
                      outer(S, X, FUN="paste", sep=sep), paste, sep=sep)
        S <- names(tmp)
        attr(S, "probs") <- as.numeric(tmp)
        i <- i+1
    }

    return(RV(outcomes = as.character(S), probs = attr(S, "probs"), fractions = fractions))
}

#' Turn a probability vector with possible outcome values in the 'names()' attribute
#' into a random variable:
#'
#' @param px A probability vector with possible outcome values in the 'names()' attribute
#' @param fractions If TRUE, return the probabilities as fractions
#' 
#' @export
as.RV <- function(px, fractions = TRUE) {
    X <- as.numeric(names(px))
    
    class(X) <- "RV"
    
    attr(X, "probs") <- px
    attr(X, "fractions") <- fractions
    attr(X, "odds") <- FALSE
    attr(X, "range") <- FALSE
    
    X
}

#' Calculate probabilities of events
#'
#' @param event A logical vector
#' @export
#' @examples
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' P(X.fair.die>3)
#' 
#' X.loaded.die <- RV(1:6, odds = c(1,1,1,1,2,4))
#' P(X.loaded.die>3)
#' P(X.loaded.die==6)
P <- function(event) { UseMethod("P") } 

#' @export
P.default <- function(event) { sum(attr(event, "probs")[event]) }

#' @export
P.RVcond <- function(event) { return(event) }

#' Expected value of a random variable
#' 
#' @param X random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' E(X.Bern)
#' 
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' E(X.fair.die)
E <- function(X) { 
    isjoint <- length(grep(",", X)) > 0
    if (isjoint) {
        val <- lapply(strsplit(outcomes(X), ","), function(test) {
            prod(as.numeric(test)) * P(X == (paste(test, collapse = ",")))
        })
        return(sum(unlist(as.numeric(val))))
    }
    return(sum(as.numeric(X)*probs(X)))
}

#' Variance of a random variable
#' 
#' @param X random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' E(X.Bern)
V <- function(X) { E((X-E(X))^2) }

#' Standard deviation of a random variable
#' 
#' @param X random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' E(X.Bern)
SD <- function(X) { sqrt(V(X)) }

#' Skewness of a random variable
#' 
#' @param X random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' SKEW(X.Bern)
SKEW <- function(X) { E((X-E(X))^3)/SD(X)^3 }

#' Kurtosis of a random variable
#'
#' @param X random variable
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' KURT(X.Bern)
KURT <- function(X) { E((X-E(X))^4)/V(X)^2 }

#' Sum of independent random variables
#' 
#' @param ... Arbitrary number of random variables
#' @param fractions If TRUE, return the probabilities as fractions
#' @export
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' X.fair.die <- RV(1:6, rep(1/6,6))
#' 
#' S5 <- SofI(X.Bern, X.Bern, X.Bern, X.Bern, X.Bern)  
#' S.mix <- SofI(X.Bern, X.fair.die)  # Independent but not IID
SofI <- function(..., fractions=attr(list(...)[[1]], "fractions")) {
    LIST <- list(...)
    S <- LIST[[1]]
    LIST <- LIST[-1]
    while(length(LIST)>0) {
        X <- LIST[[1]]
        tmp <- tapply(outer(probs(S), probs(X), FUN="*"),
                      outer(as.numeric(outcomes(S)), as.numeric(outcomes(X)), FUN="+"), sum)
        S <- as.numeric(names(tmp))  
        attr(S, "probs") <- tmp
        LIST <- LIST[-1]
    }
    
    return(RV(as.numeric(S), attr(S, "probs"), fractions = fractions))
}

#' Sum of independent identically distributed random variables
#' 
#' @param X A random variable
#' @param n The number of Xs to sum
#' @param fractions If TRUE, return the probabilities as fractions
#' @param progress If TRUE, display a progress bar
#' @export
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @examples
#' X.Bern <- RV(c(1,0), c(.5,.5))
#' 
#' S5 <- SofIID(X.Bern, 5)
#' S128 <- SofIID(X.Bern, 128)
SofIID <- function(X, n=2, progress=TRUE, fractions=attr(X, "fractions")) {
    S <- X;  i <- 2
    pb <- txtProgressBar(min = 1, max = n)
    while(i<=n) {
        tmp <- tapply(outer(probs(S), probs(X), FUN="*"),
                      outer(as.numeric(outcomes(S)), as.numeric(outcomes(X)), FUN="+"), sum)
        
        S <- as.numeric(names(tmp))  
        attr(S, "probs") <- tmp
        
        if(i%%100==0 & progress) setTxtProgressBar(pb, i)
        i <- i+1
    };
    close(pb)
    
    return(RV(as.numeric(S), attr(S, "probs"), id = attr(X, "id"), fractions = fractions))
}

#' Plot a random variable of class "RV"
#' 
#' @method plot RV
#' @param x A random variable
#' @param ... Additional arguments to be passed to the "plot" function
#' @param tol Only display outcomes with probabilities above tol
#' @param pch Either an integer specifying a symbol or a single character to be used as the default in plotting points.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param lwd The line width, a positive number, defaulting to 2.
#' @param col A specification for the default plotting color
#' @param xlab Label for the X axis
#' @param ylab Label for the Y axis
#' @export
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics points
#' @examples
#' fair.die <- RV(1:6, rep(1/6,6))
#' plot(fair.die)
plot.RV <- function(x, ..., tol=1e-10, pch=16, cex=1.2, lwd=2, col="black",
                    xlab="Possible Values",
                    ylab="Probabilities") {
    ## args <- list(...);  print(args)
    if (!is.numeric(outcomes(x))) stop("Plot method only defined with numeric outcomes.")
    
    xx <- x[probs(x) > tol]
    px <- probs(x)[probs(x) > tol]
    plot(xx, px, type="h", lwd=lwd, col=col, xlab=xlab, ylab=ylab, ...)
    abline(h=0, col="gray")
    points(xx, px, pch=pch, cex=cex, col=col)
}

#' Print a random variable of class "RV"
#' 
#' @importFrom MASS fractions
#' 
#' @method print RV
#' @author Eric Hare \email{erichare@@iastate.edu}
#' @param x A random variable
#' @param odds If TRUE, print as odds instead of probs
#' @param fractions If TRUE, print probs as fractions instead of decimals
#' @param all.outcomes If TRUE, print all outcomes rather than the first ten
#' @param digits Number of digits to print for probabilities
#' @param ... Additional arguments to be passed to the "format" function
#' @importFrom utils type.convert
#' @importFrom utils write.table
#' @export
#' @examples
#' fair.die <- RV(1:6, rep(1/6,6))
#' print(fair.die)
print.RV <- function(x, odds = attr(x, "odds"), fractions = attr(x, "fractions"), all.outcomes = FALSE, digits = 3, ...) {
    attributes(x)$class <- NULL
        
    vec <- attr(x, "probs")
    if (!fractions) vec <- round(vec, digits = digits)
    
    if (odds) vec <- vec / min(vec[vec > 0])
    if (fractions) {vec <- fractions(vec)}
    names(vec) <- x
    
    if (odds) type <- "Odds" else type <- "Probs"
    if (fractions & !odds) vec <- as.character(vec)
    if (odds) vec <- paste(round(vec, digits = digits), round(sum(as.numeric(vec)) - as.numeric(vec), digits = digits), sep = ":")
    
    df <- eval(parse(text = paste("data.frame(Outcomes = as.character(x), ", type, " = vec)", sep = "")))
    names(df)[2] <- paste(type, paste(rep(" ", nchar("Outcomes") - nchar(type) - 1), collapse = ""))
    
    df$test <- type.convert(as.character(df$Outcomes), as.is = TRUE)
    if (is.numeric(df$test)) df <- df[with(df, order(test)), ]
    df$test <- NULL
    
    old.df <- df
    
    if (attr(x, "range")) {
        cat(paste("Random variable with outcomes from", attr(x, "outcomes")[1], "to", attr(x, "outcomes")[2], "\n\n"))
    } else {
        cat(paste("Random variable with", length(x), "outcomes\n\n"))
    }
    
    if (nrow(df) > 12 & !all.outcomes) df <- df[1:12,]
    write.table(format(t(df), justify = "right", ...), col.names = FALSE, quote = FALSE)

    if (nrow(old.df) > 12) cat("\nDisplaying first 12 outcomes\n")
}

#' Normal quantile plot for RVs to answer the question how close to normal it is
#'
#' @method qqnorm RV
#' @param y A random variable
#' @param ... Additional arguments to be passed to the "plot" or "points" function
#' @param pch Either an integer specifying a symbol or a single character to be used as the default in plotting points.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.
#' @param add A logical indicating whether to add to an existing plot
#' @param xlab Label for the X axis
#' @param ylab Label for the Y axis
#' @param tol tolerance for the zero probability case
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats qnorm
#' @export
#' @examples
#' fair.die <- RV(1:6, rep(1/6,6))
#' qqnorm(fair.die)
qqnorm.RV <- function(y, ..., pch=16, cex=.5, add=FALSE, xlab="Normal Quantiles", ylab="Random Variable Quantiles", tol = 1e-10) {
    if (!is.numeric(outcomes(y))) stop("qqnorm method only defined with numeric outcomes.")
    
    ind <- which(probs(y) > tol)
    outcomes <- sort(y[ind])
    pc <- cumsum(probs(y))[ind]
    
    if (!add) {
        plot(qnorm(pc), outcomes, pch=pch, cex=cex, xlab=xlab, ylab=ylab, ...)
    } else {
        points(qnorm(pc), outcomes, pch=pch, cex=cex, ...)
    }
}

#' Marginal distributions of a joint random variable
#'
#' Extracts the marginal probability mass functions from a joint distribution.
#' @author Heike Hofmann \email{hofmann@@iastate.edu}
#' @importFrom plyr alply
#' @param X a random variable
#' @param sep parameter specifying the separator between dimensions, defaults to ","
#' @export
#' @importFrom stats xtabs
#' @importFrom utils type.convert
#' @examples
#' X <- RV(1:6, 1/6)
#' X3 <- iid(X, 3)
#' margins(X3)
margins <- function(X, sep=",") {
    dframe <- sapply(strsplit(as.character(X), split=sep, fixed=TRUE), function(x) as.matrix(x))
    
    res <- alply(dframe, .margins=1, function(x) {
        dtab <- xtabs(probs(X)~x)
        my.rv <- RV(type.convert(names(dtab)), as.numeric(dtab))
        attr(my.rv, "joint") <- X
        my.rv
    })   
    
    attributes(res) <- NULL
    for (i in 1:length(res)) attr(res[[i]], "num") <- i
    
    return(res)
}

#' Marginal distribution of a joint random variable
#'
#' Extracts the marginal probability mass functions from a joint distribution.
#' @author Eric Hare \email{erichare@@iastate.edu}
#' @param X A random variable
#' @param num Number indicating which marginal distribution to extract
#' @export
#' @examples
#' AandB <- jointRV(outcomes = list(1:3, 0:2), probs = 1:9 / sum(1:9))
#' marginal(AandB, 1)
#' marginal(AandB, 2)
marginal <- function(X, num) {
    return(margins(X)[[num]])
}

#' Tests whether the random variables X and Y are independent
#' @author Eric Hare \email{erichare@@iastate.edu}
#' @param X A random variable
#' @param Y A random variable
#' @export
#' @examples
#' AandB <- jointRV(outcomes = list(1:3, 0:2), probs = 1:9 / sum(1:9))
#' A <- marginal(AandB, 1)
#' B <- marginal(AandB, 2)
#' independent(A, B) # FALSE
#' CandD <- jointRV(outcomes = list(1:3, 0:2))
#' C <- marginal(CandD, 1)
#' D <- marginal(CandD, 2)
#' independent(C, D) # FALSE
independent <- function(X, Y) {
    firstjoint <- attr(X, "joint")
    secondjoint <- attr(Y, "joint")
    
    if (is.null(firstjoint) | is.null(secondjoint)) return(TRUE)
    
    cond1 <- all(outcomes(firstjoint) == outcomes(secondjoint))
    cond2 <- all(probs(firstjoint) == probs(secondjoint))
    cond3 <- all(as.vector(outer(probs(X), probs(Y))) == as.numeric(probs(firstjoint)))
    
    return(cond1 && cond2 && cond3)
}
