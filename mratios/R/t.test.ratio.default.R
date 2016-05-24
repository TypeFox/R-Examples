
"CIratioiter" <- 
function (nx, ny, mx, my, vx, vy, alternative = "two.sided", conf.level = 0.95, ul=1e10, ll=-1e10) 
{

    est <- mx/my

tratio <- function(rho, nx, ny, mx, my, vx, vy, alpha) {
degf <- ((vx/nx + (rho^2) * vy/ny)^2)/((vx^2)/((nx^2) * (nx - 1)) + (rho^4) * (vy^2)/((ny^2) * (ny - 1)))
crit <- qt(p = alpha, df = degf, lower.tail = FALSE)
stat <- (mx - my * rho)/sqrt(vx/nx + (rho^2) * vy/ny)
return(abs(stat)-crit)}

switch(alternative,
"two.sided"={
alpha2 <- (1-conf.level)/2

UPPER <- try(uniroot(f=function(x){tratio(rho=x, nx=nx, ny=ny, mx=mx, my=my, vx=vx, vy=vy, alpha=alpha2)},
interval=c(est,ul)), silent=TRUE)

LOWER <- try(uniroot(f=function(x){tratio(rho=x, nx=nx, ny=ny, mx=mx, my=my, vx=vx, vy=vy, alpha=alpha2)},
interval=c(ll, est)), silent=TRUE)

if(class(LOWER)=="try-error")
{lower<-NA; warning(paste("Lower confidence limit can not be found in [",ll, ",", signif(est,3),"].", LOWER))}
else{lower <- LOWER$root}

if(class(UPPER)=="try-error")
{upper <- NA; warning(paste("Upper confidence limit can not be foundin [", signif(est,4), ",", ul,"].", UPPER))}
else{upper <- UPPER$root}
},
"less"={
alpha <- (1-conf.level)
lower <- (-Inf)

UPPER <- try(uniroot(f=function(x){tratio(rho=x, nx=nx, ny=ny, mx=mx, my=my, vx=vx, vy=vy, alpha=alpha)},
lower=est, upper=ul), silent=TRUE)
if(class(UPPER)=="try-error")
{upper <- NA; warning(paste("Upper confidence limit can not be found in [", signif(est,4), ",", ul,"].", UPPER))}
else{upper <- UPPER$root}
},

"greater"={
alpha <- (1-conf.level)
upper <- Inf

LOWER <- try(uniroot(f=function(x){tratio(rho=x, nx=nx, ny=ny, mx=mx, my=my, vx=vx, vy=vy, alpha=alpha)},
lower=ll, upper=est), silent=TRUE)
if(class(LOWER)=="try-error")
{lower<-NA; warning(paste("Lower confidence limit can not be found in [",ll, ",", signif(est,3),"].", LOWER))}
else{lower <- LOWER$root}

})

return(c(lower=lower, upper=upper))

}


###################################################


t.test.ratio.default <- function (x, y, alternative = "two.sided", rho = 1, var.equal = FALSE, conf.level = 0.95, iterativeCI=FALSE, ul=1e10, ll=-1e10, ...) 
{
    addargs <- list(...)
    alternative <- match.arg(alternative, choices = c("two.sided", 
        "less", "greater"))
    if (!is.numeric(rho) | length(rho) != 1) {
        stop("Argument 'rho' must be a single numeric value")
    }
    if (!is.logical(var.equal) | length(var.equal) != 1) {
        stop("Argument'var.equal' must be either TRUE or FALSE")
    }
    if (!is.numeric(conf.level) | length(conf.level) != 1 | conf.level <= 
        0.5 | conf.level >= 1) {
        stop("Argument 'conf.level' must be a single numeric value between 0.5 and 1")
    }
    if (!is.numeric(c(x, y))) {
        stop("x, y, must be numeric vectors")
    }
    if (length(x) < 2 | length(y) < 2) {
        stop("x and y must contain at least two observations each")
    }
    mx <- mean(x)
    my <- mean(y)
    nx <- length(x)
    ny <- length(y)
    vx <- var(x)
    vy <- var(y)
    est <- mx/my

    if (sqrt(vx) < 10 * .Machine$double.eps * abs(mx)) {
        stop("data in x are essentially constant")
    }
    if (sqrt(vy) < 10 * .Machine$double.eps * abs(my)) {
        stop("data in y are essentially constant")
    }
    if (is.null(addargs$namex) || is.null(addargs$namey)) {
        namex = "x"
        namey = "y"
    }
    else {
        namex = addargs$namex
        namey = addargs$namey
    }
    if(any(c(my,mx)<0)){warning("Sample means are smaller than 0! References for this test do not consider this case explicitly!")}
    if(my<0){mxI<-(-mx); myI<-(-my)}else{mxI<-mx; myI<-my}

    if (var.equal == TRUE) {
        degf <- nx + ny - 2
        spool <- sqrt((vx * (nx - 1) + vy * (ny - 1))/degf)
        statistic <- (mxI - myI * rho)/(spool * sqrt(1/nx + (rho^2)/ny))
        if (alternative == "less") {
            p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "greater") {
            p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "two.sided") {
            p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
                lower.tail = FALSE))
            alpha <- (1 - conf.level)/2
        }
        method <- "Ratio-t-test for equal variances"
        vpool <- (vx * (nx - 1) + vy * (ny - 1))/degf
        quant <- qt(p = 1 - alpha, df = degf, lower.tail = TRUE)
        tA <- ((vpool * quant^2)/ny) - my^2
        tB <- 2 * mxI * myI
        tC <- ((vpool * quant^2)/nx) - mx^2
        if (tA >= 0) {
            warning("Confidence set unbounded.")
            upper <- NA
            lower <- NA
        }
        else {
            upper <- tB/(-2 * tA) - sqrt(((tB/2)^2) - tA * tC)/tA
            lower <- tB/(-2 * tA) + sqrt(((tB/2)^2) - tA * tC)/tA
        }
    }

    if (var.equal == FALSE & iterativeCI == FALSE) {

        degf <- max(1, ((vx/nx + (rho^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
            (nx - 1)) + (rho^4) * (vy^2)/((ny^2) * (ny - 1))) )

        stderr <- sqrt(vx/nx + (rho^2) * vy/ny)
        statistic <- (mxI - myI * rho)/stderr

        if (alternative == "less") {
            p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "greater") {
            p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "two.sided") {
            p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
                lower.tail = FALSE))
            alpha <- (1 - conf.level)/2
        }
        method <- "Ratio t-test for unequal variances"

        degfest <- max(1, ((vx/nx + (est^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
            (nx - 1)) + (est^4) * (vy^2)/((ny^2) * (ny - 1))) )

        quant <- qt(p = 1 - alpha, df = degfest, lower.tail = TRUE)

        tA <- ((vy * quant^2)/ny) - my^2
        tB <- 2 * mxI * myI
        tC <- ((vx * quant^2)/nx) - mx^2

        if (tA >= 0) {
            warning("Confidence set unbounded.")
            upper <- NA
            lower <- NA
        }
        else {
            upper <- tB/(-2 * tA) - sqrt(((tB/2)^2) - tA * tC)/tA
            lower <- tB/(-2 * tA) + sqrt(((tB/2)^2) - tA * tC)/tA
        }
    }
    
    
    if (var.equal == FALSE & iterativeCI == TRUE) {

        degf <- ((vx/nx + (rho^2) * vy/ny)^2)/((vx^2)/((nx^2) * 
            (nx - 1)) + (rho^4) * (vy^2)/((ny^2) * (ny - 1)))
        stderr <- sqrt(vx/nx + (rho^2) * vy/ny)
        statistic <- (mxI - myI * rho)/stderr
        if (alternative == "less") {
            p.value <- pt(q = statistic, df = degf, lower.tail = TRUE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "greater") {
            p.value <- pt(q = statistic, df = degf, lower.tail = FALSE)
            alpha <- (1 - conf.level)
        }
        if (alternative == "two.sided") {
            p.value <- min(1, 2 * pt(q = abs(statistic), df = degf, 
                lower.tail = FALSE))
            alpha <- (1 - conf.level)/2
        }

        method <- "Ratio t-test for unequal variances"
        
        conf.int <- CIratioiter(nx=nx, ny=ny, mx=mxI, my=myI, vx=vx, vy=vy, alternative = alternative, conf.level = conf.level, ul=ul, ll=ll) 
        lower<-conf.int[1]
        upper<-conf.int[2]
    }
    
    
    if (alternative == "two.sided") {
        conf.int <- c(lower, upper)
    }
    else {
        if (alternative == "less") {
            conf.int <- c(-Inf, upper)
        }
        else {
            if (alternative == "greater") {
                conf.int <- c(lower, Inf)
            }
        }
    }

    
    names(statistic) <- "t"
    estimate <- c(mx, my, est)
    names(estimate) <- c(paste("mean", namex), paste("mean", 
        namey), paste(namex, namey, sep = "/"))
    names(degf) <- "df"
    names(rho) <- "ratio of means"
    data.name <- paste(namex, namey, sep = " and ")
    conf.int<-as.numeric(conf.int)
    attr(conf.int, "conf.level") <- conf.level

    out <- list(statistic = statistic, parameter = degf, p.value = p.value, 
        conf.int = conf.int, estimate = estimate, null.value = rho, 
        alternative = alternative, method = method, data.name = data.name)
    class(out) <- "htest"
    return(out)
}

