


dens_estim <- function(x, point, bw, eachSide = TRUE) {
    
    N <- length(x)
    
    if (missing(bw)) 
        bw <- 1.84 * sd(x) * N^(-1/5)
    
    if (eachSide) {
        isIn_bw_left <- x >= (point - bw) & x < point
        isIn_bw_right <- x >= point & x <= (point + bw)
        
        NisIn_bw_left <- sum(isIn_bw_left, na.rm = TRUE)
        NisIn_bw_right <- sum(isIn_bw_right, na.rm = TRUE)
        
        res <- (NisIn_bw_left + NisIn_bw_right)/(2 * N * bw)
    } else {
        isIn_bw_both <- x >= (point - bw) & x <= (point + bw)
        NisIn_bw_both <- sum(isIn_bw_both, na.rm = TRUE)
        res <- NisIn_bw_both/(2 * N * bw)
    }
    res
}

dens_estim2 <- function(x, point, bw, kernel = "gaussian", ...) {
    
    
    if (missing(bw)) 
        bw <- "SJ"
    
    d <- density(x, from = point, to = point, n = 1, na.rm = TRUE, kernel = kernel, bw = bw, ...)
    d$y
}


var_estim <- function(x, y, point, bw, eachSide = TRUE) {
    
    
    N <- length(x)
    if (missing(bw)) 
        bw <- 1.84 * sd(x) * N^(-1/5)
    
    if (eachSide) {
        isIn_bw_left <- x >= (point - bw) & x < point
        isIn_bw_right <- x >= point & x <= (point + bw)
        var_inh_left <- var(y[isIn_bw_left], na.rm = TRUE)
        var_inh_right <- var(y[isIn_bw_right], na.rm = TRUE)
        res <- c(var_inh_left, var_inh_right)
    } else {
        isIn_bw <- x >= (point - bw) & x <= point + bw
        var_inh <- var(y[isIn_bw], na.rm = TRUE)
        res <- var_inh
    }
    res
}


#' @importFrom locpol locpol 
#' @importFrom locpol gaussK

### Add locpol kernel for uniform:
uniK <- function(x) ifelse(abs(x) <= 1, 1/2, 0)
attr(uniK, "RK") <- 1/2  ## Rk: kernel(u)^2
attr(uniK, "mu0K") <- 1
attr(uniK, "mu2K") <- 1/3  ## second orde rmoment of K
attr(uniK, "K4") <- NA  ## see with author!
attr(uniK, "RdK") <- NA  ## see with author!
attr(uniK, "dom") <- c(-1, 1)  ##

var_estim2 <- function(x, y, point, bw, estim = c("var", "NW", "NW_loc", "LL_kern", "LL_loc", "var_loc"), sides = c("both", "uni"), 
    kernel = c("Normal", "Uniform"), dfadj = TRUE) {
    
    sides <- match.arg(sides)
    estim <- match.arg(estim)
    kernel <- match.arg(kernel)
    N <- length(x)
    if (missing(bw)) 
        bw <- 1.84 * sd(x) * N^(-1/5)
    
    if (sides == "uni") {
        isIn_bw_left <- x >= (point - bw) & x < point
        isIn_bw_right <- x >= point & x <= (point + bw)
        var_inh_left <- var(y[isIn_bw_left], na.rm = TRUE)
        var_inh_right <- var(y[isIn_bw_right], na.rm = TRUE)
        res <- c(var_inh_left, var_inh_right)
    } else {
        if (estim == "NW") {
            ker <- switch(kernel, Uniform = "box", Normal = "normal")
            m <- ksmooth(x = x, y = y, bandwidth = bw * 2, x.points = point, kernel = ker)$y
            s <- ksmooth(x = x, y = y^2, bandwidth = bw * 2, x.points = point, kernel = ker)$y
        } else if (estim == "NW_loc") {
            ker <- switch(kernel, Uniform = uniK, Normal = gaussK)
            df_xy <- data.frame(y = y, x = x, y2 = y^2)
            # a <<- locCteSmootherC(x=x, y=y, xeval=point, bw=bw, kernel=uniK) aa <<- locCteSmootherC(x=x, y=y, xeval=point, bw=bw,
            # kernel=gaussK)
            m <- locpol(y ~ x, data = df_xy, bw = bw, xeval = point, deg = 0, kernel = ker)
            s <- locpol(y2 ~ x, data = df_xy, bw = bw, xeval = point, deg = 0, kernel = ker)
            m <- m$lpFit["y"]
            s <- s$lpFit["y2"]
        } else if (estim == "LL_kern") {
            if (kernel != "Normal") 
                warning("Kernel set to Normal for locpoly")
            m <- locpoly(x = x, y = y, bandwidth = bw, gridsize = 200)
            s <- locpoly(x = x, y = y^2, bandwidth = bw, gridsize = 200)
            m <- m$y[which.min(abs(m$x - point))]
            s <- s$y[which.min(abs(s$x - point))]
        } else if (estim == "LL_loc") {
            ker <- switch(kernel, Uniform = uniK, Normal = gaussK)
            df_xy <- data.frame(y = y, x = x, y2 = y^2)
            m <- locpol(y ~ x, data = df_xy, bw = bw, xeval = point, kernel = ker)
            s <- locpol(y2 ~ x, data = df_xy, bw = bw, xeval = point, kernel = ker)
            m <- m$lpFit["y"]
            s <- s$lpFit["y2"]
        } else {
            s <- m <- 1
        }
        sh <- s - m^2
        res <- sh
        if (estim == "var_loc") {
            ker <- switch(kernel, Uniform = uniK, Normal = gaussK)
            df_xy <- data.frame(y = y, x = x, y2 = y^2)
            m <- locpol(y ~ x, data = df_xy, bw = bw, xeval = point, kernel = ker)
            res <- m$lpFit$var
        } else if (estim == "var") {
            isIn_bw <- x >= (point - bw) & x <= (point + bw)
            var <- var(y[isIn_bw], na.rm = TRUE)
            res <- if (dfadj) 
                var * (sum(isIn_bw) - 1)/sum(isIn_bw) else var
        }
        
    }
    names(res) <- NULL
    as.numeric(res)
}


## Formula: \sqrt[ (C_2 * \sigma(x)^2 / f(x)) / ( n * h) ] Imbens & Kalyan: C_2/N*h (sigma_l^2 + \sigma_r^2)/f(x) value of
## constant: 4.8 (using boundary kernel: Triangular (value of constant: 33.6 (using boundary kernel: Triangular
## library(locpol) computeRK(equivKernel(TrianK, nu=0, deg=1, lower=0, upper=1), lower=0, upper=Inf) or:
## computeRK(equivKernel(TrianK, nu=0, deg=1, lower=-1, upper=1), lower=-Inf, upper=Inf)

all_var_low <- function(x, y, point, bw, eachSide = TRUE, return = c("se", "all")) {
    
    return <- match.arg(return)
    
    N <- length(x)
    if (missing(bw)) 
        bw <- 1.84 * sd(x) * N^(-1/5)
    
    var <- var_estim(x = x, y = y, point = point, bw = bw, eachSide = eachSide)
    dens <- dens_estim(x = x, point = point, bw = bw, eachSide = eachSide)
    
    C2 <- if (eachSide) 
        4.8 else 2/3
    const <- C2/(N * bw)
    all <- const * sum(var)/dens
    res <- sqrt(all)
    names(res) <- "se"
    if (return == "all") 
        res <- c(res, cons = const, dens = dens, var = sum(var))
    res
    
}


all_var <- function(...) all_var_low(...)

all_var.rdd_reg.np <- function(x) {
    
    bw <- getBW(x)
    dat <- getOriginalData(x)
    cutpoint <- getCutpoint(x)
    res <- all_var_low(dat$x, dat$y, point = cutpoint, bw = bw, eachSide = TRUE, return = "se")
    res
} 
