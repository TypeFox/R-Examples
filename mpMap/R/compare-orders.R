#' Compare potential orders for linkage groups 
#' 
#' Compares potential orderings on the basis of likelihood or number of expectedcrossovers. Based off of the ripple function in R/qtl. 
#' @export
#' @param cross Object of class \code{mpcross} or \code{cross}
#' @param chr Selected chromosomes
#' @param orders Orders to be compared
#' @param window Window size for order comparison
#' @param method Method for comparison, either counting the number of crossovers or likelihood value. See \code{\link[qtl]{ripple}} for further details.
#' @param error.prob See \code{\link[qtl]{ripple}} for details.
#' @param map.function See \code{\link[qtl]{ripple}} for details.
#' @param maxit See \code{\link[qtl]{ripple}} for details.
#' @param tol See \code{\link[qtl]{ripple}} for details.
#' @param sex.sp See \code{\link[qtl]{ripple}} for details.
#' @param verbose See \code{\link[qtl]{ripple}} for details.
#' @details Uses the core of the ripple function from R/qtl in order to compare orderings output from mporder. Note that \code{method="likelihood"} is substantially slower than \code{method="countxo"} and should not be used for a large number of orderings or a large window size. 
#' @return The matrix of orders with crossover counts or likelihood values appended.  Similar output to \code{ripple}. 
#' @references R/qtl
#' @seealso \code{\link[mpMap]{mporder}}, \code{\link[qtl]{ripple}}
#' @examples
#' map <- sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1)
#' compare_orders(sim.dat, chr=1, orders=rbind(1:11, c(1:3, 6:4, 7:11)))

compare_orders <- function(cross, chr, orders, window=2, method = c("countxo", "likelihood"), error.prob = 1e-04, map.function = c("haldane", "kosambi", "c-f", "morgan"), maxit = 4000, tol = 1e-06, sex.sp = TRUE, verbose = TRUE) 
{
    if (inherits(cross, "mpcross")) {
	write2cross(cross, "tmp", chr=chr)
	cross <- qtl:::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(cross, "type"))
    }

    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    if (missing(chr)) {
        chr <- names(cross$geno)[1]
        warning("chr argument not provided; assuming you want chr ", 
            chr)
    }
    else {
        if (length(chr) > 1) 
            stop("ripple only works for one chromosome at a time.")
        if (!qtl:::testchr(chr, names(cross$geno))) 
            stop("Chr ", chr, "not found.")
    }
    cross <- subset(cross, chr = chr)
    chr.name <- names(cross$geno)[1]
    if (nmar(cross)[1] < 3) {
        warning("Less than three markers.")
        return(NULL)
    }
    if (error.prob < 1e-50) 
        error.prob <- 1e-50
    if (error.prob > 1) {
        error.prob <- 1 - 1e-50
        warning("error.prob shouldn't be > 1!")
    }
    if (window < 2) {
        warning("The window argument must be > 1; using window=2.")
        window <- 2
    }
    window <- round(window)
    method <- match.arg(method)
    map.function <- match.arg(map.function)
    n.mar <- totmar(cross)

    # BEH
    if (missing(orders)) 
    {
     if (n.mar <= window) 
        orders <- qtl:::ripple.perm2(n.mar)
     else {
        temp <- qtl:::ripple.perm1(window)
        n <- nrow(temp)
        orders <- cbind(temp, matrix(rep((window + 1):n.mar, 
            n), byrow = TRUE, ncol = n.mar - window))
        for (i in 2:(n.mar - window + 1)) {
            left <- matrix(rep(1:(i - 1), n), byrow = TRUE, ncol = i - 
                1)
            if (i < n.mar - window + 1) 
                right <- matrix(rep((i + window):n.mar, n), byrow = TRUE, 
                  ncol = n.mar - window - i + 1)
            else right <- NULL
            orders <- rbind(orders, cbind(left, temp + i - 1, 
                right))
        }
        orders <- as.numeric(unlist(strsplit(unique(apply(orders, 
            1, paste, collapse = ":")), ":")))
        orders <- matrix(orders, ncol = n.mar, byrow = TRUE)
     }
    } # end of if orders missing # BEH

    n.orders <- nrow(orders)
    if (n.orders > 49) 
        print.by <- 10
    else if (n.orders > 14) 
        print.by <- 5
    else print.by <- 2
    if (method == "likelihood") {
        loglik <- 1:n.orders
        chrlen <- 1:n.orders
        m <- seq(0, by = 5, length = n.mar)
        temcross <- cross
        if (is.matrix(cross$geno[[1]]$map)) 
            temcross$geno[[1]]$map <- rbind(m, m)
        else temcross$geno[[1]]$map <- m
        for (i in 1:n.orders) {
            if (verbose && i == 1) 
                cat("  ", n.orders, "total orders\n")
            if (verbose && (i%/%print.by) * print.by == i) 
                cat("    --Order", i, "\n")
            temcross$geno[[1]]$data <- cross$geno[[1]]$data[, 
                orders[i, ]]
            newmap <- est.map(temcross, chr, error.prob, map.function, 
                m = 0, p = 0, maxit, tol, sex.sp, FALSE)
            loglik[i] <- attr(newmap[[1]], "loglik")
            chrlen[i] <- diff(range(newmap[[1]]))
        }
        loglik <- (loglik - loglik[1])/log(10)
        o <- order(loglik[-1], decreasing = TRUE) + 1
        orders <- cbind(orders, LOD = loglik, chrlen)[c(1, o), 
            ]
    }
    else {
        type <- class(cross)[1]
        if (type == "f2") {
            if (class(cross$geno[[1]]) == "A") 
                func <- "R_ripple_f2"
            else func <- "R_ripple_bc"
        }
        else if (type == "bc" || type == "riself" || type == 
            "risib" || type == "dh") 
            func <- "R_ripple_bc"
        else if (type == "4way") 
            func <- "R_ripple_4way"
        else if (type == "ri4self" || type == "ri8self" || type == 
            "ri4sib" || type == "ri8sib") 
            func <- "R_ripple_ril48"
        else stop("ripple not available for cross ", type)
        genodat <- cross$geno[[1]]$data
        genodat[is.na(genodat)] <- 0
        n.ind <- nind(cross)
        if (verbose) 
            cat("  ", n.orders, "total orders\n")
        z <- .C(func, as.integer(n.ind), as.integer(n.mar), as.integer(genodat), 
            as.integer(n.orders), as.integer(orders - 1), oblxo = as.integer(rep(0, 
                n.orders)), as.integer(print.by), PACKAGE = "qtl")
        oblxo <- z$oblxo
        o <- order(oblxo[-1]) + 1
        orders <- cbind(orders, obligXO = oblxo)[c(1, o), ]
    }
    rownames(orders) <- c("Initial", paste(1:(nrow(orders) - 
        1)))
    class(orders) <- c("ripple", "matrix")
    attr(orders, "chr") <- chr.name
    attr(orders, "window") <- window
    attr(orders, "error.prob") <- error.prob
    attr(orders, "method") <- method
    orders[, 1:n.mar] <- t(apply(orders[, 1:n.mar, drop = FALSE], 
        1, function(a) {
            n <- length(a)
            if ((1:n)[a == 1] > (1:n)[a == n]) 
                return(rev(a))
            else return(a)
        }))
    orders
}
