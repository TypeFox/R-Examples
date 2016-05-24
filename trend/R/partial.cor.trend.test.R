partial.cor.trend.test <- function(x, z, method = c("pearson", "spearman")){
##    Copyright (C) 2015, 2016  Thorsten Pohlert
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    This function computes the partial correlation trend test.
##
    method <- match.arg(method)
    na.fail(x)
    na.fail(z)
    n1 <- length(x)
    n2 <- length(z)
    if (n1 != n2) {
        stop("Error: time-series x and y must be of same length")
    }
    n <- n1
    dat <- data.frame(t=(1:n),x, z)
    rho <- cor(dat, method=method)
    a <- rho[1,2] - rho[2,3] * rho[1,3]
    b <- sqrt((1 - rho[2,3]^2)) * sqrt((1 - rho[1,3]^2))
    if (a == 0 & b == 0){
        rho.xt.z <- 0
    }
    else {
        rho.xt.z <- a / b
    }
    T <- (sqrt(n -2) * rho.xt.z) / sqrt((1 - rho.xt.z^2))
    pvalue <- 2 * (1 - pt(abs(T), df=(n-2)))


    # Ausgabe fuer Klasse 'htest'

    res <- list(statistic=NULL,
                parameter = NULL,
                estimate = NULL,
                p.value =NULL,
                statistic =NULL,
                alternative = "true correlation is not equal to 0",
                method = NULL,
                data.name = NULL,
                cor)

    DNAMEX <- deparse(substitute(x))
    DNAMEY <- deparse(substitute(z))
    partrval <- paste("r(t",DNAMEX,".", DNAMEY,")", sep="")
    res$data.name <- paste("t AND ",DNAMEX," . ",DNAMEY)
    names(rho.xt.z) <- partrval
    res$estimate <- rho.xt.z
    names(T) <- "t"
    res$statistic <- T
    df <- n -2
    names(df) <- "df"
    res$parameter <- df
    names(pvalue) <- "p-value"
    res$p.value <- pvalue
    res$cor <- rho
    dimnames(res$cor)[[1]] <- c("t", DNAMEX, DNAMEY)
    dimnames(res$cor)[[2]] <- c("t", DNAMEX, DNAMEY)

    if(method == "pearson"){
        res$method <- "Pearson's partial correlation trend test"
    }
    else {
        res$method <- "Spearman's partial correlation trend test"
    }
    class(res) <- "htest"
    return(res)
}
