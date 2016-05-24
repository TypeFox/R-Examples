#' Genomic control for over-dominant model of inheritance using VIF
#' 
#' This function estimates the corrected statistic using genomic control
#' for the over-dominant model,
#' using VIF. VIF coefficients are estimated
#' by optimizing different error functions: regress,
#' median and ks.test.
#'  
#' @param data Input vector of Chi square statistic
#' @param method Function of error to be optimized. Can be
#' "regress", "median" or "ks.test"
#' @param p Input vector of allele frequencies
#' @param index.filter Indexes for variables that will be use for analysis in data vector
#' @param n size of the sample
#' @param proportion The proportion of lowest P (Chi2) to be used when
#'   estimating the inflation factor Lambda for "regress" method only
#' @param plot If TRUE, plot of lambda will be produced
#' @param lmax The threshold for lambda for plotting (optional)
#' @param color The color of the plot
#' @param clust For developers only
#' @param vart0 For developers only
#' @param tmp For developers only
#' 
#' @return A list with elements
#' \item{Zx}{output vector corrected Chi square statistic}
#' \item{vv}{output vector of VIF}
#' \item{exeps}{output vector of exepsons (NA)}
#' \item{calrate}{output vector of calrate}
#' \item{F}{F}
#' \item{K}{K}
#' 
#' @author Yakov Tsepilov
#' 
#' @examples
#' data(ge03d2)
#' # truncate the data to make the example faster
#' ge03d2 <- ge03d2[seq(from=1,to=nids(ge03d2),by=2),seq(from=1,to=nsnps(ge03d2),by=3)]
#' qts <- mlreg(phdata(ge03d2)$dm2~1,data=ge03d2,gtmode = "overdominant")
#' chi2.1df <- results(qts)$chi2.1df
#' s <- summary(ge03d2)
#' freq <- s$Q.2
#' result <- VIFGC_ovdom(p=freq,method = "median",data=chi2.1df,n=nids(ge03d2))
#' 
#' @keywords htest
#'

VIFGC_ovdom=function (data, p, method = "regress", n, index.filter = NULL, 
    proportion = 1, clust = 0, vart0 = 0, tmp = 0,plot=TRUE,lmax=NULL,color="red") 
{
    if (is.null(index.filter)) {
        ind.function = which(!is.na(data))
    }
    else ind.function = index.filter
    ind.function = ind.function[which(!is.na(data[ind.function]))]
    if (!(method == "regress" | method == "median" | method == 
        "ks.test")) {
        print("Error. I do not know this method")
        break
    }
    Zx = data
    inf = which(!is.na(Zx))
    notinf = which(is.na(Zx))
    Clust <- function(k) {
        data.dist <- as.dist(0.5 - tmp)
        data.mds <- cmdscale(data.dist)
        return(data.mds)
    }
    ClustN <- function(k, data1.mds) {
        km <- kmeans(data1.mds, centers = k, nstart = 1000)
        i = 1
        cl = 0
        while (i <= k) {
            cl[i] <- length(which(km$cluster == i))
            i = i + 1
        }
        return(cl)
    }
    VIF <- function(p, N, F, K) {
        q = 1 - p
        Var = 2 * p * q * x^2 - 2 * F * p * q * x^2 - 4 * p^2 * 
            q^2 * x^2 + 8 * F * p^2 * q^2 * x^2 - 4 * F^2 * p^2 * 
            q^2 * x^2
        S = (1 + F) * (1 + 2 * F)
        Cov = x^2 * ((4 * ((1 - F)^3 * p^2 * q^2 + F * (1 - F) * 
            p * q))/S - (2 * (1 - F) * p * q)^2)
        if (vart0 == 1) {
            VarT0 = N * Var - N * Cov
        }
        if (vart0 == 0) {
            VarT0 = N * Var
        }
        VarT = VarT0 + Cov * K
        VIF = VarT/VarT0
        return(VIF)
    }
    GC_VIF_nlm = function(FK) {
        F = FK[1]
        K = FK[2]
        vector_vif = VIF(p, n, F, K)
        Zxl = Zx/vector_vif
        Zxl = sort(Zxl[ind.function])
        Zxl_r = Zxl[1:ntp]
        if (method == "ks.test") {
            dMedian = -log(ks.test(Zxl, "pchisq", df = 1)$p.value)
        }
        if (method == "median") {
            dMedian = abs(qchisq(0.5, 1) - median(Zxl))
        }
        if (method == "regress") {
            dMedian = sum((Zxl_r - Chi2)^2)
        }
        return(1 * dMedian)
    }
    GC_VIF = function(F, K) {
        vector_vif = VIF(p, n, F, K)
        Zxl = Zx/vector_vif
        Zxl = sort(Zxl[ind.function])
        Zxl_r = Zxl[1:ntp]
        if (method == "ks.test") {
            dMedian = -log(ks.test(Zxl, "pchisq", df = 1)$p.value)
        }
        if (method == "median") {
            dMedian = abs(qchisq(0.5, 1) - median(Zxl))
        }
        if (method == "regress") {
            dMedian = sum((Zxl_r - Chi2)^2)
        }
        return(1 * dMedian)
    }
    F = 0.5
    K = 0.2 * n[1]
    data_p <- data[ind.function]
    if (proportion > 1 || proportion <= 0) 
        stop("proportion argument should be greater then zero and less than or equal to one")
    ntp <- round(proportion * length(data_p))
    if (ntp < 1) 
        stop("no valid measurments")
    if (ntp == 1) {
        warning(paste("One measurment, Lambda = 1 returned"))
        return(list(estimate = 1, se = 999.99))
    }
    if (ntp < 10) 
        warning(paste("number of points is too small:", ntp))
    if (min(data_p) < 0) 
        stop("data argument has values <0")
    if (max(data_p) <= 1) {
        data_p <- qchisq(data_p, 1, lower.tail = FALSE)
    }
    data_p <- sort(data_p)
    ppoi <- ppoints(data_p)
    Chi2 <- sort(qchisq(1 - ppoi, 1))
    Chi2 <- Chi2[1:ntp]
    x = 1
    if (clust == 1) {
        data.mds0 <- Clust(0)
        k = 2
        cl0 <- ClustN(k, data.mds0)
        data.mds1 <- Clust(1)
        k = 2
        cl1 <- ClustN(k, data.mds1)
        if (length(cl1) > length(cl0)) {
            cl0[(length(cl0) + 1):length(cl1)] = 0
        }
        if (length(cl0) > length(cl1)) {
            cl1[(length(cl1) + 1):length(cl0)] = 0
        }
        S = sum(cl0)
        R = sum(cl1)
        K = sum((cl1 * S - cl0 * R)^2)
        K1 = K/(R * S)
        opt_2 = optimize(GC_VIF, c(0, 1), tol = 1e-04, K = K1)
        F = opt_2$minimum
        K = K1
    }
    if (clust == 0) {
        lda = median(Zx, na.rm = T)/qchisq(0.5, 1) + 0.05
        c = (lda - 1) * (n[1]/2)
        F1 = 1/(c)^0.5
        K1 = c * (1 + F1)/F1
        FK = nlm(GC_VIF_nlm, c(F = F1, K = K1))$estimate
        F = FK[1]
        K = FK[2]
    }
	#
	if (plot){
		ppp=seq(from = 0, to = 1, length = 1000)
		vv=VIF(ppp, n, F, K)
		if (!is.null(lmax)){
			ylim=0
			ylim[1]=1-(lmax*5/90)
			ylim[2]=lmax+(lmax*5/90)
			plot(ppp,vv,main="Lambda",xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color,ylim=ylim)
		} else{
			plot(ppp,vv,main="Lambda",xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color)
		}
		abline(v=min(p,na.rm=TRUE),col="green")
		abline(v=max(p,na.rm=TRUE),col="green")
	}
	#
	vector_vif = VIF(p, n, F, K)
    Zx = Zx/vector_vif
    exeps <- is.na(Zx)
    out = list()
    out$Zx <- Zx
    out$vv <- vector_vif
    out$exeps <- exeps
    out$calrate <- n/max(n)
    out$F = F
    out$K = K
    out
}