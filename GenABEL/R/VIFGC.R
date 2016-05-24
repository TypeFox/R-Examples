#' Genomic control for various model of inheritance using VIF
#' 
#' This function estimates corrected statistic using genomic control
#' for different models (recessive, dominant, additive etc.),
#' using VIF. VIF coefficients are estimated
#' by optimizing different error functions: regress,
#' median and ks.test.
#'  
#' @param data Input vector of Chi square statistic
#' @param method Function of error to be optimized. Can be
#' "regress", "median" or "ks.test"
#' @param p Input vector of allele frequencies
#' @param x Model of inheritance (0 for recessive,0.5 for additive, 1 for dominant, also it could be arbitrary)
#' @param index.filter Indexes for variables that will be use for analysis in data vector
#' @param n The size of the sample
#' @param proportion The proportion of lowest P (Chi2) to be used when
#'   estimating the inflation factor Lambda for "regress" method only
#' @param plot If TRUE, plot of lambda will be produced
#' @param type_of_plot For developers only
#' @param lmax The threshold for lambda for plotting (optional)
#' @param color The color of the plot
#' @param F The estimation of F (optional)
#' @param K The estimation of K (optional)
#' @param ladd The estimation of lambda for additive model (optional)
#' @param clust For developers only
#' @param vart0 For developers only
#' @param tmp For developers only
#' @param CA For developers only
#' @param p.table For developers only
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
#' qts <- mlreg(dm2~sex,data=ge03d2,gtmode = "dominant")
#' chi2.1df <- results(qts)$chi2.1df
#' s <- summary(ge03d2)
#' freq <- s$Q.2
#' result <- VIFGC(p=freq,x=1,method = "median",CA=FALSE,data=chi2.1df,n=nids(ge03d2))
#' 
#' @keywords htest
#'

VIFGC=function (data, p, x, method = "regress", n, index.filter = NULL,
    proportion = 1, clust = 0, vart0 = 0, tmp = 0, CA = FALSE,
    p.table = 0,plot=TRUE,lmax=NULL,color="red",F=NULL,K=NULL,type_of_plot="plot",ladd=NULL)
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
    if (CA) {
        p00 <- p.table[1, 1, ]
        p01 <- p.table[1, 2, ]
        p02 <- p.table[1, 3, ]
        n0 <- p.table[1, 4, ]
        p10 <- p.table[2, 1, ]
        p11 <- p.table[2, 2, ]
        p12 <- p.table[2, 3, ]
        n1 <- p.table[2, 4, ]
        p0 <- p.table[3, 1, ]
        p1 <- p.table[3, 2, ]
        p2 <- p.table[3, 3, ]
        n <- p.table[3, 4, ]
        Zx = 0
        exeps = 0
        Dx = (p12 - p02) + x * (p11 - p01)
        i = 1
        while (i <= length(n)) {
            vart = ((p2[i] + p1[i] * x^2) - (p2[i] + x * p1[i])^2) *
                ((1/n1[i]) + (1/n0[i]))
            if (vart == 0) {
                Zx[i] = NA
                exeps[i] = T
            }
            else {
                Zx[i] = ((Dx[i])^2)/(vart)
                exeps[i] = F
            }
            i = i + 1
        }
        inf = which(!is.na(Zx))
        notinf = which(is.na(Zx))
    }
    else {
        Zx = data
        inf = which(!is.na(Zx))
        notinf = which(is.na(Zx))
    }
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
        Var = ((2 * p * (1 - F) * q * x^2 + F * p + (1 - F) *
            p^2) - (2 * p * q * (1 - F) * x + F * p + (1 - F) *
            p^2)^2)
        S = (1 + F) * (1 + 2 * F)
        Cov1 = ((6 * F^3 * p + 11 * (1 - F) * F^2 * p^2 + (1 -
            F)^3 * p^4 + 6 * (1 - F)^2 * F * p^3)/S) - ((1 -
            F) * p^2 + F * p)^2
        Cov2 = x * (((4 * (2 * (1 - F) * F^2 * p * q + (1 - F)^3 *
            p^3 * q + 3 * (1 - F)^2 * F * p^2 * q))/S) - 2 *
            ((1 - F) * p^2 + F * p) * (2 * (1 - F) * p * q))
        Cov3 = x^2 * ((4 * ((1 - F)^3 * p^2 * q^2 + F * (1 -
            F) * p * q))/S - (2 * (1 - F) * p * q)^2)
        Cov = Cov1 + Cov2 + Cov3
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
	lambda=(median(data, na.rm = T)/qchisq(0.5, 1))
    data_p <- sort(data_p)
    ppoi <- ppoints(data_p)
    Chi2 <- sort(qchisq(1 - ppoi, 1))
    Chi2 <- Chi2[1:ntp]
    if (clust == 1) {
			F = 0.5
			K = n[1]
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
       
		if (!is.null(F) & !is.null(K)){
			FK=0
			FK[1]=F
			FK[2]=K
		} else{

			if ((median(Zx, na.rm = T)/qchisq(0.5, 1)) >= 1) {
				if (!is.null(ladd)) {
					lda=ladd
				}else{
					lda = median(Zx, na.rm = T)/qchisq(0.5, 1)
				}
				c = (lda - 1) * (n[1]/2)
				F1 = 1/(c)^0.5
				K1 = c * (1 + F1)/F1
				
			}
			else {
				F1 = 0.5
				K1 = n[1]
			}
			FK = nlm(GC_VIF_nlm, c(F = F1, K = K1))$estimate
			F = FK[1]
			K = FK[2]
		}
    }
	
	
	if (plot){
		if (type_of_plot=="plot"){
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
			abline(h=lambda,col=color)
		} else{
			ppp=seq(from = 0, to = 1, length = 1000)
			vv=VIF(ppp, n, F, K)
			lines(ppp,vv,xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color)
			abline(h=lambda,col=color)
		}
		
	}
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
