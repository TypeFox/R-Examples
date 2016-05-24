#' Polynomial genomic control
#' 
#' This function estimates the genomic controls
#' for different models and degrees of freedom,
#' using polinomial function. Polinomial coefficients are estimated
#' by optimizing different error functions: regress,
#' median, ks.test or group regress.
#' 
#' @param data Input vector of Chi square statistic
#' @param method Function of error to be optimized. Can be
#'  "regress", "median", "ks.test" or "group_regress"
#' @param p Input vector of allele frequencies
#' @param df Number of degrees of freedom
#' @param pol.d The degree of polinomial function
#' @param plot If TRUE, plot of lambda will be produced
#' @param start.corr For regress method use it only when you want to make calculations faster
#' @param index.filter Index of variables in data vector, that will be used in analysis 
#'  if zero - all variables will be used
#' @param proportion The proportion of lowest P (Chi2) to be used when
#'   estimating the inflation factor Lambda for "regress" method only
#' @param n_quiantile The number of groups for "group_regress" method
#' @param title_name The title name for plot
#' @param type_of_plot For developers only
#' @param lmax The threshold for lambda for plotting (optional)
#' @param color The color of the plot
#' 
#' @return A list with elements
#' \item{data}{Output vector corrected Chi square statistic}
#' \item{b}{Polinomial coefficients}
#' 
#' @author Yakov Tsepilov
#' 
#' @examples
#' data(ge03d2)
#' ge03d2 <- ge03d2[seq(from=1,to=nids(ge03d2),by=2),seq(from=1,to=nsnps(ge03d2),by=3)]
#' qts <- mlreg(dm2~1,data=ge03d2,gtmode = "additive")
#' chi2.1df <- results(qts)$chi2.1df
#' s <- summary(ge03d2)
#' freq <- s$Q.2
#' result=PGC(data=chi2.1df,method="median",p=freq,df=1, pol.d=2, plot=TRUE, lmax=1.1,start.corr=FALSE)
#' #"group_regress" is better to use when we have more than 50K SNPs
#' #result=PGC(data=chi2.1df,method="group_regress",p=freq,df=1, pol.d=2, plot=TRUE, start.corr=FALSE,n_quiantile=3)
#' 
#' @keywords htest
#' 

PGC=function (data, method = "group_regress", p, df, pol.d = 3, plot = TRUE, 
    index.filter = NULL, start.corr = FALSE, proportion = 1,n_quiantile=5,title_name="Lambda",type_of_plot="plot",lmax=NULL,color="red") 
{
    if (is.null(index.filter)) {
        ind.function = which(!is.na(data))
    }
    else ind.function = index.filter
    ind.function = ind.function[which(!is.na(data[ind.function]))]
	
    pol.m = function(p, b, pol.d) {
        out = 0
        for (i in (1:pol.d)) {
            out = (p^i) * b[i] + out
        }
        out = out + b[pol.d + 1]
        return(out)
    }
	
    if (!(method == "regress" | method == "median" | method == "ks.test" | method == "group_regress")) {
        print("Error. I do not know this method")
        break
    }
	
    if (start.corr) {
        data = data/(1 * median(data[ind.function], na.rm = T)/qchisq(0.5,df))
    }
    data_p <- data[ind.function]
	p_p=p[ind.function]
	lambda=median(data_p)/qchisq(0.5, df)
    
			if (proportion > 1 || proportion <= 0) stop("proportion argument should be greater then zero and less than or equal to one")
			ntp <- round(proportion * length(data_p))
			if (ntp < 1) stop("no valid measurments")
			if (ntp < 100) warning(paste("number of points is too small:", ntp))
			if (min(data_p) < 0) stop("data argument has values <0")
			
    data_p <- sort(data_p)
    ppoi <- ppoints(data_p)
	if (method=="regress"){
		Chi2 <- sort(qchisq(ppoi,df=df,lower.tail=FALSE))
		Chi2 <- Chi2[1:ntp]
	}
	if (method=="group_regress"){
		qt=quantile(p_p, probs = seq(0, 1,length=n_quiantile))
		pp=p_p
		i=1
		for (i in (2:n_quiantile)){
			pp[p_p>=qt[i-1] & p_p<=qt[i]]=i-1
		}
		tb=table(pp)
		chi2_array=array(NA,c((n_quiantile-1),max(tb)))
		for (i in (1:(n_quiantile-1))){
			ppoi <- ppoints(1:tb[i])
			ntp[i] <- round(proportion * tb[i])
			chi2_array[i,1:tb[i]] <- sort(qchisq(ppoi,df=df,lower.tail=FALSE))
		}
	}
    ppp_real = seq(from = min(p), to = max(p), length = 1000)
	ppp_abstr= seq(from = 0, to = 1, length = 1000)
	#ppp = seq(from = 0, to = 1, length = 1000)
	
    f2 = function(b) {
        Zx = data
		#if ((min(pol.m(ppp_abstr, b, pol.d))<1) & (max(pol.m(ppp_abstr, b, pol.d))>lambda*2)){
		if (min(pol.m(ppp_abstr, b, pol.d))<1){
			vv=1e5
			#vv=pol.m(p, b, pol.d)
		} else{
			if (max(pol.m(ppp_abstr, b, pol.d))>10){
				vv=1e5*max(pol.m(ppp_abstr, b, pol.d))
			#vv=pol.m(p, b, pol.d)
			}else{
				vv=pol.m(p, b, pol.d)
			}
		}
		
        Zx = Zx/vv
		if (method == "group_regress"){
			Zxl=Zx[ind.function]
			F=0
			for (i in (1:(n_quiantile-1))){
				Zxl_r=sort(Zxl[pp==i])
				Zxl_r = Zxl_r[1:ntp[i]]
				Chi2=chi2_array[i,1:ntp[i]]
				F = sum((Zxl_r - Chi2)^2)+F
			}
		}
        if (method == "ks.test") {
			Zxl = sort(Zx[ind.function])
            F = -log(ks.test(Zxl, "pchisq", df = df)$p.value)
        }
        if (method == "median") {
			 Zxl = sort(Zx[ind.function])
            F = abs(median(Zxl) - qchisq(0.5, df))
        }
        if (method == "regress") {
			Zxl = sort(Zx[ind.function])
			Zxl_r = Zxl[1:ntp]
            F = sum((Zxl_r - Chi2)^2)
        }
        return(F)
    }
    b = rep(0, pol.d + 1)
    b[pol.d + 1] = 1
    nlm = nlm(f2, b)
    b = nlm$estimate
    if (plot) {
		if (type_of_plot=="plot"){
			vv=pol.m(ppp_abstr, b, pol.d)
			if (!is.null(lmax)){
				ylim=0
				ylim[1]=1-(lmax*5/90)
				ylim[2]=lmax+(lmax*5/90)
				plot(ppp_abstr,vv,main=title_name,xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color,ylim=ylim)
			} else{
				plot(ppp_abstr,vv,main=title_name,xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color)
			}
			#ylim=0
			#ylim[1]=1-(lmax*5/90)
			#ylim[2]=lmax+(lmax*5/90)
			#plot(ppp_abstr,vv,main="Lambda",xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color,ylim=ylim)
			abline(v=min(p,na.rm=TRUE),col="green")
			abline(v=max(p,na.rm=TRUE),col="green")
		} else{
			lines(ppp_abstr,pol.m(ppp_abstr, b, pol.d),main=title_name,xlab="Frequency",ylab="lambda(freq)",typ = "l",col=color)
		}
	}
    vv = pol.m(p, b, pol.d)
	data = data/vv
    out = list()
    out$data = data
    out$b = b
    return(out)
}

