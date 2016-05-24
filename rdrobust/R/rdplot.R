### version 0.1  18Nov2013
### version 0.2  26Nov2013
### version 0.3  21Abr2014
### version 0.5  06Jun2014
### version 0.6  17Jun2014
### version 0.61 03Sep2014
### version 0.7  14Oct2014
### version 0.8  04Feb2015

rdplot = function(y, x, subset = NULL, c=0, p=4, nbins=NULL, binselect="esmv", scale=NULL, kernel = "uni", h=NULL, 
                          hide=FALSE, ci=NULL, shade=FALSE, par=NULL, title=NULL, x.label=NULL, y.label=NULL, 
                          x.lim=NULL, y.lim=NULL, col.dots=NULL, col.lines=NULL, type.dots = NULL,...) {

  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  na.ok <- complete.cases(x) & complete.cases(y)
  x <- x[na.ok]
  y <- y[na.ok]
  
  if (is.null(col.lines)) col.lines = "blue"
  if (is.null(col.dots))  col.dots  = 1
  if (is.null(type.dots)) type.dots = 20
    
	x_l = x[x<c]; x_r = x[x>=c]	
  y_l = y[x<c];	y_r = y[x>=c]
	x_min = min(x);	x_max = max(x)
	range_l = c - min(x_l)
	n_l = length(x_l)
	range_r = max(x_r) - c
	n_r = length(x_r)
	n = n_l + n_r
  meth="es"
  
  if (is.null(scale)) {
    scale = scale_l = scale_r = 1  
  } else{
    if (length(scale)==1) scale_l = scale_r = scale
    if (length(scale)==2) {
      scale_l = scale[1]
      scale_r = scale[2]
    }
  }
  
  if (!is.null(nbins)) {
    if (length(nbins)==1) nbins_l = nbins_r = nbins
    if (length(nbins)==2) {
      nbins_l = nbins[1]
      nbins_r = nbins[2]
    }
  }
  
  if (is.null(h)) {
    h_l = range_l
		h_r = range_r
  } else{
    if (length(h)==1) h_l = h_r = h
    if (length(h)==2) {
      h_l = h[1]
      h_r = h[2]
    }
  }
  k=4
  
  #####********************* ERRORS
  exit=0
	if (c<=x_min | c>=x_max){
		print("c should be set within the range of x")
		exit = 1
	}

	if (p<=0 ){
		print("p should be a positive number")
		exit = 1
	}

	if (scale<=0 |scale_l<=0 |scale_r<=0){
		print("scale should be a positive number")
		exit = 1
	}

	p_ceiling = ceiling(p)/p

	if (p_ceiling!=1) {
		print("p should be an integer number")
		exit = 1
	}

	if (exit>0) stop()
	
	###################################################################
  rp_l = matrix(NA,n_l,p+1);  rp_r = matrix(NA,n_r,p+1)
  for (j in 1:(p+1)) {
    rp_l[,j] = x_l^(j-1)
    rp_r[,j] = x_r^(j-1)
  }
  
  wh_l = rdrobust_kweight(x_l,c,h_l,kernel)
  wh_r = rdrobust_kweight(x_r,c,h_r,kernel)
  select_l = wh_l> 0
  select_r = wh_r> 0
  n_h_l=sum(select_l)
  n_h_r=sum(select_r)
  
  gamma_p1_l = qrXXinv((sqrt(wh_l)*rp_l))%*%crossprod(rp_l*wh_l, y_l)	
  gamma_p1_r = qrXXinv((sqrt(wh_r)*rp_r))%*%crossprod(rp_r*wh_r, y_r)
  
  y_hat_l = rp_l%*%gamma_p1_l
  y_hat_r = rp_r%*%gamma_p1_r
  
  X_l = Yhat_l = matrix(NA,n_l,1)
  X_r = Yhat_r = matrix(NA,n_r,1)
  X_l[select_l] = x_l[select_l]
  X_r[select_r] = x_r[select_r]
  Yhat_l[select_l]=y_hat_l[select_l]
  Yhat_r[select_r]=y_hat_r[select_r]
  
  x_sup = c(X_l , X_r)
  y_hat = c(Yhat_l , Yhat_r)	
  
  #**** Optimal Bins (using polynomial order k ***********
  rk_l = matrix(NA,n_l,(k+1))
  rk_r = matrix(NA,n_r,(k+1))
  
  for (j in 1:(k+1)) {
    rk_l[,j] = x_l^(j-1)
    rk_r[,j] = x_r^(j-1)
  }
  
  gamma_k1_l = qrXXinv(rk_l)%*%crossprod(rk_l, y_l)  
  gamma_k2_l = qrXXinv(rk_l)%*%crossprod(rk_l, y_l^2)
  gamma_k1_r = qrXXinv(rk_r)%*%crossprod(rk_r, y_r)  
  gamma_k2_r = qrXXinv(rk_r)%*%crossprod(rk_r, y_r^2)
  
                   
  #*** Bias w/sample
  mu0_k1_l = rk_l%*%gamma_k1_l
  mu0_k1_r = rk_r%*%gamma_k1_r
  mu0_k2_l = rk_l%*%gamma_k2_l
  mu0_k2_r = rk_r%*%gamma_k2_r
  drk_l = matrix(NA,n_l,k)
  drk_r = matrix(NA,n_r,k)
  for (j in 1:k) {
    drk_l[,j] = j*x_l^(j-1)
    drk_r[,j] = j*x_r^(j-1)
  }
                     
  ind_l = order(x_l); ind_r = order(x_r)
  x_i_l = x_l[ind_l] 
  y_i_l = y_l[ind_l]
  
  x_i_r = x_r[ind_r] 
  y_i_r = y_r[ind_r]
                       
  dxi_l=(x_i_l[2:length(x_i_l)]-x_i_l[1:(length(x_i_l)-1)])
  dxi_r=(x_i_r[2:length(x_i_r)]-x_i_r[1:(length(x_i_r)-1)])
  dyi_l=(y_i_l[2:length(y_i_l)]-y_i_l[1:(length(y_i_l)-1)])
  dyi_r=(y_i_r[2:length(y_i_r)]-y_i_r[1:(length(y_i_r)-1)])
                       
  x_bar_i_l = (x_i_l[2:length(x_i_l)]+x_i_l[1:(length(x_i_l)-1)])/2
  x_bar_i_r = (x_i_r[2:length(x_i_r)]+x_i_r[1:(length(x_i_r)-1)])/2
                       
  drk_i_l = matrix(NA,n_l-1,k);	rk_i_l  = matrix(NA,n_l-1,(k+1))
  drk_i_r = matrix(NA,n_r-1,k);	rk_i_r  = matrix(NA,n_r-1,(k+1))
                       
  for (j in 1:(k+1)) {
  rk_i_l[,j] = x_bar_i_l^(j-1)
  rk_i_r[,j] = x_bar_i_r^(j-1)
  }
                       
  for (j in 1:k) {
  drk_i_l[,j] = j*x_bar_i_l^(j-1)
  drk_i_r[,j] = j*x_bar_i_r^(j-1)
  }
  
  mu1_i_hat_l = drk_i_l%*%(gamma_k1_l[2:(k+1)])
  mu1_i_hat_r = drk_i_r%*%(gamma_k1_r[2:(k+1)])
  
  mu0_i_hat_l = rk_i_l%*%gamma_k1_l
  mu0_i_hat_r = rk_i_r%*%gamma_k1_r
  mu2_i_hat_l = rk_i_l%*%gamma_k2_l
  mu2_i_hat_r = rk_i_r%*%gamma_k2_r
                       
  mu0_hat_l = rk_l%*%gamma_k1_l
  mu0_hat_r = rk_r%*%gamma_k1_r
  mu2_hat_l = rk_l%*%gamma_k2_l
  mu2_hat_r = rk_r%*%gamma_k2_r
                       
  mu1_hat_l = drk_l%*%(gamma_k1_l[2:(k+1)])
  mu1_hat_r = drk_r%*%(gamma_k1_r[2:(k+1)])
                       
  mu1_i_hat_l = drk_i_l%*%(gamma_k1_l[2:(k+1)])
  mu1_i_hat_r = drk_i_r%*%(gamma_k1_r[2:(k+1)])
  
  sigma2_hat_l_bar = mu2_i_hat_l - mu0_i_hat_l^2
  sigma2_hat_r_bar = mu2_i_hat_r - mu0_i_hat_r^2
  sigma2_hat_l = mu2_hat_l - mu0_hat_l^2
  sigma2_hat_r = mu2_hat_r - mu0_hat_r^2
  
  J.fun = function(B,V) {ceiling((((2*B)/V)*n)^(1/3))}
  var_y_l = var(y_l)
  var_y_r = var(y_r)
  
  B_es_hat_dw = c( ((c-x_min)^2/(12*n))*sum(mu1_hat_l^2),((x_max-c)^2/(12*n))*sum(mu1_hat_r^2))
  V_es_hat_dw = c((0.5/(c-x_min))*sum(dxi_l*dyi_l^2),(0.5/(x_max-c))*sum(dxi_r*dyi_r^2))
  V_es_chk_dw = c((1/(c-x_min))*sum(dxi_l*sigma2_hat_l_bar),(1/(x_max-c))*sum(dxi_r*sigma2_hat_r_bar))
  J_es_hat_dw = J.fun(B_es_hat_dw, V_es_hat_dw)
  J_es_chk_dw = J.fun(B_es_hat_dw, V_es_chk_dw)
  
  B_qs_hat_dw = c((n_l^2/(24*n))*sum(dxi_l^2*mu1_i_hat_l^2), (n_r^2/(24*n))*sum(dxi_r^2*mu1_i_hat_r^2))
  V_qs_hat_dw = c((1/(2*n_l))*sum(dyi_l^2),(1/(2*n_r))*sum(dyi_r^2))
  V_qs_chk_dw = c((1/n_l)*sum(sigma2_hat_l), (1/n_r)*sum(sigma2_hat_r))
  J_qs_hat_dw = J.fun(B_qs_hat_dw, V_qs_hat_dw)
  J_qs_chk_dw = J.fun(B_qs_hat_dw, V_qs_chk_dw)
  
  J_es_hat_mv  = c(ceiling((var_y_l/V_es_hat_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_es_hat_dw[2])*(n/log(n)^2)))
  J_es_chk_mv  = c(ceiling((var_y_l/V_es_chk_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_es_chk_dw[2])*(n/log(n)^2)))
  J_qs_hat_mv  = c(ceiling((var_y_l/V_qs_hat_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_qs_hat_dw[2])*(n/log(n)^2)))
  J_qs_chk_mv  = c(ceiling((var_y_l/V_qs_chk_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_qs_chk_dw[2])*(n/log(n)^2)))
  
  #########################################################
  if (binselect=="es") {
    J_star_orig = J_es_hat_dw
    meth="es"
    binselect_type="IMSE-optimal evenly-spaced method using spacings estimators"
    J_IMSE = J_es_hat_dw
    J_MV   = J_es_hat_mv
  }
  if (binselect=="espr") {
    J_star_orig = J_es_chk_dw
    meth="es"
    binselect_type="IMSE-optimal evenly-spaced method using polynomial regression"
    J_IMSE = J_es_chk_dw
    J_MV   = J_es_chk_mv
  }
  if (binselect=="esmv" ) {
    J_star_orig = J_es_hat_mv
    meth="es"
    binselect_type="mimicking variance evenly-spaced method using spacings estimators"
    J_IMSE = J_es_hat_dw
    J_MV   = J_es_hat_mv
  }
  if (binselect=="esmvpr" ) {
    J_star_orig = J_es_chk_mv
    meth="es"
    binselect_type="mimicking variance evenly-spaced method using polynomial regression"
    J_IMSE = J_es_chk_dw
    J_MV   = J_es_chk_mv
  }
  if (binselect=="qs" ) {
    J_star_orig = J_qs_hat_dw
    meth="qs"
    binselect_type="IMSE-optimal quantile-spaced method using spacings estimators"
    J_IMSE = J_qs_hat_dw
    J_MV   = J_qs_hat_mv
  }
  if (binselect=="qspr" ) {
    J_star_orig = J_qs_chk_dw
    meth="qs"
    binselect_type="IMSE-optimal quantile-spaced method using polynomial regression"
    J_IMSE = J_qs_chk_dw
    J_MV   = J_qs_chk_mv
  }
  if (binselect=="qsmv" ) {
    J_star_orig = J_qs_hat_mv
    meth="qs"
    binselect_type="mimicking variance quantile-spaced method using spacings estimators"
    J_IMSE = J_qs_hat_dw
    J_MV   = J_qs_hat_mv
  }
  if (binselect=="qsmvpr" ) {
    J_star_orig = J_qs_chk_mv
    meth="qs"
    binselect_type="mimicking variance quantile-spaced method using polynomial regression"
    J_IMSE = J_qs_chk_dw
    J_MV   = J_qs_chk_mv
  }

  J_star_l = scale_l*J_star_orig[1]
  J_star_r = scale_r*J_star_orig[2]

  if (!is.null(nbins)) {
    J_star_l = nbins_l
    J_star_r = nbins_r
    binselect_type="manually evenly spaced"
  }
  
  scale_l = J_star_l / J_IMSE[1]
  scale_r = J_star_r / J_IMSE[2]
  
  bin_x_l = rep(0,length(x_l)); bin_x_r = rep(0,length(x_r))
  jump_l = range_l/J_star_l;jump_r = range_r/J_star_r;
  
  if (meth=="es") {
    jumps_l=seq(min(x_l),max(x_l),jump_l)
    jumps_r=seq(min(x_r),max(x_r),jump_r)
    #binselect_type="Evenly-Spaced"
  }   else if (meth=="qs") {
    jumps_l=quantile(x_l,probs=seq(0,1,1/J_star_l))
    jumps_r=quantile(x_r,probs=seq(0,1,1/J_star_r))
   # binselect_type="Quantile-Spaced"
  }
  
  for (k in 1:(J_star_l-1)) bin_x_l[x_l>=jumps_l[k] & x_l<jumps_l[k+1]] = -J_star_l+k-1 
	bin_x_l[x_l>=jumps_l[(J_star_l)]] = -1
  for (k in 1:(J_star_r-1)) bin_x_r[x_r>=jumps_r[k] & x_r<jumps_r[k+1]] = k 
	bin_x_r[x_r>=jumps_r[(J_star_r)]] = J_star_r
  
  bin_xlmean=bin_ylmean=rep(0,J_star_l)
	bin_xrmean=bin_yrmean=rep(0,J_star_r)
  	for (k in 1:(J_star_l)) {
	  bin_xlmean[k]=mean(c(jumps_l[k],jumps_l[k+1]))
	  #bin_xlmean[k]=mean(x_l[bin_x_l==-k])
	  bin_ylmean[J_star_l-k+1]=mean(y_l[bin_x_l==-k])
  	}
	for (k in 1:(J_star_r)) {
	  bin_xrmean[k]=mean(c(jumps_r[k],jumps_r[k+1]))
	  #bin_xrmean[k]=mean(x_r[bin_x_r==k])
	  bin_yrmean[k]=mean(y_r[bin_x_r==k]) 
	}
	
	bin_xlmean[J_star_l]=mean(c(jumps_l[J_star_l],c))
	bin_xrmean[J_star_r]=mean(c(jumps_r[J_star_r],max(x_r)))
  
  
  bin_x=c(bin_x_l,bin_x_r)
  bin_xmean=c(bin_xlmean,bin_xrmean)
	bin_ymean=c(bin_ylmean,bin_yrmean)
	x_sup = c(x_l, x_r)
	#y_hat = c(mu0_p1_l, mu0_p1_r)
  

	
  if (hide=="FALSE") {
  
  if (is.null(title)){
    title="RD Plot"
  } 
  
  if (is.null(x.label)){
    x.label="X axis"
  }
  
  if (is.null(y.label)){
    y.label="Y axis"
  }
  
  if (is.null(x.lim)){
    x.lim=c(min(x_l),max(x_r))
  }
  
  if (is.null(y.lim)){
    y.lim=c(min(c(y_l,y_r)),max(c(y_l,y_r)))
  }
    par=par
    if (is.null(ci)) {
      plot(bin_xmean,bin_ymean, main=title, xlab=x.label, ylab=y.label, ylim=y.lim, xlim=x.lim, col=col.dots, pch=type.dots,...)
  	  #points(x_l[order(x_l)],mu0_p1_l[order(x_l)],type="l",col=2) 
  	  #points(x_r[order(x_r)],mu0_p1_r[order(x_r)],type="l",col=2)  
      lines(x_l[order(x_l)],y_hat_l[order(x_l)],type="l",col=col.lines) 
      lines(x_r[order(x_r)],y_hat_r[order(x_r)],type="l",col=col.lines)  
      abline(v=c)
    } else {
      
         bin_ySD_l=bin_yN_l=bin_ySD_r=bin_yN_r=0
         for (j in 1:(J_star_l)) {
          bin_ySD_l[j]=sd(y_l[bin_x_l==-j])
          bin_yN_l[j]=length(y_l[bin_x_l==-j])
        }
        
        for (j in 1:(J_star_r)) {
          bin_ySD_r[j]=sd(y_r[bin_x_r==j])
          bin_yN_r[j]=length(y_r[bin_x_r==j])
        }
        
      bin_ySD_l[is.na(bin_ySD_l)]=0
      bin_ySD_r[is.na(bin_ySD_r)]=0
      bin_ySD=c(rev(bin_ySD_l),bin_ySD_r)
      bin_yN=c(rev(bin_yN_l),bin_yN_r)
      quant = -qnorm(((1-(ci/100))/2))
      cil_bin = bin_ymean - quant*bin_ySD/sqrt(bin_yN)
      cir_bin = bin_ymean + quant*bin_ySD/sqrt(bin_yN)
    
       
     if (shade==TRUE){
      plot(bin_xmean,bin_ymean, main=title, xlab=x.label, ylab=y.label, ylim=y.lim, xlim=x.lim, col=col.dots, pch=type.dots,...)
      polygon(c(rev(bin_xmean),bin_xmean),c(rev(cil_bin),cir_bin),col = "grey75")
      lines(x_l[order(x_l)],y_hat_l[order(x_l)],type="l",col=col.lines) 
      lines(x_r[order(x_r)],y_hat_r[order(x_r)],type="l",col=col.lines)  
      abline(v=c)
     } else {
       plot(bin_xmean,bin_ymean, main=title, xlab=x.label, ylab=y.label, ylim=y.lim, xlim=x.lim, col=col.dots, pch=type.dots,...)
       arrows(bin_xmean,cil_bin,bin_xmean,cir_bin,code=3,length=0.1,angle=90,col='grey')
       lines(x_l[order(x_l)],y_hat_l[order(x_l)],type="l",col=col.lines) 
       lines(x_r[order(x_r)],y_hat_r[order(x_r)],type="l",col=col.lines)  
       abline(v=c)
     }
    }
    
    }

    tabl1.str=matrix(NA,14,2)
    tabl1.str[1,]  = formatC(c(n_l,n_r),digits=0, format="f")
    tabl1.str[2,]  = formatC(c(p,p),digits=0, format="f")
    tabl1.str[3,]  = formatC(c(scale_l,scale_r),digits=0, format="f") 
    tabl1.str[4,]  = c("","")
    tabl1.str[5,]  = formatC(c(J_star_l,J_star_r),digits=0, format="f")  
    tabl1.str[6,]  = formatC(c(jump_l,jump_r),digits=4, format="f") 
    tabl1.str[7,]  = c("","")
    tabl1.str[8,]  = formatC(c(J_IMSE),digits=0, format="f")
    tabl1.str[9,]  = formatC(c(J_MV),digits=0, format="f")
    tabl1.str[10,]  = c("","")
    tabl1.str[11,]  = c("","")
    tabl1.str[12,]  = formatC(c(scale_l,scale_r),digits=4, format="f")
    tabl1.str[13,]  = formatC(c(1/(1+scale_l^3), 1/(1+scale_r^3)),digits=4, format="f")
    tabl1.str[14,] = formatC(c(scale_l^3/(1+scale_l^3), scale_r^3/(1+scale_r^3)),digits=4, format="f")

    rownames(tabl1.str)=c("Number of Obs.","Polynomial Order","Scale", "","Selected Bins","Bin Length","", "IMSE-optimal bins","Mimicking Variance bins","","Relative to IMSE-optimal:","Implied scale","WIMSE variance weight","WIMSE bias weight")
    colnames(tabl1.str)=c("Left","Right")
    
    results=matrix(NA,10,2)
    results[1,] = c(n_l,n_r)
    results[2,] = c(p,p)
    results[3,] = c(scale_l,scale_r)
    results[4,] = c(J_star_l,J_star_r)
    results[5,] = c(jump_l,jump_r)
    results[6,] = J_IMSE
    results[7,] = J_MV
    results[8,] = c(scale_l,scale_r)
    results[9,] = c(1/(1+scale_l^3), 1/(1+scale_r^3))
    results[10,] = c(scale_l^3/(1+scale_l^3), scale_r^3/(1+scale_r^3))
    rownames(results)=c("Number of Obs.","Polynomial Order","Chosen Scale","Selected bins","Bin Length","IMSE-optimal bins","Mimicking Variance bins","Implied scale","WIMSE variance weight","WIMSE bias weight")
    colnames(results)=c("Left","Right")
    
    coef=matrix(NA,p+1,2)
    coef[,1] = c(gamma_p1_l)
    coef[,2] = c(gamma_p1_r)
    colnames(coef)=c("Left","Right")
    out=list(method=binselect_type,results=results,coef=coef,tabl1.str=tabl1.str)
    
    out$call <- match.call()
    class(out) <- "rdplot"
    return(invisible(out))
}

#rdplot <- function(y,x, ...) UseMethod("rdplot")

#rdplot.default <- function(y,x, ...){
#  est <- rdplotEst(y,x, ... )
#  est$call <- match.call()
#  class(est) <- "rdplot"
#  est
#}

print.rdplot <- function(x,...){
  cat("Call:\n")
  #print(x$call)
  cat(deparse(x$call, width.cutoff=getOption("width")), sep = "\n")
  cat("\n")
  #cat(paste("Method: ",x$method))
  cat(strwrap(paste("Method: ", x$method)), sep = "\n")
  cat("\n\n")
  print(x$tabl1.str,quote=F)
}

summary.rdplot <- function(object,...) {
  TAB <- object$results
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.rdplot"
  res
}

#print.summary.rdplot <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
#}

