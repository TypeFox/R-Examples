stats <- function (x, by, quantiles=c(.25,.75))

  ## an improved descriptive statistics ... version 6/3/12
  ## Dave McArthur, David Geffin School of Medicine at UCLA   dmca <at> ucla.edu
  ## adapted from http://legacy.ncsu.edu/ST370/distance/rlab/
 	##    
  ##   'by' grouping variable presumes that x is univariate
  ##   'quantiles' can be any pair of values >0 : <1
{
  if (!missing(by)) {               
    if (!is.numeric(x)) x<-as.numeric(simplify2array(x))
	if (!is.numeric(by)) by<-as.numeric(simplify2array(by))
		 a <- as.character(by)
		    for (k in 1:length(a)) {		
			q<-a[[k]]
            q <- q[is.na(q)]<- 'missing' # work around for missing 'recode'
			a[[k]]<- q # recode(q,"NA='missing'")  #car package function	
			}
		 by<-a
      x<-.cat2list(x[order(by)],sort(by))
  }
  QA<-noquote(paste('Q:',as.character(quantiles[1]),sep=''))
  QB<-noquote(paste('Q:',as.character(quantiles[2]),sep=''))
    desc <- c("N", "mean", "SD", "robust min", "min", QA, "median",
		  QB, "max", "robust max",
		  "skew", "kurtosis","Huber mu","missing values","unique values","Shapiro p")  
  if (!is.list(x) & !is.matrix(x))
      x <- matrix(x, ncol = 1)
  if (is.list(x)) {
      nc <- length(x)
      out <- matrix(NA, ncol = nc, nrow = length(.description()))
      dimnames(out) <- list(desc, names(x))
      for (j in (1:nc)) {
        if (!is.factor(x[[j]])) {
          if (is.numeric(x[[j]])) {
              out[, j] <- .description(x[[j]],quantiles)  #list
          }
        }
      }
      return(noquote(out))
  }
  if (is.matrix(x)) {
      nc <- ncol(x)
      out <- matrix(NA, ncol = nc, nrow = length(.description()))
      dimnames(out) <- list(desc, dimnames(x)[[2]])
      for (k in (1:nc)) {
         if (!is.factor(x[[k]])) {
          out[, k] <- .description(x[, k],quantiles)  #matrix
         }
      }
      return(noquote(out))
  }
}

.description <- function (x,quantiles=c(.25,.75))
  {
	QA<-noquote(paste('Q:',as.character(quantiles[1]),sep=''))
  QB<-noquote(paste('Q:',as.character(quantiles[2]),sep=''))
    lab <- c("N", "mean", "SD", "robust min", "min", QA, "median",
		  QB, "max", "robust max",
		  "skew", "kurtosis","Huber mu","missing values","unique values","Shapiro p")
	  if (missing(x)) {
		  return(lab)
	  }
	  temp <- rep(0, length(lab))
	  xt <- x[!is.na(x)]
	  n <- length(xt)

	  if (!is.numeric(xt) || all(is.na(x))) {
		  return(c(n, rep(NA, 12), length(x), NA, NA))
	  }
	  else {

	  if (n < 4) {
		  xt <- sort(xt)
		  if (n == 0)  return(sprintf("%10.4f",length(x),rep(NA,15)))
		  if (n == 1) {
			  return(sprintf("%10.4f",c(n, xt[1], NA, rep(c(NA,xt[1]), 3), rep(NA,4), length(x) -
				  length(xt),length(unique(xt)),NA)))
		  }
		  if (n == 2) {
			  return(sprintf("%10.4f",c(n, mean(xt), sqrt(var(xt)), .rrange(xt)[1], c(xt[1], xt[1],
				  mean(xt), xt[2], xt[2]), .rrange(xt)[2], rep(NA,3), length(x) - length(xt),length(unique(xt)),NA)))
		  }
		  if (n == 3) {
			  return(sprintf("%10.4f",c(n, mean(xt), sqrt(var(xt)), .rrange(xt)[1], c(xt[1], xt[1],
				  xt[2], xt[3], xt[3]), .rrange(xt)[2], rep(NA,3), length(x) - length(xt),length(unique(xt)),NA)))
		  }
	  }
	  else {
		 if (length(unique(x[!is.na(x)]))>1) {
			 if(length(x)<5001) {    shapp<-shapiro.test(x)  }}
		 else  { shapp<-NULL              #trap shapiro's with no variability
			shapp$p<-NA
			}
		 if(length(x)>5000) { shapp<-NULL
			shapp$p<-NA    }

      quantiles<-append(quantiles[1],c(.5,quantiles[2]))
	  return(sprintf("%10.4f",c(length(xt), mean(xt), sqrt(var(xt)), .rrange(xt)[1], min(xt),
		 quantile(xt, quantiles), max(xt), .rrange(xt)[2], .skewed(xt), .kurtosis(xt),.huber.mu(xt),
		 length(x) - length(xt),length(unique(xt)),shapp$p)))     #
	  }
	}
  }

.kurtosis <- function (x)				            			# from Tom Fletcher
	{
		 n <- length (x[!(is.na(x))])
		 sd <- sqrt(var(x,na.rm=T))
		 m <- mean(x,na.rm=T)
		 (((n*(n+1))/((n-1)*(n-2)*(n-3)))*(sum(((x-m)/sd)^4, na.rm=T)))-((3*(n-1)^2)/((n-2)*(n-3)))
	}

.mad <- function (x, center = median(x), constant = 1.4826, na.rm = FALSE,
		low = FALSE, high = FALSE)
	{
		if (na.rm)
			x <- x[!is.na(x)]
		n <- length(x)
		constant * if ((low || high) && n%%2 == 0) {
			if (low && high)
				stop("'low' and 'high' cannot be both TRUE")
			n2 <- n%/%2 + as.integer(high)
			sort(abs(x - center), partial = n2)[n2]
		}
		else median(abs(x - center))
	}

.rrange <- function (x, range = 1, coef = 1.5, na.rm = TRUE)
  {                                    					# from package sfsmisc
	if (!missing(range)) {
		if (!missing(coef))
			stop("Must use either 'range' or 'coef'")
		coef <- 1.5 * range
	}
	if (!na.rm && any(is.na(x)))
		return(0 + c(NA, NA))
	boxplot.stats(x, coef = coef, do.conf = FALSE, do.out = FALSE)$stats[c(1,5)]
  }

.skewed <- function (x)								# from Tom Fletcher, University of Missouri - St. Louis
	{      												     
  	 n <- length (x[!(is.na(x))])
		 sd <- sqrt(var(x,na.rm=T))
		 m <- mean(x,na.rm=T)
		 (n/((n-1)*(n-2)))*sum(((x-m)/sd)^3, na.rm=T)
	}
	
.cat2list<-function (x, a)            			# from Institute for Mathematics Applied Geosciences
	{ 									      		           	# University Corporation for Atmospheric Research
      if(is.character(levels(a))){
      label <- levels(a)
      out <- as.list(1:length(label))
      names(out) <- levels(a)
      }
      else {a <- as.character(a)
      label <- unique(a)
      out <- as.list(1:length(label))
      names(out) <- label
      }

      for (k in 1:length(label)) {
          out[[k]] <- x[label[k] == a]
      }
      out
  }

.huber.mu <- function (x, c = 1.28, iter = 20, conv = 1e-07) # from package asbio
  {
      mu.hat <- .huber.NR(x, c, iter)
          mu.est = max(mu.hat)
          mu.est
  }
.huber.NR <- function (x, c = 1.28, iter = 20)
  {
      require(MASS)
      mu.k <- matrix(nrow = iter, ncol = 1)
      mu.k[1] <- median(x)
      for (i in 1:iter) {
         # if (..mad(x) == 0)
          #    mu.k[1] = NA
          # if (..mad(x) != 0)
          {
              A1 <- (x - mu.k[i])/.mad(x)
              A <- sum(sapply(A1, function(x) {
                  max(-c, min(c, x))
              }))
              B1 <- (A1 >= c | A1 <= -c)
              B <- length(B1[B1 == FALSE])
              mu.k[i + 1] <- mu.k[i] + ((.mad(x) * A)/B)
          }
      }
      mu.k
  }
    
    

.kurtosys <-function (x)
{ # from Tom Fletcher
  n <- length (x[!(is.na(x))])
  sd <- sqrt(var(x,na.rm=T))
  m <- mean(x,na.rm=T)
  (((n*(n+1))/((n-1)*(n-2)*(n-3)))*(sum(((x-m)/sd)^4, na.rm=T)))-((3*(n-1)^2)/((n-2)*(n-3)))
}

.mad <- function (x, center = median(x), constant = 1.4826, na.rm = FALSE, 
                  low = FALSE, high = FALSE) 
{
  if (na.rm) 
    x <- x[!is.na(x)]
  n <- length(x)
  constant * if ((low || high) && n%%2 == 0) {
    if (low && high) 
      stop("'low' and 'high' cannot be both TRUE")
    n2 <- n%/%2 + as.integer(high)
    sort(abs(x - center), partial = n2)[n2]
  }
  else median(abs(x - center))
}

.rrange <-function (x, range = 1, coef = 1.5, na.rm = TRUE)
{                                    ### robust range from package sfsmisc Zurich
  if (!missing(range)) {
    if (!missing(coef))
      stop("Must use either 'range' or 'coef'")
    coef <- 1.5 * range
  }
  if (!na.rm && any(is.na(x)))
    return(0 + c(NA, NA))
  boxplot.stats(x, coef = coef, do.conf = FALSE, do.out = FALSE)$stats[c(1,5)]
} 


# Tom Fletcher, University of Missouri - St. Louis
# http://www.umsl.edu/~fletchert/quant/DataScreen.txt
.skewed <-function (x)
{
  n <- length (x[!(is.na(x))])
  sd <- sqrt(var(x,na.rm=T))
  m <- mean(x,na.rm=T)
  (n/((n-1)*(n-2)))*sum(((x-m)/sd)^3, na.rm=T)
}


# from package asbio
# Maintainer: Ken Aho <kenaho1@gmail.com>
.ci.median<-function(x,conf=.95){
  n<-nrow(as.matrix(x))
  if(qbinom((1-conf)/2,n,0.5)==0)
    stop("CI not calculable")
  L<- qbinom((1-conf)/2,n,0.5)
  U<-n-L+1
  if(L>=U)
    stop("CI not calculable")
  order.x<-sort(x)
  res<-list()
  res$head<-paste(paste(as.character(conf*100),"%",sep=""),
                  c("Confidence Interval for Population Median"))
  res$ci<-c(median=median(x),lower=order.x[L],upper=order.x[n-L+1])
  res$ends<-c("Estimate",
              paste(as.character(c((1-conf)/2,1-((1-conf)/2))*100),"%",sep=""))
  res$coverage<-1-(2*pbinom(q=L-1,n,0.5))
  class(res)<-"ci"
  res
}
    
    
    
violins <- function (x, by, range = 1.5, h = NULL, ylim = NULL, names = NULL,
    horizontal = FALSE, col = "transparent", border = "black", lty = 1,
    lwd = 1, rectCol = "grey50", colMed = "grey80", pchMed = 19,
    at, add = FALSE, wex = 1, drawRect = TRUE, main = "", xlab="",ylab="",
    connect = c('median', 'mean', 'hubermu','deciles'), SD.or.SE = c('SD'),
    connectcol = c('lightblue','cyan','darkred','grey'),  # colors for median, mean, hubermu, deciles
    las=2, stats=FALSE, quantiles=c(.1,.9), CImed = TRUE, deciles = TRUE)
    
  ## an improved violin plot ... version 6/4/12
  ## Dave McArthur, David Geffin School of Medicine at UCLA   dmca <at> ucla.edu
  ## adapted from package 'caroline'  
  ## David M. Schruth
  ## who migrated the function out of the vioplot library
  ##
  ## Includes 95% ci.median and huber.mu M-estimator, 
  ##     both from package asbio (see Applied Statistics for Biologists)
  ##   and handles variables that contain no data or are factors
  ##
  ## Provide either a list or dataframe as input, with names or not (list need not be rectilinear).
  ##   'by' grouping variable presumes that x is univariate and length(grouping var)=length(x).
  ##   'at' can be specific list of positions or single incrementer (but is not needed initially).
  ##   'add' enables overplotting [recommendation: when add=TRUE set 'at' to a small positive increment, 
  ##         e.g., 0.15, and set border to a color other than default black.
  ##   'SD.or.SE' plots s.d. or s.e. as wider line over the usual whisker, and can remain blank.
  ##   'median' is plotted as dot, 'mean' is plotted as X, 'hubermu' is plotted as square.
  ##   'connect' and 'connectcol' are active when 2 or more variables are plotted or can be blank.
  ##   'CImed' portrays 95% confidence intervals for the median (as solid box)
  ##   'quantiles' map any pair of quantiles (as dotted box) in addition to Q1 & Q3, 
  ##         but are not shown when "c(0,0)" and arg is passed to descriptive stats when 'stats'=TRUE.
  ##   'deciles' maps deciles 0.1:0.9 (as thin lines) independently of 'quantiles' and
  ##         can be connected when 2 or more variables are plotted.

{
    options(warnings=-1)
    require(sm)
    if(is.data.frame(x)) x<-as.list.data.frame(x)  ## convert dataframe to list if needed
    if (!missing(by)) {                            ## cope with 'by' variable
    if(is.numeric(by)) x<-.cat2list(x[order(by)],sort(by))  ## if 'by' is numeric so sort it
    if(!is.numeric(by)) x <- .cat2list(x, by)               ## if 'by' is categorical
                      }
    if(is.list(x)){
    	datas <- x
    	if(length(names) == 0)
    		names <- names(x)
    }else{
    	datas <- list(x)
    }
    n <- length(datas)
    if (missing(at))
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q.1 <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    q.9 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    hubermu <- vector(mode = "numeric", length = n)
    average <- vector(mode = "numeric", length = n)
    stddevlower <- vector(mode = "numeric", length = n)
    stddevupper <- vector(mode = "numeric", length = n)
    stderrlower <- vector(mode = "numeric", length = n)
    stderrupper <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    medCI05 <- vector(mode = "list", length = n)
  	medCI95 <- vector(mode = "list", length = n)
  	decile <- matrix(NA, nrow=n,ncol=9)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h)))
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- (datas[[i]])
        data.min <- min(data,na.rm=TRUE)
        data.max <- max(data,na.rm=TRUE)
        q.1[i] <- quantile(data, quantiles[1],na.rm=TRUE)
        q1[i] <- quantile(data,.25,na.rm=TRUE)
        q3[i] <- quantile(data,.75,na.rm=TRUE)
        q.9[i] <- quantile(data, quantiles[2],na.rm=TRUE)
        med[i] <- median(data,na.rm=TRUE)
        medCI05[i] <- .ci.median(data)$ci[2]
        medCI95[i] <- .ci.median(data)$ci[3]
        hubermu[i] <- .huber.mu(data)
        average[i] <- mean(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        stddevlower[i] <- average[i]-sd(data)
        stddevupper[i] <- average[i]+sd(data)
		if (deciles) for (j in 1:9) decile[i,j] <- quantile(data,j/10)
        N <- length(data)
        stderrlower[i] <- average[i]-(sd(data)/sqrt(N))
        stderrupper[i] <- average[i]+(sd(data)/sqrt(N))
        est.xlim <- c(min(lower[i], data.min), max(upper[i],
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
            args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1)
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
        if(length(at)==1) at<- 1:n+at      ## handle 'at' as incrementer
    }
    boxwidth <- 0.05 * wex
    if (!add)
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim, las=las)
            axis(2,las=las)
            axis(1, at = at, labels = label, las=las)
            title(main,xlab=xlab,ylab=ylab)
        }
        box()
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
                c(base[[i]], rev(base[[i]])), col = col[i], border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
    			if (deciles) for (j in 1:9) 
    			rect(at[i] - boxwidth*wex, decile[i,j],
    			    at[i] + boxwidth*wex, decile[i,j], lwd=.3*lwd)
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
                  lty = lty)
                rect(at[i] - boxwidth*wex, q.1[i], at[i] + boxwidth*wex,
                  q.9[i], col = 'transparent',lty=3)
                rect(at[i] - boxwidth/3*wex, q1[i], at[i] + boxwidth/3*wex,
                  q3[i], col = rectCol)
                if(any(SD.or.SE=='SD')) lines(at[c(i+.05, i+.05)], c(stddevlower[i], 
                    stddevupper[i]), lwd=lwd*4*wex, lty=lty)
                if(any(SD.or.SE=='SE')) lines(at[c(i+.05, i+.05)], c(stderrlower[i],
                    stderrupper[i]), lwd=lwd*4*wex, lty=lty)
                points(at[i], med[i], pch = pchMed, col = colMed)
                if (CImed) rect(at[i] -boxwidth/1.6*wex, medCI05[i], at[i]+boxwidth/1.6*wex,
                  medCI95[i])
   
                points(at[i], hubermu[i], pch=12, col=colMed)
                points(at[i], average[i], pch=13, col=colMed)
            }
    			s <- seq(length(datas))
    			s <- s[-length(s)]
  if (any(connect=='median'))  segments(at[s], med[s], at[s+1], med[s+1], col= connectcol[1])
  if (any(connect=='hubermu')) segments(at[s], hubermu[s], at[s+1], hubermu[s+1], col=connectcol[2])
  if (any(connect=='mean')) segments(at[s], average[s], at[s+1], average[s+1], col=connectcol[3])
  if (deciles & any(connect=='deciles')) for (j in 1:9) segments(at[s], decile[s,j], at[s+1], 
                decile[s+1,j],lwd=.6*lwd,col=connectcol[4])
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim, las=las)
            axis(1,las=las)
            axis(2, at = at, labels = label,las=las)
        }
        box()

        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                rev(at[i] + height[[i]])), col = col[i], border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
                if (deciles) for (j in 1:9) rect(decile[i,j], 
                  at[i] - boxwidth*wex, decile[i,j], at[i] + boxwidth*wex, lwd=.5*lwd) 
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
                  lty = lty)
                rect(q.1[i], at[i] - boxwidth*wex, q.9[i], at[i] + boxwidth*wex,
                  col = 'transparent', lty=3)
                rect(q1[i], at[i] - boxwidth/3*wex, q3[i], at[i] + boxwidth/3*wex,
                  col = rectCol)
                if(any(SD.or.SE=='SD')) lines(c(stddevlower[i],  stddevupper[i]), 
                  at[c(i+.05, i+.05)], lwd=lwd*4*wex, lty=lty)
                if(any(SD.or.SE=='SE')) lines(c(stderrlower[i], stderrupper[i]),
                  at[c(i+.05, i+.05)], lwd=lwd*4*wex, lty=lty)
                if (CImed) rect(medCI05[i], at[i] -boxwidth/1.6*wex,
                  medCI95[i], at[i]+boxwidth/1.6*wex)

               # if (deciles) lines(decile[i], at[c(i-boxwidth/3*wex, i+boxwidth/3*wex)],lty=3)

                points(med[i], at[i], pch = pchMed, col = colMed)
                points(average[i], at[i], pch=13, col=colMed)
            }
                s <- seq(length(datas))
			    s <- s[-length(s)]
  if (any(connect=='median'))  segments(med[s], at[s], med[s+1], at[s+1], col=connectcol[1])
  if (any(connect=='hubermu')) segments(hubermu[s], at[s], hubermu[s+1], at[s+1], col=connectcol[2])
  if (any(connect=='mean')) segments(average[s], at[s], average[s+1], at[s+1], col=connectcol[3])
  if (deciles & any(connect=='deciles')) for (j in 1:9) segments(decile[s,j], at[s], decile[s+1,j],
          at[s+1], lwd=.6*lwd,col=connectcol[4])
        }
    }
    if(stats) {
     if (all(quantiles==c(0,0))) quantiles=c(.25,.75)
     stats(x,by,quantiles)
     }
}
