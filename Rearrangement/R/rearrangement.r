
rearrangement<- function(x, y, n=1000, stochastic=FALSE, avg=TRUE, order=1:length(x)) {
	
	if (class(x) != "list" && class(x) != "data.frame"){
		stop("x must be a list or data frame")}
	if (length(x) > 3){	stop("too many covariates")}
	d <- dim(as.array(y))
	if (length(x) != length(d)){
		stop("number of regressors is not compatable with dimensions of y")}
	if (length(order[order > length(x)]) > 0){
		stop("order must include only integers <= length(x)")} 
			
	unirearrangement <- function(x, y, n, stochastic)	{
		nx    <- (x - min(x)) / (max(x) - min(x)) 
		if (stochastic){
			uus <- sort(runif(n, 0, 1))		} 
		else { 	uus <- c(0:n) / n 	}
		fus   <- approxfun(nx, y)
		rfus  <- fus(uus) 		
		return(quantile(rfus, nx)) 	
	}
	
	index <- function(arr, dimen, num){
		dum <- NULL
		d <- dim(arr) 	
		count <- 0
		for (i in 1:length(d)){ 
			if (i == dimen){ 
				dum <- c(dum, rep(num, prod(d[-dimen]))) 
			} else if (count == 0){ 
				dum <- c(dum, rep(1:d[-dimen][1],prod(d[-dimen][-1]))) 
				count <- count + 1
			} else { 
				dum <- c(dum, rep(1:d[-dimen][-1][1],rep(d[-dimen][1],prod(d[-dimen][-1]))))}
		}
		return (array(dum, dim=c(prod(d[-dimen]),length(d))))
	}	
	
	if (length(x)==1){
		unirearrangement(x[[1]], y, n, stochastic)	} 
	else {
	
		if (avg){
			if (length(x)==2) order <- data.frame(1:2,2:1)			
			if (length(x)==3) order <- data.frame(1:3,3:1,c(1,3,2),c(2,3,1),c(3,1,2),c(2,1,3)) 
		} else{ order<-data.frame(order) 		}
		
		rr <- NULL 
		name<-c('a','b','c','d','e','f','g')
		name<-name[1:length(order)]
		count=1
		for (h in order){ 
			r <- y 			
			for (i in 1:length(h)){ 
				t <- h[-i] 				
				if (length(t)==1){
					if (t==2) t<-1 
				} else {					
					if (sum(t==c(3,2))==2 || sum(t==c(2,3))==2){ #t= 1,3 or t=3,1 or 1,2 2,1
						t <- t - 1 } 
					else{ 
						t<-ifelse(t==3,2,1)	
					}
				}					
				for (j in 1:d[h[i]]){ 				
					k <- index(y,h[i],j) 				
					p = array(r[k],dim=d[-h[i]])	
					r[k] <- rearrangement(x[-h[i]], p, n, stochastic, avg=FALSE, order=t) 
				}
			}
			rr[[name[count]]] <- r 
			count=count+1
		}
		rsum<-y
		rsum[]<-0
		rsum	<- Reduce("+",rr); 
		r<-rsum/length(order) 
		return (r)}
}

simconboot<- function(x, y, estimator, formula, B = 200, alpha = .05, sampsize = length(x), seed = 8,colInt=c(5:39)/2,...){
	set.seed(seed)
	sorteddata <- unique(sort(x))
	
	l<-match.call() 
	data = data.frame(x, y)
	names(data)<-c(l[[2]],l[[3]])
	m.taus = match("tau",names(l))
	if ( !is.na(m.taus)) {
		taus = eval(l[[m.taus]])
	}else{taus=0}	
	
	if(length(taus)>1){ 
		if(formula==0){ 	
			m.xx = match("xx",names(l))
			coldata0  = eval(l[[m.xx]])
			
			subsamplek <- function(data = data, num = length(x), sampsize = sampsize, B=B){
				temp.bootstrap<-list()
				bootstrap<-list()      
				for (sam in 1:B){
					sdata<-data[sample(num,sampsize,replace=TRUE),]
					funcsorted <- estimator(sdata[[as.character(l[[2]])]],sdata[[as.character(l[[3]])]],...)$fitted.values
					temp.bootstrap<-c(temp.bootstrap,list(funcsorted))
				}
				bootstrap <- temp.bootstrap
				return(bootstrap)
			}
			cef.0          <- estimator(x,y,...)$fitted.values		
			bootstrap      <- subsamplek(data = data, num = length(x), sampsize = sampsize, B=B)	
		}else{
			lc = attr(terms(formula),"variables");
			numReg = length(lc);
			reg = matrix(0,nrow=length(x),ncol=(numReg-2));	
			for(i in 3:numReg) { 
				j = i-2
				reg[,j]<- array(eval(lc[[i]]),c(length(x),1))
			}
			reg = cbind(1,reg); 
			
			m.colInt = match("colInt",names(l))
			coldata0  = eval(l[[m.colInt]]);
			
			iscoldata0 <- findInterval(coldata0,sorteddata);
						
			subsamplek <- function(data = data, num = length(x), sampsize = sampsize, B=B){
				temp.bootstrap<-list()
				bootstrap<-list()     
				for (sam in 1:B){					
					sdata	<- data[sample(num,sampsize,replace=TRUE),];		
					coefs	<- coef(estimator(formula,sdata,...));						
					cefs	<- t(reg%*%coefs); 
					cefs.0	<- cefs[,iscoldata0]
					temp.bootstrap<-c(temp.bootstrap,list(cefs.0))
				}
				bootstrap <- temp.bootstrap
				return(bootstrap)
			}
	
			coefs			<- coef(estimator(formula,data,...));	
			cef				<- t(reg%*%coefs);
			
			cef.0			<- cef[,iscoldata0];				
			bootstrap       <- subsamplek(data = data, num = length(x), sampsize = sampsize, B=B)	
		}
			
		mboot			<- Reduce("+",bootstrap)/B;
		m2boot			<- Reduce("+",lapply(bootstrap,function(x) x^2))/B;
		vboot			<- (B/(B-1))*(m2boot-mboot^2);

		kboot			<- lapply(bootstrap,function(x) sqrt((x-cef.0)^2/vboot));
		kmax.boot		<- sapply(kboot,max,na.rm=TRUE);
		kalpha.boot		<- quantile(kmax.boot,1-alpha/2);
		boot.se			<- sqrt(vboot);	
		u.boot  		<- cef.0 + kalpha.boot*boot.se 
		l.boot  		<- cef.0 - kalpha.boot*boot.se 
		
		conint <- list(x = x, y = y, sortedx = sorteddata, Lower = l.boot, Upper = u.boot,cef=cef.0,rowdata=taus,coldata = coldata0)	
	}
	else{ 	
		if(formula==0){
			subsamplek <- function(data = data, num = length(x), sampsize = sampsize, B=B){
				temp.bootstrap<-NULL
				bootstrap<-NULL      
				for (sam in 1:B){
					sdata<-data[sample(num,sampsize,replace=TRUE),]
					funcsorted <- estimator(sdata[[as.character(l[[2]])]],sdata[[as.character(l[[3]])]],...)$fitted.values
					temp.bootstrap<-cbind(temp.bootstrap,funcsorted)
				}
				bootstrap <- temp.bootstrap
				return(bootstrap)
			}
			cef            <- estimator(x,y,...)$fitted.values		
			bootstrap      <- subsamplek(data = data, num = length(x), sampsize = sampsize, B=B)	
		}
		else{
			subsamplek <- function(data = data, num = length(x), sampsize = sampsize, B=B){ 
				temp.bootstrap<-NULL
				bootstrap<-NULL      
				for (sam in 1:B){
					sdata			<- data[sample(num,sampsize,replace=TRUE),] 					
					func			<- approxfun(sdata[[as.character(l[[2]])]],estimator(formula,sdata,...)$fitted.values)
					temp.bootstrap	<- cbind(temp.bootstrap,func(sorteddata))
				}
				bootstrap <- temp.bootstrap
				return(bootstrap)
			}
			func            <- approxfun(x, estimator(formula,...)$fitted.values)
			cef             <- func(sorteddata)		
			bootstrap       <- subsamplek(data = data, num = length(x), sampsize = sampsize, B=B)	
		}
		
		V.boot  		<- apply(bootstrap, 1, var, na.rm = TRUE)
		nsdata			<- length(sorteddata)
		K.boot  		<- sqrt((bootstrap - matrix(cef, nsdata, B, byrow = FALSE))^2 / matrix(V.boot, nsdata, B, byrow = FALSE))
		Kmax.boot       <- apply(K.boot, 2, FUN = max, na.rm = TRUE)
		Kalpha.boot     <- quantile(Kmax.boot, 1 - alpha/2) 
		boot.se         <- sqrt(V.boot) 
		u.boot  		<- cef + Kalpha.boot*boot.se 
		l.boot  		<- cef - Kalpha.boot*boot.se 

		conint 			<- list(x = x, y = y, sortedx = sorteddata, Lower = l.boot, Upper = u.boot,cef=cef)
	}
	class(conint) <- c("conint", class(conint)) 
	return(conint)
}

rconint<- function(x, n=100, stochastic=FALSE,avg=TRUE){
	if (class(x)[1] != "conint") stop("x must be of class conint")
	if (length(x$rowdata)>1){		
		aa  = list(x$rowdata,x$coldata);
		rintl <- rearrangement(aa, x$Lower, n, stochastic,avg=TRUE)
		rintu <- rearrangement(aa, x$Upper, n, stochastic,avg=TRUE)
	}else{		
		rintl <- rearrangement(data.frame(x$sortedx), x$Lower, n, stochastic,avg=TRUE)
		rintu <- rearrangement(data.frame(x$sortedx), x$Upper, n, stochastic,avg=TRUE)
	}
	x$Lower <- rintl
	x$Upper <- rintu
	return(x)
}






lclm <- function(x, y, h, xx) 	{
	fv <- xx;
	dv <- xx;
	for (i in 1:length(xx)) {;
		z 	<- x - xx[i];
		wx	<- 1 * (abs(z) <= h);
		r	<- lm(y ~ 1, weights = wx);
		fv[i]	<- r$coef[1];
		};
	return(list(xx = xx, fitted.values = fv));
	};
	
lplm <- function(x, y, h, xx) 	{
	fv <- xx;
	dv <- xx;
	for (i in 1:length(xx)) {;
		z 	<- x - xx[i];
		wx	<- 1 * (abs(z) <= h);
		r	<- lm(y ~ z, weights = wx);
		fv[i]	<- r$coef[1];
		dv[i]	<- r$coef[2];
		};
	return(list(xx = xx, fitted.values = fv, dv = dv));
	};

lcrq2 <- function(x, y, h, xx, tau)	{
	if(length(tau)>1){
		fv <- array(0,dim=c(length(tau),length(xx)));
		for (i in 1:length(xx)) {;
			z 	<- x - xx[i];
			wx	<- 1 * (abs(z) <= h);
			r	<- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE);
			fv[,i]	<- array(r$coef[1,],dim=c(length(tau),1));
			};
	}else{
		fv <- xx;		
		for (i in 1:length(xx)) {;
			z 	<- x - xx[i];
			wx	<- 1 * (abs(z) <= h);
			r	<- rq(y ~ 1, weights = wx, tau = tau, ci = FALSE);
			fv[i]	<- r$coef[1];
			};
	}	
	return(list(xx = xx, fitted.values = fv));
	};	
	
lprq2 <- function(x, y, h, xx, tau) 	{
	if (length(tau)>1){
		fv <- array(0,dim=c(length(tau),length(xx)));
		dv <- array(0,dim=c(length(tau),length(xx)));
		for (i in 1:length(xx)) {;
			z 	<- x - xx[i];
			wx	<- 1 * (abs(z) <= h);
			r	<- rq(y ~ z, weights = wx, tau = tau, ci = FALSE);
			fv[,i]	<- array(r$coef[1,],dim=c(length(tau),1));
			dv[,i]	<- array(r$coef[2,],dim=c(length(tau),1));			
			};
	}else{
		fv <- xx
		dv <- xx
		for (i in 1:length(xx)) {;
			z 	<- x - xx[i];
			wx	<- 1 * (abs(z) <= h);
			r	<- rq(y ~ z, weights = wx, tau = tau, ci = FALSE);
			fv[i]	<- r$coef[1];
			dv[i]	<- r$coef[2];
			};	
	}
	return(list(xx = xx, fitted.values = fv, dv = dv));
	};
	
	
	
	
polygon.conint <-function(x, ...){
	if (class(x)[1] != "conint") stop("x must be of class conint")	
	polygon(c(x$sortedx,rev(x$sortedx)), c(x$Upper,rev(x$Lower)), ...)
}

points.conint <-function(x, ...){
	if (class(x)[1] != "conint") stop("x must be of class conint")
	points(x$sortedx, x$Lower, ...)
	points(x$sortedx, x$Upper, ...)
}

plot.conint <- function(x, border, col, ...){
	if (class(x)[1] != "conint") stop("x must be of class conint")
	if (missing(border)) border = NULL
	if (missing(col)) col = NA
	plot(x$x, x$y, type='n', ...)
	polygon(c(x$sortedx,rev(x$sortedx)),c(x$Upper,rev(x$Lower)),border=border,col=col)
}

lines.conint <- function(x, ...){
	if (class(x)[1] != "conint") stop("x must be of class conint")
	lines(x$sortedx, x$Lower, ...)
	lines(x$sortedx, x$Upper, ...)
}

