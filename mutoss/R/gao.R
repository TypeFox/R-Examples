# 
# 
# Author: FrankKonietschke
###############################################################################



gao<-function(formula, data, alpha = 0.05, control = NULL , silent = FALSE){
	
	dat <- model.frame(formula, data)
	if (ncol(dat) != 2) {
		stop("Specify one response and only one class variable in the formula !")
	}
	if (is.numeric(dat[, 1]) == FALSE) {
		stop("Response variable must be numeric !")
	}
	response <- dat[, 1]
	group <- as.factor(dat[, 2])
	fl <- levels(group)
	a <- nlevels(group)
	N <- length(response)
	n <- aggregate(response,list(group),FUN="length")$x
	
	if (any(n <= 2)) {
		warn <- paste("The factor level", fl[n <= 2], "has got less than two observations!")
		stop(warn)
	}
	if (is.null(control)) {
		cont <- 1
	
	}

	if(! is.null(control)){
   if (!any(fl == control)) {
		stop("The dataset doesn't contain this control group!")
	}
	cont <- which(fl == control)
}
	
	
	C<-contrMat(1:a,"Dunnett",base=cont)
	
	
	
	
# ------------- Compute the pseudo-ranks------------------ #
	
	
	#browser()
	rx <- c()
	
	for (i in 1:N){
		
		help <- expand.grid(response[i],response)
		
		help1 <- (help[,1]>help[,2])+1/2*(help[,1]== help[,2])
		help2 <- data.frame(h1=help1,h2=group)
		samples2 <- split(help2$h1, help2$h2)
		
		pseudo <- sapply(1:a, function(arg) {
					
					mean(samples2[[arg]])
				})
		
		rx[i] <-N*mean(pseudo)
	}
	new.data <-data.frame(res=rx,group=group)
	
	
# ------------------ Point estimators ---------------------#
	
	pd <- 1/N*aggregate(new.data$res,list(group), FUN="mean")$x
	
	Cpd <- C%*%pd
	
# ------------ Compute the variance estimators ----------- #
	
	v1 <- 1/N^2*aggregate(new.data$res,list(group),FUN="var")$x
	lambda <- N/n
	v11 <-c(v1*lambda)
	v2 <- diag(v1*lambda)
	
	Cv <- C%*%v2%*%t(C)
	
	
# ------------------ Test Statistics ----------------------#
	
	T <-sqrt(N)*Cpd / sqrt(c(diag(Cv))) 
	

# ------------------ Degrees of freedom--------------------#
	
	
	ncont <-which((1:a)!= cont)
	
	numerator <- c(diag(Cv))^2
	denu1<-v1[cont]^2/(n[cont]^2*(n[cont]-1))
	denu2 <- v1[ncont]^2 /(n[ncont]^2*(n[ncont]-1))
	denu <- N^2*(denu1 + denu2)
	
	df <- numerator / denu
	
	
	
#-------------------------p-Values ------------------------#
	
	pv<- c()
	
	for (h in 1:(a-1)){
		pv[h]<- min(2*pt(T[h],df[h]),2-2*pt(T[h],df[h]))
	}
	
	adj.p <- p.adjust(pv,"hochberg")
	Rejected <- (adj.p<=alpha)
	
#------------------- Build the output ---------------------# 
	vj <- which((1:a) != cont)
	vi <- rep(cont, a - 1)
	cmpid <- sapply(1:(a-1), function(arg) {
				i <- vi[arg]
				j <- vj[arg]
				paste("F", "(", fl[j], ")", "-","F","(" ,fl[i],")", sep = "")
			})

	
	result <- data.frame(Comparison=cmpid, Estimator = Cpd, df=df, Statistic = T, P.Raw=pv,P.Adj=adj.p,Rejected = Rejected )
	rownames(result)<-1:(a-1)
	
	
	output = list(Info=data.frame(Sample=fl, Size=n, Single.Effects=pd),
			Analysis=result)
	if (! silent)
	{
		cat("#----Xin Gao's (2008) Non-Parametric Multiple Test Procedure","\n") 
		cat("#----Type of Adjustment: Hochberg", "\n")
		cat("#----Level of significance", "=", alpha ,"\n")
		cat("#----The procedure compares if the distribution functions F() are equal. The FWER is strongly controlled", "\n")
		print(result)
	}
	return(output)
	
}


gao.wrapper <- function(model, data,  alpha, control) {
	control <- NULL
	result <- gao(formula=formula(model), 
			data=data, 
			alpha = alpha,control)
	
	
	

	
	pvalues <- result$Analysis$P.Adj
	estimates <- result$Analysis$Estimator
	confint <- cbind(estimates, rep(NA, length(estimates)),rep(NA,length(estimates)))
	rownames(confint)<-result$Analysis$Comparison
	rejected1 <- result$Analysis$Rejected

	
	return(list(adjPValues=pvalues,rejected=rejected1,confIntervals= confint,
					errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))
}


mutoss.gao <- function() { return(new(Class="MutossMethod",
					label="Nonparametric Multiple contrast tests",
					errorControl="FWER",
					callFunction="gao.wrapper",
					output=c("adjPValues", "rejected","confIntervals","errorControl"),
					info="<h2>Nonparametric multiple contrast tests</h2>
							<p> This function computes Xin Gao's nonparametric multiple test procedures in an unbalanced one way layout. <p> 
							<p></p>
							<h3>Reference:</h3>
							<ul>
							<li>Gao, X. et al. \"<i>Nonparametric multiple comparison procedures for unbalanced one-way factorial designs.</i>\" Journal of Statistical Planning and Inference, 77, 2574-2591, 2008.</li>
							</ul>",
					parameters=list(model=list(type="ANY"),
							hypotheses=list(type="ANY"),
							alpha=list(type="numeric")
							
							
					
					)
			)) }


