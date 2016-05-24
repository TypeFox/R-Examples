# Student - Newman - Keuls Test
# 
# Author: FrankKonietschke
###############################################################################
snk <- function(formula,data,alpha, MSE=NULL, df=NULL, silent = FALSE){
	
	dat <- model.frame(formula, data)
	if (ncol(dat) != 2) {
		stop("Specify one response and only one class variable in the formula")
	}
	if (is.numeric(dat[, 1]) == FALSE) {
		stop("Response variable must be numeric")
	}
	response <- dat[, 1]
	group <- as.factor(dat[, 2])
		fl <- levels(group)
		a <-nlevels(group)
		N <- length(response)
		samples <- split(response,group)
		n <- sapply(samples,"length")
		mm <- sapply(samples,"mean")
		vv <- sapply(samples,"var")
		if (is.null(MSE)){
		MSE <- sum((n-1)*vv)/(N-a)
	     }
		if (is.null(df)){
		df <- N-a
			}
	
		nc <- a*(a-1)/2
		order.h1 <- data.frame(Sample=fl, Size=n, Means=mm,Variance=vv)
		
		ordered <- order.h1[order(order.h1$Means,decreasing=FALSE), ]
		rownames(ordered) <- 1:a
		
		
		
  #---------------- Compute helping indices ----------#
	
	i <- 1:(a-1)
	
	h1 <- list()
	for(s in 1:(a-1)){
		
		h1[[s]]<- i[1:s]
		
	}
	vi <- unlist(h1)
	
	j <- a:2
	h2 <-list()
	
	for (s in 1:(a-1)){
		h2[[s]] <- j[s:1]
	}
	
	vj <- unlist(h2)
	
	h3 <- list()
	h4 <- list()
	
	for (s in 1:(a-1)){
		h3[[s]] <- rep(j[s],s)
		h4[[s]] <- rep(i[s],s)
	}
	Nmean <- unlist(h3)
	Step  <- unlist(h4)
	
	
	#--------Compute the Mean Differences---------#
	
	mean.difference <- sapply(1:nc,function(arg){
				i <- vi[arg]
				j <- vj[arg]
				(ordered$Means[j]-ordered$Means[i])
				
			})
	mean.difference <- round(mean.difference, 4)
	# ------- Compute the test statistics --------#
	
	T <- sapply(1:nc,function(arg){
				i<-vi[arg]
				j<-vj[arg]
				(ordered$Means[j]-ordered$Means[i])/sqrt(MSE/2*(1/ordered$Size[i] + 1/ordered$Size[j]))
				
			})


	T <- round(T, 4)
	
	# ----- Compute now the critical value -----#
	
	quantiles <- qtukey(1-alpha,Nmean,df)
	pvalues <- ptukey(T,Nmean,df,lower.tail=FALSE)
	
	#---- Calculate the rejected Hypotheses ------# 
	
	Rejected1 <- (pvalues<alpha)
	
	#------------ Names for the Output -----------# 
	
	names.ordered <- sapply(1:nc, function(arg){
				i <- vi[arg]
				j <- vj[arg]
				paste(ordered$Sample[j], "-", ordered$Sample[i], sep="")
			})
	
	# ------ Compute now the rejected statistics-----#
	
	for (s in 1:nc){
		if (Rejected1[s]==FALSE){
			Under1 <- (vj[s]>=vj)
			Under2 <- (vi[s]<=vi)
			Under3 <- Under1 * Under2
			Under4 <- which(Under3==1)
		
			Rejected1[Under4] <- FALSE
			
		}
	}
	
	
	#-----Prepare the pValues for the Output----#
	
	Out1 <- (pvalues < alpha)
	Out2 <- (Rejected1 == FALSE)
	Out3 <- Out1 * Out2
	Out4 <- (Out3 == 1)
	pvalues <- round(pvalues,4)
	pvalues[Out4] <- paste(">",alpha)
	
	variances.output <- data.frame(Overall=MSE, df=df)
	Comparison <- data.frame(Comparison=names.ordered,Diff=mean.difference, Statistic=T, Adj.P=pvalues, Rejected=Rejected1, Layer = Step)
	

      if (! silent)
	{
		cat("#----Student-Newman-Keuls (1927; 1939; 1952) rejective Multiple Test Procedure  \n\n")
		cat("#----Attention: The SNK test controls the FWER only in the WEAK sense  \n\n")	
	}
	
		 	
	result <- list(Ordered.Means = ordered, Variances=variances.output,
			SNK = Comparison)
	
	return(result)	
	
	
			}
			
			
snk.wrapper <- function(model, data,  alpha, silent=FALSE) {
				
				result <- snk(formula=formula(model), 
						data, 
						alpha = alpha)
						difference <- result$SNK$Diff#, result$SNK$Statistic, result$SNK$Layer)
						diffm<-cbind(difference,rep(NA,length(difference)),rep(NA,length(difference)))
						diffm<-matrix(diffm,nrow=length(difference))
				rownames(diffm)<-result$SNK$Comparison
				return(list(adjPValues=result$SNK$Adj.P,rejected=result$SNK$Rejected,statistics=result$SNK$Statistic,
								confIntervals= diffm,errorControl = new(Class='ErrorControl',type="FWER",alpha=alpha)))
			}
			
			
mutoss.snk <- function() { return(new(Class="MutossMethod",
								label="Student-Newman-Keuls Test",
								errorControl="FWER.weak",
								callFunction="snk.wrapper",
								output=c("adjPValues", "rejected", "confIntervals", "errorControl"),
								info="<h2>Student - Newman - Keuls rejective test procedure.\
                                       The procedure controls the FWER in the WEAK sense. </h2>\n\n\
										<p> The Newman-Keuls procedure is based on a stepwise or \
                                           layer approach to significance testing. Sample means are \
                                           ordered from the smallest to the largest. The largest \
                                           difference, which involves means that are r = p steps apart, \
                                           is tested first at alpha level of significance; if significant, \
                                           means that are r = p - 1 steps apart are tested at \alpha level \
                                           of significance and so on.\ 
										</p>\n\
										<h3>Reference:</h3>\
										<ul>\
										<li>Keuls M (1952).  \"<i> The use of the studentized range in \
                                            connection with an analysis of variance
										 </i>\" Euphytica 1 37, 112-122. </li>\n\
										</ul>",
								parameters=list(
										data=list(type="data.frame"),
										model=list(type="ANY"),
										alpha=list(type="numeric")
								)
						)) }
			
			
	
	
		
			
			

			
