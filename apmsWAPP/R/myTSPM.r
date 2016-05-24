##### Origin:
##-------------------------------------------------------------------
## 	Name: TSPM.R
##  	R code for the paper by Paul L. Auer and R.W. Doerge:
## 	"A Two-Stage Poisson Model for Testing RNA-Seq Data"
## 	Date: February 2011
## 	Contact: Paul Auer 		plivermo@fhcrc.org
##		 R.W. Doerge 		doerge@purdue.edu
##------------------------------------------------------------------------

####  !!!
####  TSPM function adapted here for the application to affinity-purification mass-spectrometry data (AP-MS) 
####  based on spectral counts



###### The TSPM function ##############################################
#######################################################################


TSPM <- function(counts, x1, x0, lib.size, alpha.wh){

## Input:
#counts:   a matrix of spectral counts (rows=proteins, columns=samples)
#x1: 		   a vector of bait/control assignement (under the alternative hypothesis)
#x0: 		   a vector of bait/control assignement (under the null hypothesis)
#lib.size: a vector of scaling factors for normalization
#alpha.wh: the significance threshold to use for deciding whether a protein is overdispersed.
#          Defaults to 0.05.


## Output:
#id:                name of the protein
#log.fold.change:		a vector containing the estimated log fold changes for each protein
#pvalues: 			a vector containing the raw p-values testing differential expression for each protein
#LRT:           a vector of Likelihood-Ratio statistics for each protein  
#dispersion:    a vector of yes/no indicating overdispersion for each protein
#padj:			a vector containing the p-values after adjusting for multiple testing using the method of Benjamini-Hochberg


######## The main loop that fits the GLMs to each protein ########################

### Initializing model parameters ####
n <- dim(counts)[1]
per.gene.disp <- NULL
LRT <- NULL
score.test <- NULL
LFC <- NULL

###### Fitting the GLMs for each protein #################
	for(i in 1:n){
		### Fit full and reduced models ###
		model.1 <- glm(as.numeric(counts[i,]) ~ x1, offset=log(lib.size), family=poisson)
		model.0 <- glm(as.numeric(counts[i,]) ~ x0, offset=log(lib.size), family=poisson)

    if (model.1$coef[2]<0) {model.1 <- model.0}     # constraint: parameter beta1>0 ! one-sided-test situation
    
		### Obtain diagonals of Hat matrix from the full model fit ###
		hats <- hatvalues(model.1)

		### Obtain Pearson overdispersion estimate ####
		per.gene.disp[i] <- sum(residuals(model.1, type="pearson")^2)/model.1$df.residual

		### Obtain Likelihood ratio statistic ####
		LRT[i] <- deviance(model.0)-deviance(model.1)

		### Obtain score test statistic ####
		score.test[i] <- 1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2

		### Obtain the estimated log fold change ###
		LFC[i] <- model.1$coef[2]
	}

## Initialize parameters for Working-Hotelling bands around the score TSs ###
qchi <- qchisq(df=1, (1:n-0.5)/n)
MSE <- 2
UL <- NULL

#### Obtain the upper boundary of the WH bands #######################################
xbar <- mean(qchi)
bottom <- sum((qchi-xbar)^2)
top <- (qchi-xbar)^2
s <- sqrt(MSE*(1/n) + (top/bottom))
W <- sqrt(2*qf(df1=1, df2=n-1, p=1-(alpha.wh/n)))
UL <- pmax(qchi + W*s,1)

###### Obtain the indices of the over-dispersed and not-over-dispersed proteins, respectively ##########

if(length(which(sort(score.test)-UL > 0))>0){     # overdispersed proteins existent
cutoff <- min(which(sort(score.test)-UL > 0))
temp <- cutoff-1 + seq(cutoff:length(score.test))
over.disp <- which(score.test %in% sort(score.test)[temp])
not.over.disp <- setdiff(1:length(score.test), over.disp)
}
else  {                                           # no overdispersed proteins
not.over.disp <- c(1:length(score.test))
over.disp <- NULL
}


###### Compute p-values ####################################
p.f <- pf(LRT[over.disp]/per.gene.disp[over.disp], df1=1, df2=model.1$df.residual, lower.tail=FALSE)
p.chi <- pchisq(LRT[not.over.disp], df=1, lower.tail=FALSE)
p <- NULL
p[over.disp] <- p.f
p[not.over.disp] <- p.chi

##### Adjust the p-values using the B-H method ####################
p.bh.f <- p.adjust(p.f, method="BH")
p.bh.chi <- p.adjust(p.chi, method="BH")
final.p.bh.tagwise <- NULL
final.p.bh.tagwise[over.disp] <- p.bh.f
final.p.bh.tagwise[not.over.disp] <- p.bh.chi

index.disp <- NULL
index.disp[over.disp] <- "yes"
index.disp[not.over.disp] <- "no"
### Output ###
data.frame(id=rownames(counts), log.fold.change=LFC, pvalues=p,
padj=final.p.bh.tagwise, LRT=LRT, dispersion=index.disp)
}
