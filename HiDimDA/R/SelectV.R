### SelectV.R  (2012-06-27)
###    
###
### Copyright 2011 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

SelectV <- function(data,grouping,Selmethod=c("ExpHC","HC","Fdr","Fair","fixedp"),
			NullDist=c("locfdr","Theoretical"),uselocfdr=c("onlyHC","always"),
			minlocfdrp=200,comvar=TRUE,Fdralpha=0.5,ExpHCalpha=0.5,HCalpha0=0.1,
			maxp=ncol(data),tol=1E-12,...)
{
   Selmethod <- match.arg(Selmethod)
   NullDist <-match.arg(NullDist)
   uselocfdr <- match.arg(uselocfdr)
   if (NullDist != "locfdr" && uselocfdr == "always") 
	stop("Error: uselocfdr argument can only be used when NullDist is set to locfdr")

   p <- ncol(data)
   if (p < minlocfdrp) NullDist <- "Theoretical"
   nk <- table(grouping)
   k <- nrow(nk)
   nk <- as.vector(nk)
   n <- sum(nk)
   if (Selmethod=="Fair" && k!=2) 
	stop("Fair method can only be used with two-group classification problems.\n")
   if (k==2) {
	tscr <- tscores(data,grouping,n,nk,comvar=comvar)
	scores <- abs(tscr$st)
   	pvalues <- 2*pt(scores,tscr$df,lower.tail=FALSE)
   }
   else {
	fscr <- fscores(data,grouping,n,nk,k)
	scores <- fscr$st
   	pvalues <- pf(scores,k-1,fscr$df,lower.tail=FALSE)
   }
   pvalues[pvalues<tol] <- tol   	# ensure that no pvalue is numerically identical to zero.
   pvalues[pvalues>1-tol] <- 1 - tol   	# ensure that no pvalue is numerically identical to one.
   if (Selmethod=="fixedp")  {
      	sortedpv <- sort(pvalues,index.return=TRUE)
	return(list(nvkpt=maxp,vkptInd=sort(sortedpv$ix[1:maxp])))
   }
   if (NullDist=="locfdr" && uselocfdr == "always") pvalues <- locfdrpval(pvalues)
   if (Selmethod=="ExpHC" || Selmethod=="Fdr")  
   {
      sortedpv <- sort(pvalues,index.return=TRUE)
	if (Selmethod=="ExpHC") usefullpv <- sortedpv$x[sortedpv$x<ExpHCalpha*1:p/(p*sum(1/1:p))]
	else usefullpv <- sortedpv$x[sortedpv$x<Fdralpha*1:p/p]
	if (length(usefullpv)==0) Fdrnvar <- 1
	else {
		maxpv <- max(usefullpv)
		Fdrnvar <- min(maxp,which(sortedpv$x==maxpv))
	}
   }
   if (Selmethod=="Fair")  {
	Stddata <- Fairstdbygrps(data,grouping,nk,n,p) 
	if (n-k>p) Fairres <- Fair(scores^2,p,nk,R=t(Stddata)%*%Stddata)
	else Fairres <- Fair(scores^2,p,nk,StdDt=Stddata)
	names(Fairres$m) <- NULL
	return(list(nvkpt=Fairres$m,vkptInd=Fairres$vkptInd))
   }
   if (NullDist=="locfdr" && uselocfdr == "onlyHC") pvalues <- locfdrpval(pvalues)
   {
 	if (Selmethod == "ExpHC" || Selmethod == "HC") {
		if (Selmethod== "ExpHC") minvar <- Fdrnvar
		else minvar <- 1
 		HCres <- HC(p,pvalues,minvkpt=minvar,alpha0=min(HCalpha0,maxp/p))
		names(HCres$nkptvar) <- NULL
		return(list(nvkpt=HCres$nkptvar,vkptInd=HCres$varkept))
	}
	else {
		if (Selmethod == "Fdr") nkptvar = Fdrnvar
		names(nkptvar) <- NULL
		return(list(nvkpt=nkptvar,vkptInd=sortedpv$ix[1:nkptvar]))
  	}
   } 
}

tscores <- function(data,grouping,n,nk,comvar)
{
  # Computes two-group t-scores

  Xbark <- apply(data,2,grpmeans,grp=grouping)
  vark <- apply(data,2,grpvar,grp=grouping)
  if (comvar==TRUE) {
	df <- n-2
	denom <- sqrt( (1/nk[1]+1/nk[2]) * ((nk[1]-1)*vark[1,]+(nk[2]-1)*vark[2,]) / df )
  }
  else {
	tmp1 <- vark[1,]/nk[1]
	tmp2 <- vark[2,]/nk[2]
	tmps <- tmp1 + tmp2
	df <- round( tmps^2/ ( tmp1^2/(nk[1]-1)+tmp2^2/(nk[2]-1) ) )
  	denom <- sqrt(tmps)
   }
   list(st=(Xbark[1,]-Xbark[2,])/denom,df=df)  
}

fscores <- function(data,grouping,n,nk,k)
{
  # Computes ANOVA f-scores

  df <- n - k
  vark <- apply(data,2,grpvar,grp=grouping)
  W <- apply((nk-1)*vark,2,sum)
  B <- (n-1)*apply(data,2,var) - W
  list(st=(B/(k-1))/(W/df),df=df)   # return(list(st=(B/(k-1))/(W/df),df=df))
}

locfdrpval <- function(pvalues)
{
	zscores <- qnorm(pvalues)	#  Note:  This is diferent from defining zscores directly from tscores
	empnull <- mylocfdr(zscores,plot=0,silently=TRUE)
	if (class(empnull)=="error1") empnull <- mylocfdr(zscores,plot=0,nulltype=2,silently=TRUE)
	if (class(empnull)=="error3") empnull <- mylocfdr(zscores,plot=0,nulltype=1,silently=TRUE)
	if (class(empnull)!="error2")  { 
		zscores <- (zscores-empnull$fp0[3,1])/empnull$fp0[3,2]
		pvalues <- pnorm(zscores)
	}
	return(pvalues)
}

HC <- function(p,pvalues,HCplus=FALSE,minvkpt=1,alpha0=0.1)
{
#    Computes Donoho and Jin's Higher Criticism threshold

#    Arguments:

#        p       -  the original number of variables 
#        pvalues -  a set of p pvalues  
#        HCplus  -  a boolean flag indicating if the HCplus version (always keep variables with pvalues below 1/p)
#                   of the HC criterion should be used 
#        minvkpt -  a minimum number of variables to be kept in the analysis
#        alpha0  -  the maximum percentage of variables to be kept in the analysis (by default 10%)

# Value:  a list with the three components:

#        threshold  -  the value of the Higher Criticism threshold (measured on the pvalue scale)
#        varkept    -  a vector with the indices of the variables to be kept in the analysis 
#        nkptvar     - the number of variables to be kept in the analysis 

        sortedpv <- sort(pvalues,index.return=TRUE)
        if (HCplus) p0 <- max(minvkpt,length(sortedpv$x[sortedpv$data0<=1/p])+1)
        else p0 <- minvkpt
        p1 <- floor(alpha0*p)
	if (p0 >= p1) nkptvar <- p0
        else  {
		unifq <- (p0:p1)/p
		HC <- p * (unifq-sortedpv$x[p0:p1]) / sqrt( (p0:p1)*(1-unifq) )
		if (max(HC)>0.) nkptvar <- which.max(HC)+p0-1
		else nkptvar <- p0
	}
  	XPind <- sort(sortedpv$ix[1:nkptvar])

        list(threshold=sortedpv$x[nkptvar],varkept=XPind,nkptvar=nkptvar) 
}

Fair <- function(T2,p,nk,R=NULL,StdDt=NULL,ivar=FALSE,blocksize=25,maxblrun=7,maxp=p)
{
#	Returns the number of variables to be used in a linear discriminant rule according to Fan and Fan, FAIR criterion

#	Arguments:

#		T2        -- vector of squared t-tests (usign different sample variances in each group) for testing
#                            mean group differences across groups.
#		p         -- total number of original variables
#		nk        -- two-dim vector with the number of observations by group 
#		R         -- empirical correlation matrix (used only when the argument StdDt is set to NULL)
#		StdDt     -- matrix (observations in rows, variables in columns) with the data standardized by groups
#		ivar      -- boolean flag indicating if independent variables are to be assumed
#		blocksize -- number of observations to be analysed simultaneously 
#                            (technical parameter used to control the algorithm speed)
#		maxblrun  -- maximum number of consecutive decreasing maximum values (in each block) of comparison criterion, 
#                            before the algorithm is stoped 
#  		maxp      -- maximum number of predictors to be used in the discriminant rule.

#	Value: A list with three components:

#		m          -- the number of variables to be used in the linear discriminant rule
#		vkptInd    -- an m-dimensional vector with the indices of the variables kept in the discriminat rule
#      		threshold  -- the threshold for variable selection (measured in standardized mean differences).


	maxeigvl <- function(m,M) { 
		indices <- srtdT2$ix[1:m]
		return(eigen(M[indices,indices],symmetric=TRUE,only.values=TRUE)$values[1])
	}  
	maxsingvl <- function(m,M) rghtsngv(M[,srtdT2$ix[1:m],drop=FALSE],nv=0)$d[1]

	n <- nk[1]+nk[2]
    	n1n2 <- nk[1]*nk[2]
	n1minusn2 <- nk[1]-nk[2]
    	T2sum <- array(dim=blocksize)
    	srtdT2 <- sort(T2,decreasing=TRUE,index.return=TRUE)

	bestval <- pT2sum <- pcrtval <- 0.
    	run <- 0
    	for (a in 1:(floor(maxp/blocksize)+1))  {
		m0 <- (a-1)*blocksize
		if (m0==maxp) break
		blocksize <- min(blocksize,maxp-m0)
		indices <- (m0+1):(m0+blocksize)
		T2sum[1] <- pT2sum + srtdT2$x[m0+1]
    		if (blocksize>1) for (b in 2:blocksize) T2sum[b] <- T2sum[b-1] + srtdT2$x[m0+b]
		pT2sum <- T2sum[blocksize]
		if (ivar==FALSE)  {
			if (!is.null(StdDt)) lambdamax <- sapply(indices,maxsingvl,M=StdDt)^2/(nrow(StdDt)-2) 
			else lambdamax <- sapply(indices,maxeigvl,M=R) 
			crtval <- n*(T2sum + indices*n1minusn2/n)^2/(lambdamax*n1n2*(indices+T2sum))
		}
		else crtval <- n*(T2sum + indices*n1minusn2/n)^2/(n1n2*(indices+T2sum))
                m1 <- which.max(crtval)
		if (crtval[m1] > bestval) {
			m <- m0+m1
			bestval <- crtval[m1]
		}
        	else if (crtval[m1] > pcrtval) run <- 0
                     else run <- run+1
                if (run > maxblrun)  break         
		pcrtval <- crtval[m1]
    }
    if (m > maxp) m <- maxp

   return(list(m=m,vkptInd=srtdT2$ix[1:m],threshold=sqrt(srtdT2$x[m])))
}

Fairstdbygrps <- function(data,grouping,nk=NULL,n=NULL,p=NULL) 
{
#  Standardizes a data set of observations divided by two groups, so that all variables have unit variances.
#  Uses a simple (unweighted) mean of the within group sample variances to estimate the unknown common variance	

	if (is.null(nk)) nk <- as.vector(table(grouping))
	if (is.null(n))  n <- nk[1] + nk[2]
	if (is.null(p))  p <- ncol(data)

	grplevels <- levels(grouping)
	vark <-  apply(data,2,grpvar,grp=grouping)
	meank <-  apply(data,2,grpmeans,grp=grouping)
	globals <- sqrt((vark[1,]+vark[2,])/2)
	for (i in 1:2) data[grouping==grplevels[i],] <- scale(data[grouping==grplevels[i],],center=meank[i,],scale=globals)
	data
}



