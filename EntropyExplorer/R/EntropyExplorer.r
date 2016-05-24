#input is a matrix
#nrows == gene count
#ncols == sample count
#
#code sample to produce such a matrix and call the functions
#data file should contain both the header and rownames
#>chick=read.table(file="Chicken_Metabolite.csv",header=TRUE, sep=',',row.names=1) #sep='\t'
#>chickm <- as.matrix(chick)		#both rownames(chick) and rownames(chickm) work fine
# chickm<-as.matrix(read.table(file="Chicken_Metabolite.csv",header=TRUE, sep='\t',row.names=1))

EntropyExplorer <- function(expm1,expm2,dmetric, otype,ntop,nperm,shift=c(0,0),padjustmethod="fdr")
{
	if(missing(expm1))
	{
		stop("expm1 is a required input.")
	}
	if(missing(expm2))
	{
		stop("expm2 is a required input.")
	}
	
	if(missing(dmetric))
	{
		stop("dmetric is a required input.")
	}
	if(missing(otype))
	{
		stop("otype is a required input.")
	}
	if(tolower(dmetric)!="de" & tolower(dmetric)!="dcv" & tolower(dmetric)!="dse")
	{
		stop("supported dmetric: de, dcv, dse; put quotes around dmetric names")
	}
	if(tolower(otype)!="v" & tolower(otype)!="pr" & tolower(otype)!="pa" & tolower(otype)!="bv" & tolower(otype)!="br" & tolower(otype)!="ba" & tolower(otype)!="vu" & tolower(otype)!="pu" & tolower(otype)!="bu")
	{
		stop("supported otype: v, pr, pa, bv, br, ba, vu, pu, bu; put quotes around otype names")
	}
	if(padjustmethod!="holm" & padjustmethod!="hochberg" & padjustmethod!="hommel" & padjustmethod!="bonferroni" & padjustmethod!="BH" & padjustmethod!="BY" & padjustmethod!="fdr" & padjustmethod!="none")
	{
		stop("supported padjustmethod: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none; put quotes around padjustmethod names")
	}
	if((tolower(dmetric)=="de" |  tolower(dmetric)=="dcv") & !missing(nperm))
	{
		stop("nperm need not be specified for differential expression analysis and differential CV analysis.")
	}
	if(tolower(dmetric)=="dse" &  (tolower(otype)=="v" | tolower(otype)=="vu") & !missing(nperm))
	{
		stop("nperm need not be specified for differential nse analysis with otype of v or vu.")
	}
	if(!missing(padjustmethod) & (tolower(otype)=="v" | tolower(otype)=="vu"))
	{
		stop("padjustmethod need not be specified for otype of v or vu.")
	}
	
	min1 <- min(expm1,na.rm=TRUE)
	min2 <- min(expm2,na.rm=TRUE)
	if(min1 <=0 | min2 <=0)
	{
		if(missing(shift))
		{
			if(min1 <=0)
			{
				cat(paste("The minimum value of matrix 1 is", as.character(min1),"\n"),fill=TRUE)
			}
			if(min2 <=0)
			{
				cat(paste("The minimum value of matrix 2 is", as.character(min2),"\n"),fill=TRUE)
			}
			stop("matrix 1 or 2 has non-positive values. EntropyExplorer expects all data values to be positive. You may rerun this function using the shift parameter to shift the data upward, or adjust the data yourself.\n")
		}
		else#if(missing(shift))
		{
			if(length(shift)!=2)
			{
				stop("The shift parameter should have length 2.")
			}
			if(min1<=0)
			{
				if(tolower(shift[[1]])=="auto")
				{
					expm1<-expm1+(abs(min1)+0.001) #use 0.001 as a default shift
				}
				else if(is.na(suppressWarnings(as.numeric(shift[[1]]))))
				{
					stop("Your shift 1 is not \"auto\" or a valid number.")
				}
				else
				{
					s1<-as.numeric(shift[[1]])
					if(s1<=abs(min1))
					{
						stop(paste("The minimum value of matrix 1 is", as.character(min1), ". Your shift 1 is not large enough.\n"))
					}
					else
					{
						expm1<-expm1+s1
					}
				}
			}#if(min1<=0)
			else
			{
				cat(paste("All values of matrix 1 are positive. Your shift 1 is ignored.\n"),fill=TRUE)
			}
			if(min2<=0)
			{
				if(tolower(shift[[2]])=="auto")
				{
					expm2<-expm2+(abs(min2)+0.001) #use 0.001 as a default shift
				}
				else if(is.na(suppressWarnings(as.numeric(shift[[2]]))))
				{
					stop("Your shift 2 is not \"auto\" or a valid number.")
				}
				else
				{
					s2<-as.numeric(shift[[2]])
					if(s2<=abs(min2))
					{
						stop(paste("The minimum value of matrix 2 is", as.character(min2), ". Your shift 2 is not large enough.\n"))
					}
					else
					{
						expm2<-expm2+s2
					}
				}
			}#if(min2<=0)
			else
			{
				cat(paste("All values of matrix 2 are positive. Your shift 2 is ignored.\n"),fill=TRUE)
			}
		}
	}#if(min1 <=0 | min2 <=0)
	else
	{
		if(!missing(shift))
		{
			cat(paste("All values of the two matrices are positive. Your shifts are ignored.\n"),fill=TRUE)
		}
	}
	
	if(tolower(dmetric)=="de")#differential expression
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffexp1(expm1,expm2, 1,ntop))
		}
		else if(tolower(otype)=="vu")#differential value; no sort
		{
			return (diffexp1(expm1,expm2, 0,ntop))
		}
		else if(tolower(otype)=="pr")#differential p-value; sort by raw p-value
		{
			return (diffexp2(expm1,expm2, 1, ntop,padjustmethod))
		}
		else if(tolower(otype)=="pa")#differential p-value; sort by adjusted p-value
		{
			return (diffexp2(expm1,expm2, 2, ntop,padjustmethod))
		}
		else if(tolower(otype)=="pu")#differential p-value; no sort
		{
			return (diffexp2(expm1,expm2, 0, ntop,padjustmethod))
		}
		else if(tolower(otype)=="br")#sort by raw p-value
		{
			return (diffexp(expm1,expm2, ntop,2,padjustmethod))
		}
		else if(tolower(otype)=="ba")#sort by adjusted p-value
		{
			return (diffexp(expm1,expm2, ntop,3,padjustmethod))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffexp(expm1,expm2, ntop,1,padjustmethod))
		}
		else if(tolower(otype)=="bu")
		{
			return (diffexp(expm1,expm2, ntop,0,padjustmethod))
		}
	}
	else if(tolower(dmetric)=="dcv")#differential CV
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffcv(expm1,expm2, 1, ntop))
		}
		else if(tolower(otype)=="vu")#differential value; no sort
		{
			return (diffcv(expm1,expm2, 0, ntop))
		}
		else if(tolower(otype)=="pr")#differential p-value; sort by raw p-value
		{
			return (diffcvpFK(expm1,expm2, 1, ntop,padjustmethod))
		}
		else if(tolower(otype)=="pa")#differential p-value; sort by adjusted p-value
		{
			return (diffcvpFK(expm1,expm2, 2, ntop,padjustmethod))
		}
		else if(tolower(otype)=="pu")#differential p-value; no sort
		{
			return (diffcvpFK(expm1,expm2, 0, ntop,padjustmethod))
		}
		else if(tolower(otype)=="br")#sort by raw p-value
		{
			return (diffcvall(expm1,expm2, ntop,2,padjustmethod))
		}
		else if(tolower(otype)=="ba")#sort by adjusted p-value
		{
			return (diffcvall(expm1,expm2, ntop,3,padjustmethod))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffcvall(expm1,expm2, ntop,1,padjustmethod))
		}
		else if(tolower(otype)=="bu")
		{
			return (diffcvall(expm1,expm2, ntop,0,padjustmethod))
		}
	}
	else if(tolower(dmetric)=="dse")#differential entropy
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffnse(expm1,expm2, 1, ntop))
		}
		else if(tolower(otype)=="vu")#differential value; no sort
		{
			return (diffnse(expm1,expm2, 0, ntop))
		}
		else if(tolower(otype)=="pr")#differential p-value; sort by raw p-value
		{
			return (diffnsepperm(expm1,expm2, 1, ntop,nperm,padjustmethod))
		}
		else if(tolower(otype)=="pa")#differential p-value; sort by adjusted p-value
		{
			return (diffnsepperm(expm1,expm2, 2, ntop,nperm,padjustmethod))
		}
		else if(tolower(otype)=="pu")#differential p-value; no sort
		{
			return (diffnsepperm(expm1,expm2, 0, ntop,nperm,padjustmethod))
		}
		else if(tolower(otype)=="br")#sort by raw p-value
		{
			return (diffseall(expm1,expm2, ntop,nperm,2,padjustmethod))
		}
		else if(tolower(otype)=="ba")#sort by adjusted p-value
		{
			return (diffseall(expm1,expm2, ntop,nperm,3,padjustmethod))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffseall(expm1,expm2, ntop,nperm,1,padjustmethod))
		}
		else if(tolower(otype)=="bu")
		{
			return (diffseall(expm1,expm2, ntop,nperm,0,padjustmethod))
		}
	}
}#EntropyExplorer

adjustlabel <- function(padjustmethod)
{
	if(padjustmethod=="none")
	{
		return ("raw p-value");
	}
	else
	{
		return (paste(padjustmethod,"p-value"));
	}
}#adjustlabel

#expm1 and expm2 are of class matrix; diffexp1 is function name
#dosort==1 then sort output rows; otherwise just output the first ntop rows of the input data
diffexp1 <- function(expm1,expm2, dosort, ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(dosort==1)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	de1 <- matrix(0,nrow=genecount,ncol=4,dimnames=list(rownames(expm1)[1:genecount],c("differential expression", "avg(expm1)", "avg(expm2)", "avg(expm1)-avg(expm2)")))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			de1[i,1] <- -1
			de1[i,2] <- NA
			de1[i,3] <- NA
			de1[i,4] <- NA
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			de1[i,1] <- -1
			de1[i,2] <- NA
			de1[i,3] <- NA
			de1[i,4] <- NA
			nacount <- nacount+1
			next
		}
		de1[i,2] <- mean(x1)
		de1[i,3] <- mean(x2)
		de1[i,4] <- de1[i,2]-de1[i,3]
		de1[i,1] <- abs(de1[i,4])
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	if(dosort==1)
	{
		de1 <- de1[ order(-de1[,1]), ,drop=FALSE]#decreasing order
	}
	if(dosort==1 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (de1[c(0:ntop), c(2:4),drop=FALSE]);
}#diffexp1

#expm1 and expm2 are of class matrix; diffexp2 is function name
#returns a two column matrix of p-values (t-test) with row-number equal to ntop
diffexp2 <- function(expm1,expm2,dosort,ntop,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(dosort>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	de1 <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1)[1:genecount],c("p-value",adjustlabel(padjustmethod))))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			de1[i,1] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			de1[i,1] <- 2
			nacount <- nacount+1
			next
		}
		if(var(x1)==0 & var(x2)==0)#t.test will complain "data are essentially constant" if both case and control have 0 variance
		{
			de1[i,1] <- 2
			nacount <- nacount+1
			next
		}
		myt=t.test(x1,x2)
		de1[i,1] <- myt[['p.value']]
	}
	de1[,2] <- p.adjust(de1[,1],padjustmethod)
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2 or identical values in both expm1 and expm2.")
	}
	if(dosort==1)
	{
		de1 <- de1[order(de1[,1]), ,drop=FALSE]#increasing order of raw p-value
	}
	else if(dosort==2)
	{
		de1 <- de1[order(de1[,2]), ,drop=FALSE]#increasing order of adjusted p-value
	}
	if(dosort>0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (de1[c(0:ntop), ,drop=FALSE]);
}#diffexp2

#sorder==3 sort by adjusted p-value; ==2 sort by raw p-value; ==1 sort by value; ==0 no sort
diffexp <- function(expm1,expm2, ntop,sorder,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(sorder>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	detwo <- matrix(0,nrow=genecount,ncol=6,dimnames=list(rownames(expm1)[1:genecount],c("differential expression","avg(expm1)", "avg(expm2)", "avg(expm1)-avg(exmp2)","p-value",adjustlabel(padjustmethod))))#col1: abs diff value; col2: value1; col3: value2; col4: signed diff value; col5: p-value; col6: adjusted p-value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			detwo[i,1] <- -1
			detwo[i,2] <- NA
			detwo[i,3] <- NA
			detwo[i,4] <- NA
			detwo[i,5] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			detwo[i,1] <- -1
			detwo[i,2] <- NA
			detwo[i,3] <- NA
			detwo[i,4] <- NA
			detwo[i,5] <- 2
			nacount <- nacount+1
			next
		}
		if(var(x1)==0 & var(x2)==0)#t.test will complain "data are essentially constant" if both case and control have 0 variance
		{
			detwo[i,2] <- mean(x1)
			detwo[i,3] <- mean(x2)
			detwo[i,4] <- detwo[i,2] - detwo[i,3]
			detwo[i,1] <- abs(detwo[i,4])
			detwo[i,5] <- NA
			nacount <- nacount+1
			next
		}
		detwo[i,2] <- mean(x1)
		detwo[i,3] <- mean(x2)
		detwo[i,4] <- detwo[i,2] - detwo[i,3]
		detwo[i,1] <- abs(detwo[i,4])
		myt=t.test(x1,x2)
		detwo[i,5] <- myt[['p.value']]
	}
	detwo[,6] <- p.adjust(detwo[,5],padjustmethod)
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2 or identical values in both expm1 and expm2.")
	}
	if(sorder==2)#sort by raw p-value
	{
		detwo<-detwo[ order(detwo[,5]), , drop=FALSE];
	}
	else if(sorder==3)#sort by adjusted p-value
	{
		detwo<-detwo[ order(detwo[,6]), , drop=FALSE];
	}
	else if(sorder==1)#sort by value
	{
		detwo<-detwo[ order(-detwo[,1]), , drop=FALSE];
	}
	if(sorder>0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (detwo[c(0:ntop), c(2:6) ,drop=FALSE]);
}#diffexp

#expm1 and expm2 are of class matrix; diffcv is function name
diffcv <- function(expm1,expm2,dosort, ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(dosort==1)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	dcv <- matrix(0,nrow=genecount,ncol=4,dimnames=list(rownames(expm1)[1:genecount],c("differential CV", "CV(expm1)", "CV(expm2)", "CV(expm1)-CV(expm2)")))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dcv[i,1] <- -1
			dcv[i,2] <- NA
			dcv[i,3] <- NA
			dcv[i,4] <- NA
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dcv[i,1] <- -1
			dcv[i,2] <- NA
			dcv[i,3] <- NA
			dcv[i,4] <- NA
			nacount <- nacount+1
			next
		}
		dcv[i,2] <- sd(x1)/mean(x1)
		dcv[i,3] <- sd(x2)/mean(x2)
		dcv[i,4] <- dcv[i,2]-dcv[i,3]
		dcv[i,1] <- abs(dcv[i,4])
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	if(dosort==1)
	{
		dcv <- dcv[order(-dcv[,1]), ,drop=FALSE]
	}
	if(dosort==1 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dcv[c(0:ntop), c(2:4),drop=FALSE]);
}#diffcv

#expm1 and expm2 are of class matrix; diffcvpFK is function name
#returns a two-column matrix of pvalues based on Fligner-Killeen (median) test
diffcvpFK <- function(expm1,expm2,dosort,ntop,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(dosort>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	pcv <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1)[1:genecount],c("p-value",adjustlabel(padjustmethod))))
	nacount<-0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			pcv[i,1] <- 2
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			pcv[i,1] <- 2
			nacount<-nacount+1
			next
		}
		result <- fligner.test(list(log(x1),log(x2)))
		pcv[i,1] <- result$p.value
	}
	pcv[,2] <- p.adjust(pcv[,1],padjustmethod)
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	if(dosort==1)#sort by raw p-value
	{
		pcv <- pcv[order(pcv[,1]), ,drop=FALSE]
	}
	else if(dosort==2)#sort by adjusted p-value
	{
		pcv <- pcv[order(pcv[,2]), ,drop=FALSE]
	}
	if(dosort>0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (pcv[c(0:ntop), ,drop=FALSE]);
}#diffcvpFK

#sorder==3 sort by adjusted p-value; ==2 sort by raw p-value; ==1 sort by value; ==0 no sort
diffcvall <- function(expm1,expm2, ntop, sorder,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(sorder>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	dcvtwo <- matrix(0,nrow=genecount,ncol=6,dimnames=list(rownames(expm1)[1:genecount],c("differential CV","CV(expm1)", "CV(expm2)", "CV(expm1)-CV(expm2)", "p-value", adjustlabel(padjustmethod))))#col1: abs diff value; col2: value1; col3: value2; col4: signed diff value; col5: p-value; col6: adjusted p-value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dcvtwo[i,1] <- -1
			dcvtwo[i,2] <- NA
			dcvtwo[i,3] <- NA
			dcvtwo[i,4] <- NA
			dcvtwo[i,5] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dcvtwo[i,1] <- -1
			dcvtwo[i,2] <- NA
			dcvtwo[i,3] <- NA
			dcvtwo[i,4] <- NA
			dcvtwo[i,5] <- 2
			nacount <- nacount+1
			next
		}
		dcvtwo[i,2] <- sd(x1)/mean(x1)
		dcvtwo[i,3] <- sd(x2)/mean(x2)
		dcvtwo[i,4] <- dcvtwo[i,2]-dcvtwo[i,3]
		dcvtwo[i,1] <- abs(dcvtwo[i,4])
		result <- fligner.test(list(log(x1),log(x2)))
		dcvtwo[i,5] <- result$p.value
	}
	dcvtwo[,6] <- p.adjust(dcvtwo[,5],padjustmethod)
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	if(sorder==2)#sort by raw p-value
	{
		dcvtwo<-dcvtwo[ order(dcvtwo[,5]), , drop=FALSE];#ascending order
	}
	else if(sorder==3)#sort by adjusted p-value
	{
		dcvtwo<-dcvtwo[ order(dcvtwo[,6]), , drop=FALSE];#ascending order
	}
	else if(sorder==1)#sort by value
	{
		dcvtwo<-dcvtwo[ order(-dcvtwo[,1]), , drop=FALSE];#decending order
	}
	if(sorder > 0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dcvtwo[c(0:ntop), c(2:6),drop=FALSE]);
}#diffcvall

#expm1 and expm2 are of class matrix; diffnse is function name
diffnse <- function(expm1,expm2,dosort,ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	if(dosort==1)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	dnse <- matrix(0,nrow=genecount,ncol=4,dimnames=list(rownames(expm1)[1:genecount],c("differential entropy","SE(expm1)", "SE(expm2)", "SE(expm1)-SE(expm2)")))
	nacount<-0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dnse[i,1] <- -1
			dnse[i,2] <- NA
			dnse[i,3] <- NA
			dnse[i,4] <- NA
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dnse[i,1] <- -1
			dnse[i,2] <- NA
			dnse[i,3] <- NA
			dnse[i,4] <- NA
			nacount<-nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		dnse[i,2] <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		dnse[i,3] <- -t/log2(length(x2))
		dnse[i,4] <- dnse[i,2]-dnse[i,3]
		dnse[i,1] <- abs(dnse[i,4])
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	if(dosort==1)
	{
		dnse <- dnse[order(-dnse[,1]), ,drop=FALSE]
	}
	if(dosort==1 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dnse[c(0:ntop), c(2:4),drop=FALSE]);
}#diffnse

#expm1 and expm2 are of class matrix; diffnsepperm is function name
#returns a two-column matrix of pvalues with row number equal to ntop
#uses the permutation test method
diffnsepperm <- function(expm1,expm2,dosort,ntop,nperm,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	if(missing(nperm))
	{
		nperm<-1000
	}
	stopifnot(ntop>0)
	stopifnot(nperm>=100)
	stopifnot(ntop<=dimv1[1])
	if(dosort>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	dnse <- matrix(0,nrow=genecount,ncol=1) #to store signed differential value
	nacount<-0
	
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dnse[i,1] <- 0
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dnse[i,1] <- 0
			nacount<-nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		nse1 <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		nse2 <- -t/log2(length(x2))
		dnse[i,1] <- nse1-nse2
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}

	pe <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1)[1:genecount],c("p-value",adjustlabel(padjustmethod)))) #to store p-value and adjusted p-value
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		x1Len <- length(x1)
		if (x1Len<4)
		{
			pe[i,1] <- 2
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		x2Len <- length(x2)
		if (x2Len<4)
		{
			pe[i,1] <- 2
			next
		}
		if(dnse[i,1]==0)
		{
			pe[i,1] <- 1
			next
		}
		xc <- c(x1,x2)
		xcLen <- x1Len + x2Len
		ecount<-0
		Ncount <- nperm
		for (j in 1:Ncount)
		{
			permvec <- sample.int(xcLen)
			xctemp <- xc[permvec] #sample(1:xcLen) produces a permutation of 1,2,...,xcLen
			x1temp <- xctemp[1:x1Len]
			s <- sum(x1temp)
			t <- 0
			for (k in 1:x1Len)
			{
				t <- t+x1temp[k]/s*log2(x1temp[k]/s)
			}
			nse1 <- -t/log2(x1Len)
			x2temp <- xctemp[(1+x1Len):xcLen]
			s <- sum(x2temp)
			t <- 0
			for (k in 1:x2Len)
			{
				t <- t+x2temp[k]/s*log2(x2temp[k]/s)
			}
			nse2 <- -t/log2(x2Len)
			if(dnse[i,1]>0)
			{
				if(nse1-nse2 >= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
			else #dnse[i,1]<0
			{
				if(nse1-nse2 <= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
		} # for j
		pe[i,1] <- ecount/Ncount
	}#for i

	pe[,2] <- p.adjust(pe[,1],padjustmethod)
	if(dosort==1)#sort by raw p-value
	{
		pe <- pe[order(pe[,1]), ,drop=FALSE]
	}
	else if(dosort==2)#sort by adjusted p-value
	{
		pe <- pe[order(pe[,2]), ,drop=FALSE]
	}
	if(dosort>0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (pe[c(0:ntop), ,drop=FALSE]);
}#diffnsepperm; expand.grid? ReferenceClasses?

diffseall <- function(expm1,expm2, ntop,nperm,sorder,padjustmethod)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	if(missing(nperm))
	{
		nperm<-1000
	}
	stopifnot(ntop>0)
	stopifnot(nperm>=100)
	stopifnot(ntop<=dimv1[1])
	if(sorder>0)
	{
		genecount <- dimv1[1]
	}
	else
	{
		genecount <- ntop
	}
	dsetwo <- matrix(0,nrow=genecount,ncol=6,dimnames=list(rownames(expm1)[1:genecount],c("differential entropy","SE(expm1)", "SE(expm2)", "SE(expm1)-SE(expm2)","p-value",adjustlabel(padjustmethod))))#col1: abs diff value; col2: value1; col3: value2; col4: signed diff value; col5: p-value; col6: adjusted p-value
	#dnse <- matrix(0,nrow=genecount,ncol=1) #to store signed differential value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dsetwo[i,1] <- -1
			dsetwo[i,2] <- NA
			dsetwo[i,3] <- NA
			dsetwo[i,4] <- NA
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dsetwo[i,1] <- -1
			dsetwo[i,2] <- NA
			dsetwo[i,3] <- NA
			dsetwo[i,4] <- NA
			nacount <- nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		dsetwo[i,2] <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		dsetwo[i,3] <- -t/log2(length(x2))
		dsetwo[i,4] <- dsetwo[i,2]-dsetwo[i,3]
		dsetwo[i,1] <- abs(dsetwo[i,4])
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in expm1 or expm2.")
	}
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		x1Len <- length(x1)
		if (x1Len<4)
		{
			dsetwo[i,5] <- 2
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		x2Len <- length(x2)
		if (x2Len<4)
		{
			dsetwo[i,5] <- 2
			next
		}
		if(dsetwo[i,4]==0)
		{
			dsetwo[i,5] <- 1
			next
		}
		xc <- c(x1,x2)
		xcLen <- x1Len + x2Len
		ecount<-0
		Ncount <- nperm
		for (j in 1:Ncount)
		{
			permvec <- sample.int(xcLen)
			xctemp <- xc[permvec] #sample(1:xcLen) produces a permutation of 1,2,...,xcLen
			x1temp <- xctemp[1:x1Len]
			s <- sum(x1temp)
			t <- 0
			for (k in 1:x1Len)
			{
				t <- t+x1temp[k]/s*log2(x1temp[k]/s)
			}
			nse1 <- -t/log2(x1Len)
			x2temp <- xctemp[(1+x1Len):xcLen]
			s <- sum(x2temp)
			t <- 0
			for (k in 1:x2Len)
			{
				t <- t+x2temp[k]/s*log2(x2temp[k]/s)
			}
			nse2 <- -t/log2(x2Len)
			if(dsetwo[i,4]>0)
			{
				if(nse1-nse2 >= dsetwo[i,4])
				{
					ecount <- ecount +1
				}
			}
			else #dsetwo[i,4]<0
			{
				if(nse1-nse2 <= dsetwo[i,4])
				{
					ecount <- ecount +1
				}
			}
		} # for j
		dsetwo[i,5] <- ecount/Ncount
	}#for i

	dsetwo[,6] <- p.adjust(dsetwo[,5],padjustmethod)
	if(sorder==2)#sort by raw p-value
	{
		dsetwo<-dsetwo[ order(dsetwo[,5]), , drop=FALSE];
	}
	if(sorder==3)#sort by adjusted p-value
	{
		dsetwo<-dsetwo[ order(dsetwo[,6]), , drop=FALSE];
	}
	else if(sorder==1)#sort by value
	{
		dsetwo<-dsetwo[ order(-dsetwo[,1]), , drop=FALSE];
	}
	if(sorder > 0 & ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dsetwo[c(0:ntop), c(2:6),drop=FALSE]);
}#diffseall
