##################################
##Name: cdf.pval
##Description: calculates the pvalue based on cdf of input data
##O/S: for R
##Date: 2/17/2010
##Author: Shannon M. Bell
##Company: Michigan State University
##notes:
#cdf.data refers to data used to make the CDF (ie a null distribution)
#sample.data is the data for which the pvalue is to be determined
#note that data should be presorted prior and cannot contain labels
#also note that the pvalue returned is the prob of that value or one more extreem (lower)
#so high pval ->1 are also rare, but 1-pval or 1-F will give you a better estimate
cdf.pval<-function(cdf.data, sample.data=NULL, ...){
	new.cdf.data<-NULL
	header<-NULL
	#First check that there are observations (numeric) 
	#this checks each column (ie phenotype) creating new dataframe
	#cdf.data is the data from the NULL ie to make the distribution
	for(i in 1:ncol(cdf.data)){
		if(is.na(median(cdf.data[,i], na.rm=T)) == 0){
			new.cdf.data<-cbind(new.cdf.data, cdf.data[,i])
			header<-c(header, colnames(cdf.data[i]))
		}
	}
	colnames(new.cdf.data)<-header
	#this creates the cdf of the null data
	sample.cdf<-apply(new.cdf.data,2,ecdf)
	#check to see of user supplied test data,
	#if not assumes data supplied is for both cdf and test
	if(is.null(sample.data)== 1){
		sample.data<-cdf.data
	}
	pden<-NULL
	#sample.cdf is attached so column names can be used to call function
	#NOTE that error handleing here is an issue, must detach if error here
	#not sure how to address issue
	attach(sample.cdf)
	for(i in 1:ncol(sample.data)){
		#for each phenotype if there is num data, the cdf function is run
		#sample.cdf is a collection of sample cdfs
		if(is.na(median(sample.data[,i], na.rm=T)) == 0){
		    temp<-do.call(paste(colnames(sample.data[i])), list(sample.data[,i]))
		    pden<-cbind(pden, temp)
		}
		#if the data is not suitable for calculate cdf then it is just pasted in
		else{
		    pden<-cbind(pden, sample.data[,i])
		}
	}
	detach(sample.cdf)
	colnames(pden)<-colnames(sample.data)
	as.data.frame(pden)
}
################
