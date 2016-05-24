##Name: rm.outliers
##Description: removes outling sample groups based on their distance from the mean of sample groups
##O/S: for R
##Date: 2/17/2010
##Author: Shannon M. Bell
##Company: Michigan State University
##notes:
#this function removes the outlier groups for Quality Control purposes
#the removal is based on whatever attribute is assigned, input data is expected to have a specific format
#first column is the 'parameter' the subsequent are all data to be evaluated
#this will remove only the data where the median value of the group is > n*MAD of the all group medians
#very important that data is ordered based on parameter and that parameter is the first column of the dataframe
#updated 11/11 to permit use of any column label as parameter, validated to ensure same result as with prior version
rm.outliers<-function(data, parameter='FLATCODE', n=3, ...){
	#first calculate the median of the each value by the grouping parameter to get initial.median
	#times is the number of observations made for each grouping parameter obtaitained by cal length,
	#this is used later on to maintain dimension
	if(! is.null(parameter)){
		form<-as.formula(paste('.~', parameter))
		initial.median<-summaryBy(form, data=data, FUN=median, na.rm=TRUE, keep.names=TRUE)
		times<-summaryBy(form, data=data, FUN=length, keep.names=TRUE)
	}
	r<-ncol(initial.median)
	#this is the median of the parameter medians
	pre.median<-apply(initial.median[,2:r], 2, median, na.rm=T)
	#this is the MAD for the parameter medians (so like the stdev)
	pre.mad<-apply(initial.median[,2:r], 2, mad, na.rm=T)
	#this is the difference between the paramter median and the median of the parameter medians
	pdiff<-sweep(initial.median[,2:r],2,FUN='-', pre.median)
	#adding back on the label for parameter to the df of differences
	fpdiff<-cbind(initial.median[,1], pdiff)
	#This goes through and looks to see if the difference b/w the parameter median and median of parameter medians
	#is beyound the tollerance n
	#for example, less than n=3*MAD or ~ 3 stddev 
	#this is eq to saying that I am considering  groups outside ~99.6% confidence interval as faulty samples.
	#by setting value for indivuduals w/in cutoff to 1 it will retain the value in the multiplication step later on
	#if missing data or median is outside the cutoff (ie 3MAD), the value will be NA (droped) in multiplication
	listing2<-NULL
	for (i in 1:nrow(pdiff)){
		calls<-NULL
		for (j in 1:ncol(pdiff)){
		#3 stdevis is ~99.6%; 95%~1.98 stdev
		#note the correction for if the value is na
			if(abs(pdiff[i,j]) <= n*pre.mad[j] && is.na(pdiff[i,j]) != 1){
				calls<-cbind(calls, 1)
			}
			else{
				calls<-cbind(calls, NA)
			}
		}
		temp<-cbind(fpdiff[i,1], as.data.frame(calls))
		listing2<-rbind(listing2, temp)
	}
	colnames(listing2)<-colnames(fpdiff)
	#repeat the rows in median info to the same dimensions as the input data based on labels
	#the output is a df with the same dimensions but just 1 and na for values
	f.listing2<-listing2[rep(1:nrow(listing2),times[,2]),]
	#this makes the values for any sample group deemed outside the cutoffs equal to NA and wont be used for further analysis
	#this is by attribute so a sample group can be outside the bounds in just one feature (e.g. assay)
	ps.flats<-data[,2:r] * f.listing2[,2:r]
	#adding back on labels
	ps.flats<-cbind(data[,1], as.data.frame(ps.flats))
	colnames(ps.flats)<-colnames(data)
	ps.flats
}


###############
