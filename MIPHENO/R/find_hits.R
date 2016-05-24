###########################################################################
##Name: find_hits
##Description: function that returns list of putative hits, using Zvalues or MIPHENO emperical pvalues
##O/S: for R x64 v 2.11.0
##Date: 09/27/2010
##Author: Shannon M. Bell
##Company: Michigan State University
#This function takes a data set either with values from MIPHENO (0-1 range) Zscore method (zero centered)
#and returns a list where >2 or 50% of the observations have an observation more extreem than a threshold
#data=a dataframe containing the data and a column named 'LOCUS' where values to sort on exist
#source=a list of loci to iterate over, if NULL will iterate over all loci
#values= the columns (start->stop) in the dataframe to process
#var.cuts=will variable cutoffs be used? If TRUE give vectors low.cut and high.cut for the low and high pvalue cutoffs
#this shoudl be the same length as the number of columns
#if var.cuts=FALSE, supply a Z, which is the absolute value Zscore that is considered a hit
#eg if Z=2 then values <=-2 and values >= 2 will be considered putative hits
#if var.cuts=FALSE, and Z=NULL you must supply a cutoff,
#which is the total from both ends that is considered a hit
#eg if cutoff=0.05 then values <=0.025 and values >= 0.975 will be considered putative hits
#note that so long is the data has zero (and values close to it) and 'normal',
#and hits are away from it, you can use the 'Z' parameter, even if you only have positive values
#and the ones you are interested in are say >=1
#updated 11/11 to include parameter 'ID' which is the Identification column....should correspond to 'source' input and defaults to LOCUS
###############################################################################
find_hits<-function(data=data, ID= 'LOCUS', source=NULL, values=list(start=11, stop=21), var.cuts=FALSE, low.cut=NULL, high.cut=NULL, cutoff=0.05, Z=NULL, ...){
    
    flats<-unique(data[[ID]])
    emp.d<-subset(data, data[[ID]] == flats[i])
    if(is.null(source)){
	source<-as.vector(unique(data[[ID]]))
    }
    ll<-length(values$start:values$stop)
    nc<-ncol(data)
    #if using variable cuts, need to get the cutoff points
    if(var.cuts==TRUE){
	l<-low.cut[1,]
	h<-high.cut[1,]
    }
    #cutoff needs to be length of phenotypes to check
    if(is.null(Z) == FALSE){
	l<-rep(-1*abs(Z), nc)
	h<-rep(abs(Z), nc)
    }
    else{
    	l<-rep(cutoff/2, nc)
	h<-rep((1-cutoff/2), nc)
    }
    hi.temp<-NULL
    hi.temp<-as.data.frame(hi.temp)
    low.temp<-NULL
    low.temp<-as.data.frame(low.temp)
    #First go through each locus
    for (i in 1:length(source)){
	temp.data<-subset(data, data[[ID]] == source[i])
	if(nrow(temp.data)>=2){
	    ph<-nrow(hi.temp)
	    pl<-nrow(low.temp)
	    #check each phenotype subseting data which is as or more extreem than cut
	    for (j in values$start:values$stop){
		#note the [j-10] is for the 10 colums with labels in the input set
		#this makes the h.cut aling with the input dataset
		#may need to be changed to make fit a given set
		htemp<-subset(temp.data, temp.data[j] >= h[j-2])
		ltemp<-subset(temp.data, temp.data[j] <= l[j-2])
		ar<-length(na.omit(temp.data[,j]))
		#if the number of samples more extreem is > 1/2 total observations
		#then data is added and the loop breaks...
		#do not itter through rest for that loci, data are added
		#note that data has to be same measure (high or low)
		if(nrow(htemp) > ar/2 & nrow(htemp) >= 2){
		    hi.temp<-rbind(hi.temp, temp.data)
		    break
		}
		if(nrow(ltemp) > ar/2 & nrow(ltemp) >= 2){
		    low.temp<-rbind(low.temp, temp.data)
		    break
		}
	    }
	    hit<-rbind(hi.temp, low.temp)
	}
    }
	#this just checks to make sure that there was data, dont need it to be 30
	#prob should change this to check rows instead....
    if(nrow(hit) >1){
        hit<-hit[do.call(order, subset(hit, hit[[ID]])),]
        unique(hit)
    }
    else{
        print('No hits found')
    }
}
###############################################################################
