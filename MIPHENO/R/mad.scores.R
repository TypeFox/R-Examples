##################################
##Name: mad.scores
##Description: calculates the mad score (zscore) based on the C2010 method published in Lu 2008
##O/S: for R
##Date: 2/17/2010
##Author: Shannon M. Bell
##Company: Michigan State University
##notes:
#made to reflect the C2010 method which calculates the a 'zscore'
#based on MAD of individual w/in flat
#note that there data for pipeline includes different ecotypes
#NOTE:data must be sorted by flatcode prior, as well as by identifier to maintain order
#n is the amount of difference before considered mutant, C2010 uses 3
#out parameter of 'Zscore' will give the MAD score of the individual, 'label' uses n and clasifies as high, low, or wt
#note: 11/11 updated function to accept any input column lable into parameter, verified consistency with prior version
mad.scores<-function(data, parameter='FLATCODE', n=3, out=c('Zscore', 'label'), ...){
	#this is just borrowing from another function, with the idea this could be expanded to include other groupings
	#calculates the median and mad for each flat
	if(! is.null(parameter)){
		form<-as.formula(paste('.~', parameter))
		initial.median<-summaryBy(form, data=data, FUN=median, na.rm=TRUE, keep.names=TRUE)
                initial.mad<-summaryBy(form, data=data, FUN=mad, na.rm=TRUE, keep.names=TRUE)
	}
   	r<-ncol(initial.median)
        flats<-unique(data[[parameter]])
        temp.cat<-NULL
        temp.cat<-as.data.frame(temp.cat)
	#this calculates the 'zscore' as per 2010 method
	#this is done by subtracting each value in a flat by the median of that flat
	#that difference is divided by the MAD for the flat
        for (i in 1:length(flats)){
		temp.d<-subset(data, data[[parameter]] == flats[i])
		temp.dif<-temp.d[,2:r]-initial.median[rep(i, nrow(temp.d)),2:r]
		temp.mad<-temp.dif/ initial.mad[rep(i, nrow(temp.dif)), 2:r]
		temp.cat<-rbind(temp.cat, temp.mad)
        }
        temp.lab<-temp.cat
        #this part gives the factor call if wanted
	temp.lab[][temp.lab[] > -n] <-'wt'
        temp.lab[][temp.cat[] <= -n] <-'low'
        temp.lab[][temp.cat[] >= n] <-'high'
        temp.lab2<-apply(temp.lab, 2, as.factor)
        cat2<-cbind(data[[parameter]], temp.cat)
        colnames(cat2)<-colnames(data)
        lab2<-cbind(data[[parameter]], temp.lab2)
        colnames(lab2)<-colnames(data)
        if(out =='Zscore'){
		cat2
        }
        else{
		lab2
        }       
}


#####################
