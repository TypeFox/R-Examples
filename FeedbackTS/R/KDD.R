 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #          KDD            # # #
  # # # # # # # # # # # # # # # # # 

.KDD.validity <- function(object){
	m=nrow(object@before.after)
	n=ncol(object@before.after)
	if(m<3){
			stop("[KDD:validity] Error: There should be at least one day before and one day after the keydays.")
		return(FALSE)
	} else {}
	if(round(m/2)==m/2){
			stop("[KDD:validity] Error: There should be the same numbers of days before and after the keydays and, consequently, the number of rows in matrix before.after should be odd.")
		return(FALSE)
	} else {}
	if(length(object@date)!=n || length(object@year)!=n ||  length(object@day)!=n){
			stop("[KDD:validity] Error: length of vectors date, year and day and number of columns of matrix before.after do not match.")
		return(FALSE)
	} else {}
	if(length(object@keyday.threshold)!=1){
			stop("[KDD:validity] Error: length of keyday.threshold should be one.")
		return(FALSE)
	} else {}
	if(sum(object@keyday.threshold>object@before.after[ceiling(m/2),])){
			stop("[KDD:validity] Error: values for keydays (central row in before.after matrix) should be larger than or equal to the value of keyday.threshold.")
		return(FALSE)
	} else {}
	return(TRUE)
}


setClass(
    Class="KDD",
    slots=c(
        before.after="matrix",
        date="character",
        year="numeric",
        day="numeric",
		keyday.threshold="numeric"
   ),
    prototype=prototype(
        before.after=matrix(),
        date=character(),
        year=numeric(),
        day=numeric(),
		keyday.threshold=numeric()
   ),
    validity=.KDD.validity
)


cat("### Constructor 1 ###\n")
kdd.from.raw.data <- function(raw.data,keyday.threshold,nb.days,col.series,col.date,na.rm=TRUE,filter=NULL){
	rain=raw.data[,col.series]
	range.days=(-nb.days):nb.days
	keyday=(1:nrow(raw.data))[rain>=keyday.threshold]
	keyday=keyday[!is.na(keyday)]
	if(length(keyday)==0){ 
		before.after=NULL
		keyday=NULL
		date=NULL
		year=NULL
	} else {
		before.after=matrix(NA,nb.days*2+1,length(keyday))	
		for(i in 1:length(keyday)){
			range.days.i=range.days[keyday[i]+range.days>=1 & 
				keyday[i]+range.days<=length(rain)]
			before.after[1+nb.days+range.days.i,i]=rain[keyday[i]+
				range.days.i]
		}
		if(na.rm){
			no.NA=apply(before.after,2,function(r) sum(is.na(r))==0)
			if(nrow(before.after)==1){
				before.after=rbind(before.after[,no.NA])
			} else {
				before.after=cbind(before.after[,no.NA])
			}
			keyday=keyday[no.NA]
		}
		if(length(filter)>0){
			for(f in 1:length(filter)){
				if(length(keyday)>0){
					if(filter[[f]]$apply.over=="keyday"){
						keep=raw.data[keyday,filter[[f]]$column]==filter[[f]]$value
						before.after=cbind(before.after[,keep])
						keyday=keyday[keep]
					}			
					if(filter[[f]]$apply.over=="range"){
						keep=NULL
						for(i in 1:length(keyday)){
							range.days.i=range.days[keyday[i]+range.days>=1 & 
								keyday[i]+range.days<=length(rain)]
							filter.values=raw.data[keyday[i]+range.days.i,
								filter[[f]]$column]
							keep=c(keep,sum(filter.values!=filter[[f]]$value)==0)
						}
						before.after=cbind(before.after[,keep])
						keyday=keyday[keep]
					}
				}
			}
		}
		year=raw.data[keyday,col.date[1]]
		date=year
		if(length(col.date)>1){
			for(j in 2:length(col.date)){
				date=paste(date,raw.data[keyday,col.date[j]],sep=".")
			}
		}
	}
	day=keyday
	if(length(keyday)==0){ 
		print("[KDD:constructor from raw data] Warning: no keyday given the applied filters.")	
	}
    return(kdd(before.after,date,year,day,keyday.threshold))
}

cat("### Constructor 2 ###\n")
kdd <- function(before.after,date,year,day,keyday.threshold){
    return(new("KDD",before.after=before.after,date=date,year=year,day=day,
		keyday.threshold=keyday.threshold))
}

cat("### Show ###\n")
.KDD.show <- function(object){
    cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ before.after: \n")		
    print(object@before.after)
    cat("\n ~ date: ",object@date)
    cat("\n ~ year: ",object@year)
    cat("\n ~ day: ",object@day)
    cat("\n ~ keyday.threshold: ",object@keyday.threshold)
    cat("\n")
    return(invisible())
}
setMethod(f="show",signature="KDD",definition=.KDD.show)


cat("### Getteur ###\n")
.KDD.get <- function(x,i,j,drop){
    switch(EXPR=i,
       "before.after"={return(x@before.after)},
       "date"={return(x@date)},
       "year"={return(x@year)},
       "day"={return(x@day)},
       "keyday.threshold"={return(x@keyday.threshold)},
       stop("[KDD:get] ",i," is not a 'KDD' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="KDD",definition=.KDD.get)


cat("### Setteur ###\n"	)
.KDD.set <- function(x,i,j,value){
    switch(EXPR=i,
       "before.after"={x@before.after=value},
       "date"={x@date=value},
       "year"={x@year=value},
       "day"={x@day=value},
       "keyday.threshold"={
       		if(x@keyday.threshold>value){
       			cat("[KDD:set] Warning: the value of keyday.threshold has been modified but, by decreasing this value, you should maybe include other keydays in your KDD object.")
       			x@keyday.threshold=value
       		} else {
       			if(x@keyday.threshold<value){
					m=nrow(x@before.after)
					rows2delete=(value>x@before.after[ceiling(m/2),])
       				if(sum(rows2delete)>0){
		       			cat("[KDD:set] Warning: the value of keyday.threshold has been increased and your KDD object has been reduced accordingly (i.e. the numbers of keydays has been decreased).")
		       			x@before.after=cbind(x@before.after[,!rows2delete])
		       			x@date=x@date[!rows2delete]
		       			x@year=x@year[!rows2delete]
		       			x@day=x@day[!rows2delete]
   		    		} else {}
   		    		x@keyday.threshold=value
   		    	}
       		}	
       	},
       stop("[KDD:set] Error:",i," is not a 'KDD' slot")
    )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="KDD",definition=.KDD.set)

.KDD.summary=function(object){
	cat("Keyday data object:\n")
	total.nb.keydays=length(object@date)
	cat("Data period: ", object@year[1],"-",object@year[length(object@year)],"\n")
	cat("Keyday threshold: ", object@keyday.threshold,"\n")
	cat("Total number of keydays: ", total.nb.keydays,"\n")
	kdd.values=object@before.after[ceiling(nrow(object@before.after)/2),]
	cat("Mean and standard deviation of rainfalls at keydays: ", 
			round(mean(kdd.values),digits=2)," and ",round(sd(kdd.values),digits=2),"\n")
	return(invisible())
}
setMethod(f="summary",signature="KDD",definition=.KDD.summary)

.KDD.names=function(x){
	slotNames(x)
}
setMethod(f="names",signature="KDD",definition=.KDD.names)
