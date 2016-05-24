 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #    KDD.yearly.average   # # #
  # # # # # # # # # # # # # # # # # 


.KDD.yearly.average.validity <- function(object){
	m=nrow(object@before.after)
	n=ncol(object@before.after)
	if(m<3){
		stop("[KDD.yearly.average:validity] Error: There should be at least one day before and one day after the keydays.")
	} else {}
	if(round(m/2)==m/2){
		stop("[KDD.yearly.average:validity] Error: There should be the same numbers of days before and after the keydays and, consequently, the number of rows in matrix before.after should be odd.")
	} else {}
	if(length(object@yearly.nb.keydays)!=n || length(object@year)!=n){
		stop("[KDD.yearly.average:validity] Error: lengths of vectors year and yearly.nb.keydays, and number of columns of matrix before.after do not match.")
	} else {}
	if(length(object@keyday.threshold)!=1){
		stop("[KDD.yearly.average:validity] Error: length of keyday.threshold should be one.")
	} else {}
	if(sum(object@keyday.threshold>object@before.after[ceiling(m/2),])){
		stop("[KDD.yearly.average:validity] Error: values for keydays (central row in before.after matrix) should be larger then or equal to the value of keyday.threshold.")
	} else {}
	return(TRUE)
}


## definition of a new class heriting from class KDD
setClass(Class="KDD.yearly.average",
	slots=c(before.after="matrix",
        	year="numeric",
			keyday.threshold="numeric",
			yearly.nb.keydays="numeric"),
	validity=.KDD.yearly.average.validity
) 

cat("### Constructor 3 ###\n")
kdd.yearly.average <- function(object){
	year=NULL
	before.after=NULL
	yearly.nb.keydays=NULL
	for(i in unique(object@year)){
		year=c(year,i)
		before.after=cbind(before.after,
			apply(cbind(object@before.after[,object@year==i]),1,mean))
		yearly.nb.keydays=c(yearly.nb.keydays,sum(object@year==i))
	}
    return(new("KDD.yearly.average",before.after=before.after,year=year,
    	keyday.threshold=object@keyday.threshold,
    	yearly.nb.keydays=yearly.nb.keydays))
}

cat("### Show ###\n")
.KDD.yearly.average.show <- function(object){
    cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ before.after: \n")		
    print(object@before.after)
    cat("\n ~ year: ",object@year)
    cat("\n ~ keyday.threshold: ",object@keyday.threshold)
    cat("\n ~ yearly.nb.keydays: ",object@yearly.nb.keydays)
    cat("\n")
    return(invisible())
}
setMethod(f="show",signature="KDD.yearly.average",definition=.KDD.yearly.average.show)


cat("### Getteur ###\n")
.KDD.yearly.average.get <- function(x,i,j,drop){
    switch(EXPR=i,
       "before.after"={return(x@before.after)},
       "year"={return(x@year)},
       "keyday.threshold"={return(x@keyday.threshold)},
       "yearly.nb.keydays"={return(x@yearly.nb.keydays)},
       stop("[KDD.yearly.average:get] ",i," is not a 'KDD.yearly.average' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="KDD.yearly.average",definition=.KDD.yearly.average.get)


.KDD.yearly.average.summary=function(object){
	cat("Yearly average of a keyday data object:\n")
	total.nb.keydays=sum(object@yearly.nb.keydays)
	cat("Data period: ", object@year[1],"-",object@year[length(object@year)],"\n")
	cat("Keyday threshold: ", object@keyday.threshold,"\n")
	cat("Total number of keydays: ", total.nb.keydays,"\n")
	kdd.values=object@before.after[ceiling(nrow(object@before.after)/2),]
	cat("Unweighted mean and standard deviation of rainfalls at keydays: ", 
			round(mean(kdd.values),digits=2)," and ",round(sd(kdd.values),digits=2),"\n")
	return(invisible())
}
setMethod(f="summary",signature="KDD.yearly.average",definition=.KDD.yearly.average.summary)

.KDD.yearly.average.names=function(x){
	slotNames(x)
}
setMethod(f="names",signature="KDD.yearly.average",definition=.KDD.names)
