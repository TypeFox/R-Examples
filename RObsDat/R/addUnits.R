addUnits <- function(Name, Type, Abbreviation){
	stopifnot(length(Name)==length(Type))
	stopifnot(length(Name)==length(Abbreviation))
	for(i in seq(along=Name)){
		#check existing
		if(NROW(IgetUnits(getOption("odm.handler"), Name=Name[i], Type=Type[i], Abbreviation=Abbreviation[i], exact=TRUE))>0){
			warning(paste("Skiping existing entry:", Name[i]))
			next
		}
		warning(paste("Extending Units table which should not be necessary. Please propose new term to CUASHI at http://his.cuahsi.org/mastercvreg/", sep=""))
		IaddUnits(getOption("odm.handler"), Name=Name[i], Type=Type[i], Abbreviation=Abbreviation[i])
	}
}
