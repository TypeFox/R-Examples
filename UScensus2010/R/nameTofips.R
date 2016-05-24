nameTofips <-
function(name,state,type=c("county","msa"),statefips=FALSE){
	
nameTofips.aux<-function(name,state,type,statefips){	

if(!is.null(state)){
########## Guarntee State is in the right form
############Check state
state<-check.state(state,statefips)

if(is.null(state)){
	stop("Not a State! \n")
	}
}

	
##########Load and bind data
	data("MSAfips",envir = parent.frame())
	data("MSAnames",envir = parent.frame())
	data("countyfips",envir = parent.frame())
	assign("temp",countyfips)
	assign("countyfips",temp)
	assign("temp",MSAfips)
	assign("MSAfips",temp)
	assign("temp",MSAnames)
	assign("MSAnames",temp)
##########Load and bind data

###Put state in the right format
state.full<-unique(countyfips$statename)
state.ab<-unique(countyfips$acronym)
state<-state.ab[which(state==state.full)]
###Put state in the right format

########### SimpleCap from help(tolower); capitalizes first letter in each word
SimpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}
########### SimpleCap from help(tolower); capitalizes first letter in each word





msaTofips<-function(msaname,state){
		msaname<-SimpleCap(tolower(msaname))
		link<-MSAnames$name
		state<-toupper(state)
		msafips<-MSAnames$msa.cmsa.fips[which(regexpr(msaname,link)>0 & regexpr(state,link)>0)]
		msafips
}

countyTofips<-function(county,state){
		county<-tolower(tolower(county))
		link<-countyfips$countyname
		link2<-countyfips$acronym
		state<-tolower(state)
		cfips<-countyfips$fips[which(regexpr(county,link)>0 & regexpr(state,link2)>0)]
		cfips
}

if(type=="county"){
	out<-countyTofips(name,state)
	return(out)
	}

if(type=="msa"){
	out<-unique(msaTofips(name,state))
	return(out)
	}	
	
stop("Not Available!")	
}
nameTofips.aux(name,state,type,statefips)
}

