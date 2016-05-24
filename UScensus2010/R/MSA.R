
MSA <-
function(msafips=NULL,msaname=NULL,state=NULL,statefips=FALSE,level=c("tract","blk","blkgrp"),proj=NULL){
MSA.aux <-
function(msafips=NULL,msaname=NULL,state=NULL,statefips=FALSE,level=c("tract","blk","blkgrp"),proj=NULL){

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

########### SimpleCap from help(tolower); capitalizes first letter in each word
SimpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}
########### SimpleCap from help(tolower); capitalizes first letter in each word

if(!is.null(state)){
########## Guarntee State is in the right form
############Check state
state<-check.state(state,statefips)

if(is.null(state)){
	stop("Not a State! \n")
	}
############Check state


state.full<-unique(countyfips$statename)
state.ab<-unique(countyfips$acronym)

state<-state.ab[which(state==state.full)]

########## Guarntee State is in the right form
}
	
msa.fips<-function(msafips,level){	
	m<-MSAfips$msa.cmsa.fips%in%msafips
	state<-MSAfips$fips.state[m]
	fips<-MSAfips$fips.county[m]

if(length(unique(state))==1){
	out<-county(fips=unique(fips),state=state[1],level=level,statefips=TRUE)
	return(out)
}else{
	#fips<-unique(fips)
	state.u<-unique(state)
	
 out<-county(fips=fips[which(state==state.u[1])],state=state.u[1],level=level,statefips=TRUE)
	for(i in 2:length(state.u)){
		out<-spRbind(out,county(fips=fips[which(state==state.u[i])],state=state.u[i],level=level,statefips=TRUE))
		}
		return(out)
	}
}

if(!is.null(msafips)){
	#####MSA fips code provided
	out<-msa.fips(msafips,level)
	return(out)
	
}else if(!is.null(msaname)){
	######### MSA name (e.g. "portland" and state="or")
	if(!is.null(state)){
		msaname<-SimpleCap(tolower(msaname))
		link<-MSAnames$name
		state<-toupper(state)
		msafips<-MSAnames$msa.cmsa.fips[which(regexpr(msaname,link)>0 & regexpr(state,link)>0)]
		out<-msa.fips(unique(msafips),level)

		}else{
		#######Case full MSA name plus state
		msafips<-MSAnames$msa.cmsa.fip[MSAnames$name%in%msaname]
		out<-msa.fips(msafips,level)
	}
}

##Check proj
if(is.null(proj)==FALSE){
	require(rgdal)
	out<-spTransform(out,proj)
}
##Check proj
out
	
}
out<-MSA.aux(msafips,msaname,state,statefips,level,proj)
}