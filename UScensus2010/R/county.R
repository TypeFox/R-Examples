county <-
function(fips=NULL,name=NULL,state,level=c("tract","blk","blkgrp"),statefips=FALSE,sp.object=NULL,proj=NULL){

county.aux <-
function(fips=NULL,name=NULL,state,level=c("tract","blk","blkgrp"),statefips=FALSE,sp.object=NULL,proj=NULL){

data("countyfips",envir = parent.frame())
assign("temp",countyfips)
assign("countyfips",temp)
############Check state
state<-check.state(state,statefips)

if(is.null(state)){
	cat("Not a State! \n")
	return()
	}
############Check state
	
	
########### Function to pullout counties
build.county.shape<-function(fips,state,sp.object){		
require(paste("UScensus2010",level,sep=""),character.only=TRUE)

if(is.null(sp.object)==FALSE){
	temp<-sp.object
	}else{
	x<-paste(state,".",level,"10",sep="")
	data(list=x,envir = parent.frame())
	temp<-get(x)
}



if(level=="blk"){
	out<-temp[which(substr(temp$fips,1,3)%in%fips==TRUE),]
	}else{
	out<-temp[which(temp$county%in%fips==TRUE),]
}

return(out)
}
########### Function to pullout counties


########### Function to handle county names
check.county.name<-function(name){
	
	            name<-tolower(name)
                link <- countyfips$countyname
                # state
                fips <- countyfips$fips[which(regexpr(name,link) > 0 & state==countyfips$statename)]
                
                if(identical(fips, character(0)))
                	fips<-NULL
                
                return(fips)
	}
########### Function to handle county names


########### Check on sp object

if(is.null(sp.object)==FALSE & class(sp.object)[1]!="SpatialPolygonsDataFrame"){
		cat("Not a SpatialPolygonsDataFrame object! \n")
		return()
		}



###########


######Case 1
if(!is.null(fips)){
	if(sum(fips%in%substr(countyfips$fips,3,5))==0){
			cat("Not a valid county FIPS!")
			return()
		}
	
	out<-build.county.shape(fips,state,sp.object)
	}else{
		
		if(length(name)>1){
			fips<-vector(length=length(name))
			for(i in 1:length(name)){
				if(is.null(check.county.name(name[i]))==FALSE){
					fips[i]<-check.county.name(name[i])
				}else{
						cat("Not a valid county name!")
						return()
				}
			}
			}else{
				fips<-check.county.name(name)
				if(is.null(fips)){
					cat("Not a valid county name!")
					return()
					}
			}
out<-build.county.shape(substr(fips,3,5),state,sp.object)
}
		

		



## Check proj
if(is.null(proj)==FALSE){
	require(rgdal)
	out<-spTransform(out,proj)
}
##check proj
out
}

out<-county.aux(fips,name,state,level,statefips,sp.object,proj)
}