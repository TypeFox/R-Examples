check_init_length <- function(object,width,jump,type){

GL <- new.env()

if(type==1){

GL$count <- 1

init <- sapply(object@region.data@biallelic.matrix,function(x){

        if(length(x)==0){GL$count<- GL$count + 1;return(NaN)}
	repeatlength <- ceiling( (dim(x)[2]-width+1)/jump)
	if(repeatlength<=0){GL$count<-GL$count+1;return(NaN)}  
        GL$count <- GL$count + 1
        return(repeatlength)
        })

return(init)

}


if(type==2){
GL$count <- 1
init  <- sapply(object@region.data@biallelic.matrix,function(x){

        if(length(x)==0){GL$count  <- GL$count + 1;return(NaN)}
        repeatlength <- ceiling( (object@n.sites[GL$count]-width+1)/jump)
        if(repeatlength<=0){GL$count<-GL$count+1;return(NaN)}  
        GL$count <- GL$count + 1
        return(repeatlength)
        })
return(init)
}

}#END of function
