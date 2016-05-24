##' getwikicoords function
##'
##' A function to return the lat/lon coordinates of towns in the UK from Wikipedia. Does not always work. Sometimes
##' the county has to be specified too.
##'
##' @param place character, ther name of the town
##' @param county character, the county it is in
##' @param rmslash remove slash from place name. Not normally used.
##' @return The lat/lon coordinates from Wikipedia
##' @export

getwikicoords <- function(place,county=NULL,rmslash=TRUE){
    cat("getting:",place,"\n")
    if(rmslash & length(grep("/",place))>=1){
        place <- strsplit(place,"/")[[1]][1]
    }
    if(substr(place,nchar(place),nchar(place))==" "){
        while(substr(place,nchar(place),nchar(place))==" "){
            place <- substr(place,1,nchar(place)-1)
        }
    }
    place <- sub(" ","_",place)
    webadd <- paste("http://en.wikipedia.org/wiki/",place,sep="")
    x <- try(readLines(webadd),silent=TRUE) 
    idx <- grep("class=\"geo\"",x)
    if(inherits(x,"try-error") | length(idx)==0){
        place <- sub("_"," ",place)
        place <- paste(place,county,sep=", ")
        place <- sub(" ","_",place)
        
        webadd <- paste("http://en.wikipedia.org/wiki/",place,sep="")
        x <- try(readLines(webadd),silent=TRUE)
        if(inherits(x,"try-error")){
            return(c(NA,NA))
        }
    }  
    idx <- grep("class=\"geo\"",x)
    if(length(idx)==0){
        return(c(NA,NA))      
    }    
    str <- x[idx[1]]
    coords <- getstrbetween(str,1,startmark="class=\"geo\">",endmark="<")
    coords <- as.numeric(unlist(strsplit(coords,";")))
    return(coords)
}
