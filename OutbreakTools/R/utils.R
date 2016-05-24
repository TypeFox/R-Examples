

####################
## extract_string ##
####################
.extract.string <- function(v.string,my.split,position,from="first"){
	if(!from%in%c("first","last")){stop("incorrect \"from\" argument")}
	res<-sapply(v.string,function(x){
        tmp<-strsplit(as.character(x),split=my.split,fixed=T)[[1]]
        if(from=="last"){
        	position=length(tmp)-position+1
        }
        return(tmp[position])})
	return(res)
}


## #####################
## ## .findDateFormat ##
## #####################
## .findDateFormat <- function(x){
##     x <- as.character(x)[1L]
##     symb <- unlist(strsplit(gsub("[[:digit:]]","", x),""))[1L]
##     temp <- unlist(strsplit(x, paste("[",symb,"]",sep="")))

##     if(nchar(temp)[1]==4){
##         return(paste("%Y","%m","%d", sep=symb))
##     }

##     if(nchar(temp)[3]==4){
##         return(paste("%d","%m","%Y", sep=symb))
##     }

##     return("")
## }

###################
## .process.Date ##
###################
.process.Date <- function(x, format=NULL){
    if(inherits(x, "Date")) return(x)

    if(is.null(format)){

        bshape <- 0
        if(is.factor(x)) x <- as.character(x)

        ## find the first non-NA value ##
        x1 <- x[!is.na(x)][1]

        ## escape if all values are NA ##
        if(is.na(x1)){
            return(as.Date(rep(NA, length(x))))
        }
        symb <- unlist(strsplit(gsub("[[:digit:]]","", x1),""))[1L]
        temp <- unlist(strsplit(x1, paste("[",symb,"]",sep="")))

        if(nchar(temp)[2]>2){
            ## months are written in letters
            ## needs to make sure the locale is fine
            bshape <- 1
            lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")

            on.exit(Sys.setlocale("LC_TIME", lct))
        }

        if(nchar(temp)[1]==4){
            if(bshape)
                format <- paste("%Y","%b","%d", sep=symb)
            else
                format <- paste("%Y","%m","%d", sep=symb)
        }
        else if(!is.na(temp[3]) && nchar(temp)[3]==4){
            if(bshape)
                format <- paste("%d","%b","%Y", sep=symb)
            else
                format <- paste("%d","%m","%Y", sep=symb)
        }
        else{
            format <- ""
            warning(paste("date is provided in an ambiguous format\n"))
            return(as.Date(rep(NA, length(x))))
        }
    }

    return(as.Date(x,format=format))
} # end .process.Date





####################
## .inlineSummary ##
####################
## FUNCTION TO GET INLINE SUMMARY FOR ONE VECTOR ##
.inlineSummary <- function(x, ...){
    cat("class: ", class(x), ",  ", sep="")
    if(!all(is.na(x))){
    if(inherits(x,"Date")) {
        cat("mean: ", format(mean(x, na.rm=TRUE), ...),
            ", range: [", format(min(x, na.rm=TRUE), ...), ";", format(max(x, na.rm=TRUE), ...), "],  ",
            sum(is.na(x)), " NAs", sep="")
    } else if(is.numeric(x)){
        cat("mean: ", format(mean(x, na.rm=TRUE), ...), ",  sd:", format(sd(x, na.rm=TRUE), ...),
            ", range: [", format(min(x, na.rm=TRUE), ...), ";", format(max(x, na.rm=TRUE), ...), "],  ",
            sum(is.na(x)), " NAs", sep="")
    } else if(is.character(x) || is.factor(x)){
        cat(length(unique(x)), " unique values,  frequency range: [",
            min(table(x)), ";", max(table(x)), "],  ",
            sum(is.na(x)), " NAs", sep="")
    }
    }else{
      if(inherits(x,"Date")) {
        cat("mean: ", format(mean(x), ...),
            ", range: [", format(min(x), ...), ";", format(max(x), ...), "],  ",
            sum(is.na(x)), " NAs", sep="")
      } else if(is.numeric(x)){
        cat("mean: ", format(mean(x), ...), ",  sd:", format(sd(x), ...),
            ", range: [", format(min(x), ...), ";", format(max(x), ...), "],  ",
            sum(is.na(x)), " NAs", sep="")
      } else if(is.character(x) || is.factor(x)){
        cat(length(unique(x)), " unique values,  frequency range: [",
            min(table(x)), ";", max(table(x)), "],  ",
            sum(is.na(x)), " NAs", sep="")
      }
    }
    cat("\n")
}





#########################
## .restoreNumericType ##
#########################
## RESTORE NUMERIC TYPE TO A VECTOR IF NEEDED ##
.restoreNumericType <- function(x){
    if(all(is.na(x))) return(x)
    x.nona <- x[!is.na(x)] 

    ## find numeric types ##
    temp <- sub("-{0,1}[[:digit:].]*e{0,1}-{0,1}[[:digit:].]*","", x.nona)
    if(all(temp=="")) return(as.numeric(as.character(x)))  
    return(x)
}





##################################
## .retrieveLabelsFromDataframe ##
##################################
## FUNCTION TO MAKE NEW UNIQUE LABELS FROM FIELDS MATCHING 'NAME' IN A DATA.FRAME ##
.retrieveLabelsFromDataframe <- function(x, sep=" ", unique=TRUE){
    ## look for fields 'name', generate unique ID ##
    temp <- grep("name", names(x), ignore.case=TRUE)
    if(length(temp)>0){
        name.fields <- names(x)[temp]

        ## check order ##
        ## 'first' name first
        temp <- grep("first", name.fields, ignore.case=TRUE)
        if(length(temp)>0) name.fields <- c(name.fields[temp],name.fields[-temp])

        ## 'last' or 'sur' names last
        temp <- grep("last", name.fields, ignore.case=TRUE)
        if(length(temp)>0) name.fields <- c(name.fields[-temp],name.fields[temp])
        temp <- grep("sur", name.fields, ignore.case=TRUE)
        if(length(temp)>0) name.fields <- c(name.fields[-temp],name.fields[temp])

        ## generate new unique names ##
        out <- apply(x[,name.fields],1,paste, collapse=sep)
        if(unique) {
            out <- make.unique(out)
        }
        return(out)
    } else
    return(NULL)
} # end .retrieveLabelsFromDataframe





