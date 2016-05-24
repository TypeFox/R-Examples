setClass(
		Class = "NG_data",
		representation = representation(
				name = "character",
				shortnames = "character",
				group = "OptionalCharNum",
				labels = "character",
				data = "data.frame"
		),
		validity = function(object){
			goThrough <- TRUE
			
			n <- dim(object@data)[1]
			nvar <- dim(object@data)[2]
			
			if(length(object@shortnames) != 0){
				if(length(object@shortnames) != nvar){
					cat("[NG_data:validation] shortnames has wrong dimension.\n")
					goThrough <- FALSE
				}
				if(all(grepl(' ', object@shortnames, fixed=TRUE))) {
					cat("[NG_data:validation] shortnames can not contain white spaces.\n")
					goThrough <- FALSE
				}
			}else {
				if(all(grepl(' ', names(object@data), fixed=TRUE))) {
					cat("[NG_data:validation] names(data) contains at least one white space. Define shortnames without white spaces.\n")
					goThrough <- FALSE
				}
			}
			if(length(object@group) != 0){
				if(length(object@group) != n){
					cat("[NG_data:validation] group has wrong dimension.\n")
					goThrough <- FALSE
				}
			}
			if(length(object@labels) != 0){
				if(length(object@labels) != n){
					cat("[NG_data:validation] labels has wrong dimension.\n")
					goThrough <- FALSE
				}
			}
			return(goThrough)
		}
)



## create NG_data object
ng_data <- function(name, data, shortnames = character(0), group = numeric(0), labels = character(0)){
	if(length(group)!=0) {
		if(is(group,"factor")) {
			group <- as.numeric(group)
		}
	}
	if(length(labels)!=0) {
		if(is(labels,"factor")) {
			labels <- as.character(labels)
		}
	} else {
          labels <- as.character(1:dim(data)[1])
        }
	
	name = gsub(' ','',name, fixed=TRUE)
	
	## navGraph can only handle numeric values within the data fram
	isNumeric <- sapply(data,function(var){is.numeric(var)})
	if(sum(isNumeric == FALSE) != 0) {
		ii <- which(isNumeric == FALSE)
		nam <- names(data)
		for(i in ii) {
			 cat(paste('Variable: ',nam[i], ' is not numeric!\n',sep=''))	
		}
		stop("[ng_data] data argument must be a data.frame with solely numeric variables.\n")
	}
	
	new("NG_data", data = data, name = name, shortnames = shortnames, group = group, labels= labels)
}



## Printvector and Tableelement are important for the show method
printvector <- function(x){
	n <- length(x)
	if(n>1){
		txt <- paste(paste(x[1:(n-1)], collapse = ', '), ', ',x[n], '.', sep = '')
	}else{
		txt <- paste(x[1],'.', sep = '')
	}
	return(txt)
}


tableelement <- function(x, indent, n){
	txt <- paste(
			paste(rep(' ',indent), collapse = ''),
			x,
			paste(rep(' ', max(n - nchar(x),0)), collapse = ''),
			sep = ''
	)
	return(txt)
}


## show Method
setMethod("show","NG_data",
		function(object){
			cat("object from NG_data class.\n")
			cat(paste("  name:", object@name,'\n'))
			
			dimData <- dim(object@data)
			
			cat(paste("  data:", dimData[1] ,'x', dimData[2]),'\n')
			
			maxChar <- max(nchar(c(names(object),'Variable Names')))
			nam <- names(object)
			
			if(length(object@shortnames) !=0){
				snam <- object@shortnames
				cat(paste(tableelement("Variable Names", 4, maxChar),' | ', 'Short Names\n', sep = ''))
			}else{
				snam <- rep('',length(nam))
				cat(paste(tableelement("Variable Names", 4, maxChar),' | ', 'No Short Names\n', sep = ''))
			}
			cat(paste('    ',paste(rep('-',2*maxChar),collapse = ''),'\n', sep = ''))
			
			for(i in 1:dimData[2]){
				cat(paste(tableelement(nam[i], 4, maxChar),' | ', snam[i] ,'\n', sep = ''))
			}
			
			
			if(length(object@group) == 0){
				cat("  group:  No group variable defined.\n")
			}else{
				cat(paste("  group: ", length(unique(object@group)), 'groups.\n'))
			}
			
			if(length(object@labels) == 0){
				cat("  label: labels weren't defined.\n")
			}else{
                          if(length(unique(object@labels))>25) {
                            cat(paste("  labels: ",	printvector(unique(object@labels)[1:25]),"..\n", sep = ''))
                          } else { 
                            cat(paste("  labels:",	printvector(unique(object@labels)),"\n"))
                          }
                        }
		})


### names
setMethod("names","NG_data",
		function(x){
			return(names(x@data))
		})


### names<-
setReplaceMethod(
		f = "names","NG_data",
		function(x,value){
			names(x@data) <- value
			return(x)
		})


## shortnames
setMethod("shortnames","NG_data",
		function(x){
			return(x@shortnames)
		})


## shortnames<-
setReplaceMethod(
		f = "shortnames","NG_data",
		function(x,value){
			x@shortnames <- value
			validObject(x)
			return(x)
		})

## acessor funtion []
#setMethod(f = '[',
#		signature = "NG_data",
#		definition = function(x,i,j,drop){
#				return(x[,match(j,x@shortnames)])				
#			}else{
#				callNextMethod()
#		}
#)

setMethod(f = "ng_get",
		signature = "NG_data",
		definition = function(obj, what=NULL, ...){
			possibleOptions <- c("name","data","group","labels")
			
			if(is.null(what)){
				cat("Get what? Possible options are: ")
				cat(paste(possibleOptions, collapse = ", "))
				cat("\n")
			}else{
				if(any(is.na(match(what,possibleOptions)))){
					stop(paste("[ng_get] object",what,"is not defined."))
				}
				
				if(what == "name"){
					return(obj@name)
				}else if(what == "data"){
					return(obj@data)
				}else if(what == "group"){
					return(obj@group)
				}else if(what == "labels"){
					return(obj@labels)
				}
			}
		})

setMethod(f = "ng_set",
		signature = "NG_data",
		definition = function(object){
		possibleOptions <- c("name","data","group","labels")
		cat(paste("Replace what? Possible options are:", paste(possibleOptions, collapse = ", "),'\n'))
		cat('Use ng_set<- to set a value.\n')
	})

setReplaceMethod(
		f = "ng_set","NG_data",
		function(object,what="",value){
						
			possibleOptions <- c("name","data","group","labels")
			tmp <- object
			if(!(what %in% possibleOptions)) {
				stop(paste("Replace what? Possible options are:", paste(possibleOptions, collapse = ", "),'\n'))
			}	
			
			slot(object, what, check = TRUE) <- value 
			
			if(validObject(object)){
				return(object)
			}else {
				cat("[ng_set] assignment is wrong\n")
				return(tmp)
			} 
		})
