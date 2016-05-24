### input fac for factors in character
###       condition the matrix column: fac
###       condition the matrix row: well position (col first)
###       position the matrix of (wellRow,wellCol)

setClass(Class="Screen",
    representation=representation(fac="character",condition="matrix",position="matrix"),
    validity=function(object){if(length(object@fac)!=ncol(object@condition)){
       stop ("[Screen: validation] the number of factors does not correspond")
      }else{}
     return(TRUE)}
)

screen<-function(fac,condition,position) {
    colnames(condition)=fac
    rownames(condition)=paste(position[,1],position[,2])
    return(new(Class='Screen',fac=fac,condition=condition,position=position))
}

setGeneric("screenCsv",function(object,fileName){standardGeneric("screenCsv")})

setMethod(f="screenCsv",
	signature=c(object='Screen',fileName='character'),
	def=function(object,fileName){
		data=data.frame(object@position,object@condition,row.names=NULL)
		write.csv(data,fileName)
		return(invisible())
	}
)

setGeneric("getCondition",function(object){standardGeneric("getCondition")})

setMethod(f="getCondition",
        signature=c(object='Screen'),
        def=function(object){
                return(object@condition)
        }
)


