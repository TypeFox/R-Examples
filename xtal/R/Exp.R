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



### model response cureve of crystal score
#source('Screen.R')

setClass(Class='Exp',
    representation=representation(screen='Screen',score='numeric'),
    validity=function(object){if(nrow(getCondition(object@screen))!=length(object@score)){
        stop ("[Screen: validation] the number of scores does not correspond")
       }else{}
       return(TRUE)}
)

setGeneric('getOptimal',function(zga){standardGeneric('getOptimal')})

setMethod(f='getOptimal',
    signature='Exp',
    def=function(zga){
       data=data.frame(score=zga@score)
       condition=getCondition(zga@screen)
       for (i in 1:ncol(condition)) {
          data=cbind(data,as.numeric(condition[,i]))
       }
       colnames(data)[2:ncol(data)]<-colnames(condition)
       mod=loess(data[,1]~data[,2]+data[,3]+data[,4],degree=2,model=TRUE)
       opt=which(predict(mod)==max(predict(mod)))
       return(condition[opt,])
    }
)



