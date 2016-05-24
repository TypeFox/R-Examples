### input fac for factors in character
###       condition the matrix

setClass(Class="Design",
    representation=representation(volume="numeric",stock="data.frame",portion="array","VIRTUAL")
)

setGeneric("writeTecan",function(object,fileName,source,destination,liquidType){standardGeneric("writeTecan")})

setGeneric("design2Screen",function(object){standardGeneric("design2Screen")})

setMethod(
	f="design2Screen",
	signature='Design',
	def=function(object){
		stock<-object@stock
		fac<-colnames(stock)
		portion<-object@portion
		col<-rep(1:dim(portion)[1],dim(portion)[2])
		row<-rep(LETTERS[1:dim(portion)[2]],each=dim(portion)[1])
		position=cbind(row,col)
		dim(portion)<-c(dim(portion)[1]*dim(portion)[2],dim(portion)[3])
		condition<-matrix(nrow=dim(portion)[1],ncol=length(fac))
		colnames(condition)<-fac
		rownames(condition)<-paste(row,col)
		weightSum<-function(x1,x2) sum(x1*x2)
		for (i in seq(along=fac)) {
			condition[,fac[i]]<-apply(portion,1,weightSum,x2=stock[,fac[i]])
		}
		screen(fac=fac,condition=condition,position=position)
	}

)