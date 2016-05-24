inmap <- function(source,layer='',attribute,type='',label='',col='',symbol='')
{
#read data
#read data
if(layer=="gisobject"||is.null(layer)||layer=="gisdata"||layer==""){
object <- source

}
else{

object <- readOGR(source, layer)
}

#map chart Region
slices=object[[attribute]]
main=label
xlab=toupper(attribute)
#default label
if(label==""){
main=paste("Map of", attribute)
}

#define colours
colour=col
if(col==""){
colour="black"
}

#construct map
if(type=="points"){
if(symbol==""){
symbol="o"
}
scot_SP_LL <- SpatialPointsDataFrame(coordinates(object), proj4string = CRS(as.character(NA)), data = as(object, "data.frame")[c(attribute)])
 plot(bubble(scot_SP_LL, main=main,attribute,col=colour,pch = symbol,scales = list(draw = TRUE)))
#convert to points
}
else{
plot(spplot(object[attribute],main=main,scales = list(draw = TRUE)))
print(spplot(object,attribute,main=main,as.table = TRUE,scales = list(draw = FALSE)))
print(spplot(object,attribute,main=main,as.table = TRUE,scales = list(draw = TRUE)))


}

data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label,colour=col)
}

#generic function
map<-function(source,layer='',attribute,type='',label="",col='',symbol='') UseMethod ("map")

#default function
map.default<-function(source,layer='',attribute,type='',label="",col='',symbol='')
{

factor<-inmap(source,layer,attribute,type,label,col,symbol)

factor$call<-match.call()
class(factor)<-"map"
factor
}

#summary function
summary.map<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Summary  \n")
print(summary(x$call))


}

#print function
print.map<-function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\nCall: Data\n")
print(x$data)

cat("\n Table  \n")
print(table(x$table))

}

#plot function
plot.map<-function(x,...)
{
map(source=x$source,layer=x$layer,attribute=x$attribute,type=x$type,col=x$colour,label=x$label)
}
#attribute=("F_AREA","Sex")
#attribute=c("F_AREA")

#graph=map("E:/geokings/advancegis/R-process","KwabibiremAtiwa",attribute=attribute,type="",label="Plant popu",col="red")

