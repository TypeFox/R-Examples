inpie <- function(source,layer='',attribute,type='',label='',col='')
{
#read data
if(layer=="gisobject"||is.null(layer)||layer=="gisdata"||layer==""){
object <- source

}
else{

file=c("[.]csv$","[.]txt$","[.]tab$","[.]dat$")
i=1
nonspatial=0
while (i <= length(file)) {
log=grep(file[i],source)
if(length(log)>0)
{
file.type=file[i]
thisfile=source
nonspatial=nonspatial+1
}
else{
log=grep(file[i],layer)
if(length(log)>0){
file.type=file[i]
thisfile=paste(source,layer,sep = "/")
nonspatial=nonspatial+1
}
}
i=i+1
}
if(nonspatial<1){
if(layer==""||layer=="nofile"||layer=="notfile"||layer=="nonfile"){
object=source
}
else{
#check input data
if(source==""){
print("Full Folder path is needed:source")
return
}
if(layer==""){
print("Map name is needed: layer")
return
}

object <- readOGR(source, layer)
}
}
else{

log=grep('[.]csv$',file.type)
if(length(log)>0)
{
file.type=file[i]
object=read.csv(thisfile,header=TRUE)

}
else{
object=read.table(thisfile,header=TRUE)
}

}
}
#pie chart Region
slices=table(toupper(object[[attribute]]))
pct <- round(slices/sum(slices)*100)
lbl=names(slices)
lbls <- paste(lbl, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="")
main=label

#default label
if(label==""){
main=paste("Pie Chart of", attribute)
}

#label only names without %
if(type=="simple"){
lbls <- names(slices)
}

#define colours
colour=rainbow(length(lbls))
if(col=="gray"||col=="black"||col=="white"){
colour=gray(seq(0.4,1.0,length=10))
}



#construct pie
if(type=="map"||type=="spatial"){
plot(spplot(object[attribute],scales = list(draw = TRUE)))
}
else{

pie(slices,labels = lbls, col=colour,  
main=main)
}
data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label)
}

#generic function
piemap<-function(source,layer='',attribute,type='',label="",col='') UseMethod ("piemap")

#default function
piemap.default<-function(source,layer='',attribute,type='',label="",col='')
{

factor<-inpie(source,layer,attribute,type,label,col)

factor$call<-match.call()

class(factor)<-"piemap"
factor
}

#summary function
summary.piemap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Table ")
print (x$attribute)
print(x$table)

}

#print function
print.piemap<-function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\nCall: Data\n")
print(x$data)

cat("\n Table ")
print (x$attribute)
print(x$table)

}

#plot function
plot.piemap<-function(x,...)
{
piemap(x$source,x$layer,x$attribute,x$type,x$label)
}

