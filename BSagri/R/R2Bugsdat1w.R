`R2Bugsdat1w` <-
function(formula, data)
{
if(class(data)!="data.frame")
 {stop("argument 'data' must be of class data.frame")}

if(length(formula[[2]])!=1)
 {stop("The left hand side of 'formula' must be a single variable name")}

if(length(formula[[3]])!=1)
 {stop("The right hand side of 'formula' must be a single variable name")}

mf<-model.frame(formula=formula, data=data)

print(mf)

resp<-mf[,1]

if(!class(resp) %in% c("numeric", "integer"))
 {stop("Ther response variable in data must be integer or numeric")}

# create the bugsdat, appropriate for the model

mm<-as.data.frame(model.matrix(object=formula, data=mf))

print(mm)

Y<-resp

pnames<-colnames(mm)

ngroups<-length(pnames)

lnames<-paste("X", 1:ngroups, sep="")

bugsdat<-as.list(mm)
names(bugsdat)<-lnames

bugsdat$Y<-resp
bugsdat$P<-ngroups
bugsdat$N<-length(resp)

ini<-list(
beta=rep(1,ngroups),
r=1)

splitlist<-split(mf[,1],f=mf[,2])
ni<-unlist(lapply(splitlist, length))

out<-list(bugsdat=bugsdat,
parameters="muvec",
inits=ini,
data=data,
names=list(pnames=pnames, ni=ni),
Intercept=TRUE
)
class(out)<-"R2Bugsdat1w"
return(out)
}

