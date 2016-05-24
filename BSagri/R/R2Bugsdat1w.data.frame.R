`R2Bugsdat1w.data.frame` <-
function(data, response, treatment, Intercept=FALSE)
{
if(class(data)!="data.frame")
 {stop("argument 'data' must be of class data.frame")}

if(!is.character(response) | !is.character(treatment))
 {stop("Arguments 'response' and 'treatment' must be character strings")}

if(length(response)!=1 | length(treatment)!=1)
 {stop("Arguments 'response' and 'treatment' must be single character strings")}

dnames<-names(data)

if(!response %in% dnames)
 {stop("response could not be found in data")}

if(!treatment %in% dnames)
 {stop("treatment could not be found in data")}

# formula without intercept

if(!Intercept)
 {form<-as.formula(paste(response, paste(0, treatment, sep="+"), sep="~"))}
 else
  {form<-as.formula(paste(response,  treatment, sep="~"))}

mf<-model.frame(formula=form, data=data)

resp<-mf[,1]

if(!class(resp) %in% c("numeric", "integer"))
 {stop("Ther response variable in data must be integer or numeric")}

# create the bugsdat, appropriate for the model

mm<-as.data.frame(model.matrix(object=form, data=mf))

pnames<-colnames(mm)

ngroups<-length(pnames)

lnames<-paste("X", 1:ngroups, sep="")

bugsdat<-as.list(mm)
names(bugsdat)<-lnames

bugsdat$Y<-resp
bugsdat$P<-ngroups
bugsdat$N<-length(resp)

ini<-list(
beta=rep(0,ngroups),
r=1)


splitlist<-split(mf[,1],f=mf[,2])
ni<-unlist(lapply(splitlist, length))

out<-list(bugsdat=bugsdat,
parameters="muvec",
inits=ini,
data=data,
names=list(pnames=pnames, ni=ni),
Intercept=Intercept
)
class(out)<-"R2Bugsdat1w"
return(out)
}

