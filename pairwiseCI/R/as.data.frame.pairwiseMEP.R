`as.data.frame.pairwiseMEP` <-
function(x,row.names = NULL, optional = FALSE, whichep=NULL, ...)
{
args<-list(...)
args$row.names<-row.names
args$optional<-optional

CIs<-x$conf.int

NP<-length(CIs)

if(is.null(whichep))
{
# all in one data.frame
outd<-NULL
for(i in 1:NP)
{
outd<-rbind(outd,CIs[[i]])
}
args$x<-outd

out<-do.call("as.data.frame", args)
}
# end all in one data.frame

else{

if(!is.numeric(whichep) &!is.character(whichep) &!is.integer(whichep))
 {stop("Argument 'whichep' must be either numeric or integer! ")}

if( length(table(whichep)) < length(whichep) )
 {warning("Some elements in 'whichep' are doubled!")}

if(is.numeric(whichep))
{

if(max(whichep)>length(CIs))
 {stop("Maximal number in 'whichep' is larger than the number of elements in 'x$conf.int'!")}

outd<-NULL
for(i in 1:length(whichep))
{
outd<-rbind(outd,CIs[[whichep[i]]])
}
args$x<-outd

out<-do.call("as.data.frame", args)
}

if(is.character(whichep))
{

if(any(!whichep %in% names(CIs) ))
 {stop("At least one element of 'whichep' can not be found in 'x$conf.int'!")}

outd<-NULL
for(i in 1:length(whichep))
{
outd<-rbind(outd, CIs[[whichep[i]]])
}
args$x<-outd

out<-do.call("as.data.frame", args)
}
}
return(out)
}

