
`pairwiseMEP` <-
function(x,...){UseMethod("pairwiseMEP")}


`pairwiseMEP.data.frame` <-
function(x, ep, f,
 control=NULL, conf.level=0.95,
 alternative=c("two.sided","less","greater"),
 method="Param.diff",
 ...)
{

PCargs<-list(...)

data<-x

if(!is.data.frame(data)) {stop("Argument x must be of class 'data.frame'!")}
if(!is.character(ep)) {stop("Argument ep must be a character vector!")}
if(!is.character(f)) {stop("Argument f must be a single character string!")}
if(!is.character(method)) {stop("Argument method must be a character vector!")}
if(length(f)!=1) {stop("Argument f must be a single character string!")}
if(any(method%in%c("Prop.diff", "Prop.ratio","Prop.or"))){stop("Methods for binary data are currently not available for multiple endpoints!")}

tabep<-table(ep)
if(any(tabep!=1)){
wg1 <- which(tabep>1)
warning(paste("The value (s) ", paste(names(tabep)[wg1], collapse=", ") ," occur(s) more than once in the 'ep'!" , sep="" ))
}

NP<-length(ep)

if(length(method)!=1)
 {
 if(length(method)!=NP)
  {stop("Argument method must be either a single character string or a vector of character strings equal to the length of the vector in argument 'ep'")}
 }
else
 {
  method<-rep(method, length.out=NP)
 }



alternative<-match.arg(alternative)

nam<-names(data)

 if(any(!ep %in% nam))
 {
 wepnn<-which(!ep %in% nam)
 stop(paste("Variable(s) ",paste(ep[wepnn], collapse=", ")," could not be found in 'x'!", sep=""))
 }

 if(all(nam!=f))
 {stop(paste("Variable f (",f,") could not be found in 'x'!"))}

CLI <- conf.level

outl<-list()

 for(i in 1:NP)
{
form<-as.formula(paste(ep[i],"~",f))

PCargs$formula<-form
PCargs$data<-data 
PCargs$control<-control
PCargs$conf.level<-CLI
PCargs$alternative<-alternative
PCargs$method<-method[i]
temp<-do.call("pairwiseCI", args=PCargs)

dat<-as.data.frame(temp)
ni<-nrow(dat)
mout<-rep(method[i], ni)
epout<-rep(ep[i], ni)

outl[[i]] <- cbind(dat, "method"=mout, "response"=epout)

}

names(outl)<-ep

out<-list(
 conf.int=outl,
 data=data,
 ep=ep,
 f=f,
 control=control,
 conf.level=conf.level,
 alternative=alternative,
 method=method
)

class(out)<-"pairwiseMEP"

return(out)
}

