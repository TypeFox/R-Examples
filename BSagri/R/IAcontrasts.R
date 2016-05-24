
##############################

c2compnames<-function(cmat, ntype="aggr")
{

if(is.null(colnames(cmat)))
 {colnames(cmat)<-1:ncol(cmat)}

if(!is.matrix(cmat))
 {stop("cmat must be a matrix")}

if(!is.numeric(cmat))
 {stop("cmat must be a numeric matric")}

if(length(ntype)!=1)
 {stop("ntype must be a single character string!")}

if(!ntype %in% c("aggr", "sequ"))
 {stop("ntype must be one of 'aggr', 'sequ'")}


if(ntype=="aggr")
{
cnames<-colnames(cmat)

rnames<-character(length=nrow(cmat))

for(i in 1:nrow(cmat))
{
si<-sign(cmat[i,])

wsip<-si==1
wsin<-si==(-1)
rnames[i]<-paste( "(", paste(cnames[wsip], collapse="+"), ")-(", paste(cnames[wsin], collapse="+"), ")", collapse="", sep="" )
}
}


if(ntype=="sequ")
{
cnames<-colnames(cmat)

rnames<-character(length=nrow(cmat))

for(i in 1:nrow(cmat))
{
si<-sign(cmat[i,])


wsin0<-si!=0
wsip<-si[wsin0]==1
wsin<-si[wsin0]==(-1)

nam<-cnames[wsin0]
sic<-character(length=length(nam))

sic[wsip]<-"+"
sic[wsin]<-"-"

rnames[i]<-paste(paste(sic, nam, sep=""), collapse="")

}
}

rownames(cmat)<-rnames

return(cmat)

}

###########################

IAcontrasts<-function(type, k)
{

if(!all(type %in% c("Dunnett", "Tukey", "Sequence", "Identity")))
 {stop("all elements of type must be one of 'Dunnett','Tukey' or 'Sequence'")}

if ( any(c(length(k),length(type))!=2))
 {stop("k and type must be vectors of length 2")}

if(!is.numeric(k) & !is.integer(k))
 {stop("k must be an integer vector")}

n1<-rep(3,k[1])
names(n1)<-as.character(1:k[1])

n2<-rep(3,k[2])
names(n2)<-as.character(1:k[2])


if(type[1]!="Identity"){CM1 <- contrMat(n=n1, type=type[1])}
else{CM1<-diag(rep(1,k[1]))}

if(type[2]!="Identity"){CM2 <- contrMat(n=n2, type=type[2])}
else{CM2 <-diag(rep(1,k[2]))}

out <- kronecker(CM1, CM2)

cnames <- paste( rep(colnames(CM1), each=length(colnames(CM2))),
 rep(colnames(CM2), times=length(colnames(CM1))), sep="")

colnames(out)<-cnames

return(out)

}

############################

IAcontrastsCMAT<-function(CMAT1, CMAT2)
{

out <- kronecker(CMAT1, CMAT2)

cnames <- paste( rep(colnames(CMAT1), each=length(colnames(CMAT2))),
 rep(colnames(CMAT2), times=length(colnames(CMAT1))), sep="")

colnames(out)<-cnames

return(out)

}






