Resdata <-
function(enter=list)
{
prepdata = function(g,nomVar=NULL, nomInd=NULL){
nligne<-nrow(g)
ncolon<-( ncol(g) )/2
if( is.null(nomInd) )
{
nomInd<- paste("Row", 1:nligne, sep = "")
}
g = t(g)
g=data.frame(g)
x <- matrix(2:(nrow(g)+1), nrow = nrow(g), ncol = 1 )
colnames(x) <- paste("obs")
x=data.frame(x)
g <- cbind(x , g)
v=x%%2
colnames(v) <- paste("obs1")
v=data.frame(v)
g <- cbind(v, g)
tmin <- sqldf("select * from g where obs1 = 0")
tmax <- sqldf("select * from g where obs1 = 1")
tmin = t(tmin)
tmax = t(tmax)
tmin=tmin[-c(1:2),]
tmax=tmax[-c(1:2),]
rownames(tmin)<-nomInd
rownames(tmax)<-nomInd
colnames(tmin)<-paste("Var", 1:ncolon, sep = 'Min')
colnames(tmax)<-paste("Var", 1:ncolon, sep = 'Max')

return(list(tmin=tmin,tmax=tmax))
}
enter<-as.list(enter)
   long<-length(enter)
  resmin<-vector('list',long)
for(i in 1:long)
  resmin[[i]]<-prepdata(enter[[i]])$tmin

resmax<-vector('list',long)
for(i in 1:long)
  resmax[[i]]<-prepdata(enter[[i]])$tmax
return(list(tablemin=resmin,tablemax=resmax))
}
