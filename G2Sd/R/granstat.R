granstat <-
function(x,web_interface=FALSE,statistic="all",aggr=TRUE,modes=FALSE){

if(web_interface==TRUE){
  .G2Sd_web()}
  
if(web_interface==FALSE){  
  
x <- .grancompat(x)
            
STAT=c("arithmetic","geometric","folk.ward","all")
statistic <- pmatch(statistic, STAT)

if (sum(as.numeric(row.names(x)))>45) 
{
um=as.numeric(gsub("<40","0",row.names(x)))
phi=rep(0,length(um))
for (i in 1:length(um))
if (um[i]>10) phi[i]=-log2(um[i]/1000) else phi[i]=10
}

if (sum(as.numeric(row.names(x)))<=45) 
{
phi=as.numeric(row.names(x))
um=rep(0,length(phi))
for (i in 1:length(phi)) um[i]=(1/2^phi[i])*1000
}

if (modes==TRUE) mod=.mode.sedim(x,um) else mod=NULL

index=.index.sedim(x,phi,um)
textur=.texture.sedim(x,um)
descript=.sedim.descript(x,um)

if (statistic==1)
{
arith=.moment.arith(x,um)

if (aggr==TRUE)
{
result=data.frame(rbind(arith,mod,index,textur,descript))
result
}
if (aggr!=TRUE)
{
result=new.env()
result$stat=arith
result$mode=mod
result$index=index
result$sedim=rbind(textur,descript)
result=as.list(result)
}
}

if (statistic==2)
{
geom=.moment.geom(x,phi)
if (aggr==TRUE)
{
result=data.frame(rbind(geom,mod,index,textur,descript))
result
}
if (aggr!=TRUE)
{
result=new.env()
result$stat=geom
result$mode=mod
result$index=index
result$sedim=rbind(textur,descript)
result=as.list(result)
}
}


if (statistic==3)
{
fowa=.fowa.stat(x,phi,um)
if (aggr==TRUE)
{
result=data.frame(rbind(fowa,mod,index,textur,descript))
result
}
if (aggr!=TRUE)
{
result=new.env()
result$stat=fowa
result$mode=mod
result$index=index
result$sedim=rbind(textur,descript)
result=as.list(result)
}
}


if (statistic==4)
{
arith=.moment.arith(x,um)
geom=.moment.geom(x,phi)
fowa=.fowa.stat(x,phi,um)
if (aggr==TRUE)
{
result=data.frame(rbind(arith,geom,fowa,mod,index,textur,descript))
result

}
if (aggr!=TRUE)
{
result=new.env()
result$stat=rbind(arith,geom,fowa)
result$mode=mod
result$index=index
result$sedim=rbind(textur,descript)
result=as.list(result)
}
}
result
}
}
