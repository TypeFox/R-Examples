colorsofdata2<-function(dendat,pcf,lst,paletti=NULL,
clusterlevel=NULL,nodes=NULL)
{
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs

n<-dim(dendat)[1]
d<-dim(dendat)[2]
rnum<-length(pcf$value)

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
    
if (is.null(paletti))
paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

# links from node to color
if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
if (!is.null(clusterlevel))
col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
if (!is.null(nodes))
col<-colobary.nodes(lst$parent,nodes,paletti)

# links from rec to node (invert the links in infopointer)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

# find links from dendat to rec
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
value<-matrix(0,n,1)
for (i in 1:n){
    j<-1
    while (j<=rnum){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
            den2pcf[i]<-j
            pcf2den[j]<-i
            value[i]<-pcf$value[j]
         }
         j<-j+1
    }
}

datcol<-matrix("white",n,1)
for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
}

or<-order(value,decreasing=FALSE)
return(list(datacolo=datcol,ord=or))
}



