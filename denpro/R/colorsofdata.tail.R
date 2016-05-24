colorsofdata.tail<-function(dendat,lst,paletti=NULL)
{
# links from dendat to node to color
# "lst$infopointer" gives links from nodes to data

n<-dim(dendat)[1]
d<-dim(dendat)[2]
    
if (is.null(paletti))
paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

# links from node to color
col<-colobary(lst$parent,paletti)

# links from dendat to node (invert the links in infopointer)
nodefinder<-matrix(0,n,1)
for (i in 1:n) nodefinder[lst$infopointer[i]]<-i

datcol<-matrix("white",n,1)
for (i in 1:n){
    tok<-nodefinder[i]
    datcol[i]<-col[tok]
}

return(datacolo=datcol)
}



