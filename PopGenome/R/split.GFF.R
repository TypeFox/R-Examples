split.GFF <- function(GFF,pos1,pos2){

# find the GFF subset for a given chunk
GFF.pos   <- GFF[,4]
regionX   <- GFF.pos<=pos2                 #<- .Call("find_windowC",GFF.pos,pos1,pos2,1)

GFF.pos   <- GFF[,5]
regionY   <- GFF.pos>=pos1                 #<- .Call("find_windowC",GFF.pos,pos1,pos2,1)

#region    <- c(regionX,regionY)
 region    <- regionX & regionY
#----------------------------------

if(sum(region)!=0){

 #region   <- sort(region)
 #start    <- region[1]
 #end      <- region[length(region)]

}else{return(NULL)}

return(GFF[region,,drop=FALSE])
#return(GFF[start:end,,drop=FALSE])
}
