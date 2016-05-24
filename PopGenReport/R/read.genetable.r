
#must provide a filename and if pop, ind and in 2 or two columns per allele
read.genetable <- function(filename, pop=NULL, ind=NULL,lat=NULL, long=NULL, x=NULL, y=NULL, other.min=NULL, other.max=NULL,oneColPerAll,NA.char=NA,  sep=NULL, ncode=NULL ,ploidy=2 )
{
gfile <-read.csv(filename)

popnr=NULL
pops=NULL
indnr=NULL
inds=NULL

latlong=NULL

if ("pop" %in% colnames(gfile))            
  {
  pops <- as.factor(gfile$pop)
  popnr <- which(colnames(gfile)==tolower("pop"))
  }
if ("ind" %in% colnames(gfile)) 
  {
  inds <- as.factor(gfile$ind)
  indnr <- which(colnames(gfile)==tolower("ind"))
  }





genes <- gfile
rem=NULL

if (!is.null(other.min) & !is.null(other.max)) rem <- other.min:other.max         
if (!is.null(popnr)) rem <- c(rem, popnr)
if (!is.null(indnr)) rem <- c(rem, indnr)
if (!is.null(lat))  rem <- c(rem, lat)
if (!is.null(long)) rem <- c(rem, long)
if (!is.null(x)) rem <- c(rem, x)
if (!is.null(y)) rem <- c(rem, y)


if (!is.null(rem)) genes <- gfile[, -c(rem)]   else genes <- gfile


if (oneColPerAll==TRUE)
{
res <-data.frame(allele= rep(NA,dim(genes)[1]) )

for (i in seq(1,dim(genes)[2]-1,ploidy))
  {
  dummy <-genes[,i]
  for (ii in (i+1):(i+ploidy-1))
  {
  dummy<- paste(dummy,genes[,ii], sep="/")
  
  }
  
  res[,ceiling(i/ploidy)] <-  dummy
  
 

    colnames(res)[ceiling(i/ploidy)] <-  colnames(genes)[i]
  }

  
sep="/"

genes <- res
} 



  
  
  df <- df2genind(genes, pop=as.character(pops), ind.names=as.character(inds), NA.char=NA.char,     ncode=ncode,loc.names=colnames(genes), sep=sep, ploidy=ploidy)
  
  if (!is.null(lat) & !is.null(long))
  df@other$latlong <- data.frame(lat=gfile[,lat], long=gfile[,long])
  if (!is.null(x) & !is.null(y))
  df@other$xy <- data.frame(x=gfile[,x], y=gfile[,y])
  
  
  if (!is.null(other.min) & !is.null(other.max)) 
  {
  df@other$data  <- data.frame(gfile[,c(other.min:other.max)])
  }

return(df)



}
