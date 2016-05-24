#' @title Write ped file plink analysis
#' @description write a ped and a map file suitable for analysis with 
#' \href{plink}{http://pngu.mgh.harvard.edu/~purcell/plink/} 
#' @usage write.ped(dat, ilab = NULL, pop = NULL, 
#'         fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL,sex=NULL) 
#' @param dat a hierfstat data frame
#' @param ilab individal labels
#' @param pop population id
#' @param fname filename for ped file
#' @param na.str character string to use for missing values
#' @param f.id father id. default to unknown
#' @param m.id mother id. default to unknown
#' @param loc.pos the loci position default to unknown
#' @param sex the individual sex. default to unknown
#' @return a map file containing the loci positions 
#' @return a ped file containing genotypes etc... 
#' @references \href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/}{Purcell etal (2007) PLINK:} A Tool Set for Whole-Genome Association 
#'  and Population-Based Linkage Analyses 81:559-575
#' @export  
####################################################################################
write.ped<-function (dat, ilab = NULL, pop = NULL, fname = "dat",na.str="0",f.id=NULL,m.id=NULL,loc.pos=NULL,sex=NULL) 
{
    dum <- getal.b(dat[, -1])
    dum[is.na(dum)] <- na.str
    nind <- dim(dum)[1]
    nloc <- dim(dum)[2]
	ddum<-matrix(numeric(nind*nloc*2),nrow=nind)
	al2<-(1:nloc)*2
	al1<-al2-1
	ddum[,al1]<-dum[,,1]
	ddum[,al2]<-dum[,,2]
  if (is.null(pop)) popid <- dat[, 1] else popid<-pop
	if (is.null(ilab)) ind.id<-1:nind else ind.id<-ilab
	if(is.null(f.id))  f.id<-rep(0,nind) 
	if(is.null(m.id)) m.id<-rep(0,nind)
	if(is.null(sex)) sex<-rep(0,nind)
    locnames <- paste("L", names(dat)[-1], sep = "")
	mapf<-cbind(0,locnames,0,0)
    write.table(mapf, paste(fname,".map",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	datn<-data.frame(fam.id=popid,ind.id=ilab,f.id=f.id,m.id=m.id,sex=sex,pheno=rep(0,nind),ddum)
	write.table(datn,paste(fname,".ped",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}
