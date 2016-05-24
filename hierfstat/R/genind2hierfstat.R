########################################################
#' @title Converts genind objects from adegenet into a hierfstat data frame
#' @description  Converts genind objects from adegenet into a hierfstat data frame
#' @usage genind2hierfstat(dat,pop=NULL)
#' @param dat a genind object
#' @param pop a vector containing the population to which each individual belongs. 
#' If pop=NULL, pop taken from slot pop of the genind object
#'  
#' @return a data frame with nloci+1 columns and ninds rows.  The first column
#' contains the population identifier, the following the genotypes at each locus
#' 
#' @examples
#' 
#' \dontrun{
#' library(adegenet)
#' data(nancycats)
#' genind2hierfstat(nancycats)
#' basic.stats(nancycats)
#' genet.dist(nancycats)
#' data(H3N2)
#' basic.stats(genind2hierfstat(H3N2,pop=rep(1,dim(H3N2@@tab)[1])),diploid=FALSE)
#' }
#' 
#' @export
##########################################################
genind2hierfstat<-function(dat,pop=NULL){
  if (!is.genind(dat)) stop("dat must be a genind object. Exiting")

  if(is.null(pop)){
    if (is.null(adegenet::pop(dat))){
      stop("population factor must be defined")
    } else {
      pop <- adegenet::pop(dat)
    }
  }
  
  if (dat@type!="codom") stop("data type must be codominant. Exiting")  
  ploid<-unique(dat@ploidy)
  if (length(ploid)!=1) stop("data must contain only diploids or only haploids. Exiting")
  if (ploid>2L) stop("Data must come from diploids or haploids. Exiting")
  if (ploid==2L) diploid<-TRUE #diploid
  if (ploid==1L) diploid<-FALSE #haploid
  nucleotides<-c("A","C","G","T")
  alleles.name<-toupper(names(table(unlist(dat@all.names))))
  nuc<-FALSE
  if(identical(alleles.name,nucleotides)) nuc<-TRUE
  

  x<-genind2df(dat,sep="",usepop=FALSE)
  #to catch alleles encoded with letters, e.g. H3N2
  if (length(grep("[A-Z]",alleles.name))==0) x<-as.matrix(data.frame(lapply(x,as.integer)))
  else {
    if (nuc){
      
      tmp<-lapply(x,function(a) gsub("A","1",a))
      tmp<-lapply(tmp,function(a) gsub("C","2",a))
      tmp<-lapply(tmp,function(a) gsub("G","3",a))
      tmp<-lapply(tmp,function(a) gsub("T","4",a))
      tmp<-lapply(tmp,function(a) gsub("a","1",a))
      tmp<-lapply(tmp,function(a) gsub("c","2",a))
      tmp<-lapply(tmp,function(a) gsub("g","3",a))
      tmp<-lapply(tmp,function(a) gsub("t","4",a))
      tmp<-lapply(tmp,as.numeric)
      x<-as.matrix(data.frame(tmp))
    }
    else (stop("alleles must be encoded as integers or nucleotides. Exiting"))
  }
  x<-data.frame(pop=pop,x)
  return(x)
}
  