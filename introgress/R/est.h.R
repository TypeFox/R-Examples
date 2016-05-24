est.h <-
function(introgress.data=NULL, loci.data=NULL, ind.touse=NULL,
                fixed=FALSE, p1.allele=NULL, p2.allele=NULL){
  ## This makes sure that the genotype matrix, locus data, and parental
  ## population data were supplied
  if (is.null(dim(loci.data))==TRUE)
    stop("Locus information was not supplied")
  else if (is.null(introgress.data)==TRUE)
    stop("The count was not supplied")
  if (fixed==FALSE & is.list(introgress.data)==FALSE)
    stop("introgress.data must be a list if fixed=FALSE")
  ## lets user know est.h is working
  cat("est.h is working; this may take a few minutes", fill=TRUE)
  if (is.list(introgress.data)==TRUE) admix.gen<-as.matrix(introgress.data$Admix.gen)
  ## sets up files if fixed==TRUE
  if (fixed==TRUE & (sum(loci.data[,2]=="D") + sum(loci.data[,2]=="d")) > 0)
    stop("dominant data can not be modeled as fixed")
  if (fixed==TRUE & is.null(admix.gen)==TRUE){
    if (is.null(p1.allele)==TRUE | is.null(p2.allele)==TRUE)
      stop("parental alleles must be provided if fixed==TRUE")
    admix.gen<-array(dim=dim(introgress.data))
    for (i in 1:dim(admix.gen)[1]){
      if (loci.data[i,2]=="C" | loci.data[i,2]=="c"){
        for (j in 1:dim(admix.gen)[2]){
          if (is.na(introgress.data[i,j])==TRUE) admix.gen[i,j]=NA
          else {    
            if (as.numeric(introgress.data[i,j])==2)
              admix.gen[i,j]<-as.character(paste(p1.allele,"/",p1.allele))
            else if (as.numeric(introgress.data[i,j])==1)
              admix.gen[i,j]<-as.character(paste(p1.allele,"/",p2.allele))
            else if (as.numeric(introgress.data[i,j])==0)
              admix.gen[i,j]<-as.character(paste(p2.allele,"/",p2.allele))
          }
        }
      }
      else if (loci.data[i,2]=="H" | loci.data[i,]=="h"){
        for (j in 1:dim(admix.gen)[2]){
          if (is.na(introgress.data[i,j])==TRUE) admix.gen[i,j]<-NA
          else{
            if (as.numeric(introgress.data[i,j])==1) admix.gen[i,j]<-p1.allele
            else if (as.numeric(introgress.data[i,j])==0) admix.gen[i,j]<-p2.allele
          }  
        }
      }
    }
    p1.freq<-cbind(rep(1,dim(admix.gen)[1]),rep(0,dim(admix.gen)[1]))
    p2.freq<-cbind(rep(0,dim(admix.gen)[1]),rep(1,dim(admix.gen)[1]))
    alleles<-cbind(rep(p1.allele,dim(admix.gen)[1]),rep(p2.allele,dim(admix.gen)[1]))
    introgress.data<-list(NULL,introgress.data,NULL,p1.freq,p2.freq,alleles)
    names(introgress.data)<-c("Individual.data","Count.matrix","Combos.to.use",
                              "Parental1.allele.freq","Parental2.allele.freq","Alleles")
  }
  if (fixed==TRUE & is.null(admix.gen)==FALSE){
    if (is.null(p1.allele)==TRUE | is.null(p2.allele)==TRUE)
      stop("parental alleles must be provided if fixed==TRUE")
    p1.freq<-cbind(rep(1,dim(admix.gen)[1]),rep(0,dim(admix.gen)[1]))
    p2.freq<-cbind(rep(0,dim(admix.gen)[1]),rep(1,dim(admix.gen)[1]))
    alleles<-cbind(rep(p1.allele,dim(admix.gen)[1]),rep(p2.allele,dim(admix.gen)[1]))
    introgress.data[[4]]<-p1.freq
    introgress.data[[5]]<-p2.freq
    introgress.data[[6]]<-alleles
  }
  ## subset by individual
  if (is.null(ind.touse)==FALSE) {
    if (is.character(ind.touse)==TRUE & is.character(colnames(introgress.data[[2]]))==FALSE){
      stop ("individual names were not supplied for subsetting")
    }
    admix.gen<-admix.gen[,ind.touse]
  }
  hi<-data.frame(lower=numeric(ncol(admix.gen)),
                 h=numeric(ncol(admix.gen)),
                 upper=numeric(ncol(admix.gen)))
  for(i in 1:ncol(admix.gen)){
    hi[i, ] <- h.func(geno=admix.gen[,i],
                      locustype=loci.data[,"type"],
                      r=introgress.data$Parental2.allele.freq,
                      s=introgress.data$Parental1.allele.freq,
                      alleles=introgress.data$Alleles)
  }
  return(zapsmall(hi))
}

