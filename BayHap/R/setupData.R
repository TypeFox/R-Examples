`setupData` <-
function(data,col.snps=NULL,snps.name=NULL,sep){
   
   
   if (!is.data.frame(data)) stop("Data is not a data.frame object. \n")
   if (any(apply(data,2,is.factor))) stop("Some data columns are factors. \n")
   if (any(names(data)=="haplotypes")) stop('The variable name "haplotypes" is forbidden. \n')
   if (length(names(data))!=length(unique(names(data)))) stop("There are variables names repeated.\n")
   if (((!is.null(col.snps))& length(col.snps)==1)||((!is.null(snps.name)) & length(snps.name)==1)){ 
         stop("There are an unique SNP, is not possible to perform a haplotype analysis. \n")
   }
   if(!is.null(col.snps)){
         snps.name<-names(data)[col.snps]        
   }
   is.snp(data[,names(data)%in%snps.name],sep)
  
   haplotypes<-rep(0,nrow(data))
   data.res<-cbind(data[,names(data)%in%snps.name],haplotypes,data[,!(names(data)%in%snps.name)])
   names(data.res)<-c(snps.name,"haplotypes",names(data)[!(names(data)%in%snps.name)])
   data.res
}

