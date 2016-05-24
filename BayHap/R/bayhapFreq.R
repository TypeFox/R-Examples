`bayhapFreq` <-
function(data, na.snp.action = "omit", snps.name = NULL,
                 col.snps = NULL, sep = "/", total.iter.haplo = 10000,
                 devhaplo = 0.1, dist = 1, lag.number = 10, sign =
                 0.05, file = FALSE, verbose = 2)




{
       burn.in<-5000
       type.pheno<-0
       freqmin<-0
       devlogit<-0
       type.model<-0
       varcont<-0
       ncovars<-0
       covariates.vec<-0
       covars.name<-"a"
       inters<-0
       intcovars<-0
       prior<-0
       mu<-0
       varp<-0
       burn.in.pheno<-100
       #1000
       total.iter<-100
       #total.iter.haplo 
       covar<-FALSE
    

    
    if (!is.data.frame(data)) stop("Data is not a data.frame object. \n")
    if (sum(apply(data,2,is.factor))>0) stop("Some data columns are factors. \n")
    if ((freqmin>1)||(freqmin<0)) stop("freqmin is negative or higher than 1 \n")
    if (devhaplo<0) stop("devhaplo value is negative \n")
    if (lag.number>total.iter.haplo) stop("lag number must to be lower than number of iters \n")
    if (sign>1) stop("sign value is out of the interval (0,1) \n")
    if ((length(col.snps)!=0 & length(col.snps)==1)||(length(snps.name)!=0 & length(snps.name)==1)){ 
         stop("There are an unique SNP, is not possible to compute haplotype frequencies. \n")
    }
    
    
    if (length(snps.name)!=0) geno<-data[,names(data)%in%snps.name]
    if (length(col.snps)!=0){
        geno<-data[,col.snps]
        snps.name<-names(data)[col.snps]
    }
    if ((length(snps.name)==0)&(length(col.snps)==0)) snps.name<-names(data)
      
    if (na.snp.action=="omit"){
        geno<-geno[(apply(apply(geno,2,is.na),1,sum)==0),]  
        miss<-0
    }else{
        if (na.snp.action=="keep"){
    geno<-geno[apply(geno,1,FUN=function(x) sum(is.na(x))!=ncol(geno)),]
    miss<-ifelse(sum(apply(apply(as.data.frame(geno),2,is.na),1,sum)==0)==ncol(geno),0,1)
        }else{
           stop("na.snp.action value is not correct. \n")
        }
    } 
    if (nrow(geno)==0) 
        stop("After apply na.geno.action there are not enough rows with non missing values to perform the analyse.\n")

   
    is.snp(geno,sep) 
   

    nindiv<-nrow(geno)
    nlocus<-ncol(geno)
    pheno<-rep(0,nindiv)
    time<-rep(0,nindiv)

    mat<-NULL

    alleles_mat<-NULL

    alleles<-NULL

    
    dades<-geno 
    numindivexp<-nindiv
    weights<-rep(1,nindiv)
  
    j<-1

    i<-1

    for (j in 1:nlocus)

        {

                geno<-genotype(dades[,j],sep=sep)

                a<-allele(geno) 

                alleles<- sort(attr(a,"allele.names"))

                if (length(alleles)==1) alleles<-c(alleles,alleles)

                alleles_mat<-rbind(alleles_mat,alleles)

                a[is.na(a)]<-"8"

                a[a==alleles[1]]<-"0"

                a[a==alleles[2]]<-"1"

                mat<-cbind(mat,a)

        } 

    mat<-apply(mat,2,as.numeric)

    
    cat("Computing haplotype pairs. \n")

    pairs<-haplo.pairs(mat)
    n.pairs.indiv<-pairs[[1]]    
    pairs<-as.numeric(unlist(strsplit(unlist(pairs[2]),split="-")))
    pairs.m<-matrix(pairs,byrow=TRUE,ncol=2)
    codi.haplo.cromo1<-pairs.m[,1]
    codi.haplo.cromo2<-pairs.m[,2]


       

    genotypes.vec<-c(t(mat))

    genotypes.vec<-as.numeric(genotypes.vec)

  
    id.outhaplo<-paste(gsub("[./: ]","",format(Sys.time(), "%d/%m/%y %H:%M:%OS3")),paste(sample(LETTERS)[1:10],collapse=""),sep="")

    res.c <- .C("mcmc",numindiv=as.integer(nindiv),

               numindivexp=as.integer(numindivexp),

               weights=as.double(weights),

               numlocus=as.integer(nlocus),

               burn.in=as.integer(burn.in),

               burn.in.haplo=as.integer(total.iter.haplo),

               burn.in.pheno=as.integer(burn.in.pheno),

               total.iter=as.integer(total.iter),

               genotypes=as.integer(0),

               phenotype=as.integer(pheno),

               max_haplo=as.integer(1000),

               distribution=as.integer(dist),

               type.pheno=as.integer(type.pheno),

               freqmin=as.double(freqmin),

               devhaplo=as.double(devhaplo),

               devlogit=as.double(devlogit),

               lag.number=as.integer(lag.number),

               time_var=as.double(time),

               mem_error=as.integer(0),

               res_freq=as.double(numeric(2^nlocus)),

               res_se_freq=as.double(numeric(2^nlocus)),

               res_coeff=as.double(numeric(2^nlocus+ncovars+2^nlocus*ncovars+1)),

               res_se_coeff=as.double(numeric(2^nlocus+ncovars+2^nlocus*ncovars+1)),

               haplos_int=as.double(numeric(2^nlocus)),

               llista_haplos=as.integer(numeric((2^nlocus)*nlocus)),

               n_coeff=as.integer(0),

               type_model=as.integer(type.model),

               var_cont=as.double(varcont),

               num_covar=as.integer(ncovars),

               covariates=as.double(covariates.vec),

               interactions=as.integer(inters),

               int_covars=as.integer(intcovars),

               prior=as.integer(prior),

               mu=as.double(mu),

               varp=as.double(varp),
               hap1_R=as.integer(codi.haplo.cromo1),
               hap2_R=as.integer(codi.haplo.cromo2),
               n_pairs=as.integer(n.pairs.indiv),
               verbose=as.integer(verbose),
               versem_max=as.double(0),
               onlyfreq=as.integer(1),
               id_outhaplo=as.character(id.outhaplo)

               )


        if (res.c$mem_error==1) stop("Memory error due to number of subjects.")

        if (res.c$mem_error==2) stop("Memory error due to number of locus.") 
 
        unlink(paste("outhaplo",id.outhaplo,".txt",sep=""))
        unlink("a.txt")
        unlink("b.txt")
        mat.iters<-import(paste("outhaplo2",id.outhaplo,".txt",sep=""))  
        if (file==FALSE)  unlink(paste("outhaplo2",id.outhaplo,".txt",sep="")) 
        
        res.c$llista_haplos<-haplo.list.all(nlocus)
        r<-list(res.c,sign,sep,snps.name,alleles_mat,mat.iters)                                

class(r)<-"freq"
return(r)

}

