`bayhapReg` <-
function(formula = formula, data, family = "gaussian", t.model
                 = "additive", sep = "/", na.snp.action = "omit",
                 freqmin = 0.01, burn.in.haplo = 10000, burn.in.pheno =
                 10000, total.iter = 10000, devhaplo = 0.1, dist = 1,
                 lag.number = 10, sign = 0.05, file = FALSE, prior.val
                 = haplo.prior(), verbose = 2)
{
    if (class(formula)!="formula") stop("Formula is not correct. \n")
    if (!is.data.frame(data)) stop("Data is not a data.frame object. \n")
    if (length(names(data))!=length(unique(names(data)))) stop("Some columns of data object have the same name.\n")
    if (any(apply(data,2,is.factor))) stop("Some data columns are factors. \n")
    if ((freqmin>1)||(freqmin<0)) stop("Freqmin value must be non negative and lower than 1. \n")
    if (devhaplo<0) stop("devhaplo values must to be positive \n")
    if (lag.number>total.iter) stop("The total number of iterations must to be higher than lag number. \n")
    if (sign>1) stop("sign value is out of the interval (0,1) \n")
    if (!(any(names(data)=="haplotypes"))) stop("Execute setupData function before running bayhap.reg. \n") 
    
    devlogit<-0.1
    burn.in<-5000

    #checking the family 
   
    if (family[1]=="gaussian"){
             type.pheno<-2        
    }else{
         if (family[1]=="weibull"){
             type.pheno<-1
         }else{
             if (family[1]=="binomial"){
                type.pheno<-0
             }else{
                stop("The family name is not correct. Please check it.\n")
             }
         }
    }

    #checking type model
    
    if (t.model=="additive"){
              type.model<-0
    }else{
         if (t.model=="dominant"){
              type.model<-1
         }else{
              if (t.model=="recessive"){
                  type.model<-2
              }else{
                  stop("The model type is not correct. Please check it.\n")
              }
         }
    } 

    #names 

    i<-1
    names.tot<-NULL
    for (i in 1:length(attr(terms(formula),"variables"))){
       names.tot<-c(names.tot,as.character(attr(terms(formula),"variables")[[i]]))
    }
    d<-1:ncol(data)
    n.snps<-min(d["haplotypes"==names(data)])-1
    snps.name<-names(data)[1:n.snps]
    names.tot<-c(names.tot,snps.name)
    names.tot<-names.tot[names.tot!="list"&names.tot!="survival"]
   
    if (length(unique(names.tot))!=length(names.tot)) stop("There are repeated variable names, please check it.\n")
   
    data.form<-data[,names(data)%in%names.tot]
    col.na<-as.numeric(apply(apply(data.form,1,is.na),1,sum))
    names.tot.na<-names(data.form)[col.na==nrow(data.form)]
    
    resp.name<-as.character(attr(terms(formula),"variables")[[2]])
    if (type.pheno==1) resp.name<-resp.name[resp.name!="survival"]
    covars.name<-names(data.form)[names(data.form)%in%names.tot[!(names.tot %in% c(resp.name,snps.name,"haplotypes"))]]
    
    #missings

    if (any(resp.name%in%names.tot.na)) 
       stop("Response variable values are all missing. The analysis can not be performed.\n")
    if (any(snps.name%in%names.tot.na)) 
       stop("All values are missing for at least one SNP.\n")
    if (any(covars.name%in%names.tot.na)) 
       stop("All values are missing for at least one covariate. \n")     
    
    if (na.snp.action=="omit"){
        data.form<-data.form[(apply(apply(data.form[,1:length(snps.name)],2,is.na),1,sum)==0),]  
        miss<-0
    }else{
        if (na.snp.action=="keep"){
    data.form<-data.form[apply(data.form[,1:length(snps.name)],1,FUN=function(x) sum(is.na(x)))!=length(snps.name),]
    miss<-ifelse(sum(apply(apply(as.data.frame(data.form[,1:length(snps.name)]),2,is.na),1,sum)==0)==length(snps.name),0,1)
        }else{
           stop("na.snp.action value is not correct. \n")
        }
    } 
    if (nrow(data.form)==0) 
        stop("After apply na.geno.action, there are not enough rows with non missing values to perform the analysis.\n")

    geno<-data.form[,(names(data.form)%in%snps.name)] 
    is.snp(geno,sep) 

    #phenotype variables    

    if (type.pheno==0){
       pheno <- model.response(model.frame(formula,data=data.form,na.action=NULL), "any")
       if (is.factor(pheno)) pheno<-as.character(levels(pheno))[pheno]
       if (length(unique(pheno))>2) stop("Response variable takes more than two different values for a logistic regression!\n")
       pheno.fact<-as.factor(pheno)
       pheno[pheno.fact==levels(pheno.fact)[1]]<-0
       pheno[pheno.fact!=levels(pheno.fact)[1]]<-1
       pheno<-as.numeric(pheno)
       time<-rep(0,nrow(data.form))
    }else{
       if (type.pheno==1){
          pheno<-data.form[,names(data.form)==resp.name[2]]
          if (is.factor(pheno)) pheno<-as.character(levels(pheno))[pheno]
          if (length(unique(pheno))>2) stop("Censure variable takes more than two different values.\n")
          pheno.fact<-as.factor(pheno)
          pheno[pheno.fact==levels(pheno.fact)[1]]<-0
          pheno[pheno.fact!=levels(pheno.fact)[1]]<-1
          pheno<-as.numeric(pheno)
          time<-data.form[,names(data.form)==resp.name[1]]
          if (sum(is.na(time))>0) stop("At least one value is missing in time variable.")
       }else{
          if (type.pheno==2){
             pheno<-model.response(model.frame(formula,data=data.form,na.action=NULL), "any")
             if (is.factor(pheno)) pheno<-as.character(levels(pheno))[pheno]
             if (length(unique(pheno))<5) {
                 warning("Response variable takes few than 5 different values for a linear model.\n")
             }
             time<-pheno
          }
       }
       
    }

 
    nindiv<-nrow(data.form)
    nlocus<-length(snps.name)
   
    #covariates
    if (length(covars.name)>0){
        covars<-as.data.frame(data.form[,names(data.form)%in%covars.name])
        names(covars)<-covars.name
        covars.val<-unlist(apply(covars,2,FUN=function(o) length(unique(o))))
        if (sum(covars.val==1)>0) stop("At least one covariate has all its values identical. Please, remove it.")
        if (sum(is.na(covars))>0) stop("At least one covariate has missing values.") 
        covars.cat<-covars[,covars.val<6]
        covars.cat<-as.data.frame(covars.cat)
        names(covars.cat)<-names(covars)[covars.val<6]
        covars.cat.name<-names(covars.cat)
     
        covars.no.cat<-covars[,!(names(covars)%in%covars.cat.name)]
        covars.no.cat<-as.data.frame(covars.no.cat)
        names(covars.no.cat)<-names(covars)[covars.val>5]
        covars.no.cat.name<-names(covars.no.cat)

        if (ncol(covars.cat)!=0){
            dummies.list<-apply(covars.cat,2,make.dummies)
            dummies.mat<-dummies.list[[1]]
            i<-2
            while (i<=length(dummies.list)) {
           dummies.mat<-cbind(dummies.mat,dummies.list[[i]])
          i<-i+1
        }
        dummies<-as.data.frame(dummies.mat,nrow=nindiv)
            i<-1
            n<-NULL
            for (i in 1:length(dummies.list)){       
                 n<-c(n,paste(names(dummies.list[i]),names(dummies.list[[i]]),sep="."))       
            }
            names(dummies)<-n
            covars.c<-cbind(covars.no.cat,dummies)
         
       }else{
            covars.c<-covars
       }
       ncovars<-ncol(covars)
       ncovars.c<-ncol(covars.c)
       covars.name<-c(covars.no.cat.name,covars.cat.name)
       covars.c.name<-names(covars.c)
       covar<-TRUE    
    }else{
       covar<-FALSE  
       ncovars<-0
       ncovars.c<-0
       covars.c<-"a"
       covars.c.name<-names(covars.c)
    }
       
    
     #interactions
    
    inters.name<-attr(terms(formula),"term.labels")[attr(terms(formula),"order")==2]
    if (length(inters.name)>0) {
       list.names<-strsplit(colnames(attr(terms(formula),"factors")),split=":")
       list.names.inter<-lapply(list.names,length)==2
       num.inters<-sum(list.names.inter)
       if (sum(colnames(attr(terms(formula),"factors"))[colnames(attr(terms(formula),"factors"))!="haplotypes"][!list.names.inter]%in%unlist(list.names[list.names.inter]))!=num.inters) stop("All interaction covariates must be main effect variables. \n")
 
       inters<-length(inters.name)
       u<-unlist(strsplit(inters.name,split=":",fixed=TRUE))
       vars.int<-u[u!="haplotypes"]
       intcovars<-NULL
       i<-1
       for (i in 1:length(covars.name)){
            if (covars.name[i]%in%vars.int){
                if (covars.val[names(covars.val)==covars.name[i]]>5){
                    intcovars<-c(intcovars,1)
                }else{
                    intcovars<-c(intcovars,rep(1,(covars.val[names(covars.val)==covars.name[i]])-1))
                }                
            }else{
                if (covars.val[names(covars.val)==covars.name[i]]>5){
                    intcovars<-c(intcovars,0)
                }else{
                    intcovars<-c(intcovars,rep(0,(covars.val[names(covars.val)==covars.name[i]])-1))
                }
            }
       }
      
    }else{

      intcovars<-0
      vars.int<-NULL
    }
    inters<-sum(intcovars)

    #priors

    if (!is.null(dim(prior.val))){
        prior<-1   
        mu<-prior.val$mean.haplo
        varp<-prior.val$sd.haplo
        mu<-as.double(levels(mu))[mu]
        varp<-as.double(levels(varp))[varp]
        mu.interc<-mu[prior.val$Freq==max(prior.val$Freq)]
        varp.interc<-varp[prior.val$Freq==max(prior.val$Freq)]
        mu<-c(mu.interc,mu[prior.val$Freq!=max(prior.val$Freq)])
        varp<-c(varp.interc,varp[prior.val$Freq!=max(prior.val$Freq)])
        mu[is.na(mu)]<--9
        varp[is.na(varp)]<--9
    }else{
        prior<-0
        mu<-rep(-9,2^nlocus)
        varp<-rep(-9,2^nlocus)
    }
  
   c(data,data.form,geno,pheno,nindiv,nlocus,covars.c,covar,ncovars,covars.name,inters.name,intcovars,time)
    
    #preparing genotypes 

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
    if ((verbose==1)||(verbose==2)) cat("Computing haplotype pairs... \n")

    #codifying haplotypes


    pairs<-haplo.pairs(mat)
    n.pairs.indiv<-pairs[[1]]    
    pairs<-as.numeric(unlist(strsplit(unlist(pairs[2]),split="-")))
    pairs.m<-matrix(pairs,byrow=TRUE,ncol=2)
    codi.haplo.cromo1<-pairs.m[,1]
    codi.haplo.cromo2<-pairs.m[,2] 

    genotypes.vec<-c(t(mat))

    genotypes.vec<-as.numeric(genotypes.vec)

   
    if (covar==TRUE){

       covariates.vec<-c(t(covars.c))

       covariates.vec<-as.numeric(covariates.vec)

       covarsm<-as.matrix(covars.c)
       ncovars<-ncovars.c

      

    }

    if (covar==FALSE){

       covariates.vec<-0

       covars.name<-"a"

    }


    varcont<-var(time)

    

 id.outhaplo<-paste(gsub("[./: ]","",format(Sys.time(), "%d/%m/%y %H:%M:%OS3")),paste(sample(LETTERS)    
 [1:10],collapse=""),sep="")

 #id.outhaplo<-unlist(strsplit(id.outhaplo,""))

 #id.outhaplo<-as.numeric(gsub("[./: ]","",format(Sys.time(), "%d/%m/%y %H:%M:%OS3")))
 # call to C 
  
    res <- .C("mcmc",numindiv=as.integer(nindiv),

               numindivexp=as.integer(numindivexp),

               weights=as.double(weights),

               numlocus=as.integer(nlocus),

               burn.in=as.integer(burn.in),

               burn.in.haplo=as.integer(burn.in.haplo),

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

               onlyfreq=as.integer(0),

               id_outhaplo=as.character(id.outhaplo)               

               )
      

        if (res$mem_error==1) stop("Memory error due to number of subjects.")

        if (res$mem_error==2) stop("Memory error due to number of locus.") 
        
        #id.outhaplo<-paste(id.outhaplo,collapse="")

        mat.iters<-import(paste("outhaplo",id.outhaplo,".txt",sep=""))
        unlink("a.txt")
        unlink("b.txt")
        if (file==FALSE)  unlink(paste("outhaplo",id.outhaplo,".txt",sep="")) 
        res$llista_haplos<-haplo.list.all(nlocus)
        res.end<-list(res,sign,alleles_mat,snps.name,covars.c.name,vars.int,formula,family[1],mat.iters)
        class(res.end)<-"reg"


return(res.end)
}

