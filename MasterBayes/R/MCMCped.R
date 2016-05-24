"MCMCped" <-
function(PdP=PdataPed(),
         GdP=GdataPed(),
         sP=startPed(),
         tP=tunePed(),
         pP=priorPed(),
         mm.tol=999,
         nitt=13000,
         thin=10,
         burnin=3000,
         write_postG=FALSE,
         write_postA=FALSE,  
         write_postP="MARGINAL",
         checkP = FALSE,
         jointP = TRUE,
         DSapprox = FALSE,
         verbose=TRUE
){

    if(is.null(PdP$id)==FALSE){                    # if phenotypic data exists then use id's from PdP as reference
      unique_id<-as.character(unique(PdP$id))     
      PdP$id<-match(PdP$id, unique_id)             # convert phenotypic id's to numeric
      if(length(PdP$USsire)>1){
        if(length(PdP$USsire)==length(PdP$offspring)){
          PdP$USsire<-as.character(PdP$USsire[which(PdP$offspring==1)])
        }else{
          stop("USsire in PdP should either be TRUE/FALSE or the same length as id")
        }
      }
      if(length(PdP$USdam)>1){
        if(length(PdP$USdam)==length(PdP$offspring)){
          PdP$USdam<-as.character(PdP$USdam[which(PdP$offspring==1)])
        }else{
          stop("USdam in PdP should either be TRUE/FALSE or the same length as id")
        }
      }
    }else{             
      unique_id<-as.character(unique(GdP$id))
    }    

    nind<-length(unique_id)                        # number of sampled individuals

############################### re-arrange genetic data to match phenotype data ####################################
 
    if(is.null(GdP$G)==FALSE){        
      if(is.null(sP$A)==FALSE & (GdP$marker.type=="MSC" || GdP$marker.type=="MSW")){
        for(i in 1:length(GdP$G)){
          sP$A[[i]]<-sP$A[[i]][order(sP$A[[i]], decreasing=T)]
          if(any((allele.names(GdP$G[[i]])%in%names(sP$A[[i]]))==FALSE)){
            stop("some alleles in GdP$G are not in sP$A")
          }
          GdP$G[[i]]<-genotype(GdP$G[[i]], alleles=names(sP$A[[i]]), reorder="no")
        }
      }
      grouped_by_id<-order(match(GdP$id, unique_id))        
      GdP$id<-GdP$id[grouped_by_id]                                   # reorder genotype id's
      GdP$id<-match(GdP$id, unique_id)                                # convert genetic id's to numeric
      GdP$G<-lapply(GdP$G, function(x){x[grouped_by_id]})             # reorder genotype records
      GdP$categories<-GdP$categories[grouped_by_id]                   # reorder error catgeories
      GdP$categories<-match(GdP$categories, unique(GdP$categories))   # convert error cat's to numeric
    }

    nsamp<-length(GdP$id)
    maxrep<-max(0,table(as.numeric(GdP$id)))
    ncat<-max(1,length(unique(GdP$categories)))

############################### build design matrices ############################################################

    if(is.null(GdP$G)==FALSE){    
       if(sP$estG & GdP$marker.type=="MSW"){
          GdP$marker.type<-"MSC" # speed calculation up is genotypes are going to be calculated anyway 
          X.list<-getXlist(PdP, GdP, A=sP$A, E1=sP$E1, E2=sP$E2, mm.tol=mm.tol) 
          GdP$marker.type<-"MSW"
       }else{
          X.list<-getXlist(PdP, GdP, A=sP$A, E1=sP$E1, E2=sP$E2, mm.tol=mm.tol)
       }
    }else{
      X.list<-getXlist(PdP, GdP, A=sP$A, E1=sP$E1, E2=sP$E2, mm.tol=mm.tol) 
    }
    
 
    noff<-length(X.list$X)	
    ndam<-c(unlist(lapply(X.list$X,function(x){length(x$restdam.id)})))	
    nsire<-c(unlist(lapply(X.list$X,function(x){length(x$restsire.id)})))
    ntdam<-c(unlist(lapply(X.list$X,function(x){length(x$dam.id)})))	
    ntsire<-c(unlist(lapply(X.list$X,function(x){length(x$sire.id)})))

	
    nbeta<-c(ncol(X.list$X[[1]]$XDus), 
             ncol(X.list$X[[1]]$XDs),
             ncol(X.list$X[[1]]$XSus),
             ncol(X.list$X[[1]]$XSs),
             ncol(X.list$X[[1]]$XDSus),
             ncol(X.list$X[[1]]$XDSs))

    nlinked<-length(c(grep("linked", colnames(X.list$X[[1]]$XDus)), grep("linked", colnames(X.list$X[[1]]$XDs))))

    # number of linked parameters

    if(length(nbeta)==0){nbeta<-rep(0,6)}

    off_id<-as.numeric(names(X.list$X))	                                
    dam_id<-c(unlist(lapply(X.list$X,function(x){x$restdam.id})))
    sire_id<-c(unlist(lapply(X.list$X,function(x){x$restsire.id})))

################################ make sure evrything is OK ####################################################

  sP<-ErrorCheck(sP, tP, pP, PdP=PdP, GdP=GdP, unique_id=unique_id, nbeta=nbeta, nlinked=nlinked)

  if(sP$estG==FALSE & write_postG==TRUE){write_postG<-FALSE}

################################ get starting/tuning/prior parameterisation if not given ############################

  sPtP<-getsPandtP(sP, tP, PdP=PdP, GdP=GdP, X.list=X.list, nbeta=nbeta, unique_id=unique_id, checkP=checkP)
  sP<-sPtP$sP
  tP<-sPtP$tP

################################ write some C++ stuff #########################################################

  nloci<-length(sP$A)                      # number of loci
  l_name<-names(sP$A)                      # loci names
  if(nloci!=0){
    nall<-unlist(lapply(GdP$G, nallele))    # number of alleles per locus
  }else{
    nall<-0
  }
  maxall<-max(0,nall)                      # number of alleles at the most polymorhic locus

######################### empty vectors to which posterior distributions are written to #######################

  post<-list(beta=NULL, USdam=NULL, USsire=NULL, E1=NULL, E2=NULL, G=NULL, A=NULL, P=NULL)
  post<-getPost(post, sP, X.list, nitt, thin, burnin, write_postG, write_postP, write_postA, unique_id, GdP$marker.type)

############################## Format for C++ ########################################################

    if(is.null(PdP$USdam)==FALSE){
      if(length(PdP$USdam)==1){
        if(PdP$USdam==TRUE){
          PdP$USdam<-rep(1, sum(PdP$offspring))
          nusd<-1
        }else{
          PdP$USdam<--999
          nusd<-0
        }
      }else{
        PdP$USdam<-match(PdP$USdam, unique(PdP$USdam))
        nusd<-length(unique(PdP$USdam))
      }
    }else{
      nusd<-0
      PdP$USdam<--999
    }

   if(is.null(PdP$USsire)==FALSE){
    if(length(PdP$USsire)==1){
      if(PdP$USsire==TRUE){
        PdP$USsire<-rep(1, sum(PdP$offspring))
        nuss<-1
      }else{
        PdP$USsire<--999
        nuss<-0
      }
    }else{
      PdP$USsire<-match(PdP$USsire, unique(PdP$USsire))
      nuss<-length(unique(PdP$USsire))
    }
   }else{
     nuss<-0
     PdP$USsire<--999
   }

   X_design_betaDus<-c(unlist(lapply(X.list$X,function(x){(x$XDus)})))	
   X_design_betaSus<-c(unlist(lapply(X.list$X,function(x){(x$XSus)})))
   X_design_betaDSus<-c(unlist(lapply(X.list$X,function(x){(x$XDSus)})))
   X_design_betaDs<-c(unlist(lapply(X.list$X,function(x){(x$XDs)})))
   X_design_betaSs<-c(unlist(lapply(X.list$X,function(x){(x$XSs)})))
   X_design_betaDSs<-c(unlist(lapply(X.list$X,function(x){(x$XDSs)})))

   if(nbeta[1]>0){
     X_design_betaDus[which(is.na(X_design_betaDus)==TRUE)]<-0
   }
   if(nbeta[3]>0){
     X_design_betaSus[which(is.na(X_design_betaSus)==TRUE)]<-0
   }
   if(nbeta[5]>0){
   X_design_betaDSus[which(is.na(X_design_betaDSus)==TRUE)]<-0
   }

   if(sP$estG==TRUE || sP$estE1==TRUE || sP$estE2==TRUE || sP$estA==TRUE){            # if geneotype data is present convert to c++ format
     GdP$G<-GtoC(GdP$G, biallelic=(GdP$marker.type!="MSC" & GdP$marker.type!="MSW"))
   }else{
     GdP$G<-0  
     GdP$id<-0
     GdP$categories<-0
   }

#return(list(G=sP$G, ped=sP$ped))

   if(length(sP$G)!=0){              # if sP$G is specified then convert to c++ format
       sP$G<-GtoC(sP$G, biallelic=(GdP$marker.type!="MSC" & GdP$marker.type!="MSW"))
    }else{
      sP$estG<-0  
      sP$estA<-0
      sP$estE1<-0
      sP$estE2<-0
   }
    	
   sP$ped[,2][which(is.na(sP$ped[,2])==TRUE)]<-nind+1
   sP$ped[,3][which(is.na(sP$ped[,3])==TRUE)]<-nind+1+nusd
   tPus<-NULL

   if(((nusd+nuss)*(sP$estUSdam+(sP$estUSsire==TRUE | sP$estUSsire=="USdam")))>0){
     tPus<-sqrt(diag(nusd+nuss)*c(tP$USdam, tP$USsire))
   }

   if(is.null(sP$USdam) & is.null(sP$USsire)){sP$USdam<--999}
   pPUSmu<-c(pP$USdam$mu, pP$USsire$mu)
   pPUSsigma<-c(pP$USdam$sigma, pP$USsire$sigma)

   if(FALSE%in%(as.integer(pPUSmu)%in%999) & TRUE%in%(as.integer(pPUSmu)%in%999)){    
     pPUSsigma<-pPUSsigma[-which(as.integer(pPUSmu)==999)]
     pPUSmu<-pPUSmu[-which(as.integer(pPUSmu)==999)]
   }

   if(pP$beta$mu[1]!=999 & length(pP$beta$mu)>1){
     if(length(pP$beta$mu)!=length(unique(X.list$beta_map))){stop("pP$beta$mu is the wrong length")}
     if(dim(pP$beta$sigma)[1]!=length(unique(X.list$beta_map))){stop("pP$beta$sigma is the wrong dimension")}
     if(dim(pP$beta$sigma)[2]!=length(unique(X.list$beta_map))){stop("pP$beta$sigma is the wrong dimension")}
   }

# get linked parameters

  if(sum(nbeta)>0){
    beta_map<-X.list$beta_map-1
    nunique_beta<-length(unique(beta_map))
  }else{
    nunique_beta<-0
    beta_map<--999
  }  

  mtype.numeric<-sum(c("MSC", "AFLP", "MSW", "SNP")%in%GdP$marker.type*c(1:4))

  rel.mate<-unlist(lapply(PdP$formula, function(x){length(grep("relational.*MATE",x))==TRUE}))  # another hack - relational=MATE variables

  DSapprox<-as.numeric(DSapprox)

  if(any(rel.mate)){
    if(DSapprox==1){
      if(length(grep("Female", PdP$formula[[match(TRUE, rel.mate)]]))>0){
        DSapprox<-2   # mate=male
      }else{
        DSapprox<-1   # mate=female
      }
    }
  }else{
    DSapprox<-0    # no mate variables 
  }

  estimating<-c(sP$estP,sP$estG,sP$estA,sP$estE1, sP$estE2, sP$estbeta, sP$estUSdam==TRUE, sP$estUSsire==TRUE, GdP$perlocus, mtype.numeric, sP$estUSsire=="USdam", checkP, jointP,DSapprox)

  store_post<-c(write_postG,write_postA,write_postP=="JOINT",verbose)

 Merge4C<-function(X.list){
 Merge<-matrix(unlist(lapply(X.list$X, function(x){x$mergeN})), length(X.list$X[[1]]$mergeN), length(X.list$X))
 Mmat<-list()
 for(i in 1:(length(X.list$X[[1]]$mergeN)/2)){
     Mmat[[i]]<-Merge[(2*i-1):(2*i),]
 } 
# Mmat<-unlist(lapply(Mmat, t))
 unlist(Mmat)
}

 if(length(X.list$merge)>0){
   MergeN<-Merge4C(X.list)
   MergeV<-X.list$merge-1
   MergeUS<-X.list$mergeUS-1
   nMerge<-length(X.list$mergeUS)
 }else{
   MergeN<--999
   MergeV<--999
   MergeUS<--999
   nMerge<-0
 }
################# Error Check ##########################################################################

output<-.C("MCMCped",
        as.integer(c(nind, noff)),		# number of individuals sampled
        as.integer(ndam), 	
	as.integer(nsire), 
        as.integer(ntdam), 	
	as.integer(ntsire), 	
        as.integer(nsamp),              # number of samples 
        as.integer(nloci),		# number of loci
        as.integer(nall),               # number of alleles per locus
        as.integer(maxall),             # number of alleles at most polymorphic locus
        as.integer(maxrep),             # maximum number of repeat samples per individual
        as.integer(ncat),               # number of categories 
        as.integer(nusd),
        as.integer(nuss),
	as.integer(nbeta), 
        as.integer(beta_map),
        as.integer(nunique_beta),
        as.double(MergeN),
        as.integer(MergeV),
        as.integer(MergeUS),
        as.integer(nMerge),
 	as.integer(nitt),		# number of itterations
	as.integer(thin),		# thinning interval
	as.integer(burnin),     	# burn in
        as.integer(GdP$id-1),           # numeric id relating samples to individuals
	as.integer(GdP$G),              # observed genotypes	
	as.integer(off_id-1), 
	as.integer(dam_id-1),
        as.integer(sire_id-1),
        as.double(X_design_betaDus),         # design matrices for dam variables
        as.double(X_design_betaSus),         # design matrices for sire variables 
        as.double(X_design_betaDSus),        # design matrices for dam:sire interactions
	as.double(X_design_betaDs),          # design matrices for dam variables
        as.double(X_design_betaSs),          # design matrices for sire variables 
        as.double(X_design_betaDSs),         # design matrices for dam:sire interactions/
	as.double(unlist(sP$A)),	     # starting allele frequencies
	as.double(sP$E1),	             # starting values of E1 and E2
	as.double(sP$E2),	             # starting values of E1 and E2
	as.double(c(sP$beta)),	             # starting vector of beta
	as.double(c(sP$USdam, sP$USsire)),   # starting vector of US numbers
        as.integer(sP$G),                    # starting true genotypes  
	as.integer(as.numeric(sP$ped[,2])-1),   # starting vector of dams
	as.integer(as.numeric(sP$ped[,3])-1),   # starting vector of sires
	as.double(post$A),	                # posterior distribution of allele frequencies
	as.double(post$E1),	                # posterior distribution of E1
	as.double(post$E2),	                # posterior distribution of E2
	as.double(post$beta),	                # posterior distribution of beta
	as.double(c(post$USdam, post$USsire)),	# posterior distribution of USs'
        as.integer(unlist(post$G)),     # posterior distribution of true genotypes
	as.integer(post$P),	        # posterior distribution of Pedigree
	as.double(c(pP$E1)),            # prior distribution for E1 across categories; beta(a,b)
	as.double(c(pP$E2)),            # prior distribution for E2 across categories; beta(a,b)
	as.double(c(pP$beta$mu)),       # prior distribution for MVN(mu, sigma)
	as.double(c(solve(pP$beta$sigma))),
        as.double(c(sum(log(eigen(pP$beta$sigma)$values)))),
	as.double(pPUSmu),              # prior distribution for MVN(mu, sigma)
	as.double(pPUSsigma),
        as.double(c(tP$E1)),                  # standard deviation of normal candidate generating function for E1
        as.double(c(tP$E2)),                  # standard deviation of normal candidate generating function for E1
        as.double(c(tP$beta)),                # standard deviation of normal candidate generating function for beta
        as.double(tPus),                      # standard deviation for proposal US
        as.integer(GdP$categories-1),
        as.integer(PdP$USdam-1),
        as.integer(PdP$USsire-1),
        as.integer(estimating),
        as.integer(store_post))


  post$A<-output[[43]]
  post$E1<-output[[44]]
  post$E2<-output[[45]]
  post$beta<-output[[46]]

  if(sP$estUSdam==TRUE | sP$estUSsire==TRUE){
    UStmp<-t(matrix(output[[47]], (sP$estUSsire==TRUE)*nuss+sP$estUSdam*nusd, ceiling((nitt-burnin)/thin)))
    if(sP$estUSdam==TRUE){
      post$USdam<-mcmc(as.matrix(UStmp[,1:nusd]))
      if(nusd>1){
        colnames(post$USdam)<-paste("USdam.", unique(PdP$USdam), sep="")
      }else{
        colnames(post$USdam)<-"USdam"
      }
    }
    if(sP$estUSsire==TRUE | sP$estUSsire=="USdam"){
      post$USsire<-mcmc(as.matrix(UStmp[,(1:nuss)+nusd*sP$estUSdam]))
      if(nuss>1){
        colnames(post$USsire)<-paste("USsire.", unique(PdP$USsire), sep="")
      }else{
        colnames(post$USsire)<-"USsire"
      }
    }
  }

  post$G<-output[[48]]
  post$P<-output[[49]]

  getPost(post, sP, X.list, nitt, thin, burnin, write_postG, write_postP, write_postA, unique_id, GdP$marker.type)

}

