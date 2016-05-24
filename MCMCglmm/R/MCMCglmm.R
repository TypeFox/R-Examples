"MCMCglmm"<-function(fixed, random=NULL, rcov=~units, family="gaussian", mev=NULL, data, start=NULL, prior=NULL, tune=NULL, pedigree=NULL, nodes="ALL",scale=TRUE, nitt=13000, thin=10, burnin=3000, pr=FALSE, pl=FALSE, verbose=TRUE, DIC=TRUE, singular.ok=FALSE, saveX=TRUE, saveZ=TRUE, saveXL=TRUE, slice=FALSE, ginverse=NULL){

    orig.na.action<-options("na.action")[[1]]
    options("na.action"="na.pass")	

    if(missing(data)){stop("data argument is missing")}
    if(missing(fixed)){stop("fixed is missing")}
    if(class(fixed)!="formula"){stop("fixed should be a formula")}
    if(class(rcov)!="formula"){stop("rcov should be a formula")}
    if(class(random)!="formula" & class(random)!="NULL"){stop("random should be a formula")}

    reserved.names<-c("units", "MCMC_y", "MCMC_y.additional","MCMC_liab","MCMC_meta", "MCMC_mev", "MCMC_family.names", "MCMC_error.term", "MCMC_dummy")
    family.types<-c("gaussian", "poisson", "multinomial", "notyet_weibull", "exponential", "cengaussian", "cenpoisson", "notyet_cenweibull", "cenexponential",  "notyet_zigaussian", "zipoisson", "notyet_ziweibull", "notyet_ziexponential", "ordinal", "hupoisson", "ztpoisson", "geometric", "zapoisson", "zibinomial", "threshold")

    if(any(names(data)%in%reserved.names)){
      stop(paste(names(data)[which(names(data)%in%reserved.names)], " is a reserved variable please rename it"))
    }

    covu<-0

    if(is.null(prior)==FALSE & any(names(prior)%in%c("R", "G", "B")==FALSE)){stop("prior list should contain elements R, G, and/or B only")}

    if(!is.null(prior)){ # reformat old-style prior 
      if(!is.list(prior$R[[1]])){
        prior$R<-list(R1=prior$R)
      }
    }
    if(!is.null(prior$R[[1]]$covu)){
      if(prior$R[[1]]$covu){
        covu<-1
        prior$G<-c(prior$G, list(prior$R[[1]]))
        prior$G[[length(prior$G)]]$covu<-NULL
        prior$G[[length(prior$G)]]$fix<-NULL
      }
    }
    if(!is.null(start$R)){
      if(!is.list(start$R)){
        start$R<-list(R1=start$R)
      }
      if(covu!=0){
        start$G<-c(start$G,list(start$R[[1]]))
      }
    }

    if(is.null(prior)==FALSE & any(names(prior)%in%c("R", "G", "B")==FALSE)){stop("prior list should contain elements R, G, and/or B only")}
    if(((is.null(start$G) & is.null(random)==FALSE) & is.null(start$R)==FALSE) | (is.null(start$R) & is.null(start$G)==FALSE)){stop("need both or neither starting R and G structures")}
    if(((is.null(prior$G) & is.null(random)==FALSE) & is.null(prior$R)==FALSE) | (is.null(prior$R) & is.null(prior$G)==FALSE)){stop("either both or neither R and G structures need a prior")}

    if(is.null(start$QUASI)){
      QUASI=TRUE
    }else{
      QUASI<-start$QUASI
      if(is.logical(QUASI)==FALSE){stop("start$QUASI should be TRUE or FALSE")}
    }

    original.fixed<-fixed                                                                            # original model specification
    original.random<-random 
    original.rcov<-rcov
    original.family<-family

    response.names<-names(get_all_vars(as.formula(paste(as.character(fixed)[2], "~1")), data))       # response variable names 
    data[,response.names]<-model.frame(as.formula(paste(as.character(fixed)[2], "~1")), data)[[1]]

#########################################################################
# for ginverse analyses form A and augment with missing nodes if needed #
#########################################################################

	if(is.null(pedigree)==FALSE){  # back compatibility if pedigree is directly passed
           if(is.null(ginverse)){
             ginverse<-list(animal=inverseA(pedigree=pedigree, scale=scale, nodes=nodes)$Ainv)
           }else{
             if("animal"%in%names(ginverse)){stop("animal ginverse appears but pedigree has been passed")}
             ginverse$animal<-inverseA(pedigree=pedigree, scale=scale, nodes=nodes)$Ainv
           }
        }

	if(is.null(ginverse)==FALSE){
	  for(i in 1:length(ginverse)){           
            if(is.null(rownames(ginverse[[i]]))){stop(paste(names(ginverse)[i], "ginverse must have non-null rownames"))}
            if(names(ginverse)[i]%in%names(data)==FALSE){stop(paste(names(ginverse)[i], "does not appear in data"))}
            if(any(na.omit(unique(data[,names(ginverse)[i]]))%in%rownames(ginverse[[i]])==FALSE)){stop(paste("some levels of", names(ginverse)[i], "do not have a row entry in ginverse"))}
            if(any(duplicated(rownames(ginverse[[i]])))){stop(paste("rownames of ", names(ginverse)[i], "ginverse must be unique"))}
            if(determinant(ginverse[[i]])$sign<=0){stop(paste(names(ginverse)[i], "ginverse is not positive definite"))}
	    data[,names(ginverse)[i]]<-factor(data[,names(ginverse)[i]], levels=rownames(ginverse[[i]]))                        
          }
	}

    MVasUV=FALSE  # Multivaraite as Univariate
    if(is.null(family)){
       if(is.null(data$family)){
         stop("no family specified")
       }else{
         if(is.null(data$trait)){
           stop("data.frame needs to have a column indexing traits if fitting multivaraite models as univariate")
         }else{
           if(is.factor(data$trait)==FALSE){
             stop("trait must be a factor")
           }
           if(any(tapply(data$family, data$trait, function(x){length(unique(x))})!=1)){
             stop("all data from the same trait must come from the same distribution")
           }
           rterm.family<-data$family[match(levels(data$trait), data$trait)]
         }           
         if(length(unique(data$trait))==1){
           family.names<-as.character(data$family[1])
         }else{     
           family.names<-as.character(data$family)   
           MVasUV=TRUE
         }
         if(length(grep("cen|multinomial|zi|hu|za", family.names))>0){ 
           stop("For setting up multi-trait models as univariate the responses cannot come from distributions that require more than one data column or have more than one liability (i.e. censored, multinomial, zero-inflated, categorical with k>2): set it up as multivariate using cbind(...)")
         }
       }
    }else{
       if(is.null(data$family)==FALSE){
          stop("family column exists in data and in family argument to MCMCglmm: specify family=NULL")
       }else{
          family.names<-family 
       }                                                                            # response distribution
    }

    # zero-infalted and multinomial>2 need to be preserved in the same R-structure even if idh 

    diagR<-1  # should the residual structure be diagonal even if us is used?

    if(length(grep("hu|zi|za|multinomial[3-9]|[1-9][0-9]", family.names))>0){ 
      if(length(grep("idh\\(trait|us\\(trait|trait:units|units:trait|corg\\(trait|corgh\\(trait|cors\\(trait|idv\\(trait|ante.*\\(trait|sub\\(trait", rcov))==0){
        stop("For error structures involving multinomial data with more than 2 categories, or zero-infalted/altered/hurdle models")
      }else{
        if(length(grep("idh\\(", rcov))>0){
          diagR<-2
          rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
        }
        if(length(grep("trait:units|units:trait", rcov))>0){
          rcov=~us(trait):units
          diagR<-3
        }
      }	  
    }

    if(any(grepl("path\\(", fixed))){ 
      if(length(grep("idh\\(", rcov))>0){
        diagR<-2
        rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
      } 
    }

	
    nS<-dim(data)[1]                               # number of subjects
    y.additional<-matrix(NA, nS,0)                 # matrix equal in dimension to y holding the additional parameters of the distribution (n, upper interval etc.)
    nt<-1                                          # number of traits (to be iterated because y may change dimension with multinomial/categorical/censored distributions)

    if(MVasUV){
      mfac<-rep(0, nlevels(data$trait))
      y.additional<-matrix(NA, nS,1)      # matrix equal in dimension to y holding the additional parameters of the distribution (n, upper interval etc.)
      if(any(family.names=="categorical")){
        ncat<-tapply(data[,response.names][which(family.names=="categorical")], as.character(data$trait[which(family.names=="categorical")]), function(x){length(unique(x))})
        for(i in 1:length(ncat)){
           trait.i<-unique(data$trait[which(family.names=="categorical")])[i]
           data[,response.names][which(data$trait==trait.i)]<-as.numeric(as.factor(data[,response.names][which(data$trait==trait.i)]))-1
        }
        if(any(ncat>2)){
           stop("For setting up multi-trait models as univariate the responses cannot come from distributions that require more than one data column or have more than one liability (i.e. censored, multinomial, zero-inflated, categorical with k>2): set it up as multivariate using cbind(...)")
         }
        y.additional[which(family.names=="categorical")]<-1
        family.names[which(family.names=="categorical")]<-"multinomial"
      }
      if(any(family.names=="ordinal" |  family.names=="threshold")){
        ordinal.traits<-unique(as.numeric(data$trait)[which(family.names=="ordinal" |  family.names=="threshold")])
        for(i in 1:length(ordinal.traits)){
           trait.i<-levels(data$trait)[ordinal.traits[i]]
           data[,response.names][which(data$trait==trait.i)]<-as.numeric(as.factor(data[,response.names][which(data$trait==trait.i)]))
        }
        if(length(ordinal.traits)>2){slice<-FALSE}
        mfac[ordinal.traits]<-0:(length(ordinal.traits)-1)
        ncutpoints<-tapply(data[,response.names][which(family.names=="ordinal" | family.names=="threshold")], data$trait[which(family.names=="ordinal" |  family.names=="threshold")], function(x){length(unique(na.omit(x)))})+1
        stcutpoints<-tapply(data[,response.names][which(family.names=="ordinal" | family.names=="threshold")], data$trait[which(family.names=="ordinal" |  family.names=="threshold")], function(x){qnorm(cumsum(c(0, prop.table(table(x)))))})
        ordinal.names<-levels(data$trait)[ordinal.traits]
        if(any(is.na(ncutpoints))){
          ordinal.names<-ordinal.names[-which(is.na(ncutpoints))]
          stcutpoints<-stcutpoints[-which(is.na(ncutpoints))]
          ncutpoints<-ncutpoints[-which(is.na(ncutpoints))]
        }
      }else{
        ncutpoints<-c()
      }
    }else{
      if(any(names(data)=="trait")){
        stop("trait is a reserved variable please remove or rename this column in the data.frame")
      }
      mfac<-c()                 # stores additional information for each R-structure term. For multinomial models it stores the number of k-2 categories
                                # for zero inflated models it indicates whether the R-structure term is for the Poisson part (0) or zero infaltion part (1)
                                # for ordinal responses it indicates whether its the i^th ordinal response. For example R-structure terms for data A:E with 
                                # multinomialA, multinomialA,  multinomialA, multinomialB, multinomialB, ordinalC, zero-inflationD, ordinalE: mfac=c(1,1,0,1,0,1,2) other 
                                # traits all have 0
      ncutpoints<-c()           # number of cutpoints for ordinal variables = k+1
      stcutpoints<-list()         
      ordinal.names<-c()
      rterm.family<-c()

      for(i in 1:length(family)){

        rterm.family<-c(rterm.family, family[i])

        dist.preffix<-substr(family[i],1,2)                  
        if(nt>length(response.names)){stop("family is the wrong length")}
        if(any(dist.preffix%in%c("ce", "mu", "ca", "tr", "zi", "hu", "za"))){

######################
# categorical traits #
######################
				
          if(dist.preffix=="ca"){
            if(all(is.na(data[[response.names[nt]]]))){
              cont<-matrix(NA, nrow(data),1)
              new.names<-paste(response.names[nt], ".", 1, sep="")
            }else{
              cont<-as.matrix(model.matrix(~ as.factor(data[[response.names[nt]]]))[,-1])                            # form new J-1 variable   
              new.names<-paste(response.names[nt], ".", levels(as.factor(data[[response.names[nt]]]))[-1], sep="")   # give new variables names
            }
            colnames(cont)<-new.names  
            nJ<-dim(cont)[2]
            rterm.family[length(rterm.family)]<-paste("multinomial", nJ+1, sep="") 
            if(nJ>1){
              slice<-FALSE
              rterm.family<-c(rterm.family, rep(rterm.family[length(rterm.family)], nJ-1))  
            }  
            if(length(grep("idh\\(trait|us\\(trait|trait:units|units:trait|corg\\(trait|corgh\\(trait|cors\\(trait|idv\\(trait|ante.*\\(trait", rcov))==0 & nJ>1){
             stop("For error structures involving catgeorical data with more than 2 categories pleasue use trait:units or variance.function(trait):units.")
            }else{
              if(length(grep("idh\\(", rcov))>0  & nJ>1){
                diagR<-2
                rcov=as.formula(paste(gsub("idh\\(", "us\\(", rcov), collapse=""))
              }
              if(length(grep("trait:units|units:trait", rcov))>0  & nJ>1){
                rcov=~us(trait):units
                diagR<-3
              }
            }	  
            mfac<-c(mfac, rep(nJ-1,nJ))           
            data<-data[,-which(names(data)==response.names[nt]), drop=FALSE]    # remove original variable
	    row.names(data)<-row.names(cont)
            data<-cbind(data, cont)     # add new variables to data.frame
            ones<-rep(1, length(response.names))
            ones[which(response.names==response.names[nt])]<-nJ
            response.names<-rep(response.names, ones)
            family.names<-rep(family.names, ones)
            family.names[which(response.names==response.names[nt])]<-"multinomial"
            response.names[which(response.names==response.names[nt])]<-new.names
            nt<-nt+nJ
            y.additional<-cbind(y.additional, matrix(1,nS,nJ))
          }

######################
# multinomial traits #
######################

         if(dist.preffix=="mu"){
           nJ<-as.numeric(substr(family[i],12,nchar(family[i])))-1
           if(nJ>1){
             slice<-FALSE
             rterm.family<-c(rterm.family, rep(rterm.family[length(rterm.family)], nJ-1))  
           }                                                                           # number of J-1 categories
 	   if(nJ<1){stop("Multinomial must have at least 2 categories")}	 
           mfac<-c(mfac, rep(nJ-1,nJ))  
           if(all(data[,match(response.names[0:nJ+nt], names(data))]%%1==0, na.rm=T)==FALSE | all(data[,match(response.names[0:nJ+nt], names(data))]>=0, na.rm=T)==FALSE){
             stop("multinomial data must be positive integers")
           }
           y.additional<-cbind(y.additional, matrix(rowSums(data[,match(response.names[0:nJ+nt], names(data))]), nS,nJ))             # get n of the multinomial
           if(any(na.omit(y.additional)>1)){slice<-FALSE}                                                                           # number of J-1 categories
	   if(any(is.na(y.additional[,dim(y.additional)[2]]) & apply(data[,match(response.names[0:nJ+nt], names(data))], 1, function(x){any(is.na(x)==FALSE)}))){
             stop("multinomial responses must be either completely observed or completely missing")
           }	 
           data<-data[,-which(names(data)==response.names[nt+nJ]),drop=FALSE]                                                        # remove first category
           response.names<-response.names[-(nt+nJ)]
           family.names[nt]<-"multinomial"
           ones<-rep(1,length(family.names))
           ones[nt]<-nJ
           family.names<-rep(family.names, ones)
           nt<-nt+nJ
         }

			
###################
# censored traits #
###################

         if(dist.preffix=="ce"){       
           mfac<-c(mfac, 0)  
           if(any(data[,which(names(data)==response.names[nt+1])]<data[,which(names(data)==response.names[nt])], na.rm=T)){stop("for censored traits left censoring point must be less than right censoring point")}
           y.additional<-cbind(y.additional, data[,which(names(data)==response.names[nt+1])])                               # get upper interval
           if(family.names[nt]=="cenpoisson"){
	     if(all(data[,response.names[0:1+nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[0:1+nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be positive integers")}
             data[[response.names[nt]]][which(data[[response.names[nt]]]!=data[[response.names[nt+1]]])]<-data[[response.names[nt]]][which(data[[response.names[nt]]]!=data[[response.names[nt+1]]])]-1
             }
	   if(family.names[nt]=="cenexponential"){
	     if(any(data[,response.names[0:1+nt]]<0, na.rm=T)){stop("Exponential data must be positive")}  
	   }
	   data<-data[,-which(names(data)==response.names[nt+1]),drop=FALSE]                                                # remove upper interval from the response
           response.names<-response.names[-(nt+1)]
           nt<-nt+1
         }

###################
# truncated traits #
###################
		  
	 if(dist.preffix=="tr"){   
           mfac<-c(mfac, 0)  
	   y.additional<-cbind(y.additional, data[,which(names(data)==response.names[nt+1])])                     # get upper interval
	   data<-data[,-which(names(data)==response.names[nt+1]),drop=FALSE]                                      # remove upper interval from the response
	   response.names<-response.names[-(nt+1)]
	   nt<-nt+1
	 }

########################
# zero-infalted traits #
########################

	 if(dist.preffix=="zi" || dist.preffix=="hu" || dist.preffix=="za"){

           rterm.family<-c(rterm.family, rterm.family[length(rterm.family)])
 
           if(grepl("poisson", family[i])){
	     y.additional<-cbind(y.additional, rep(1,nS), rep(0,nS))
	     if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){
               stop("Poisson data must be non-negative integers")
             }
	     cont<-as.matrix(as.numeric(data[,which(names(data)==response.names[nt])]==0))
           }
           if(grepl("binomial", family[i])){
	     y.additional<-cbind(y.additional, rowSums(data[,match(response.names[0:1+nt], names(data))]), rep(0,nS))
             if(all(data[,match(response.names[0:1+nt], names(data))]%%1==0, na.rm=T)==FALSE | all(data[,match(response.names[0:1+nt], names(data))]>=0, na.rm=T)==FALSE){
               stop("binomial data must be non-negative integers")
             } 
	     cont<-as.matrix(as.numeric(data[,which(names(data)==response.names[nt])]==0))
             response.names<-response.names[-(nt+1)]
           }
           colnames(cont)<-paste(dist.preffix, response.names[nt], sep="_")
	   data<-cbind(data, cont)
	   mfac<-c(mfac, c(0,1))
	   ones<-rep(1, length(response.names))
	   ones[which(response.names==response.names[nt])]<-2
	   response.names<-rep(response.names, ones)
	   family.names<-rep(family.names, ones)
	   family.names[which(response.names==response.names[nt])]<-family.names[nt]
	   response.names[which(response.names==response.names[nt])]<-c(response.names[nt], paste(dist.preffix, response.names[nt], sep="_"))
	   nt<-nt+2			  
	 }

###############################################
# gaussian/poisson/exponential/ordinal traits #
###############################################

        }else{
          mfac<-c(mfac, 0)  
          if(family.names[nt]=="poisson"){ 
	    if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){stop("Poisson data must be non-negative integers")}
          }
          if(family.names[nt]=="ztpoisson"){ 
	    if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=1, na.rm=T)==FALSE){stop("Zero-truncated Poisson data must be positive integers")}
          }
	  if(family.names[nt]=="exponential"){ 
	    if(any(data[,response.names[nt]]<0, na.rm=T)){stop("Exponential data must be positive")}
	  }	
	  if(family.names[nt]=="geometric"){ 
	    if(all(data[,response.names[nt]]%%1==0, na.rm=T)==FALSE | all(data[,response.names[nt]]>=0, na.rm=T)==FALSE){stop("Geometric data must be non-negative integers")}
	  }	
	  if(family.names[nt]=="ordinal" | family.names[nt]=="threshold"){
            mfac[length(mfac)]<-length(ncutpoints)  
            data[,response.names[nt]]<-as.numeric(as.factor(data[,response.names[nt]]))
            if(all(is.na(data[,response.names[nt]]))){ 
              ncutpoints<-c(ncutpoints, 3)
              stcutpoints[[length(stcutpoints)+1]]<-c(-Inf, 0, Inf)
            }else{
              if(max(data[,response.names[nt]], na.rm=T)>2){slice<-FALSE}  
              ncutpoints<-c(ncutpoints, max(data[,response.names[nt]], na.rm=T)+1)
              stcutpoints[[length(stcutpoints)+1]]<-qnorm(cumsum(c(0, prop.table(table(data[,response.names[nt]])))))
              ordinal.names<-c(ordinal.names, response.names[nt]) 	
            }    
	  }
          y.additional<-cbind(y.additional,matrix(NA,nS,1))     
          nt<-nt+1
        }
      }	
      nt<-nt-1
    }

    if(sum((family.names%in%family.types)==FALSE)!=0){stop(paste(unique(family.names[which((family.names%in%family.types)==FALSE)]), "not a supported distribution"))}

###**************************************########################

    if(MVasUV){
      data<-reshape(data, varying=response.names, v.names="MCMC_y", direction="long", idvar="units")       # reshape the data into long format 
      data$MCMC_family.names<-family.names
    }else{
      data<-reshape(data, varying=response.names, v.names="MCMC_y", direction="long", timevar="trait", idvar="units")       # reshape the data into long format 
      data$trait<-factor(response.names[data$trait], response.names)
      if(length(response.names)!=length(family.names)){stop("family must have the same length as the number of responses")}
      data$MCMC_family.names<-rep(family.names, each=nS)     
    }
    data$units<-as.factor(data$units)
    data$MCMC_y.additional<-c(y.additional) 

######################################################
# for (random) meta-analysis add weights/model terms #
###################################################### 	

    if(is.null(mev)==FALSE){
      if(any(dim(mev)!=dim(y.additional))){stop("mev has to be the same dimension as y")}
      data$MCMC_mev<-mev	
      data$MCMC_meta<-factor(1:dim(data)[1], levels=1:dim(data)[1])
      if(is.null(random)){
        random = ~us(sqrt(MCMC_mev)):MCMC_meta
        if(is.null(prior$R)==FALSE){
          prior$G<-list(G1=list(V=as.matrix(1), nu=1, fix=1))
        }
      }else{
        random<-as.formula(paste(paste(as.character(random), collapse=""),"+us(sqrt(MCMC_mev)):MCMC_meta", sep=""))
        if(is.null(start$G)==FALSE){
          start$G[[length(start$G)+1]]<-as.matrix(1)
        }
        if(is.null(prior$G)==FALSE){
          prior$G[[length(prior$G)+1]]<-list(V=as.matrix(1), nu=1, fix=1)
        }
      } 
    }

    rmodel.terms<-split.direct.sum(as.character(random)[2])

    ngstructures<-length(rmodel.terms)
    rmodel.terms<-c(rmodel.terms, split.direct.sum(as.character(rcov)[2])) 

    nfl<-c()                                                       # number of fixed levels the random term is structured by
    nrl<-c()                                                       # number of random levels
    nrt<-c()                                                       # number of new terms per original term

    variance.names<-c()
    GR<-list()
    GRprior<-list()
    Aterm<-c()                                                     # 1 if associated with ginverse
    update<-c()                                                    # variance structure type
    trait.ordering<-c()

    if(is.null(start$R)){
      NOstartG=TRUE
    }else{
      NOstartG=FALSE
    }
    if(is.null(prior$R)){
      NOpriorG=TRUE
    }else{
      NOpriorG=FALSE
    }

##########################
# Build Z ################
##########################

    nr<-1

    if(!NOpriorG){
      if(length(prior$G)!=ngstructures){stop("prior$G has the wrong number of structures")}
    }
    if(!NOstartG){
      if(length(start$G)!=ngstructures){stop("start$G has the wrong number of structures")}
    }

    for(r in 1:length(rmodel.terms)){

       if(r==(ngstructures+1)){nG<-nr-1}  # number of (new) G structures
       if(r<=ngstructures){
         Zlist<-buildZ(rmodel.terms[r], data=data, nginverse=names(ginverse))
       }else{
         Zlist<-buildZ(rmodel.terms[r], data=data)
       }

       update<-c(update, Zlist$vtype)

       nfl<-c(nfl,Zlist$nfl)          
       nrl<-c(nrl,Zlist$nrl)
       nrt<-c(nrt,length(Zlist$nrl))

       Aterm<-c(Aterm, Zlist$Aterm)
       variance.names<-c(variance.names, Zlist$vnames)

       if(r<=ngstructures){

         if(covu!=0 & r==ngstructures){
           covu<-Zlist$nfl
           if(!NOpriorG){
             prior$G[[r]]$V<-prior$G[[r]]$V[1:covu, 1:covu]
             prior$G[[r]]$alpha.V<-prior$G[[r]]$alpha.V[1:covu, 1:covu]
             prior$G[[r]]$alpha.mu=prior$G[[r]]$alpha.mu[1:covu]
           }
           if(!NOstartG){
             start$G[[r]]<-start$G[[r]][1:covu, 1:covu]
           }
         }

         GRtmp<-priorformat(if(NOpriorG){NULL}else{prior$G[[r]]},if(NOstartG){NULL}else{start$G[[r]]}, Zlist$nfl, meta=any(grepl("MCMC_meta", rmodel.terms[r])), diagR=0, vtype=Zlist$vtype)

         if(r==1){
           Z<-Zlist$Z
         }else{
           Z<-cBind(Z, Zlist$Z)     
         }
       }else{
         trait.ordering<-c(trait.ordering, Zlist$trait.ordering)
         if(r==(ngstructures+1)){
           ZR<-Zlist$Z
           Zlist$nfl<-Zlist$nfl+covu             
         }else{
           ZR<-cBind(ZR, Zlist$Z)     
         }

         GRtmp<-priorformat(if(NOpriorG){NULL}else{prior$R[[r-ngstructures]]},if(NOstartG){NULL}else{start$R[[r-ngstructures]]}, Zlist$nfl, meta=any(grepl("MCMC_meta", rmodel.terms[r])), diagR=diagR, vtype=Zlist$vtype)

         if(r==(ngstructures+1) & covu>0){
           beta_rr<-t(solve(GRtmp$start[[1]]$start[1:covu, 1:covu, drop=FALSE], GRtmp$start[[1]]$start[1:covu, (covu+1):Zlist$nfl, drop=FALSE]))+1e-16
           if(ngstructures==1){
             ustart<-0
           }else{
             ustart<-sum(nfl[1:(ngstructures-1)]*nrl[1:(ngstructures-1)])
           }
           if((nfl[length(nfl)-1]*nrl[length(nrl)-1])!=(covu*Zlist$nrl)){
             stop("number of levels in G and R structure do match for cov=TRUE")
           }
           Z[,ustart+1:(covu*Zlist$nrl)]<-Z[,ustart+1:(covu*Zlist$nrl)]+Zlist$Z%*%kronecker(beta_rr, Diagonal(Zlist$nrl))
         }
       }

       for(i in 1:length(GRtmp$prior)){
         GRprior[[nr]]<-GRtmp$prior[[i]]  
         GR[[nr]]<-GRtmp$start[[i]]$start   
         nr<-nr+1
       }
     }

   nR<-nr-nG-1  # number of R structures

   if(nR>1 & diagR!=1){stop("sorry - block-diagonal R structures not yet implemented for responses involving multinomial data with more than 2 categories or zero-infalted/altered/hurdle models")}

     missing<-which(colSums(ZR)==0)
     nadded<-length(missing)

     if(nadded>0){
       ZRaug<-Diagonal(ncol(ZR),as.numeric(colSums(ZR)==0))[missing,]
       ZR<-rBind(ZR, ZRaug)
     }

     if(any(colSums(ZR)>1)){stop("R-structure miss-specified: each residual must be unique to a data point")}
     if(any(rowSums(ZR)==0)){stop("R-structure miss-specified: each data point must have a residual")}
     if(any(ZR@x!=1)){stop("R-structure miss-specified: random regressions not permitted")}
 
     data$MCMC_dummy<-rep(0,dim(data)[1])
                                      
     if(nadded>0){
       data<-rbind(data,data[dim(data)[1]+1:nadded,])       
       data$MCMC_dummy[dim(data)[1]-(nadded-1):0]<-1 
       data$MCMC_family.names[dim(data)[1]-(nadded-1):0]<-"gaussian"  
       if(ngstructures!=0){
         Z<-rBind(Z, as(matrix(0,nadded,ncol(Z)), "sparseMatrix"))  
       }
     }
 
     ordering<-ZR@i+1

    data<-data[ordering,]         

    data$MCMC_error.term<-rep(1:sum(nfl[nG+1:nR]), rep(nrl[nG+1:nR], nfl[nG+1:nR]))

    if(ngstructures==0){
      Z<-as(matrix(0,1,0), "sparseMatrix")
    }else{                                                     # rearrange data to match R-structure
      Z<-Z[ordering,]                                                     
    }

    mfac<-mfac[trait.ordering]

    if(any(rterm.family=="threshold" | rterm.family=="ordinal")){

      nthordinal<-cumsum(rterm.family=="threshold" | rterm.family=="ordinal")*(rterm.family=="threshold" | rterm.family=="ordinal") 

      if(any(tapply(data$MCMC_family.names[which(data$MCMC_dummy==0)], data$MCMC_error.term[which(data$MCMC_dummy==0)], function(x){length(unique(x))>1 & any(x=="threshold" | x=="ordinal")}))){
        stop("threshold/ordinal traits must have sperate residual variance from other traits")
        #stop("threshold/ordinal traits can't have common cutpoints because they differ in their number of categories")
      }

      ncutpoints<-ncutpoints[nthordinal[trait.ordering]]

      stcutpoints<-stcutpoints[nthordinal[trait.ordering]]
      stcutpoints<-lapply(stcutpoints, function(x){
         x<-x-x[2]
         x[1]<--1e+64
         x[length(x)]<-1e+64
         x})
      stcutpoints<-unlist(stcutpoints)
      ordinal.names<-ordinal.names[nthordinal[trait.ordering]]
    }

    rterm.family<-rterm.family[trait.ordering]

    if(is.null(tune)){
      AMtune=c(rep(FALSE, nG), rep(TRUE, nR))
      tune<-as.list(1:nR)
      for(i in 1:nR){
        tune[[i]] = diag(nfl[nG+i])
      }
    }else{
      AMtune=rep(FALSE, nR+nG)
      if(nR>1){
        tune<-sapply(diag(tune), as.matrix, simplify=FALSE)
      }else{
        tune<-list(tune)
      }
      for(i in 1:nR){
        tune[[i]]<-as.matrix(tune[[i]])
        if(dim(tune[[i]])[1]!= dim(tune[[i]])[2] |  dim(tune[[i]])[2]!= nfl[nG+i]){stop(paste("proposal distribution ", i, " is the wrong dimension"))}
        if(is.positive.definite(tune[[i]])==FALSE){stop(paste("proposal distribution ", i, " is not positive definite"))}
      }
    }	

############################
# Build Fixed Effect Model #
############################

   fixed.text<-paste(deparse(fixed[[3]]), collapse="")

   if(grepl("path\\(", fixed.text)){
     if(nadded>0){stop("observations in residual blocks must be complete")}
     path.terms<-close.bracket("path\\(", fixed.text)
     path.terms<-as.formula(paste("~", paste(apply(path.terms,1,function(x){substr(fixed.text, x[1], x[2])}), collapse="+"), "-1"))
     fixed.text<-gsub("(\\+|\\-) path\\(.*\\)", "", fixed.text)  # remove path analytic terms
   }else{
     path.terms<-NULL
   }

   fixed<-as.formula(paste("~",fixed.text))

   ffterms<-names(get_all_vars(fixed, data=data))

   if(any(!ffterms%in%names(data))){
    stop(paste(ffterms[which(!ffterms%in%names(data) & !ffterms%in%reserved.names)], "not in data"))
   }

   X<-model.matrix(fixed,data)

   if(nadded>0){
     X[which(data$MCMC_dummy==1),]<-0
   }

   sir<-any(grepl("sir\\(", dimnames(X)[[2]]))

   if(sir & !is.null(path.terms)){"path and sir functions cannot be used together: fit using sir function only"}

   if(sir | !is.null(path.terms)){  # Do structural parameters exist?
     if(sir){
       if(any(family.names!="gaussian")){stop("currently simultaneous/recursive models can only be fitted to Gaussian data unless path() is used")}
       L<-X[,grep("sir\\(", dimnames(X)[[2]]),drop=FALSE]
       L<-as(L, "sparseMatrix")
       if(any(is.na(data$MCMC_y) & rowSums(L)!=0)){
         stop("currently missing observations cannot be augmented if they depend on simultaneous/recursive terms unless path() is used")
       }
       pL<-unique(cumsum(((duplicated(attr(X, "assign"))==FALSE)+(grepl("sir\\(", dimnames(X)[[2]])==FALSE))>0)[grep("sir\\(", dimnames(X)[[2]])])
       # position of structural parameters
       Lnames<-grep("sir\\(", attr(terms(fixed), "term.labels"))
       Lnames<-attr(terms(fixed), "term.labels")[Lnames]
       X<-X[,-grep("sir\\(", dimnames(X)[[2]]),drop=FALSE]
     }else{
       if(nR>1){stop("path models can only have one R-term")}
       if(any(family.names=="threshold")){stop("path models are not implemented for threshold responses")}
       if(sum(ncutpoints)>0){stop("path models are not implemented for ordinal responses with >2 categories")}
       L<-as(model.matrix(path.terms), "sparseMatrix")
       Lnames<-attr(terms(path.terms), "term.labels")
     }
     nL<-dim(L)[2]/dim(L)[1]  # number of structural parameters
     pL<-FALSE
     warning("priors for sir/path parameters not implemented")
   }else{
     L<-Matrix(0,0,0)
     nL<-0
     Lnames=NULL
   }

   if(any(is.na(X))){stop("missing values in the fixed predictors")}

   if(singular.ok==FALSE){
     if(all(is.na(data$MCMC_y))){stop("all data are missing. Use singular.ok=TRUE to sample these effects, but use an informative prior!")}
     sing.rm<-lm(replace(data$MCMC_y, which(data$MCMC_y%in%c(-Inf, Inf)), 0)~X-1, subset=(is.na(data$MCMC_y)==FALSE))
     sing.rm<-which(is.na(sing.rm$coef))
     if(length(sing.rm)){
       warning("some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!")
       X<-X[,-sing.rm]
     }
   }	
   X<-as(X, "sparseMatrix")
   if(any(apply(X,1, function(x){all(x==0)}))){
     X[,1][which(apply(X,1, function(x){all(x==0)}))]<-1e-18 
   }

#####################################
# Check prior for the fixed effects #
#####################################

   if(is.null(prior$B)){
      prior$B=list(V=diag(dim(X)[2])*1e+10, mu=matrix(0,dim(X)[2],1))
      prior$L=list(V=diag(nL)*1e+10, mu=matrix(0,nL,1))
   }else{
      if(is.matrix(prior$B$mu)==FALSE){
	prior$B$mu<-matrix(prior$B$mu,  length(prior$B$mu), 1)
      }
      if(is.matrix(prior$B$V)==FALSE){
	prior$B$V<-as.matrix(prior$B$V)
      }	
      if(is.positive.definite(prior$B$V)==FALSE){stop("fixed effect V prior is not positive definite")}
      if((dim(X)[2]+nL)!=dim(prior$B$mu)[1]){stop("fixed effect mu prior is the wrong dimension")}
      if((dim(X)[2]+nL)!=dim(prior$B$V)[1] | (dim(X)[2]+nL)!=dim(prior$B$V)[2]){stop("fixed effect V prior is the wrong dimension")}

      if(nL>0){
        if(any(prior$B$V[pL,-pL, drop=FALSE]!=0)){stop("sorry - sir effects and classic fixed effects have to be independent a priori")}
        prior$L$mu<-prior$B$mu[pL,1, drop=FALSE]
        prior$L$V<-prior$B$V[pL,pL, drop=FALSE]
        prior$B$mu<-prior$B$mu[-pL,1, drop=FALSE]
        prior$B$V<-prior$B$V[-pL,-pL, drop=FALSE]
      }else{
        prior$L=list(V=diag(0)*1e+10, mu=matrix(0,nL,1))
      }   
   }	   

################
# Missing data #
################

    if(any(substr(data$MCMC_family.names, 1,3)=="cen" & data$MCMC_y==data$MCMC_y.additional)){      # replace liabilities of family of censored variables with y=y.additional          
      cen_areknown<-which(substr(data$MCMC_family.names, 1,3)=="cen" & data$MCMC_y==data$MCMC_y.additional)                          
      data$MCMC_family.names[cen_areknown]<-substr(data$MCMC_family.names[cen_areknown], 4, nchar(as.character(data$MCMC_family.names[cen_areknown])))
    }

    cnt<-1
    mvtype<-c()
    proposal<-c()
    for(i in 1:nR){
       mp<-matrix(data$MCMC_y[cnt+1:(nfl[i+nG]*nrl[i+nG])-1], nrl[i+nG], nfl[i+nG])
       fp<-matrix(match(data$MCMC_family.names, family.types)[cnt+1:(nfl[i+nG]*nrl[i+nG])-1], nrl[i+nG], nfl[i+nG])
       proposal<-c(proposal, AMtune[i+nG]*((fp==11)*(mp==0 & is.na(mp)==FALSE))[,1])
       missing.pattern<-fp*(is.na(mp)==FALSE) # 0 for missing data 1 for gaussian >1 for other
       mvtype_tmp<-apply(missing.pattern, 1,function(x){all(x==1 | x==0 | x==20)})-2 # -2 if observed non-gaussian non-threshold present, -1 otherwise 

       if(nfl[i+nG]==1 & slice){    # if univariate 
          mvtype_tmp[which(missing.pattern==14)]<-0  # ordinal
          if(max(mp, na.rm=T)==1){                   # binary
            mvtype_tmp[which(missing.pattern==3)]<-0
          }    
       }   
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){all(x==1 | x==20)}))]<-0  
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){all(x==1)}))]<-2          
       mvtype_tmp[which(apply(missing.pattern, 1,function(x){all(x==0)}))]<-1  
  
       if(all(mvtype_tmp>c(-1))){AMtune[i+nG]=FALSE}  # If everything can be gibbsed/sliced do not tune
       mvtype<-c(mvtype,mvtype_tmp) 
       cnt<-cnt+nfl[i+nG]*nrl[i+nG]

       # missing data codes: 2 complete gaussian (ignore);
       #                     1 completely missing (unconditional Gibbs)
       #                     0 univariate slice sampling or complete threshold/gaussian response
       #                    -1 partial missing with observed being Gaussian (conditional Gibbs)
       #                    -2 partial or fully observed with non-gaussian (MH)  currently -1 & -2 are both MHed. 

    }

    if(any(mvtype==0) & !is.null(path.terms)){stop("slice samlping not possible with path/sir models")}

    if(is.null(start$liab)){
       data$MCMC_liab<-rnorm(length(data$MCMC_y)) 
       if(QUASI==TRUE){ 
        for(i in unique(data$MCMC_error.term)){
          trait_set<-which(!is.na(data$MCMC_y) & data$MCMC_dummy==0 & data$MCMC_error.term==i)
          missing_set<-which(is.na(data$MCMC_y) & data$MCMC_dummy==0 & data$MCMC_error.term==i)
          dummy_set<-which(is.na(data$MCMC_y) & data$MCMC_dummy==1 & data$MCMC_error.term==i)
          if(length(trait_set)<2){
            if(length(trait_set)==0){
              warning(paste("all observations are missing for error term ", i, ": liabilities sampled from Norm(0,1)", sep=""))
            }
            mu<-0
            v<-1
          }else{ 
            data_tmp<-data[trait_set,]          
            family_set<-data_tmp$MCMC_family.names[1]

            if(family_set=="poisson" | family_set=="ztpoisson"){
              mu<-mean(data_tmp$MCMC_y)
              v<-abs(log(((var(data_tmp$MCMC_y)-mu)/(mu^2))+1))
              mu<-log(mu)-0.5*v
            }
            if(family_set=="cenpoisson"){
              mu<-mean(data_tmp$MCMC_y+1)
              v<-abs(log(abs(((var(data_tmp$MCMC_y+1)-mu)/(mu^2))+1)))
              mu<-log(mu)-0.5*v
            }          
            if(family_set=="multinomial"){
              if(length(table(data_tmp$MCMC_y))>2){
                 m1<-summary(glm(cbind(MCMC_y, MCMC_y.additional)~1, family="quasibinomial", data=data_tmp))
                 v<-abs(((as.numeric(m1$dispersion[1])-0.5)/2)^2)
                 mu<-as.numeric(m1$coef[1])
              }else{
                 v<-1
                 mu<-0
              }
            }
            if(family_set=="exponential" | family_set=="cenexponential"){
              if(any(data_tmp$MCMC_y==0)){
                data_tmp$MCMC_y[which(data_tmp$MCMC_y==0)]<-1e-6
              }
              m1<-summary(glm(MCMC_y~1, family="Gamma", data=data_tmp))
              v<-abs((as.numeric(m1$dispersion[1])-1)/2)
              mu<-as.numeric(m1$coef[1])
            }
            if(family_set=="geometric"){
              mu<-1/(mean(data_tmp$MCMC_y)+1)
              mu<-log(mu)-log(1-mu)
              v<-mu^2
            }
            if(family_set=="ordinal" | family_set=="threshold"){
              v<-1
              mu<-qnorm(cumsum(c(0,table(data_tmp$MCMC_y)/length(data_tmp$MCMC_y))), 0, sqrt(1+(data_tmp$MCMC_family.names[1]=="ordinal")))[2]   
            }
            if(family_set=="cengaussian"){ 
              v<-var(apply(cbind(data_tmp$MCMC_y,data_tmp$MCMC_y.additional),1,function(x){mean(x[which(abs(x)!=Inf)])}))
              mu<-mean(apply(cbind(data_tmp$MCMC_y,data_tmp$MCMC_y.additional),1,function(x){mean(x[which(abs(x)!=Inf)])}))
            }
            if(family_set=="gaussian"){ 
              v<-var(data_tmp$MCMC_y)
              mu<-mean(data_tmp$MCMC_y)
            }
            if(family_set=="zipoisson" | family_set=="hupoisson" | family_set=="zapoisson"){
              if(max(data_tmp$MCMC_y)==1){
                mu<-mean(data_tmp$MCMC_y==1)
                if(family_set!="zapoisson"){
                   mu<-log(mu/(1-mu))
                }else{
                   mu<-log(-log(mu))
                }
                v<-diag(GRprior[[nG+nR]]$V)[length(diag(GRprior[[nG+nR]]$V))]
              }else{
                data_tmp<-data_tmp[-which(data_tmp$MCMC_y==0),]  
                mu<-mean(data_tmp$MCMC_y)
                v<-abs(log(((var(data_tmp$MCMC_y)-mu)/(mu^2))+1))
                mu<-log(mu)-0.5*v
              }
            }
            if(family_set=="zibinomial"){
              if(max(data_tmp$MCMC_y)==1){
                mu<-mean(data_tmp$MCMC_y==1)
                mu<-log(mu/(1-mu))
                v<-diag(GRprior[[nG+nR]]$V)[length(diag(GRprior[[nG+nR]]$V))]
              }else{
                 data_tmp<-data_tmp[-which(data_tmp$MCMC_y==0),]  
                 m1<-summary(glm(cbind(MCMC_y, MCMC_y.additional)~1, family="quasibinomial", data=data_tmp))
                 v<-abs(((as.numeric(m1$dispersion[1])-0.5)/2)^2)
                 mu<-as.numeric(m1$coef[1])
              }
            }
          }
          if(is.na(v) | is.na(mu)){
            warning(paste("good starting values not obtained for error term", i, "liabilities: using Norm(0,1)"))
            if(is.na(v)){v<-1}
            if(is.na(mu)){mu<-0}
          }
          if(length(missing_set)>0){
            data$MCMC_liab[missing_set]<-rnorm(length(missing_set),mu, sqrt(v))
          }
          if(length(dummy_set)>0){
            data$MCMC_liab[dummy_set]<-rnorm(length(dummy_set), 0, sqrt(v))
          }
          if(length(trait_set)>0){
            l_tmp<-sort(rnorm(length(trait_set), mu, sqrt(v)))
            mix<-sample(1:length(l_tmp),max(1,as.integer(length(l_tmp)/2)))
            l_tmp[mix]<-l_tmp[rev(mix)] # mix up 50% of latent variables to stop complete separation
            data$MCMC_liab[trait_set][order(data$MCMC_y[trait_set])]<-l_tmp
          }
        }
      }
    }else{
      if(length(c(start$liab))!=length(data$MCMC_y)){stop("liabilities must have the same dimensions as the response")}
      if(any(is.na(start$liab))){stop("starting liabilities must not contain missing vlaues")}
      if(any(data$MCMC_family.names=="cengaussian" & (start$liab>data$MCMC_y.additional | start$liab<data$MCMC_y))){
       stop("starting liabilities for censored gaussian data must lie between consoring points")
      }
      data$MCMC_liab<-c(start$liab)
    }

    if(any(data$MCMC_family.names=="cengaussian" & (data$MCMC_liab>data$MCMC_y.additional | data$MCMC_liab<data$MCMC_y), na.rm=T)){
      outside<-which(data$MCMC_family.names=="cengaussian" & (data$MCMC_liab>data$MCMC_y.additional | data$MCMC_liab<data$MCMC_y))
       # liabilities for censored data have to lie in between censoring points.
       data$MCMC_liab[outside]<-mapply(data$MCMC_y[outside], data$MCMC_y.additional[outside], FUN=function(x,y){if(x==c(-Inf) | y==Inf){if(x==-Inf){y}else{x}}else{runif(1, x, y)}})
       if(any(data$MCMC_liab[outside]=="Inf")){
          data$MCMC_liab[outside][which(data$MCMC_liab[outside]==Inf)]<-0
       }
    }

    observed<-as.numeric(is.na(data$MCMC_y)==FALSE)
    data$MCMC_y[is.na(data$MCMC_y)]<-0
    data$MCMC_y.additional[is.na(data$MCMC_y.additional)]<-0
    data$MCMC_y[which(data$MCMC_y==-Inf | data$MCMC_y==Inf)]<-sign(data$MCMC_y[which(data$MCMC_y==-Inf | data$MCMC_y==Inf)])*1e+32
    data$MCMC_y.additional[which(data$MCMC_y.additional==-Inf | data$MCMC_y.additional==Inf)]<-sign(data$MCMC_y.additional[which(data$MCMC_y.additional==-Inf | data$MCMC_y.additional==Inf)])*1e+32
	
    if(any(data$MCMC_family.names=="gaussian")){                                       # replace liabilities of ovserved gaussian data with data                                    
      data$MCMC_liab[which(data$MCMC_family.names=="gaussian" & observed)]<-data$MCMC_y[which(data$MCMC_family.names=="gaussian" & observed)]
    }

    split<-unlist(lapply(GRprior, function(x){x$fix}))
    if(any(split[1:nG]>nfl[1:nG] | (split[1:nG]<1 & split[1:nG]!=0))){stop("fix term in priorG must be at least one less than the dimension of V")}
    if(split[nG+1]>nfl[nG+1] && covu==0){stop("fix term in priorR must be at least one less than the dimension of V")}
    if(nR>1){
     if(any(split[nG+2:nR]>nfl[nG+2:nR])){stop("fix term in priorR must be at least one less than the dimension of V")}
    }
    if(any(split>1 & grepl("corg|corgh", update))){stop("sorry, fix terms cannot yet be used in conjunction with corg/corh structures")}
    if(any(split>1 & grepl("ante", update))){stop("sorry, fix terms cannot yet be used in conjunction with antedependence structures")}

    if(diagR==2){  # need to reform priors such that the marginal distribution of us is equal to distribution of idh 
      GRprior[[nG+nR]]$V<-GRprior[[nG+nR]]$V*GRprior[[nG+nR]]$n/(GRprior[[nG+nR]]$n+(dim(GRprior[[nG+nR]]$V)[1]+1))
      GRprior[[nG+nR]]$n<-GRprior[[nG+nR]]$n+(dim(GRprior[[nG+nR]]$V)[1]+1)
    }	
    GRinv<-unlist(lapply(GR, function(x){c(solve(x))}))
    GRvpP<-lapply(GRprior, function(x){(x$V)*(x$n)})
    if(any(update=="corg")){  # correlation matrices get I which is removed from Gtmp.
      for(i in which(update=="corg")){
	GRvpP[[i]]<-diag(nrow(GRprior[[i]]$V))  	
      }
    }	 
    if(any(update=="corgh")){  # correlation matrices get Diag(V) which is removed from Gtmp.
      for(i in which(update=="corgh")){
	GRvpP[[i]]<-diag(diag(GRprior[[i]]$V))  	
      }
    }
    if(any(update=="cors")){  # correlation sub-matrix gets I which is removed from Gtmp.
      for(i in which(update=="cors")){
	GRvpP[[i]][GRprior[[i]]$fix:nrow(GRprior[[i]]$V),GRprior[[i]]$fix:nrow(GRprior[[i]]$V)]<-diag(length(GRprior[[i]]$fix:nrow(GRprior[[i]]$V))) 	
      }
    }	

    GRvpP<-unlist(GRvpP)
    if(diagR==3){  # need to add aditional prior nu because of way trait:units are updated
      GRprior[[nG+nR]]$n<-GRprior[[nG+nR]]$n+(nrl[nG+nR]+1)*(dim(GRprior[[nG+nR]]$V)[1]-1)
    }
    GRnpP<-unlist(lapply(GRprior, function(x){c(x$n)}))      

    nanteP<-as.numeric(gsub("[a-z]", "", update))*grepl("ante", update)
    if(any(is.na(nanteP))){
      nanteP[which(is.na(nanteP))]<-0
    }
    if(any(!nanteP<nfl)){
      stop("order of the ante-dependence structure must be less than the dimensions of the matrix")
    }
    nanteP<-c(nanteP, unlist(lapply(GRprior, function(x){length(x$beta.mu)})))
    nanteP<-c(nanteP, grepl("ante.*v$", update))

    anteBmupP<-unlist(lapply(GRprior, function(x){c(x$beta.mu)}))      
    anteBvpP<-unlist(lapply(GRprior, function(x){if(!is.null(x$beta.V)){c(solve(x$beta.V))}}))      
      
    BvpP<-c(solve(prior$B$V), sum(prior$B$V!=0)==dim(prior$B$V)[1])
    BmupP<-c(prior$B$mu)
    if(nL>0){
      LvpP<-c(solve(prior$L$V), sum(prior$L$V!=0)==dim(prior$L$V)[1])
      LmupP<-c(prior$L$mu)
    }else{
      LvpP<-1
      LmupP<-0
    }

    AmupP<-unlist(lapply(GRprior, function(x){x$alpha.mu}))
    PXterms<-unlist(lapply(GRprior, function(x){all(x$alpha.V==0)==FALSE}))
    AVpP<-list2bdiag(lapply(GRprior, function(x){if(all(x$alpha.V==0)==FALSE){solve(x$alpha.V)}else{x$alpha.V-999}}))
    if(any(PXterms)){
      PXlevels<-which(apply(AVpP, 1, function(x){all(x!=-999)}))
      AVpP<-as(AVpP[PXlevels,PXlevels,drop=FALSE], "sparseMatrix")    
      AmupP<-AmupP[PXlevels]
    }else{ 
      PXterms<-rep(0, sum(unlist(lapply(GRprior, function(x){dim(x$V)[1]}))))
      AmupP<-1
      AVpP<-as(diag(1), "sparseMatrix")
    }

    update[which(update=="idh" | update=="idv" | update=="us")]<-1
    update[which(update==1 & split>1)]<-2
    update[grep("corg|corgh", update)]<-3
    update[which(split==1)]<-0
    update[grep("cors", update)]<-4
    update[grep("ante", update)]<-5
    update[grep("sub", update)]<-6

    update<-c(update, covu)

    if(covu>0){
      update<-c(update, update[nG+1])
      update[0:1+nG]<-0
    }else{
      update<-c(update, 0)
    }

    # update codes: 0 fixed - do not sample;
    #             : 1 unstructured 
    #             : 2 block diagonal constrained 
    #             : 3 correlation
    #             : 4 unstructured with correlation sub-matrix
    #             : 5 ante-dependence structure
    #             : 6 us + identity direct sum 

    # diagR codes:  1 normal
    #            :  2 idh(trait):units specification turned into us(trait):units to keep latent variables together
    #            :  3 trait:units specification turned into us(trait):units to keep latent variables together

    nordinal<-length(ncutpoints)
    if(nordinal==0){
      ncutpoints<-1
      stcutpoints<-1
    }  # no cutpoints need to be estimated if cutpoints=3 (standard binary)
    ncutpoints_store<-sum((ncutpoints-3)*(ncutpoints>3))

    data$MCMC_family.names<-match(data$MCMC_family.names, family.types)     # add measurement error variances and y.additional 

    if(nitt%%1!=0){stop("nitt must be integer")}
    if(thin%%1!=0){stop("thin must be integer")}
    if(burnin%%1!=0){stop("burnin must be integer")}

    nkeep<-ceiling((nitt-burnin)/thin)

    if(nkeep<1){stop("burnin is equal to, or greater than number of iterations")}
    if(nG==0){pr<-FALSE}

    Loc<-1:((sum((nfl*nrl)[1:nG])*pr+dim(X)[2]+nL*0)*nkeep)
    lambda<-1:(nL*nkeep)

    if(ncutpoints_store>0){
      CP<-1:(ncutpoints_store*nkeep)
    }else{
      CP<-1
    }
    dbar<-1:(2+nkeep)

    if(pl==TRUE){
      PLiab<-1:(length(data$MCMC_y)*nkeep)
    }else{
      PLiab<-1
    }

    Var<-1:((length(GRinv)-covu^2)*nkeep)

    if(all(Aterm==0)){
      ginverse<-list(A=as(diag(1), "sparseMatrix"))
    }

    if(DIC==TRUE){
      if(((sir | !is.null(path.terms)) & any(mvtype!=2)) | nordinal>1){
         DIC<-FALSE
      }
    }

	output<-.C("MCMCglmm",
        as.double(data$MCMC_y),   
        as.double(data$MCMC_y.additional),
        as.double(data$MCMC_liab), 
        as.integer(mvtype),   
        as.integer(length(data$MCMC_y)),
        as.integer(c(X@Dim,Z@Dim, L@Dim, unlist(lapply(ginverse, function(x){x@Dim[1]})))),  
        as.integer(c(length(X@x),length(Z@x),length(L@x),unlist(lapply(ginverse, function(x){length(x@x)})))),       
        as.integer(X@i),         
        as.integer(X@p),        
        as.double(X@x),         	  
        as.integer(Z@i),         
        as.integer(Z@p),          
        as.double(Z@x),               
        as.integer(unlist(lapply(ginverse, function(x){x@i}))),         
        as.integer(unlist(lapply(ginverse, function(x){x@p}))),       
        as.double(unlist(lapply(ginverse, function(x){x@x}))),                 
    	as.integer(Aterm-1),	  	  	  	  
        as.integer(nfl),
        as.integer(nrl),
        as.integer(update),
        as.integer(split-1),
        as.integer(c(nG,nR)), 
        as.double(GRinv), 
        as.double(GRvpP),            
        as.double(GRnpP),            
        as.integer(nitt),
        as.integer(thin),
        as.integer(burnin),
        as.integer(c(pr, pl)),
        as.double(Loc),
        as.double(Var),
        as.double(PLiab),
        as.integer(data$MCMC_family.names),
        as.double(c(unlist(tune))),
        as.integer(verbose),
        as.double(BvpP),
        as.double(BmupP),
        as.integer(mfac), 
	as.integer(observed),
        as.integer(diagR-1),
        as.integer(AMtune),
	as.integer(DIC),
        as.double(dbar),	  
        as.integer(proposal),
        as.integer(ncutpoints),
        as.integer(nordinal),
        as.double(stcutpoints),
        as.double(CP),
        as.double(AmupP),
        as.integer(AVpP@i),         
        as.integer(AVpP@p),        
        as.double(AVpP@x),    
        as.integer(length(AVpP@x)), 
        as.integer(PXterms),
        as.integer(L@i),    # Gianola & Sorensen's Lambda         
        as.integer(L@p),        
        as.double(L@x),   
        as.double(lambda),        
        as.double(LvpP),    # prior for structural parameters
        as.double(LmupP),
        as.double(nanteP),    # prior antedpendence betas  
        as.double(anteBvpP),   
        as.double(anteBmupP)             
        )

        Sol<-t(matrix(output[[30]], sum((nfl*nrl)[1:nG])*pr+dim(X)[2], nkeep))

        if(pr){     
          colnames(Sol)<-c(colnames(X), colnames(Z))
        }else{
          colnames(Sol)<-c(colnames(X))
        }
        colnames(Sol)<-gsub("MCMC_", "", colnames(Sol))

        if(nL>0){       
           lambda<-t(matrix(output[[58]], nL, nkeep))  
           colnames(lambda)<-Lnames     
           lambda<-mcmc(lambda, start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin)
        }else{
           lambda<-NULL
        }

       if(ncutpoints_store!=0){
         CP<-mcmc(t(matrix(output[[48]],ncutpoints_store, nkeep)), start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin)
         colnames(CP)<-c(paste("cutpoint.trait", rep(ordinal.names, ncutpoints-3), ".", unlist(sapply(ncutpoints-3, function(x){if(x!=0){1:x}})), sep=""))
         colnames(CP)<-gsub("MCMC_", "", colnames(CP))
       }else{
         CP<-NULL
       }


       VCV<-t(matrix(output[[31]], length(GRinv)-covu^2, nkeep))

       if(covu){
         if(nG==1){
           ustart<-0
         }else{
           ustart<-sum(nfl[1:(nG-1)]^2)
         }
         gnames<-gsub(".*MCMCsplit","", variance.names[ustart+seq(1, nfl[nG]^2, nfl[nG])])
         rnames<-gsub(".*MCMCsplit","", variance.names[ustart+nfl[nG]^2+seq(1, nfl[nG+1]^2, nfl[nG+1])])
         grnames<-apply(expand.grid(c(gnames, rnames), c(gnames, rnames)), 1, paste, collapse=":")
         if(nG==1){
           if(nR>1){
             variance.names<-c(grnames, variance.names[(sum(nfl[1:(nG+1)]^2)+1):length(variance.names)])
           }else{
             variance.names<-grnames
           }
         }else{
           if(nR>1){
             variance.names<-c(variance.names[1:ustart], grnames, variance.names[(sum(nfl[1:(nG+1)]^2)+1):length(variance.names)])
           }else{
             variance.names<-c(variance.names[1:ustart], grnames)
           }
         }
       }

        colnames(VCV)<-variance.names
        colnames(VCV)<-gsub("MCMC_", "", colnames(VCV))
        colnames(VCV)<-gsub("MCMCsplit", ":", colnames(VCV))

        if(covu>0){
          nfl[nG+1]<-nfl[nG+1]+covu
          nfl<-nfl[-nG]
          nrl<-nrl[-nG]
          nrt<-nrt[-ngstructures]
          nG<-nG-1
          ngstructures<-ngstructures-1
        }

        if(diagR==2){  # idh structures that were held as us 
          VCV<-VCV[,-c(((dim(VCV)[2]-nfl[nG+1]^2+1):dim(VCV)[2])[-diag(matrix(1:(nfl[nG+1]^2),nfl[nG+1],nfl[nG+1]))]),drop=FALSE]          
          colnames(VCV)[(dim(VCV)[2]-nfl[nG+1]+1):dim(VCV)[2]]<-sapply(colnames(VCV)[(dim(VCV)[2]-nfl[nG+1]+1):dim(VCV)[2]], function(x){substr(x, gregexpr(":", x)[[1]][ceiling(length(gregexpr(":", x)[[1]])/2)]+1, nchar(x))})
          nfl<-c(nfl,rep(1,nfl[nG+1]-1))
          nrl<-c(nrl,rep(nrl[nG+1],nfl[nG+1]-1))
          nrt[ngstructures+1]<-nfl[nG+1]
          nR<-nfl[nG+1]
          nfl[nG+1]<-1        
        }
        if(diagR==3){ # trait:unit structures that were held as us 
          VCV<-VCV[,-((dim(VCV)[2]-nfl[nG+1]^2+2):dim(VCV)[2]),drop=FALSE]
          colnames(VCV)[dim(VCV)[2]]<-"trait:units" 
          nrl[nG+1]<-nrl[nG+1]*nfl[nG+1]
          nfl[nG+1]<-1
          nrt[ngstructures+1]<-1
          nR<-1
        }

        if(DIC){
         deviance<-mcmc(-2*output[[43]][1:nkeep], start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin)
         DIC<--4*output[[43]][nkeep+1]+2*output[[43]][nkeep+2]
        }else{
         deviance<-NULL
         DIC<-NULL
        }

        dummy.data<-which(data$MCMC_dummy[order(ordering)]==1)

        if(pl==TRUE){
          Liab<-t(matrix(output[[32]], length(data$MCMC_y), nkeep))
          Liab<-Liab[,order(ordering),drop=FALSE]
          if(length(dummy.data)>0){
            Liab<-Liab[,-dummy.data,drop=FALSE]
          }
          Liab<-mcmc(Liab, start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin)
        }else{
          Liab<-NULL
        }

        nF<-dim(X)[2]

        if(saveX==FALSE){
          X<-NULL
        }else{
          X<-X[order(ordering),,drop=FALSE]
          if(length(dummy.data)>0){
            X<-X[-dummy.data,,drop=FALSE]
          }
        }

        if(saveZ==FALSE | (nG==0 & covu==0)){
          Z<-NULL
        }else{
          Z<-Z[order(ordering),,drop=FALSE]
          if(length(dummy.data)>0){
            Z<-Z[-dummy.data,,drop=FALSE]
          }
        }

        if(saveZ){
          if(length(dummy.data)>0){
            ZR<-ZR[-dummy.data,,drop=FALSE]
          }
        }

        if(saveXL==FALSE | nL==0){
          L<-NULL
        }else{
           if(!is.null(path.terms)){
             L<-kronecker(as.matrix(L), Diagonal(nrl[nG+1]))
           }
           Lordering<-rep(order(ordering),nL)+rep(((1:nL)-1)*length(ordering),each=length(ordering))
           L<-L[order(ordering),Lordering,drop=FALSE]
           Ldummy.data<-rep(dummy.data,nL)+rep(((1:nL)-1)*length(dummy.data),each=length(dummy.data))
           if(length(dummy.data)>0){
             L<-L[-dummy.data,-Ldummy.data,drop=FALSE]
           }
        }

        if(all(Aterm==0)){
          ginverse=NULL
        }

        error.term<-data$MCMC_error.term[order(ordering)]
        family<-family.types[data$MCMC_family.names[order(ordering)]]
        y.additional<-data$MCMC_y.additional[order(ordering)]

        if(any(family=="multinomial")){
           family[which(family=="multinomial")]<-paste("multinomial", y.additional[which(family=="multinomial")], sep="")
        }

        if(length(dummy.data)>0){
          error.term<-error.term[-dummy.data]
          family<-family[-dummy.data]
        }

    	options("na.action"=orig.na.action)
        if(nG==0){
          Gnfl<-NULL
          Gnrl<-NULL
          Gnat<-NULL
          Gnrt<-NULL
          Rnfl<-nfl
          Rnrl<-nrl
          Rnrt<-nrt
        }else{
          Gnfl<-nfl[1:nG]
          Gnrl<-nrl[1:nG]
          Gnat<-Aterm[1:nG]
          Gnrt<-nrt[1:ngstructures]
          Rnfl<-nfl[nG+1:nR]
          Rnrl<-nrl[nG+1:nR]
          Rnrt<-nrt[-c(1:ngstructures)]
        }

        Tune<-as.list(1:nR)
        for(i in 1:nR){
          if(covu>0 & i==1){
            Tune[[i]]<-matrix(output[[34]][1:((Rnfl[1]-covu)^2)],Rnfl[1]-covu,Rnfl[1]-covu)
          }else{
            Tune[[i]]<-matrix(output[[34]][sum(Rnfl[1:i]^2)-Rnfl[i]^2+1:(Rnfl[i]^2)],Rnfl[i],Rnfl[i])
          }
        }

        output<-list(
            Sol=mcmc(Sol, start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin),
            Lambda = lambda,
            VCV=mcmc(VCV, start=burnin+1, end=burnin+1+(nkeep-1)*thin, thin=thin),
            CP=CP,
            Liab=Liab,
            Fixed=list(formula=original.fixed, nfl=nF, nll=nL),
            Random=list(formula = original.random, nfl=Gnfl, nrl=Gnrl, nat=Gnat, nrt=Gnrt),
            Residual=list(formula = original.rcov, nfl=Rnfl, nrl=Rnrl, nrt=Rnrt, family=rterm.family, original.family=original.family),
            Deviance=deviance,
            DIC=DIC,
            X=X,
            Z=Z,
            ZR=ZR,
            XL=L,
            ginverse=ginverse,
            error.term=error.term,
            family=family, 
            Tune=list2bdiag(Tune)
        )

	class(output)<-c("MCMCglmm")
        output
}

