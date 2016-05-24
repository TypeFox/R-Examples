"additivePenal" <-
function (formula, data, correlation=FALSE, recurrentAG=FALSE, cross.validation=FALSE, n.knots, kappa,
             maxit=350, hazard="Splines", nb.int, LIMparam=1e-4, LIMlogl=1e-4, LIMderiv=1e-3, print.times=TRUE)
{

##### hazard specification ######

haztemp <- hazard
hazard <- strsplit(hazard,split="-")
hazard <- unlist(hazard)  
if(!(length(hazard) %in% c(1,2))){stop("Please check and revise the hazard argument according to the format specified in the help.")}

### longueur hazard = 1

if((all.equal(length(hazard),1)==T)==T){
   if(!(hazard %in% c("Weibull","Piecewise","Splines"))){
	stop("Only 'Weibull', 'Splines' or 'Piecewise' hazard can be specified in hazard argument.")
   }else{
	typeof <- switch(hazard,"Splines"=0,"Piecewise"=1,"Weibull"=2)
	if(typeof %in% c(0,2)){
### Splines
		if (!(missing(nb.int))){
			stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int' arguments must be deleted.")
		}
		if (typeof == 0){
			size1 <- 100
			size2 <- 100
			equidistant <- 2
			nbintervR <- 0
		}
### Weibull
		if (typeof == 2){
			equidistant <- 2
			nbintervR <- 0
			size1 <- 100
		}
	}else{
		stop ("The hazard argument is incorrectly specified.Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
	}
   } 
}else{
### Picewise
      typeof <- 1
#### longueur hazard > 1
      if(!("Piecewise" %in% hazard)){
	stop("Only 'Piecewise' hazard can be specified in hazard argument in this case")
      }
      if(!(all(hazard %in% c("Piecewise","per","equi")))){
		stop ("The hazard argument is incorrectly specified.Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
      }else{
		if (!(haztemp %in% c("Piecewise-per","Piecewise-equi"))){
			stop ("The hazard argument is incorrectly specified. Please refer to the help file of frailtypack.")
		}
		equidistant <- switch(haztemp,"Piecewise-per"=0,"Piecewise-equi"=1)
      }
 }


#AD:
	if (missing(formula))stop("The argument formula must be specified in any model")
	if(class(formula)!="formula")stop("The argument formula must be a formula")
#AD:  
	if(typeof == 0){ 
		if (missing(n.knots)) stop("number of knots are required")
#AD:	 
		n.knots.temp <- n.knots	
#AD
		if (n.knots<4) n.knots<-4
		if (n.knots>20) n.knots<-20
		
		if (missing(kappa))stop("smoothing parameter (kappa) is required")
		
	}else{
		if (!(missing(n.knots)) || !(missing(kappa)) || !(missing(cross.validation))){
			stop("When parametric hazard function is specified, 'kappa', 'n.knots' and 'cross.validation' arguments must be deleted.")	
		}
		n.knots <- 0
		kappa <- 0
		crossVal <- 0
	}


    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$correlation <- m$n.knots <- m$recurrentAG <- m$cross.validation <- m$kappa <- m$maxit <- m$hazard <- m$nb.int <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$print.times <- m$... <- NULL
    special <- c("strata", "cluster", "slope")
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)   
    ord <- attr(Terms, "order")
  #  if (length(ord) & any(ord != 1)) 
   #     stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    cluster <- attr(Terms, "specials")$cluster

    # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
    classofY <- attr(model.extract(m, "response"),"class")
    # attention le package pec rajoute un element dans l'attribut "class" des objets de survie
    if (length(classofY)>1) classofY <- classofY[2]

    typeofY <- attr(model.extract(m, "response"),"type")

#Al : tri du jeu de donnees par cluster croissant
	if (length(cluster)){
		tempc <- untangle.specials(Terms, "cluster", 1:10)
		ord <- attr(Terms, "order")[tempc$terms]
		if (any(ord > 1))stop("Cluster can not be used in an interaction")
		m <- m[order(m[,tempc$vars]),] # soit que des nombres, soit des caracteres
		cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
		uni.cluster <- unique(cluster)
	}
#Al

    if (NROW(m) == 0) 
        stop("No (non-missing) observations")
   
    Y <- model.extract(m, "response")
    if (classofY != "Surv")
        stop("Response must be a survival object")
    ll <- attr(Terms, "term.labels")

    mt <- attr(m, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, m, contrasts)


    strats <- attr(Terms, "specials")$strata
    cluster <- attr(Terms, "specials")$cluster
    slope <- attr(Terms, "specials")$slope

	if (length(cluster)){
		ll <- ll[-grep("cluster",ll)]
	}
	if (length(slope)){
		ll <- ll[-grep("slope",ll)]
	}
	if (length(strats)){
		ll <- ll[-grep("strata",ll)]
	}
#=========================================================>
	ind.place <- grep("factor",ll)
	vecteur <- NULL
	vecteur <- c(vecteur,ll[ind.place])
	mat.factor <- matrix(vecteur,ncol=1,nrow=length(vecteur))

 # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
  vec.factor <-apply(mat.factor,MARGIN=1,FUN=function(x){
  if (length(grep("factor",x))>0){
    if(length(grep(":",x))>0){
      if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
        
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
        pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos4 <- length(unlist(strsplit(x,split="")))
        return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
        pos4 <- length(unlist(strsplit(x,split="")))-1
        return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else{#both factors
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
        pos4 <- length(unlist(strsplit(x,split="")))-1
        return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
      }
    }else{
      pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      pos2 <- length(unlist(strsplit(x,split="")))-1
      return(substr(x,start=pos1,stop=pos2))}
  }else{
    return(x)
  }})

#=========================================================>
#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
#	if(length(vec.factor) > 0){
#		vect.fact <- attr(X,"dimnames")[[2]]
#
#		vect.fact <- vect.fact[grep("factor",vect.fact)]
#		vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
#		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
#		pos2 <- length(unlist(strsplit(x,split="")))-2
#		return(substr(x,start=pos1,stop=pos2))})		
#		occur <- rep(0,length(vec.factor))
	

#		for(i in 1:length(vec.factor)){
#			occur[i] = sum(vec.factor[i] == vect.fact)
#		}
#	}
  if(length(vec.factor) > 0){
  vect.fact <- attr(X,"dimnames")[[2]]
   vect.fact <- vect.fact[grep(paste(vec.factor,collapse="|"),vect.fact)]
 
  occur <- rep(0,length(vec.factor))
  
  interaction<-as.vector(apply(matrix(vect.fact,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
  which.interaction <- which(interaction==1)
  
  for(i in 1:length(vec.factor)){
      if(length(grep(":",unlist(strsplit(vec.factor[i],split=""))))>0){
        pos <- grep(":",unlist(strsplit(vec.factor[i],split="")))
        length.grep <- 0
        for(j in 1:length(vect.fact)){
         if(j%in%which.interaction){
          if(length(grep(substr(vec.factor[i],start=1,stop=pos-1),vect.fact[j]))>0 && length(grep(substr(vec.factor[i],start=pos+1,stop=length(unlist(strsplit(vec.factor[i],split="")))),vect.fact[j]))>0){
              length.grep <- length.grep + 1
              which <- i}
        }}
       occur[i] <- length.grep
      }else{
       if(length(vect.fact[-which.interaction])>0){occur[i] <- length(grep(vec.factor[i],vect.fact[-which.interaction]))
      }else{occur[i] <- length(grep(vec.factor[i],vect.fact))}
    }
  }
  }
#=========================================================>    

    dropx <- NULL

    if (length(cluster)) {
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) 
            stop("Cluster can not be used in an interaction")
        cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
        dropx <- tempc$terms
        uni.cluster<-unique(cluster)
    }
    else
     {
      stop("grouping variable is needed")   
     }

    if(length(uni.cluster)==1) 
     {
      stop("grouping variable must have more than 1 level")   
     }
 
 

    if (length(strats)){
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
        uni.strat<-length(unique(strats))

        if (missing(kappa)) stop("smoothing parameter (kappa) is required")
        if ((typeof==0) & (length(kappa)!=uni.strat)) stop("wrong length of argument 'kappa' for the current stratification")

        if (uni.strat!=2) stop("maximum number of strata is 2")

    }
    else
      {
        uni.strat <- 1
        if ((typeof==0) & (length(kappa)!=uni.strat)) stop("wrong length of argument 'kappa' for the current stratification")
        strats <- rep(1,nrow(data))
        kappa <- c(kappa,0)
      }


    if (length(slope))
     {
        temps <- untangle.specials(Terms, "slope", 1:10)
        dropx <- c(dropx, temps$terms)
        ord <- attr(Terms, "order")[temps$terms]
        if (any(ord > 1)) 
            stop("'slope' can not be used in an interaction")
        aux<-temps$vars
        nnOK<-gsub(")","",gsub("slope\\(","",aux))
     }
    else
     {
      stop("interacting (slope) variable is needed")   
     }

    
    #type <- attr(Y, "type")
    type <- typeofY
    if (type != "right" && type != "counting") 
        stop(paste("Cox model doesn't support \"", type, "\" survival data", 
            sep = ""))

    if (type != "counting" && recurrentAG)
       stop("recurrentAG needs counting process formulation")
    

    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
    Xlevels <- .getXlevels(newTerms, m)
    contr.save <- attr(X, 'contrasts')

	if(length(vec.factor) > 0){
#========================================>
		position <- unlist(assign,use.names=F)
	}
#========================================>
    if (ncol(X) == 1) 
      {
         X<-X-1
         noVar<-1 
      }
    else
      {
         X <- X[, -1, drop = FALSE]
         noVar<-0
      }     


    nn<-dimnames(X)[[2]]
#AD
    if(!(nnOK %in% ll)) stop("covariate between 'slope()' missing in the terms formula.")



    varInt<-c(1:length(nn))[nn==nnOK]


    nvar<-ncol(X) 

    var<-matrix(c(X),nrow=nrow(X),ncol=nvar)

    n<-nrow(X)    

    if (type=="right")
      {
        tt0 <- rep(0,n)
        tt1 <- Y[,1]
        cens <- Y[,2]
      }
    else
      {
        tt0 <- Y[,1]
        tt1 <- Y[,2]
        cens <- Y[,3]
      }

    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0

    AG<-ifelse(recurrentAG,1,0)
    if (typeof == 0){ 
 	crossVal<-ifelse(cross.validation,0,1)
    }

#=======================================>
#======= Construction du vecteur des indicatrice
	if(length(vec.factor) > 0){
		if(length(vec.factor) > 0){
			k <- 0
			for(i in 1:length(vec.factor)){
				ind.place[i] <- ind.place[i]+k
					k <- k + occur[i]-1
			}
		}
	}
#==================================

        if(equidistant %in% c(0,1)){
		if (missing(nb.int)) stop("Time interval 'nb.int' is required")
		if (class(nb.int) != "numeric") stop("The argument 'nb.int' must be a numeric")	
		if (nb.int < 1) stop("Number of time interval 'nb.int' must be between 1 and 20")
		if (nb.int > 20){
			 nb.int <-20
			indic.nb.int <- 1 # equals 1 for nb.int > 20
		}else{
			indic.nb.int <- 0 # equals 1 for nb.int < 20
		}
		nbintervR <- nb.int
		size1 <- 3*nbintervR
	}
	if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0
	np <- switch(as.character(typeof),
		"0"=(as.integer(uni.strat) * (as.integer(n.knots) + 2) + as.integer(nvar) + as.integer(correlation) + 2 * 1),
		"1"=(as.integer(uni.strat) * nbintervR + as.integer(nvar) + as.integer(correlation) + 2 * 1),
		"2"=(as.integer(uni.strat) * 2 + nvar + as.integer(correlation) + 2 * 1))	
 

		xSu1 <- rep(0,100)
		xSu2 <- rep(0,100)
		if (typeof==0){
			mt1 <- size1	
		}else{
			mt1 <- 100
		}
		size2 <- mt1
			
			if (print.times){
				ptm<-proc.time()
				cat("\n")
				cat("Be patient. The program is computing ... \n")
			}

			ans <- .Fortran("additive",
				as.integer(n),
				as.integer(length(uni.cluster)),
				as.integer(uni.strat),
				as.integer(n.knots),
				as.double(kappa[1]),
				as.double(kappa[2]),
				as.double(tt0),
				as.double(tt1),
				as.integer(cens),
				as.integer(cluster),
				as.integer(nvar),
				as.integer(strats),
				as.double(var),
				as.integer(varInt),
				as.integer(AG),
				as.integer(noVar),
				as.integer(maxit),
				as.integer(crossVal),
				as.integer(correlation),
				as.integer(np),
				b=as.double(rep(0,np)),
				coef=as.double(rep(0,nvar)),
				varcoef=as.double(rep(0,nvar)),
				varcoef2=as.double(rep(0,nvar)),
				rho=as.double(0),
				cov=as.double(0),
				varcov=as.double(0),
				varSigma2=as.double(c(0,0)),
				varTau2=as.double(c(0,0)),
				ni=as.integer(0),
				loglikpen=as.double(0),
				LCV=as.double(rep(0,2)),
				k0=as.double(c(0,0)),
				x1=as.double(rep(0,size1)),
				lam1=as.double(matrix(0,nrow=size1,ncol=3)),
                                xSu1=as.double(xSu1),
				surv1=as.double(matrix(0,nrow=size2,ncol=3)),
				x2=as.double(rep(0,size1)),
				lam2=as.double(matrix(0,nrow=size1,ncol=3)),
                                xSu2=as.double(xSu2),
				surv2=as.double(matrix(0,nrow=size2,ncol=3)),
				as.integer(typeof),
				as.integer(equidistant),
				as.integer(nbintervR),
				as.integer(size1),
				ier=as.integer(0),
				ddl=as.double(0),
				istop=as.integer(0),
				shape.weib=as.double(rep(0,2)),
				scale.weib=as.double(rep(0,2)),
				as.integer(mt1),
				trunc=as.integer(0),
				zi=as.double(rep(0,(n.knots+6))),
				time=as.double(rep(0,(nbintervR+1))),
				
				martingale.res=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.pred=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.pred2=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.var=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.var2=as.double(rep(0,as.integer(length(uni.cluster)))),
				frailty.cov=as.double(rep(0,as.integer(length(uni.cluster)))),
				linear.pred=as.double(rep(0,n)),
				EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
				PACKAGE = "frailtypack")  # 62
		

    if (ans$trunc == 1){
    	stop("'addivePenal' can not deal with left truncation.")
    }
    if (ans$istop == 4){
	 warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
     }

    if (ans$istop == 2){
         warning("Model did not converge.")
    }
    if (ans$istop == 3){
         warning("Matrix non-positive definite.")
    }

    flush.console()
    if (print.times){
        cost<-proc.time()-ptm
        cat("The program took", round(cost[3],2), "seconds \n")
    }

    if (noVar==1) nvar<-0

    fit <- NULL
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- length(uni.cluster)
    fit$n.events <- sum(cens)  #coded 0: censure 1:event

    if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans$loglikpen
    }else{
        fit$logLik <- ans$loglikpen
    }
#AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
#AD:
    fit$type <- type

    fit$b <- ans$b

    if (noVar==1) {
      fit$coef <- NULL
    }
    else
     {
       fit$coef <- ans$coef
       names(fit$coef) <- factor.names(colnames(X))
     }

    fit$varH <- ans$varcoef
    fit$varHIH <- ans$varcoef2

    fit$cov <- ans$cov
    fit$varcov <- ans$varcov

    fit$sigma2<-ans$b[np-nvar-1]^2
    fit$varSigma2<-ans$varSigma2
    fit$tau2<-ans$b[np-nvar]^2
    fit$varTau2<-ans$varTau2
 
    fit$rho<-ans$rho
    fit$cov<-ans$cov
    fit$varcov<-ans$varcov

    fit$formula <- formula(Terms)

    fit$n.strat <- uni.strat
    fit$n.knots<-n.knots
    if (uni.strat > 1) fit$kappa <- ans$k0
    else fit$kappa <- ans$k0[1]
    fit$n.iter <- ans$ni
    fit$cross.Val<-cross.validation
    fit$correlation<-correlation

    fit$x <- cbind(ans$x1,ans$x2)
    fit$lam <- array(c(ans$lam1,ans$lam2), dim=c(size1,3,2))
    fit$xSu <- cbind(ans$xSu1,ans$xSu2)
    fit$surv <- array(c(ans$surv1,ans$surv2), dim=c(size2,3,2))

    fit$npar <- np
    fit$type <- type
    fit$AG <- recurrentAG

#AD:
    fit$noVar <- noVar
    fit$nvar <- nvar
#AD:
    fit$typeof <- typeof
    fit$istop <- ans$istop
    if (typeof == 0){
    	fit$DoF <- ans$ddl
    	fit$n.knots.temp <- n.knots.temp
	fit$zi <- ans$zi
    }
    if(typeof == 1){
	fit$time <- ans$time
	fit$nbintervR <- nbintervR
    }
    fit$indic.nb.int <- indic.nb.int
    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib
##
	fit$martingale.res <- ans$martingale.res
	fit$frailty.pred <- ans$frailty.pred
	fit$frailty.pred2 <- ans$frailty.pred2
# 	fit$frailty.var <- ans$frailty.var
# 	fit$frailty.var2 <- ans$frailty.var2
# 	fit$frailty.cov <- ans$frailty.cov
	fit$linear.pred <- ans$linear.pred
##
    fit$EPS <- ans$EPS

#AD
    if(ans$ier==2000)
        stop("The cross validation procedure cannot be finished. Try to change 
          either the number of knots or the seed for kappa parameter")

#========================= Test de Wald pour shared
	
	if(length(vec.factor) > 0){
		Beta <- ans$coef
		VarBeta <- diag(ans$varcoef)
		nfactor <- length(vec.factor)
		p.wald <- rep(0,nfactor)

		fit$global_chisq <- waldtest(N=nvar,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta)
		fit$dof_chisq <- occur
		fit$global_chisq.test <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factor)){
			p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occur[i]), 3)
		}
		fit$p.global_chisq <- p.wald
		fit$names.factor <- vec.factor 

		
	}else{
		fit$global_chisq.test <- 0
	}
	
#===============================================	


if (length(Xlevels) >0)fit$Xlevels <- Xlevels
    fit$contrasts <- contr.save

    class(fit) <- "additivePenal"
    fit
}


