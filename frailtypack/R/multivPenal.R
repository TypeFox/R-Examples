
"multivPenal" <-
 function(formula, formula.Event2, formula.terminalEvent, data, initialize=TRUE, recurrentAG=FALSE,
  n.knots, kappa, maxit=350, hazard="Splines", nb.int, print.times=TRUE)
{

## pour l'utilisateur theta1 est theta, et theta2 est eta (variances des frailties)

NN <- colnames(get_all_vars(update(formula,"~1"),data))
## evt 1
EVENT1 <- NN[length(NN)]
# t1 t0
NOM <- NN[-length(NN)]
# evt 2
TT <- untangle.specials(terms(formula, c("event2")), "event2", 1:10)$vars 
start <- grep("\\(",unlist(strsplit(TT,"")))
stop <- grep("\\)",unlist(strsplit(TT,"")))
EVENT2 <- substr(TT,start=start+1,stop=stop-1)
# cluster
TT <- untangle.specials(terms(formula, c("cluster")), "cluster", 1:10)$vars 
start <- grep("\\(",unlist(strsplit(TT,"")))
stop <- grep("\\)",unlist(strsplit(TT,"")))
CLUSTER <- substr(TT,start=start+1,stop=stop-1)


Tab <- transfo.table(data,EVENT1,EVENT2,NOM,CLUSTER,recurrentAG)

# n.knots = c(n.knots1, n.knots2, n.knots3)
# kappa = c(kappa1, kappa2, kappa3)
# nb.int = c(nb.int1, nb.int2, nb.int3)

call <- match.call()
mCall <- match.call()
mCall$formula <- mCall$formula.terminalEvent <- mCall$formula.Event2 <- mCall$recurrentAG <-
mCall$initialize <- mCall$n.knots <- mCall$kappa <- mCall$maxit <- 
mCall$hazard <- mCall$nb.int <- mCall$print.times <- NULL




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
	### Splines (equidistant par defaut)
		if (typeof == 0){
			if (!missing(nb.int)){
				stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int' argument must be deleted.")
			}
			size1 <- 100
			size2 <- 100
			size3 <- 100
			equidistant <- 1
			nbintervR <- 0
			nbintervDC <- 0
			nbintervR2 <- 0
		}
		### Weibull
		if (typeof == 2){
			if (!missing(nb.int)){
				stop("When the hazard function equals 'Splines' or 'Weibull', 'nb.int' argument must be deleted.")
			}
			size1 <- 100
			size2 <- 100
			size3 <- 100
			equidistant <- 2
			nbintervR <- 0
			nbintervDC <- 0
			nbintervR2 <- 0
		}
		if (typeof == 1){
			stop ("The hazard argument is incorrectly specified. Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
		}
	}
}else{
	#### longueur hazard > 1
	if(all(!(c("Splines","Piecewise") %in% hazard))){
		stop("Only 'Splines' and 'Piecewise' hazard can be specified in hazard argument in this case")
	}
	### Splines percentile
	if ("Splines" %in% hazard){
		typeof <- 0
		if(!(all(hazard %in% c("Splines","per")))){
			stop ("The hazard argument is incorrectly specified. Only 'per' is allowed with 'Splines'. Please refer to the help file of frailtypack.")
		}else{
			size1 <- 100
			size2 <- 100
			size3 <- 100
			equidistant <- 0
			nbintervR <- 0
			nbintervDC <- 0
			nbintervR2 <- 0
		}
	}
	### Piecewise (per or equi)
	if ("Piecewise" %in% hazard){
		typeof <- 1
		if(!(all(hazard %in% c("Piecewise","per","equi")))){
			stop ("The hazard argument is incorrectly specified. Type of hazard are required ('per' or 'equi'). Please refer to the help file of frailtypack.")
		}else{
			if (!(haztemp %in% c("Piecewise-per","Piecewise-equi"))){
				stop ("The hazard argument is incorrectly specified. Please refer to the help file of frailtypack.")
			}
			equidistant <- switch(haztemp,"Piecewise-per"=0,"Piecewise-equi"=1)
		}
	}
}

#AD:
	if (missing(formula))stop("The argument formula must be specified in any model")
	if(class(formula)!="formula")stop("The argument formula must be a formula")
	
	if(typeof == 0){
#AD:   
		if (missing(n.knots))stop("number of knots are required")   
		if (length(n.knots)!= 3)stop("length of knots must be 3")
#AD:	 
		# permettre al'utilisateur de rentrer les n.knots dans cet ordre : loco, meta, deces
		n.knots.temp <- n.knots
		n.knots[2] <- n.knots.temp[3]
		n.knots[3] <- n.knots.temp[2]

		n.knots.temp <- n.knots	
#AD
		n.knots[n.knots<4] <- 4
		n.knots[n.knots>20] <- 20

		if (missing(kappa))stop("smoothing parameter (kappa1) is required")
		if(class(kappa)!="numeric")stop("The argument kappa must be a numeric")
		if(length(kappa)==1)stop("length of smoothing parameter (kappa) must greater than 2")

		# permettre a l'utilisateur de rentrer les kappa dans cet ordre : loco, meta, deces
		kappa.temp <- kappa
		kappa[2] <- kappa.temp[3]
		kappa[3] <- kappa.temp[2]

		if(length(kappa)==2){
			kappa1 <- kappa[1]
			kappa3 <- kappa[2] 
			indic.Kappa2 <- 1
			stop("smoothing parameter for death is required for the multiv model")
		}else{
			kappa1 <- kappa[1]
			kappa2 <- kappa[2] 	
			kappa3 <- kappa[3] 
			indic.Kappa2 <- 0
		}
 
	}else{
		if (!(missing(n.knots)) || !(missing(kappa))){
			stop("When parametric hazard is not 'Splines' function is specified, 'kappa' and 'n.knots' arguments must be deleted.")
		}
		n.knots <- 0
		kappa1 <- 0
		kappa2 <- 0
		kappa3 <- 0
		crossVal <- 0
		
	}

#==========================================================================================>
#===============================>          formula         <================================
#==========================================================================================>
	m <- mCall
	m$data <- Tab$dataR
	
	special <- c("strata", "cluster", "terminal","event2")	
	Terms <- terms(formula, special)
	Terms.formula <- Terms	
	ord <- attr(Terms, "order") ## longueur de ord=nbre de var.expli	
#	if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
	m$formula <- Terms	
	m[[1]] <- as.name("model.frame") 
		
	m.formula <- eval(m)

	if (NROW(m.formula) == 0)stop("No (non-missing) observations") 
	Y <- model.extract(m.formula, "response") 
	if (!inherits(Y, "Surv"))stop("Response must be a survival object") 	
	ll <- attr(Terms, "term.labels")	
	X <- if (!is.empty.model(attr(m.formula,"terms")))model.matrix(attr(m.formula,"terms"),m.formula,contrasts) 

	ind.place <- attr(X,"assign")[duplicated(attr(X,"assign"))]

	vec.factor <- NULL
	vec.factor <- c(vec.factor,ll[ind.place])

#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
	strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
	cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()

	if (length(cluster)){
		ll_tmp <- ll[grep("cluster",ll)]
		ll <- ll[-grep("cluster",ll)]

		pos1 <- grep("r",unlist(strsplit(ll_tmp,split="")))[1]+2
		pos2 <- length(unlist(strsplit(ll_tmp,split="")))-1
		Names.cluster <- substr(ll_tmp,start=pos1,stop=pos2)
	}
	if (length(strats)){
		ll <- ll[-grep("strata",ll)]
	}

#   plus besoin de as.factor() pour afficher le test de Wald global
	mat.factor <- matrix(vec.factor,ncol=1,nrow=length(vec.factor))
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


 ind.place <- grep(paste(vec.factor,collapse="|"),ll[-c(grep("terminal",ll),grep("event2",ll))])


if(length(vec.factor) > 0){
  vect.fact <- attr(X,"dimnames")[[2]]
  
  #vect.fact <- vect.fact[grep("factor",vect.fact)]
  vect.fact <- vect.fact[grep(paste(vec.factor,collapse="|"),vect.fact)]
  
  #   	vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
  # 		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
  # 		pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
  # 		return(substr(x,start=pos1,stop=pos2))})
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

#==========================================================================================>
#=============================>         end  formula         <==============================
#==========================================================================================>


#==========================================================================================>
#===============================>          terminalEvent        <===========================
#==========================================================================================>

	terminalEvent <- attr(Terms.formula, "specials")$terminal 
	dropx <- NULL

	if (length(cluster)){
		tempc <- untangle.specials(Terms.formula, "cluster", 1:10)
		ord <- attr(Terms.formula, "order")[tempc$terms]
	#	if (any(ord > 1))stop("Cluster can not be used in an interaction")
		
		
		cluster <- strata(m.formula[, tempc$vars], shortlabel = TRUE)
		dropx <- tempc$terms
		uni.cluster<-unique(cluster)
	}else{
		stop("grouping variable is needed")   
	}

	if(length(uni.cluster)==1){ 
		stop("grouping variable must have more than 1 level")   
	}

	if (length(strats)){

		temp <- untangle.specials(Terms.formula, "strata", 1)
		dropx <- c(dropx, temp$terms)
		if (length(temp$vars) == 1)strata.keep <- m.formula[[temp$vars]]
		else strata.keep <- strata(m.formula[, temp$vars], shortlabel = TRUE)
		strats <- as.numeric(strata.keep)
		uni.strat<-length(unique(strats))

		if (uni.strat!=2)stop("maximum number of strata is 2")
	}else{
		uni.strat<-1
		strats <- rep(1,nrow(data))
		kappa2 <- 0 
		kappa3 <- 0
	}
	
#AD: indicator of terminal()
	ind.terminal <- length(terminalEvent)
#AD:
	if (length(terminalEvent)){
	
		tempterm <- untangle.specials(Terms.formula, "terminal", 1:10) 
#ici on comme terme tempterm$vars qui est le nom dans l'appel(ex;"terminal(death)"
#et tempterm$terms qui est la position de la variable dans l'appel, ici elle vient a la position 6
	
		ord <- attr(Terms.formula, "order")[tempterm$terms] # ord[6]=1 ici dans notre exemple
		
	#	if (any(ord > 1))stop("Terminal can not be used in an interaction")
		dropx <- c(dropx,tempterm$terms) # vecteur de position
		terminal <- strata(m.formula[, tempterm$vars], shortlabel = TRUE)
		terminal <- as.numeric(as.character(terminal))

	}
#==========================================================================================>
#===============================>          end terminalEvent        <========================
#==========================================================================================>

#==========================================================================================>
#===============================>          formula.Event2        <==========================
#==========================================================================================>
	if (!missing(formula.Event2)){
		tempterm2 <- untangle.specials(Terms.formula, "event2", 1:10) 
		dropx <- c(dropx,tempterm2$terms) # vecteur de position
	}

	formula2Event <- attr(Terms.formula, "specials")$event2
	ind.event2 <- length(formula2Event)
	if (length(formula2Event)) event2 <- Tab$dataM[,EVENT2]	


#==========================================================================================>
#==========================================================================================>


	type <- attr(Y, "type")

	if (type != "right" && type != "counting"){ 
		stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
	}

	if (type != "counting" && recurrentAG){
		stop("recurrentAG needs counting process formulation")
	}
	
	#drop contient les position liees au fonction() ic ex:cluster(id) et terminal(death)

	if (length(dropx)){ 
		newTerms <- Terms.formula[-dropx]
	}else{
		newTerms <- Terms.formula
	}

#newTerm vaut Terms - les variables dont les position sont dans drop

	X.formula <- model.matrix(newTerms, m.formula)
	
	assign <- lapply(attrassign(X.formula, newTerms)[-1], function(x) x - 1)


	if(length(vec.factor) > 0){
		position <- unlist(assign,use.names=F)
	}

	if (ncol(X.formula) == 1){ 
		X.formula <- X.formula-1
		noVar1 <- 1 
	}else{
		X.formula <- X.formula[, -1, drop = FALSE]
		noVar1 <- 0
	}
	Names.formula <- colnames(X.formula)
	
	nvar.formula <- ncol(X.formula) 
	var.formula  <- matrix(X.formula,nrow=nrow(X.formula),ncol=nvar.formula)
	n.formula <- nrow(X.formula)

	if (type=="right"){
		tt0.formula <- rep(0,n.formula)
		tt1.formula <- Y[,1]
		cens.formula <- Y[,2]

		tt0.formula2 <- rep(0,dim(Tab$dataM)[1])
		if(length(NOM)==2){
			tt1.formula2 <- Tab$dataM[,NOM[2]]
		}else{
			tt1.formula2 <- Tab$dataM[,NOM]
		}
	}else{
		tt0.formula <- Y[,1]
		tt1.formula <- Y[,2]
		cens.formula <- Y[,3]

		tt0.formula2 <- Tab$dataM[,NOM[1]]
		tt1.formula2 <- Tab$dataM[,NOM[2]]

	}

	if (min(cens.formula)==0) cens.data<-1
	if (min(cens.formula)==1 && max(cens.formula)==1) cens.data<-0

	AG<-ifelse(recurrentAG,1,0)

	if (typeof == 0){
		crossVal<-1 # on fait de la cross.validation pour les shared quel que soit la situation
	}


	flush.console()
	if(print.times){
		ptm<-proc.time()
		cat("\n")
		cat("Be patient. The program is computing ... \n")
	}


#======= Construction du vecteur des indicatrice
	if(length(vec.factor) > 0){
#		ind.place <- ind.place -1
		k <- 0
		for(i in 1:length(vec.factor)){
			ind.place[i] <- ind.place[i]+k
				k <- k + occur[i]-1
		}
	}

# 	if(Frailty =="FALSE"){
# 		stop("For the multiv frailty models, 'Frailty' must be equal to 'TRUE' ")
# 	}


	if (!recurrentAG)
	{
		tt1.death<-aggregate(tt1.formula,by=list(cluster),FUN=sum)[,2]
		tt0.death<-rep(0,length(tt1.death))
	}else{
		tt1.death<-aggregate(tt1.formula,by=list(cluster),FUN=function(x) x[length(x)])[,2]
		tt0.death<-rep(0,length(tt1.death))
	}


	if (!missing(formula.terminalEvent)){
		Terms2 <- terms(formula.terminalEvent, special)
		ord2 <- attr(Terms2, "order")
	
#		if (length(ord2) & any(ord2 != 1)){ 
#			stop("Interaction terms are not valid for terminal event formula")
#		}
	}

	if (!missing(formula.Event2)){

		Terms.formula2 <- terms(formula.Event2, special)
		
		ord.formula2 <- attr(Terms.formula, "order")
	
#		if (length(ord.formula2) & any(ord.formula2 != 1)){ 
#			stop("Interaction terms are not valid for event formula.Event2")
#		}
	}

#AD:Joint model needs "terminal()"
	if (ind.terminal){
		terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
	}else{
		stop(" multiv frailty model miss specified")
	}


# terminalEvent might be 0-1 

	if (!(all(terminalEvent%in%c(1,0)))){ 
		stop("terminal must contain a variable coded 0-1 and a non-factor variable")
	}

 	if (!(all(event2%in%c(1,0)))){ 
 		stop("event2 must contain a variable coded 0-1 and a non-factor variable")
 	}

#===========================================================================>
#================                 formulaTerminal
	if (!missing(formula.terminalEvent)){
		m2 <- mCall 
		m2$formula <- Terms2
		m2[[1]] <- as.name("model.frame")
		m2 <- eval(m2)
		match.noNA<-dimnames(m2)[[1]]%in%dimnames(m.formula)[[1]]

		m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA

		if (!missing(formula.terminalEvent))newTerms2<-Terms2
		
		Xdc <- model.matrix(newTerms2, m2)
		lldc <- attr(newTerms2,"term.labels")
		#ind.placedc <- grep("factor",lldc)
		ind.placedc <- attr(Xdc,"assign")[duplicated(attr(Xdc,"assign"))]
		
		vec.factordc <- NULL
		vec.factordc <- c(vec.factordc,lldc[ind.placedc])

		mat.factordc <- matrix(vec.factordc,ncol=1,nrow=length(vec.factordc))
		# Fonction servant a prendre les termes entre "as.factor"
		# Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
		vec.factordc <-apply(mat.factordc,MARGIN=1,FUN=function(x){
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
		

		if(length(vec.factordc) > 0){
		  vect.factdc <- attr(Xdc,"dimnames")[[2]]
		  
		  #vect.fact <- vect.fact[grep("factor",vect.fact)]
		  vect.factdc <- vect.factdc[grep(paste(vec.factordc,collapse="|"),vect.factdc)]
		 
		  #   	vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
		  # 		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
		  # 		pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
		  # 		return(substr(x,start=pos1,stop=pos2))})
		  occurdc <- rep(0,length(vec.factordc))
		  
		  interactiondc<-as.vector(apply(matrix(vect.factdc,nrow=length(vect.factdc)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
		  which.interactiondc <- which(interactiondc==1)
		
      for(i in 1:length(vec.factordc)){
		  
		    if(length(grep(":",unlist(strsplit(vec.factordc[i],split=""))))>0){
		      
		      
		      pos <- grep(":",unlist(strsplit(vec.factordc[i],split="")))
		      length.grep <- 0
		      for(j in 1:length(vect.factdc)){
		        if(j%in%which.interactiondc){
		        
		          if(length(grep(substr(vec.factordc[i],start=1,stop=pos-1),vect.factdc[j]))>0 && length(grep(substr(vec.factordc[i],start=pos+1,stop=length(unlist(strsplit(vec.factordc[i],split="")))),vect.factdc[j]))>0){
		            length.grep <- length.grep + 1
		            which <- i}
		        }}
		      occurdc[i] <- length.grep
		      
		    }else{
		      
		      
		      if(length(vect.factdc[-which.interactiondc])>0){
            occurdc[i] <- length(grep(vec.factordc[i],vect.factdc[-which.interactiondc]))
		      }else{
            occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))}
		    }
		  }
		}

   
#=========================================================>
		assign <- lapply(attrassign(Xdc, newTerms2)[-1], function(x) x - 1)
#========================================>
		if(length(vec.factordc) > 0){
			positiondc <- unlist(assign,use.names=F)
		}
#========================================>
		if (ncol(Xdc) == 1)
		{
			Xdc<-Xdc-1
			noVar2 <- 1
		}else{
			Xdc <- Xdc[, -1, drop = FALSE]
			noVar2 <- 0
		}
		Names.terminal <- colnames(Xdc)

		vardc.temp<-matrix(c(Xdc),nrow=nrow(Xdc),ncol=ncol(Xdc))

# 		if(is.null(nrow(m2)))
# 		{
# 			if (length(m2) != nrow(m.formula)){
# 				stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
# 			}
# 		}else{
# 			
# 			if (nrow(m2) != nrow(m.formula)){
# 				stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
# 			}
# 
# 		}

		if (!is.null(ncol(vardc.temp))){
			vardc<-aggregate(vardc.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]
			if (ncol(vardc.temp)>1){
			
				for (i in 2:ncol(vardc.temp)){
					vardc.i<-aggregate(vardc.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
					vardc<-cbind(vardc,vardc.i)
				}
			}
		}else{
			vardc<-aggregate(vardc.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]
		}

		if(is.matrix(vardc)){
			vardc <- matrix(vardc,nrow=nrow(vardc),ncol=ncol(vardc))
			nvardc <- ncol(vardc)
		}else{
			nvardc <- 1
			vardc <- matrix(vardc,nrow=length(vardc),ncol=1)
		}
	}else{
		noVar2 <- 1
		vardc<-0
		nvardc <- 0
	}



#===========================================================================>
#================                 formula.Event2
# 
# 	varformula2 <- model.matrix(formula.Event2,data=Tab$dataM)
# 	
# 	if (!missing(formula.Event2)){
# 		mformula2 <- mCall
# 		mformula2$data <- Tab$dataM
# 		mformula2$formula <- formula.Event2
# 		mformula2[[1]] <- as.name("model.frame")
# 		mformula2 <- eval(mformula2)
# 		match.noNA<-dimnames(mformula2)[[1]]%in%dimnames(m.formula)[[1]]
# 	
# 		mformula2<-mformula2[match.noNA, ,drop=FALSE]
# 			
 #		if (!missing(formula.Event2))newTerms.formula2<-Terms.formula2

#=========================================================>
# 		if (!missing(formula.Event2)){
# 			llformula2 <- attr(newTerms.formula2,"term.labels")
# 			ind.placeformula2 <- grep("factor",llformula2)
# 			vecteur.formula2 <- NULL
# 			vecteur.formula2 <- c(vecteur.formula2,llformula2[ind.placeformula2])
# 			mat.factor.formula2 <- matrix(vecteur.formula2,ncol=1,nrow=length(vecteur.formula2))
# 
# # Fonction servant a prendre les termes entre "as.factor"

# 			vec.factor.formula2 <-apply(mat.factor.formula2,MARGIN=1,FUN=function(x){
# 			pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
# 			pos2 <- length(unlist(strsplit(x,split="")))-1
# 			return(substr(x,start=pos1,stop=pos2))})
# 		}	
# 	
#=========================================================>
# 		if (!missing(formula.Event2)){
			
# #=========================================================>
# # On determine le nombre de categorie pour chaque var categorielle
# 			if(length(vec.factor.formula2) > 0){
# 				vect.fact.formula2 <- attr(Xformula2,"dimnames")[[2]]
# 				occurformula2 <- rep(0,length(vec.factor.formula2))
# 	
# 				for(i in 1:length(vec.factor.formula2)){
# 					occurformula2[i] = length(grep(vec.factor.formula2[i],vect.fact.formula2))
# 				}
# 			}
#=========================================================>
# 			assign.formula2 <- lapply(attrassign(Xformula2, newTerms.formula2)[-1], function(x) x - 1)
# #========================================>
# 			if(length(vec.factor.formula2) > 0){
# 				position.formula2 <- unlist(assign.formula2,use.names=F)
# 			}
# #========================================>
# 			if (ncol(Xformula2) == 1) 
# 			{
# 				Xformula2<-Xformula2-1
# 				noVar3 <- 1 
# 			}else{
# 				Xformula2 <- Xformula2[, -1, drop = FALSE]
# 				noVar3 <- 0 
# 			} 
# 
# 			Names.formula2 <- colnames(Xformula2)
# 
# 			varformula2<-matrix(c(Xformula2),nrow=nrow(Xformula2),ncol=ncol(Xformula2))
# 	
# # 			if(is.null(nrow(mformula2)))
# # 			{
# # 				if (length(mformula2) != nrow(m.formula)){
# 					stop(" There are missing values in the covariates modelling the event2 event. \n Prepare data only with complete cases")
# 				}
# 			}else{
# 				
# 				if (nrow(mformula2) != nrow(m.formula)){
# 					stop(" There are missing values in the covariates modelling the event2 event. \n Prepare data only with complete cases")
# 				}
# # 	
# # 			}
# 			nvarformula2 <-ncol(Xformula2)
# 		}else{
# 			noVar3 <- 1
# 			varformula2 <-0
# 			nvarformula2 <- 0
# 		}
# 	}else{
# 		noVar3 <- 1
# 		varformula2 <-0
# 		nvarformula2 <- 0
# 	}
# 
	varformula2 <- model.matrix(formula.Event2,data=Tab$dataM)
	
	Names.formula2 <- colnames(varformula2)

	ind.placeformula2 <- attr(varformula2,"assign")[duplicated(attr(varformula2,"assign"))]
	#ind.placeformula2 <- grep("factor",colnames(varformula2))
  ind.placeformula2.temp <- ind.placeformula2
  for(i in 1:length(ind.placeformula2.temp)){
    ind.placeformula2[i] <- which(attr(varformula2,"assign")==ind.placeformula2.temp[i])[1]-1
    
  }


	vec.factorevent2 <- NULL

	vec.factorevent2 <- c(vec.factorevent2,Names.formula2[ind.placeformula2+1])


	mat.factor.formula2 <- matrix(vec.factorevent2,ncol=1,nrow=length(vec.factorevent2))
# Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions

  vec.factorevent2 <-apply(mat.factor.formula2,MARGIN=1,FUN=function(x){
  if (length(grep("factor",x))>0){
    if(length(grep(":",x))>0){
      if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
        
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
        pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos4 <- length(unlist(strsplit(x,split="")))
        return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
       pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
        pos4 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
        return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
      }else{#both factors
        pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
        pos2 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
        pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
        pos4 <- grep("\\)",unlist(strsplit(x,split="")))[2]-1
        return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
      }
    }else{
      pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      pos2 <- grep("\\)",unlist(strsplit(x,split="")))[1]-1
      return(substr(x,start=pos1,stop=pos2))}
  }else{
    return(x)
  }})


  if(length(vec.factorevent2) > 0){
  vect.fact.formula2 <- attr(varformula2,"dimnames")[[2]]
  
  #vect.fact <- vect.fact[grep("factor",vect.fact)]
  vect.fact.formula2 <- vect.fact.formula2[grep(paste(vec.factorevent2,collapse="|"),vect.fact.formula2)]
  
  #   	vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
  # 		pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
  # 		pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
  # 		return(substr(x,start=pos1,stop=pos2))})
  occurformula2 <- rep(0,length(vec.factorevent2))
  
  if(length(vect.fact.formula2)>0){

  interaction.formula2<-as.vector(apply(matrix(vect.fact.formula2,nrow=length(vect.fact.formula2)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
  which.interaction.formula2 <- which(interaction.formula2==1)
  
  for(i in 1:length(vec.factorevent2)){
    
    if(length(grep(":",unlist(strsplit(vec.factorevent2[i],split=""))))>0){
      
      
      pos <- grep(":",unlist(strsplit(vec.factorevent2[i],split="")))
      length.grep <- 0
      for(j in 1:length(vect.fact.formula2)){
        if(j%in%which.interaction.formula2){
          
          if(length(grep(substr(vec.factorevent2[i],start=1,stop=pos-1),vect.fact.formula2[j]))>0 && length(grep(substr(vec.factorevent2[i],start=pos+1,stop=length(unlist(strsplit(vec.factorevent2[i],split="")))),vect.fact.formula2[j]))>0){
            length.grep <- length.grep + 1
            which <- i}
        }}
      occurformula2[i] <- length.grep
      
    }else{
      
      
      if(length(vect.fact.formula2[-which.interaction.formula2])>0){occurformula2[i] <- length(grep(vec.factorevent2[i],vect.fact.formula2[-which.interaction.formula2]))
      }else{occurformula2[i] <- length(grep(vec.factorevent2[i],vect.fact.formula2))}
    }
  }
}
}
  if(is.na(vec.factorevent2))vec.factorevent2 <- NULL
	if (ncol(varformula2) == 1)
	{
		varformula2<-varformula2-1
		noVar3 <- 1
	}else{
		varformula2 <- varformula2[, -1, drop = FALSE]
		noVar3 <- 0
	}
	Names.formula2 <- colnames(varformula2)

	nvarformula2 <-ncol(varformula2)

# 	print("---------- Pris en compte de var expli -----------")
# 	print(" Var Recurrentes ")
# 	cat("noVar1",noVar1,"\n")
# 	cat("nvar.formula",nvar.formula,"\n")
# 	print(head(var.formula))
# 
# 	print("---------")
# 	print(" Var Death ")
# 	cat("noVar2",noVar2,"\n")
# 	cat("nvardc",nvardc,"\n")
# 
# 	print(head(vardc))
# 
# 	print("---------")
# 	print(" Var formula.Event2 ")
# 	cat("noVar3",noVar3,"\n")
# 	cat("nvarformula2",nvarformula2,"\n")
# 	print(head(varformula2))	

	nvarRec <- nvar.formula
	nvarEnd <- nvardc
	nvarRec2 <- nvarformula2
# 	cat("# Rec : ",nvarRec,"\n")
# 	cat("# Meta : ",nvarRec2,"\n")
# 	cat("# Dc : ",nvarEnd,"\n")

	if (!missing(formula.terminalEvent)){	
# 		print("---- formula.terminalEvent -----")
# 		print(vec.factordc)
# 		print(ind.placedc)
# 		print("----")
		if(length(vec.factordc) > 0){
			k <- 0
			for(i in 1:length(vec.factordc)){
				ind.placedc[i] <- ind.placedc[i]+k
					k <- k + occurdc[i]-1
			}
		}
	}

#	if (!missing(formula.Event2)){	
# # 		print("---- formula.Event2 -----")
# # 		print(vec.factor.formula2)
# # 		print(ind.placeformula2)
# # 		print("----")

 #   if(length(vec.factorevent2) > 0){
#			k <- 0
#			for(i in 1:length(vec.factorevent2)){
#				ind.placeformula2[i] <- ind.placeformula2[i]+k
#					k <- k + occurformula2[i]-1
#			}
#		}
#	}

# nombre total de variable
	if (sum(as.double(var.formula))==0) nvarRec <- 0
	if (sum(as.double(varformula2))==0) nvarEnd <- 0
	if (sum(as.double(vardc))==0) nvarRec2 <- 0

	nvar <- nvarRec + nvarEnd + nvarRec2
	
# ... end preparing data 
#--------- parametres
	effet <- 1
	indic_alpha <- 4
	indic_eta  <- 3	
	indic_rho <- 0
	indic_a1 <- 1
	indic_a2 <- 2
#--------- parametres
	nst = 3
	
        if ((equidistant %in% c(0,1)) & (typeof == 1)){
		# permettre a l'utilisateur de rentrer les nb.int dans cet ordre : loco, meta, deces
		nb.int.temp <- nb.int
		nb.int[2] <- nb.int.temp[3]
		nb.int[3] <- nb.int.temp[2]
		if (missing(nb.int)) stop("Time interval 'nb.int' is required")
		if (class(nb.int) != "numeric") stop("The argument 'nb.int' must be a numeric")
		if (length(nb.int) != 3) stop("The argument 'nb.int' must be a numeric vector of length 3")
		if (nb.int[1] < 1) stop("Number of Time interval 'nb.int[1]' must be between 1 and 20")
		if (nb.int[2] < 1) stop("Number of Time interval 'nb.int[2]' must be between 1 and 20")
		if (nb.int[3] < 1) stop("Number of Time interval 'nb.int[3]' must be between 1 and 20")

		if (nb.int[1] > 20){
			 nb.int[1] <-20
			indic.nb.int1 <- 1 # equals 1 for nb.int1 > 20
		}else{
			indic.nb.int1 <- 0 # equals 1 for nb.int1 < 20
		}
		
		if (nb.int[2] > 20){
			 nb.int[2] <-20
			indic.nb.int2 <- 1 # equals 1 for nb.int1 > 20
		}else{
			indic.nb.int2 <- 0 # equals 1 for nb.int1 < 20
		}
		
		if (nb.int[3] > 20){
			 nb.int[3] <-20
			indic.nb.int3 <- 1 # equals 1 for nb.int1 > 20
		}else{
			indic.nb.int3 <- 0 # equals 1 for nb.int1 < 20
		}

		nbintervR <- nb.int[1]
		size1 <- 3*nbintervR

		nbintervDC <- nb.int[2]
		size2 <- 3*nbintervDC

		nbintervR2 <- nb.int[3]
		size3 <- 3*nbintervR2
	
	}
	if ((typeof == 0) | (typeof == 2)){
		indic.nb.int1 <- 0
		indic.nb.int2 <- 0
		indic.nb.int3 <- 0
	}

	np <- switch(as.character(typeof),
		"0"=(sum(n.knots) + 6 + nvar + effet + indic_alpha),
		"1"=(nbintervR + nbintervDC + nbintervR2 + nvar + effet + indic_alpha),
		"2"=(2*nst + nvar + effet + indic_alpha))

# 	print("=============================")
# 	print(cens.formula)
# 	print("=============================")
# 	print(terminal)
# 	print("=============================")
	if (all(all.equal(as.numeric(cens.formula),terminal)==T)){
		stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
	}

	if (all(all.equal(as.numeric(cens.formula),event2)==T)){
		stop("'Recurrent event' variable and 'event2 event' variable need to be different")
	}

	xSu1 <- rep(0,100)
	xSu2 <- rep(0,100)
	xSu3 <- rep(0,100)

	if (typeof==0){
		mt11 <- size1
		mt12 <- size2
		mt13 <- size2
	}else{
		mt11 <- 100
		mt12 <- 100
		mt13 <- 100
	}

# 	print(" -------------------------------------------------------------")
# 	print(" -------------------------------------------------------------")
# 	print("------------------ Appel de multive Yassin -------------------")
# 	print(" -------------------------------------------------------------")
# 	print(" -------------------------------------------------------------")

	noVar <- c(noVar1,noVar2,noVar3)

# 	cat("nsujet0 :",as.integer(n.formula),"\n")
# 	cat("ng0 :",length(uni.cluster),"\n")
# 	cat("nz0 :",n.knots,"\n")
# 	cat("k0 :",kappa,"\n")

# 	print("=============                        =============")
# 	print("============= Verication du programme =============")

# 	print(n.knots)
	zi <- 0
	zidc <- 0
	zi2 <- 0

	if (typeof==0){
		zi <- rep(0,(n.knots[1] + 6))
		zidc <- rep(0,(n.knots[2] + 6))
		zi2 <- rep(0,(n.knots[3] + 6))
	}else{n.knots[1:3] <- 0} #AK 04/11/2015
  
	time <- 0
	timedc <- 0
	time2 <- 0

	if(typeof!=0)kappa <- rep(0,3)

	if (typeof==1){
		time <- rep(0,(nbintervR + 1))
		timedc <- rep(0,(nbintervDC + 1))
		time2 <- rep(0,(nbintervR2 + 1))
	}

	n <- n.formula
	n2 <- dim(Tab$dataM)[1]
	ng <- length(uni.cluster)

	nobsEvent <- c(n,n2,ng)
	nbvar <- c(nvarRec,nvarEnd,nvarRec2)

	maxIteration <- c(maxit,maxit)
	nbIntervEvent <- c(nbintervR,nbintervDC,nbintervR2)
	mtEvent <- c(size1,size2,size3)
	mt1Event <- c(mt11,mt12,mt13)
	ResMartingaleEvent <- matrix(0,nrow=ng,ncol=3)
	frailtyEstimates <- matrix(0,nrow=ng,ncol=5)

	#Recurrent
#	print(length(tt0.formula))
#	print(length(tt1.formula))
#	print(length(cens.formula))
#	print(length(Tab$dataR[,CLUSTER]))
	#Meta
#	print(length(tt0.formula2))
#	print(length(tt1.formula2))
#	print(length(event2))	
#	print(length(Tab$dataM[,CLUSTER]))
	
	#Mdeces
# 	print(length(tt0.death))
# 	print(length(tt1.death))
# 	print(length(terminalEvent))
# 	print(length(Tab$dataM[,CLUSTER]))
# 	
# 	print(dim(var.formula))
# 	print(dim(varformula2))
# 	print(dim(vardc))


	ans <- .Fortran("jointMultiv",
		as.integer(nobsEvent),
                as.integer(n.knots),
                k0=as.double(kappa),
                as.double(tt0.formula),
                as.double(tt1.formula),
                as.double(tt0.formula2),
                as.double(tt1.formula2),
                as.integer(cens.formula),
                as.integer(event2),
                as.integer(Tab$dataR[,CLUSTER]),
                as.integer(Tab$dataM[,CLUSTER]),
                as.integer(uni.cluster),
                as.double(tt0.death),
                as.double(tt1.death),
                as.integer(terminalEvent),
                as.integer(nbvar),
                as.double(var.formula),
                as.double(varformula2),
                as.double(vardc),
                as.integer(noVar),
                as.integer(maxIteration),
                as.integer(initialize),
                np=as.integer(np),
                b=as.double(rep(0,np)),
                H=as.double(matrix(0,nrow=np,ncol=np)),
                HIH=as.double(matrix(0,nrow=np,ncol=np)),
                loglik=as.double(0),
		LCV=as.double(rep(0,2)),
		critCV=as.integer(rep(0,5)),
                x1=as.double(rep(0,size1)),
                lam=as.double(matrix(0,nrow=size1,ncol=3)),
                xSu1=as.double(xSu1),
                surv=as.double(matrix(0,nrow=mt11,ncol=3)),
                x2=as.double(rep(0,size2)),
                lam2=as.double(matrix(0,nrow=size2,ncol=3)),
		xSu2=as.double(xSu2),
                surv2=as.double(matrix(0,nrow=mt12,ncol=3)),
                x3=as.double(rep(0,size3)),
                lam3=as.double(matrix(0,nrow=size3,ncol=3)),
		xSu3=as.double(xSu3),
                surv3=as.double(matrix(0,nrow=mt13,ncol=3)),
		as.integer(typeof),
		as.integer(equidistant),
		as.integer(nbIntervEvent),
		as.integer(mtEvent),
                ni=as.integer(0),
                cptEvent=as.integer(rep(0,3)),
		shape.weib=as.double(rep(0,3)),
		scale.weib=as.double(rep(0,3)),
		as.integer(mt1Event),
		as.integer(!crossVal),
		as.integer(recurrentAG),
		ResMartingaleEvent=as.double(ResMartingaleEvent),
		frailtyEstimates=as.double(frailtyEstimates),
		linear.pred=as.double(rep(0,n)),
		lineardc.pred=as.double(rep(0,ng)),
		linear2.pred=as.double(rep(0,n2)),
		zi=as.double(zi),
		zidc=as.double(zidc),
		zi2=as.double(zi2),
		time=as.double(time),
		timedc=as.double(timedc),
		time2=as.double(time2),
                PACKAGE = "frailtypack")


	if (ans$critCV[2] == 4){
		warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
	}
	
	if (ans$critCV[2] == 2){
		warning("Model did not converge. Change the 'maxit' parameter")
	}
	if (ans$critCV[2] == 3){
		warning("Matrix non-positive definite.")
	}

#AD:
    if (all(noVar==1)) nvar<-0
#AD:
    np <- ans$np
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- n
    fit$groups <- ng
    fit$n.events <- ans$cptEvent[1]
    fit$n.deaths <- ans$cptEvent[2]
    fit$n.events2 <- ans$cptEvent[3]
    fit$AG <- recurrentAG

    if(as.character(typeof)=="0"){
        fit$logLikPenal <- ans$loglik
    }else{
        fit$logLik <- ans$loglik
    }
#AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]
#AD: 
    fit$theta1 <- ans$b[np - nvar - indic_alpha]^2
    fit$theta2 <- ans$b[np - nvar - indic_eta]^2

    fit$alpha1 <- ans$b[np - nvar - indic_a1]
    fit$alpha2 <- ans$b[np - nvar - indic_a2]

    fit$rho <- 2*exp(ans$b[np - nvar - indic_rho])/(exp(ans$b[np - nvar - indic_rho])+1) - 1

    fit$npar <- np

#AD:
	Names <- NULL 
	if(noVar[1]!=1) Names <- c(Names,factor.names(Names.formula))
	if(noVar[2]!=1) Names <- c(Names,factor.names(Names.terminal))
	if(noVar[3]!=1) Names <- c(Names,factor.names(Names.formula2))

	coef <- NULL
	if (all(noVar==1)){
		fit$coef <- NULL
	}else{	
		if(noVar[1]!=1) coef <- c(coef,ans$b[(np - nvar + 1):(np - nvar + nvarRec)])
		if(noVar[2]!=1) coef <- c(coef,ans$b[(np - nvar + nvarRec + 1):(np - nvar + nvarRec + nvarEnd)])
		if(noVar[3]!=1) coef <- c(coef,ans$b[(np - nvar + nvarRec + nvarEnd + 1):np])
		fit$coef <- coef
		names(fit$coef) <- Names
     }

#AD:
    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)


    fit$nvar<-c(nvarRec,nvarEnd,nvarRec2)


    fit$varH <- temp1[(np - (nvar) - 1):np, (np - (nvar) - 1):np]
    fit$varHIH <- temp2[(np - (nvar) - 1):np, (np - (nvar) - 1):np]



    fit$formula <- formula(Terms.formula)
    fit$formula.terminalEvent <- formula(Terms2)
    fit$formula.Event2 <- formula(formula.Event2)

    fit$x1 <- ans$x1
    fit$lam1 <- matrix(ans$lam, nrow = size1, ncol = 3)
    fit$xSu1 <- ans$xSu1
    fit$surv1 <- matrix(ans$surv, nrow = mt11, ncol = 3)


    if (missing(formula.Event2)){
    	fit$x2 <- NA
    }else{
        fit$x2 <- ans$x3
    }
    fit$lam2 <- matrix(ans$lam3, nrow = size2, ncol = 3)
    fit$xSu2 <- ans$xSu3
    fit$surv2 <- matrix(ans$surv3, nrow = mt12, ncol = 3)

    if (missing(formula.terminalEvent)){
    	fit$xEnd <- NA
    }else{
        fit$xEnd <- ans$x2
    }
    fit$lamEnd <- matrix(ans$lam2, nrow = size3, ncol = 3)
    fit$xSuEnd <- ans$xSu2
    fit$survEnd <- matrix(ans$surv2, nrow = mt13, ncol = 3)


    fit$type <- type
    fit$n.iter <- ans$ni
    fit$typeof <- typeof

    if (typeof == 0){
		fit$n.knots<-n.knots
		fit$kappa <- ans$k0
		fit$n.knots.temp <- n.knots.temp
		fit$zi <- ans$zi
    }


    if(typeof == 1){
		fit$time <- ans$time
		fit$timedc <- ans$timedc
		fit$time2 <- ans$time2
    }


    fit$noVar <- noVar
    fit$nbintervR <- nbintervR
    fit$nbintervDC <- nbintervDC
    fit$nbintervR2 <- nbintervR2

    fit$nvarRec <- nvarRec
    fit$nvarEnd <- nvarEnd
    fit$nvarRec2 <- nvarRec2

    fit$istop <- ans$critCV[2]
    fit$indic.nb.int1 <- indic.nb.int1
    fit$indic.nb.int2 <- indic.nb.int2
    fit$indic.nb.int3 <- indic.nb.int3

    fit$shape.weib <- ans$shape.weib
    fit$scale.weib <- ans$scale.weib


#    if (Frailty){
	ResMartingaleEvent <- matrix(ans$ResMartingaleEvent,nrow=ng,ncol=3)
	frailtyEstimates <- matrix(ans$frailtyEstimates,nrow=ng,ncol=5)

	fit$martingale.res <- ResMartingaleEvent[,1]
	fit$martingaledeath.res <- ResMartingaleEvent[,2]
	fit$martingale2.res <- ResMartingaleEvent[,3]
	fit$frailty.pred <- frailtyEstimates[,1]
	fit$frailty2.pred <- frailtyEstimates[,2]
	fit$frailty.var <- frailtyEstimates[,3]
	fit$frailty2.var <- frailtyEstimates[,4]
	fit$frailty.corr <- frailtyEstimates[,5]
	
	fit$linear.pred <- ans$linear.pred
	fit$lineardeath.pred <- ans$lineardc.pred
	fit$linear2.pred <- ans$linear2.pred
#    }
 


#================================> For the reccurrent
#========================= Test de Wald
	ntot <- nvarEnd + nvarRec + nvarRec2
	Beta <- ans$b[(np - nvar + 1):np]
	VarBeta <- diag(diag(fit$varH)[-c(1,2)])
	if(length(vec.factor) > 0){
		nfactor <- length(vec.factor)
		p.wald <- rep(0,nfactor)
		
		fit$global_chisq <- waldtest(N=nvarRec,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarEnd,Ntot=ntot)

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

#================================> For the death
#========================= Test de Wald
	if (!missing(formula.terminalEvent)){
		if(length(vec.factordc) > 0){
			nfactor <- length(vec.factordc)
			p.walddc <- rep(0,nfactor)
			fit$global_chisq_d <- waldtest(N=nvarEnd,nfact=nfactor,place=ind.placedc,modality=occurdc,b=Beta,Varb=VarBeta,Lfirts=nvarRec,Ntot=ntot)
		 fit$dof_chisq_d <- occurdc
			fit$global_chisq.test_d <- 1
	# Calcul de pvalue globale
			for(i in 1:length(vec.factordc)){
				p.walddc[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurdc[i]), 3)
			}
			fit$p.global_chisq_d <- p.walddc
			fit$names.factordc <- vec.factordc
		}else{
			fit$global_chisq.test_d <- 0
		}
	}else{
		fit$global_chisq.test_d <- 0
	}


#================================> For the reccurrent 2
#========================= Test de Wald
	ntmp <- nvarRec+nvarEnd

	if(length(vec.factorevent2) > 0){
		nfactor <- length(vec.factorevent2)
		p.wald2 <- rep(0,nfactor)
			fit$global_chisq2 <- waldtest(N=nvarRec2,nfact=nfactor,place=ind.placeformula2,modality=occurformula2,b=Beta,Varb=VarBeta,Lfirts=ntmp,Ntot=ntot)
	 fit$dof_chisq2 <- occurformula2
		fit$global_chisq.test2 <- 1
# Calcul de pvalue globale
		for(i in 1:length(vec.factorevent2)){
			p.wald2[i] <- signif(1 - pchisq(fit$global_chisq2[i], occurformula2[i]), 3)
		}
		fit$p.global_chisq2 <- p.wald2
		fit$names.factor2 <- vec.factorevent2

		
	}else{
		fit$global_chisq.test2 <- 0
	}


	class(fit) <- "multivPenal"


	if(print.times){
		cost<-proc.time()-ptm
		cat("The program took", round(cost[3],2), "seconds \n")
	}
 fit

}



