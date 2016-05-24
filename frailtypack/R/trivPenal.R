

"trivPenal" <-
  function (formula, formula.terminalEvent, formula.LongitudinalData, data,  data.Longi, random, id, intercept = TRUE, link="Random-effects",
            left.censoring=FALSE, recurrentAG=FALSE, n.knots, kappa,
            maxit=350, hazard="Splines", init.B,
            init.Random, init.Eta, init.Alpha, method.GH = "Standard", n.nodes, LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, print.times=TRUE)
  {

    m3 <- match.call() # longitudinal
    m3$formula <- m3$formula.terminalEvent <- m3$data <- m3$recurrentAG <- m3$random <- m3$id <- m3$link <- m3$n.knots <- m3$kappa <- m3$maxit <- m3$hazard <- m3$init.B <- m3$LIMparam <- m3$LIMlogl <- m3$LIMderiv <- m3$print.times <- m3$left.censoring <- m3$init.Random <- m3$init.Eta <- m3$init.Alpha <- m3$method.GH <- m3$intercept <- m3$n.nodes <- m3$... <- NULL
    Names.data.Longi <- m3$data.Longi

    m2 <- match.call() #terminal
    m2$formula <- m2$formula.terminalEvent <- m2$formula.LongitudinalData <- m2$data.Longi <- m2$recurrentAG <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard  <-  m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$init.Alpha <- m2$method.GH <- m2$intercept <- m2$n.nodes <- m2$... <- NULL
    Names.data.Terminal <- m2$data


    #### Frailty distribution specification ####
    if (!(all(random %in% c("1",names(data.Longi))))) { stop("Random effects can be only related to variables from the longitudinal data or the intercept (1)") }
    if (!(id %in% c(names(data.Longi))) || !(id %in% c(1,names(data)))) { stop("Identification for individuals can be only related to variables from both data set") }

    #### Link function specification ####
    if(!(link %in% c("Random-effects","Current-level"))){
      stop("Only 'Random-effects' or 'Current-level' link function can be specified in link argument.")}

    ### Left-censoring
    if(!is.null(left.censoring) && left.censoring!=FALSE){
      if(!is.numeric(left.censoring))stop("If you want to include left-censored longitudinal outcome you must give the threshold value as the argument of 'left.censoring'")
    }
    ### Intercept
    if(!is.logical(intercept))stop("The argument 'intercept' must be logical")
    ### Gauss-Hermite method
    if(all(!(c("Standard","Pseudo-adaptive","HRMSYM") %in% method.GH))){
      stop("Only 'Standard', 'Pseudo-adaptive' and 'HRMSYM' hazard can be specified as a method for the Gaussian quadrature")
    }
    GH <- switch(method.GH,"Standard"=0,"Pseudo-adaptive"=1,"HRMSYM"=2)

    if(!missing(n.nodes) ){
     if(!n.nodes%in%c(5,7,9,12,15,20,32)) stop("Number of points used in the numerical integration must be chosen from following: 5, 7, 9, 12, 15, 20, 32")
     if(n.nodes%in%c(5,7,9,12,15,20,32) && GH==2) warning("Using HRMSYM algorithm the number of points cannot be chosen")
     }else{
      n.nodes <- 9
     }

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
         
          size1 <- 100
          size2 <- 100
          equidistant <- 1

        }
        ### Weibull
        if (typeof == 2){
         
          size1 <- 100
          size2 <- 100
          equidistant <- 2

        }

      }
    }else{
      #### longueur hazard > 1
      if(all(!(c("Splines") %in% hazard))){
        stop("Only 'Splines'  hazard can be specified in hazard argument in this case")
      }
      ### Splines percentile
      if ("Splines" %in% hazard){
        typeof <- 0
        if(!(all(hazard %in% c("Splines","per")))){
          stop ("The hazard argument is incorrectly specified. Only 'per' is allowed with 'Splines'. Please refer to the help file of frailtypack.")
        }else{
          size1 <- 100
          size2 <- 100
          equidistant <- 0

        }
      }
    }


    if (missing(formula))stop("The argument formula must be specified in every model")
    if (missing(formula.LongitudinalData))stop("The argument formula.LongitudinalData must be specified in every model") #AK
    if (missing(formula.terminalEvent))stop("The argument formula.terminalEvent must be specified in every model") #AK

    if(class(formula)!="formula")stop("The argument formula must be a formula")

    if(typeof == 0){
      if (missing(n.knots))stop("number of knots are required")
      n.knots.temp <- n.knots
      if (n.knots<4) n.knots<-4
      if (n.knots>20) n.knots<-20
      if (missing(kappa))stop("smoothing parameter (kappa) is required")

    }else{
      if (!(missing(n.knots)) || !(missing(kappa)) ){
        stop("When parametric hazard function is specified, 'kappa', 'n.knots' arguments must be deleted.")
      }
      n.knots <- 0
      kappa <- 0

    }
    call <- match.call()

    m <- match.call(expand.dots = FALSE) # recurrent events

    m$formula.LongitudinalData <- m$formula.terminalEvent <- m$recurrentAG <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard   <-  m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$init.Alpha <- m$method.GH <- m$intercept <- m$n.nodes <- m$... <- NULL


    special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

    Terms <- if (missing(data)){
      terms(formula, special)
    }else{
      terms(formula, special, data = data)
    }

    ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli

  #  if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
    #si pas vide tous si il ya au moins un qui vaut 1 on arrete

    m$formula <- Terms



    m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait

    #model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) +
    #as.factor(charlson) + sex + chemo + terminal(death), data = readmission)

    m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument

    cluster <- attr(Terms, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()


    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

    # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
    classofR <- attr(model.extract(m, "response"),"class")
    # attention le package pec rajoute un element dans l'attribut "class" des objets de survie
    if (length(classofR)>1) classofR <- classofR[2]

    typeofR <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..
  
    #Al : tri du jeu de donnees par cluster croissant
    if (length(cluster)){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms]
    #  if (any(ord > 1))stop("Cluster can not be used in an interaction")
      m <- m[order(m[,tempc$vars]),] # soit que des nombres, soit des caracteres
      ordre <- as.integer(row.names(m)) # recupere l'ordre du data set
      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
      uni.cluster <- unique(cluster)
    }

    # verification de la sutructure nested si besoin
    if (length(subcluster)) stop("'subcluster' is not an allowed option")


    if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

    R <- model.extract(m, "response") # objet de type Surv =Time


      if (classofR != "Surv") stop("Response in the 'formula' must be a survival object")


    ll <- attr(Terms, "term.labels")#liste des variables explicatives

    #cluster(id) as.factor(dukes) as.factor(charlson) sex chemo terminal(death)

    #=========================================================>

    mt <- attr(m, "terms") #m devient de class "formula" et "terms"

    X <- if (!is.empty.model(mt))model.matrix(mt, m, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables

    ind.place <- unique(attr(X,"assign")[duplicated(attr(X,"assign"))]) ### unique : changement au 25/09/2014


    vec.factor <- NULL
    vec.factor <- c(vec.factor,ll[ind.place])

    #=========================================================>
    # On determine le nombre de categorie pour chaque var categorielle
    strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
    cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
    num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
    vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()

    #booleen pour savoir si au moins une var depend du tps
    if (is.null(vartimedep)) timedep <- 0
    else timedep <- 1

    if (timedep==1) stop("The option 'timedep' is not allowed in this model.")
    
    if(!is.null(num.id)) stop("num.id is not an allowed option")


    subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

    if (length(subcluster))stop("'subcluster' is not an allowed option")
    if (length(strats))stop("Stratified analysis is not allowed")

    if (length(subcluster)){
      ll <- ll[-grep("subcluster",ll)]
    }
    if (length(cluster)){
      ll_tmp <- ll[grep("cluster",ll)]
      ll <- ll[-grep("cluster",ll)]
      
      pos1 <- grep("r",unlist(strsplit(ll_tmp,split="")))[1]+2
      pos2 <- length(unlist(strsplit(ll_tmp,split="")))-1
      Names.cluster <- substr(ll_tmp,start=pos1,stop=pos2) # nom du cluster
    }
    if (length(strats)){
      ll <- ll[-grep("strata",ll)]
    }
    
    #   plus besoin de as.factor() pour afficher le test de Wald global
    if (length(grep("strata",vec.factor))) vec.factor <- vec.factor[-grep("strata",vec.factor)]
    if (length(grep("cluster",vec.factor))) vec.factor <- vec.factor[-grep("cluster",vec.factor)]
    if (length(grep("subcluster",vec.factor))) vec.factor <- vec.factor[-grep("subcluster",vec.factor)]
    if (length(grep("num.id",vec.factor))) vec.factor <- vec.factor[-grep("num.id",vec.factor)]
    

  
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
  
  
  if(length(grep("terminal",ll))>0){ind.place <- grep(paste(vec.factor,collapse="|"),ll[-grep("terminal",ll)])
  }else{ind.place <- grep(paste(vec.factor,collapse="|"),ll)}
  
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
  


    #=========================================================>

    terminalEvent <- attr(Terms, "specials")$terminal #nbre de var qui sont en fonction de terminal()

    dropx <- NULL

    if (length(cluster) ){
      tempc <- untangle.specials(Terms, "cluster", 1:10)
      ord <- attr(Terms, "order")[tempc$terms]
      if (any(ord > 1))stop("Cluster can not be used in an interaction")

      cluster <- strata(m[, tempc$vars], shortlabel = TRUE)

      dropx <- tempc$terms
      uni.cluster<-unique(cluster)

    }else if (!length(cluster)){

      stop("grouping variable is needed in the formula")

    }

    if (length(num.id))stop("'num.id' is not an allowed option")

    if(length(uni.cluster)==1){
      stop("grouping variable must have more than 1 level")
    }


    if (length(subcluster))stop("subcluster is not an allowed option")



    if (length(strats))stop("Stratified analysis is not allowed")

    if (typeof==0 && (length(kappa)!=2)) stop("wrong length of argument 'kappa'. The length should be 2")

    #AD: indicator of terminal()
    ind.terminal <- length(terminalEvent)

    #AD:
    if (length(terminalEvent)){

      tempterm <- untangle.specials(Terms, "terminal", 1:10)
      #ici on comme terme tempterm$vars qui est le nom dans l'appel(ex;"terminal(death)"
      #et tempterm$terms qui est la position de la variable dans l'appel, ici elle vient a la position 6

      ord <- attr(Terms, "order")[tempterm$terms] # ord[6]=1 ici dans notre exemple

      if (any(ord > 1))stop("Terminal can not be used in an interaction")
      dropx <- c(dropx,tempterm$terms) # vecteur de position
      terminal <- strata(m[, tempterm$vars], shortlabel = TRUE)
      terminal <- as.numeric(as.character(terminal))

    }

    #type <- attr(Y, "type")
    type <- typeofR

    if (type != "right" && type != "counting") { # Cox supporte desormais la censure par intervalle
      stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
    }

    if (type != "counting" && recurrentAG) {
      stop("recurrentAG needs counting process formulation")
    }

    if (length(dropx)){
      newTerms <- Terms[-dropx]
    }else{
      newTerms <- Terms
    }

    #newTerm vaut Terms - les variables dont les position sont dans drop

    X <- model.matrix(newTerms, m)

    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
    Xlevels <- .getXlevels(newTerms, m)
    contr.save <- attr(X, 'contrasts')


    # assigne donne la position pour chaque variables
    #ncol(X) : nombre de variable sans sans les fonction speciaux comme terminal()...+id
    if(length(vec.factor) > 0){
      #========================================>
      position <- unlist(assign,use.names=F)
    }

    #========================================>

    if (ncol(X) == 1){
      X<-X-1
      noVarR <- 1
    }else{
      X <- X[, -1, drop = FALSE]
      noVarR <- 0
    }
    # on enleve ensuite la premiere colonne correspondant a id


    nvar<-ncol(X) #nvar==1 correspond a 2 situations:

    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0

    #varnotdep <- colnames(X)[-grep("timedep",colnames(X))]
    #vardep <- colnames(X)[grep("timedep",colnames(X))]
    #vardep <- apply(matrix(vardep,ncol=1,nrow=length(vardep)),1,timedep.names)

    #if (length(intersect(varnotdep,vardep)) != 0) {
    #  stop("A variable is both used as a constant and time-varying effect covariate")
    #}

    #nvartimedep <- length(vardep)

    #filtretps <- rep(0,nvar)
    #filtretps[grep("timedep",colnames(X))] <- 1

    var<-matrix(c(X),nrow=nrow(X),ncol=nvar) #matrix sans id et sans partie ex terminal(death)

    nsujet<-nrow(X)

    if (type=="right"){
      tt0 <- rep(0,nsujet)
      tt1 <- R[,1]
      cens <- R[,2]
    } else {
      tt0 <- R[,1]
      tt1 <- R[,2]
      cens <- R[,3]
    }

    if (min(cens)==0) cens.data<-1
    if (min(cens)==1 && max(cens)==1) cens.data<-0

    AG<-ifelse(recurrentAG,1,0)


    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factor) > 0){
      #         ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factor)){
        ind.place[i] <- ind.place[i]+k
        k <- k + occur[i]-1
      }
    }



    #========= Longitudinal Data preparation =========================

    TermsY <- if (missing(data.Longi)){
      terms(formula.LongitudinalData, special)
    }else{
      terms(formula.LongitudinalData, special, data = data.Longi)
    }

    ord <- attr(TermsY, "order") # longueur de ord=nbre de var.expli

    #si pas vide tous si il ya au moins un qui vaut 1 on arrete

    m3$formula <- TermsY


    m3[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait


    if (NROW(m3) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

    llY <- attr(TermsY, "term.labels")#liste des variables explicatives


    #=========================================================>

    name.Y <- as.character(attr(TermsY, "variables")[[2]])
    Y <- data.Longi[,which(names(data.Longi)==name.Y)]
 
    
  # on identifie les variables explicatives facteurs avec nombre de niveau plus que 2

  ind.placeY <- which(llY%in%names(which(lapply(data.Longi[,which(names(data.Longi)%in%llY)],function(x) length(levels(x)))>2)))
 
  defined.factor <- llY[grep("factor",llY)]
  
  
  vec.factorY.tmp <- NULL
  if(length(defined.factor)>0){
    mat.factorY.tmp <- matrix(defined.factor,ncol=1,nrow=length(defined.factor))
    
    # Fonction servant a prendre les termes entre "as.factor"
    vec.factorY.tmp <-apply(mat.factorY.tmp,MARGIN=1,FUN=function(x){
      if (length(grep("factor",x))>0){
        if(length(grep(":",x))>0){
          if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
            
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos4 <- length(unlist(strsplit(x,split="")))
            if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
          }else if(grep("\\(",unlist(strsplit(x,split="")))[1]>grep(":",unlist(strsplit(x,split="")))[1] && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[1]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
            }else{#both factors
            pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
            pos2 <- grep(":",unlist(strsplit(x,split="")))[1]-2
            pos3 <- grep("\\(",unlist(strsplit(x,split="")))[2]+1
            pos4 <- length(unlist(strsplit(x,split="")))-1
            if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2 || length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos3,stop=pos4))])))>2)return(paste(substr(x,start=pos1,stop=pos2),":",substr(x,start=pos3,stop=pos4),sep=""))
            else return(NaN)
            }
        }else{
          pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
          pos2 <- length(unlist(strsplit(x,split="")))-1
          if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==substr(x,start=pos1,stop=pos2))])))>2)return(substr(x,start=pos1,stop=pos2))
          else return(NaN)
        }
      }else{
        return(x)
      }})
    
    vec.factorY.tmp <- vec.factorY.tmp[which(vec.factorY.tmp!="NaN")]
    
    if(length(vec.factorY.tmp)>0){
      for(i in 1:length(vec.factorY.tmp)){
        if(length(grep(":",vec.factorY.tmp[i]))==0){
          if(length(levels(as.factor(data.Longi[,which(names(data.Longi)==vec.factorY.tmp[i])])))>2)ind.placeY <- c(ind.placeY,which(llY%in%paste("as.factor(",vec.factorY.tmp[i],")",sep="")))
        }
        
      }}
    
  }

  ind.placeY <- sort(ind.placeY)
  #=========================================================>
  # On determine le nombre de categorie pour chaque var categorielle
  stratsY <- attr(TermsY, "specials")$strata #nbre de var qui sont en fonction de strata()
  clusterY <- attr(TermsY, "specials")$cluster #nbre de var qui sont en fonction de cluster()
  num.idY <- attr(TermsY, "specials")$num.id #nbre de var qui sont en fonction de patkey()
  vartimedepY <- attr(TermsY, "specials")$timedep #nbre de var en fonction de timedep()
  
  #booleen pour savoir si au moins une var depend du tps
  if (is.null(vartimedep)) timedepY <- 0
  else timedepY <- 1
  
  
  if (timedepY==1) stop("The option 'timedep' is not allowed in this model.")
  
  
  subclusterY <- attr(TermsY, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()
  
  
  if (length(clusterY))stop("Only the argument 'id' can represent the clusters")
  
  Names.clusterY <- id # nom du cluster
  
  if (length(num.idY))stop("'num.id' for the interval censoring is an allowed option")
  
  if (length(stratsY))stop("Stratified analysis is not an allowed option yet")
  
  
  
  # which_n<-which(names(data.Longi)%in%random)
  #  data.Longi[which(data.Longi[,which_n]==0),which_n]<- 0.01
  
  mat.factorY2 <- matrix(llY,ncol=1,nrow=length(llY))
  
  # Fonction servant a prendre les termes entre "as.factor"
  llY2 <-apply(mat.factorY2,MARGIN=1,FUN=function(x){
    if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
      pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      pos2 <- length(unlist(strsplit(x,split="")))-1
      x<-substr(x,start=pos1,stop=pos2)
      return(paste(x,levels(as.factor(data.Longi[,which(names(data.Longi)==x)]))[2],sep=""))
    }else{
      return(x)
    }})
  
  # Fonction servant a prendre les termes entre "as.factor" - without the name of the level
  llY3 <-apply(mat.factorY2,MARGIN=1,FUN=function(x){
    if (length(grep("factor",x))>0  && length(grep(":",x))==0 && unlist(strsplit(x,split=""))[length(unlist(strsplit(x,split="")))]==")"){
      pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
      pos2 <- length(unlist(strsplit(x,split="")))-1
      return(substr(x,start=pos1,stop=pos2))
    }else{
      return(x)
    }})
  
  llY.real.names <- llY3  
  llY3 <- llY3[!llY2%in%llY]
  
  
  if(is.factor(data.Longi[,names(data.Longi)==llY.real.names[1]]))X_L<- as.numeric(data.Longi[,names(data.Longi)==llY.real.names[1]])-1
  else X_L<- data.Longi[,names(data.Longi)==llY.real.names[1]]
  
  
  
  
  if(length(llY)>1){
    for(i in 2:length(llY.real.names)){
      
      if(is.factor(data.Longi[,names(data.Longi)==llY.real.names[i]]))X_L<- cbind(X_L,as.numeric(data.Longi[,names(data.Longi)==llY.real.names[i]])-1)
      else X_L<- cbind(X_L,data.Longi[,names(data.Longi)==llY.real.names[i]])
    }}
  
  #X_L<- data.Longi[,names(data.Longi)%in%(llY)]
  
  llY.fin <- llY.real.names
  llY <- llY.real.names
  
  if(sum(ord)>length(ord)){
    
    for(i in 1:length(ord)){
      if(ord[i]>1){
        
        name_v1 <- strsplit(as.character(llY[i]),":")[[1]][1]
        name_v2 <- strsplit(as.character(llY[i]),":")[[1]][2]
        
        if(length(grep("factor",name_v1))>0){name_v1<-substring(name_v1,11,nchar(name_v1)-1)
                                             v1 <- as.factor(data.Longi[,names(data.Longi)==name_v1])}
        else{v1 <- data.Longi[,names(data.Longi)==name_v1]}
        if(length(grep("factor",name_v2))>0){name_v2<-substring(name_v2,11,nchar(name_v2)-1)
                                             v2 <- as.factor(data.Longi[,names(data.Longi)==name_v2])}
        else{v2 <- data.Longi[,names(data.Longi)==name_v2]}
        
        llY[i] <- paste(name_v1,":",name_v2,sep="")
        #   if(is.factor(v1) && length(levels(v1))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
        #   if(is.factor(v2) && length(levels(v2))>2)stop("Interactions not allowed for factors with 3 or more levels (yet)")
        if(is.factor(v1) && !is.factor(v2)){
          
          dummy <- model.matrix( ~ v1 - 1)
          # if(length(levels(v1)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
          for(j in 2:length(levels(v1))){
            X_L <- cbind(X_L,dummy[,j]*v2)
            if(i>1 && i<length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llY.fin[(i+1+j-2):length(llY.fin)])
            else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""))
            else llY.fin <- c(paste(name_v1,".",levels(v1)[j],":",name_v2,sep=""),llY.fin[(2+j-2):length(llY.fin)])
          }
          
        }else if(!is.factor(v1) && is.factor(v2)){
          
          dummy <- model.matrix( ~ v2 - 1)
          #  if(length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
          for(j in 2:length(levels(v2))){
            
            X_L <- cbind(X_L,dummy[,j]*v1)
            
            if(i>1 && i<length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llY.fin[(i+1+j-2):length(llY.fin)])
            else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2)],paste(name_v1,":",name_v2,levels(v2)[j],sep=""))
            else llY.fin <- c(paste(name_v1,":",name_v2,levels(v2)[j],sep=""),llY.fin[(2+j-2):length(llY.fin)])
          }
        }else if(is.factor(v1) && is.factor(v2)){
          
          
          dummy1 <- model.matrix( ~ v1 - 1)
          dummy2 <- model.matrix( ~ v2 - 1)
          #   if(length(levels(v1)>2) || length(levels(v2)>2))vec.factorY <- c(vec.factorY,paste(name_v1,":",name_v2,sep=""))
          for(j in 2:length(levels(v1))){
            for(k in 2:length(levels(v2))){
              
              X_L <- cbind(X_L,dummy1[,j]*dummy2[,k])
              if(i>1 && i<length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llY.fin[(i+1+j-2+k-2):length(llY.fin)])
              else if(i==length(llY.fin))llY.fin <- c(llY.fin[1:(i-1+j-2+k-2)],paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""))
              else llY.fin <- c(paste(name_v1,levels(v1)[j],":",name_v2,levels(v2)[k],sep=""),llY.fin[(2+j-2+k-2):length(llY.fin)])
            }
          } 
        }else{
          
          X_L <- cbind(X_L,v1*v2)
        }
        
      }
    }
  }
  
  
  
  if(length(grep(":",llY))>0){
    for(i in 1:length(grep(":",llY))){
      if(length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llY[grep(":",llY)[i]],":")[[1]])[1]]))>2 || length(levels(data.Longi[,which(names(data.Longi)%in%strsplit(llY[grep(":",llY)[i]],":")[[1]])[2]]))>2){
        ind.placeY <- c(ind.placeY,grep(":",llY)[i])
        #     vec.factorY <- c(vec.factorY,llY[grep(":",llY)[i]])
      }
    }
  }
  
  vec.factorY <- NULL

  if(length(vec.factorY.tmp)>0)vec.factorY <- c(llY[ind.placeY],vec.factorY.tmp)
  else vec.factorY <- c(vec.factorY,llY[ind.placeY])
  
  vec.factorY <- unique(vec.factorY)
  
  

  
  mat.factorY <- matrix(vec.factorY,ncol=1,nrow=length(vec.factorY))
  # Fonction servant a prendre les termes entre "as.factor" et (AK 04/11/2015) interactions
  vec.factorY <-apply(mat.factorY,MARGIN=1,FUN=function(x){
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
  
  
  for(i in 1:length(llY.fin)){
    
    if(sum(names(data.Longi)==llY.fin[i])>0){
      if(is.factor(data.Longi[,names(data.Longi)==llY.fin[i]]) && length(levels(data.Longi[,names(data.Longi)==llY.fin[i]]))==2){
        llY.fin[i] <- paste(llY.fin[i],levels(data.Longi[,names(data.Longi)==llY.fin[i]])[2],sep="")}
    }
  }
  
  #  llY <- llY.fin
  if(dim(X_L)[2]!=length(llY.fin))stop("The variables in the longitudinal part must be in the data.Longi")
  X_L <- as.data.frame(X_L)
  names(X_L) <- llY.fin
  
  Intercept <- rep(1,dim(X_L)[1])
  
  if(intercept){
    X_L <- cbind(Intercept,X_L)
    ind.placeY <- ind.placeY+1
  }

  
  X_Lall<- X_L
  "%+%"<- function(x,y) paste(x,y,sep="")
  if(length(vec.factorY) > 0){
    for(i in 1:length(vec.factorY)){
      if(length(grep(":",vec.factorY[i]))==0){
        
        factor.spot <- which(names(X_L)==vec.factorY[i])
       
        X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1],X_L[(factor.spot+1):ncol(X_L)])
      } }
    
 
	
	
	
    vect.factY<-names(X_L)[which(!(names(X_L)%in%llY))]
    if(intercept) vect.factY <- vect.factY[-1]
    
    
    #               vect.fact <-apply(matrix(vect.fact,ncol=1,nrow=length(vect.fact)),MARGIN=1,FUN=function(x){
    #               pos1 <- grep("r",unlist(strsplit(x,split="")))[1]+2
    #               pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
    #               return(substr(x,start=pos1,stop=pos2))})
    
    occurY <- rep(0,length(vec.factorY))
    
    #         for(i in 1:length(vec.factorY)){
    #                #occur[i] <- sum(vec.factor[i] == vect.fact)
    #               occurY[i] <- length(grep(vec.factorY[i],vect.factY))
    #      }
    
    
    
    interaction<-as.vector(apply(matrix(vect.factY,nrow=length(vect.factY)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
    which.interaction <- which(interaction==1)
   
    for(i in 1:length(vec.factorY)){
      
      if(length(grep(":",unlist(strsplit(vec.factorY[i],split=""))))>0){
        
        
        pos <- grep(":",unlist(strsplit(vec.factorY[i],split="")))
        length.grep <- 0
        for(j in 1:length(vect.factY)){
          if(j%in%which.interaction){
            
            if(length(grep(substr(vec.factorY[i],start=1,stop=pos-1),vect.factY[j]))>0 && length(grep(substr(vec.factorY[i],start=pos+1,stop=length(unlist(strsplit(vec.factorY[i],split="")))),vect.factY[j]))>0){
              length.grep <- length.grep + 1
              which <- i}
          }}
        occurY[i] <- length.grep
        
      }else{
        
        
        if(length(vect.factY[-which.interaction])>0){occurY[i] <- length(grep(vec.factorY[i],vect.factY[-which.interaction]))
        }else{occurY[i] <- length(grep(vec.factorY[i],vect.factY))}
      }
    }
  }
  

    if (ncol(X_L) == 0){
      noVarY <- 1
    }else{
      noVarY <- 0
    }



    #=========================================================>

    clusterY <- data.Longi$id
    max_rep <- max(table(clusterY))
    uni.clusterY<-as.factor(unique(clusterY))


    if(is.null(id))     stop("grouping variable is needed")

    if(is.null(random)) stop("variable for random effects is needed")


    if(length(uni.clusterY)==1){
      stop("grouping variable must have more than 1 level in the longitudinal part")
    }


    if (length(subcluster))stop("'Subcluster' is not allowed")



    if (typeof==0 && missing(kappa)) stop("smoothing parameter (kappa) is required")




    #newTerm vaut Terms - les variables dont les position sont dans drop






    #========================================>


    nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:

    # au cas ou on a aucune var explicative dans la partie rec, mais X=0
    # cas ou on a 1seul var explicative, ici X est en general different de 0

    #varnotdepY <- colnames(X_L)[-grep("timedep",colnames(X_L))]
    #vardepY <- colnames(X_L)[grep("timedep",colnames(X_L))]
    #vardepY <- apply(matrix(vardepY,ncol=1,nrow=length(vardepY)),1,timedep.names)

    #if (length(intersect(varnotdepY,vardepY)) != 0) {
    #  stop("A variable is both used as a constant and time-varying effect covariate")
    #}

    #nvartimedepY <- length(vardepY)

    #filtretpsY <- rep(0,nvarY)
    #filtretpsY[grep("timedep",colnames(X_L))] <- 1

    varY <- as.matrix(sapply(X_L, as.numeric))

    nsujety<-nrow(X_L)




    #=======================================>
    #======= Construction du vecteur des indicatrice
    if(length(vec.factorY) > 0){
      #         ind.place <- ind.place -1
      k <- 0
      for(i in 1:length(vec.factorY)){
        ind.placeY[i] <- ind.placeY[i]+k
        k <- k + occurY[i]-1
      }
    }

    # Random effects

    if(link=="Random-effects") link0 <- 1
    if(link=="Current-level") link0 <- 2

    nRE <- length(random)

    ne_re <- nRE*(nRE+1)/2

    matzy <- NULL
    names.matzy <- NULL
    if("1"%in%random){
      names.matzy<-c("Intercept",random[-which(random=="1")])
    }else{
      names.matzy<-random
    }

    matzy <- data.matrix(X_Lall[,which(names(X_Lall)%in%names.matzy)])
    if(!intercept && 1%in%random) matzy <- as.matrix(cbind(rep(1,nsujety),matzy))

    if(link0==1){netadc <- ncol(matzy)
                 netar <- ncol(matzy)}
    if(link0==2){netadc <- 1
                 netar <- 1}


    #== Left-censoring ==
    cag <- c(0,0)
    if(!is.null(left.censoring) && left.censoring!=FALSE){
      if(left.censoring<min(Y))stop("The threshold for the left censoring cannot be smaller than the minimal value of the longitudinal outcome")
      cag[1] <- 1
      cag[2] <- left.censoring
      n.censored <- length(which(Y<=left.censoring))
      prop.censored <- n.censored/nsujety
    }
    #============= pseudo-adaptive Gauss Hermite ==============
    #m <- lme(measuret ~ time+interact+treatment, data = data, random = ~ 1| idd)
    inn <-paste("pdSymm(form=~",random,")",sep="")
    rand <- list(eval(parse(text=inn)))
    names(rand) <- id

    m_lme<-lme(formula.LongitudinalData,data = data.Longi, random = rand)

    b_lme <-as.matrix(ranef(m_lme))

    formula_lme <- formula(m_lme$modelStruct$reStruct[[1]])
    model_lme <- model.frame(terms(formula_lme), data = data.Longi)
    Terms_lme <- attr(model_lme, "terms")
    Z_lme <- model.matrix(formula_lme, model_lme)
    #id <- as.vector(unclass(m2$groups[[1]]))
    # Cholesky matrices for GH
    #B_lme <- lapply(pdMatrix(m_lme$modelStruct$reStruct), "*",  m_lme$sigma^2)[[1]]


    #invBi_chol <- lapply(Bi_lme, function (x) solve(chol(solve(x))))


    ### model - frailtyPenal
#    m_penal <- match.call() #terminal
#    m_penal$formula.LongitudinalData <- m_penal$data.Longi <-   m_penal$random <-  m_penal$id <-  m_penal$link <-   m_penal$LIMparam <-  m_penal$LIMlogl <-  m_penal$LIMderiv <-  m_penal$print.times <-  m_penal$left.censoring <-  m_penal$init.Random <-  m_penal$init.Eta <-  m_penal$init.Alpha <- m_penal$init.B <- m_penal$method.GH <-  m_penal$intercept <-  m_penal$... <- NULL
#  print(paste(m_penal$formula))
#    print(paste(m_penal$formula)[[2]])

#formu<-paste("frailtyPenal(",paste(m_penal$formula)[[2]],"~",paste(m_penal$formula)[[3]],",formula.terminalEvent=~",paste(m_penal$formula.terminalEvent)[[2]],",data=data, recurrentAG=",
#        m_penal$recurrentAG,", hazard='Weibull',RandDist='LogN',print.times=FALSE)" )

 #  f<-eval(parse(text=formu))


#u_paGH <-cbind(b_lme,f$frailty.pred)

  B_lme <-pdMatrix(m_lme$modelStruct$reStruct)[[1]]
 # B_paGH <- cbind(B_lme,rep(0,dim(B_lme)[1]))
 # B_paGH <- rbind(B_paGH,rep(0,dim(B_paGH)[2]))
#B_paGH[dim(B_paGH)[1],dim(B_paGH)[2]] <- f$sigma2
#print(B_paGH)

  Bi_tmp  <- vector("list", length(uni.cluster) )
invBi_chol <- matrix(rep(0,ne_re*length(uni.cluster)),nrow=length(uni.cluster) )

invB_lme  <- solve(B_lme)
invB_chol<- solve(chol(solve(B_lme)))
invB_cholDet <- det(invB_chol)
for (i in 1:length(uni.cluster )) {
  Zi_lme <- Z_lme[clusterY  == i, , drop = FALSE]

  Bi_tmp[[i]] <- chol(solve(crossprod(Zi_lme) / m_lme$sigma^2 + invB_lme))
  for(j in 1:(nRE)){
    for(k in 1:(nRE)){
      if (k>=j)invBi_chol[i,j+k*(k-1)/2] <-  Bi_tmp[[i]][j,k]
      #  else     inv.chol.VC(j,k)=matv[k+j*(j-1)/2]
    }}
}
#Vs[[i]] <- solve(crossprod(Z.i) / m2$sigma^2 + inv.VC)
#}
invBi_cholDet <- sapply(Bi_tmp,  det)

   #==============================================
    #=== preparation survival data ===============

    if (!recurrentAG)
    {

    tt1T<-aggregate(tt1,by=list(cluster),FUN=sum)[,2]
    tt0T<-rep(0,length(tt1T))
    clusterT <- 0
    ligneT <- 0
    tempT <- 0

    }else{
      tt1T<-aggregate(tt1,by=list(cluster),FUN=function(x) x[length(x)])[,2]
      tt0T<-rep(0,length(tt1T))
      clusterT <- 0
      ligneT <- 0
      tempT <- 0

    }

    TermsT <- terms(formula.terminalEvent, special, data = data)


    #AD:
    if (!missing(formula.terminalEvent)){
      ord2 <- attr(TermsT, "order")

   #   if (length(ord2) & any(ord2 != 1)){
   #     stop("Interaction terms are not valid for terminal event formula")
   #   }
    }

    #AD:Joint model needs "terminal()"
    if (!ind.terminal)stop(" Joint frailty model miss specified ")

    terminalEvent<-aggregate(terminal,by=list(cluster),FUN=function(x) x[length(x)])[,2]
    # terminalEvent might be 0-1
    if (!all(terminalEvent%in%c(1,0)))  stop("terminal must contain a variable coded 0-1 and a non-factor variable")

    m2$formula <- TermsT

    m2[[1]] <- as.name("model.frame")
    m2 <- eval(m2, sys.parent()) #ici il prend les colonne associe au var.exp, si rien il prend tout


    match.noNA<-dimnames(m2)[[1]]%in%dimnames(m)[[1]]#masque logique pour ne pas avoir de NA

    m2<-m2[match.noNA, ,drop=FALSE]#m2 inchanger si pas de NA

    if (!missing(formula.terminalEvent))newTermsT<-TermsT

    #=========================================================>
    if (!missing(formula.terminalEvent)){
      X_T <- model.matrix(newTermsT, m2)

      llT <- attr(newTermsT,"term.labels")
      #ind.placedc <- grep("factor",lldc)
      ind.placeT <- unique(attr(X_T,"assign")[duplicated(attr(X_T,"assign"))])#changement unique le 26/09/2014
      
      vec.factorT <- NULL
      vec.factorT <- c(vec.factorT,llT[ind.placeT])
      
      mat.factorT <- matrix(vec.factorT,ncol=1,nrow=length(vec.factorT))
      # Fonction servant a prendre les termes entre "as.factor"
      vec.factorT <-apply(mat.factorT,MARGIN=1,FUN=function(x){
        if (length(grep("factor",x))>0){
          if(length(grep(":",x))>0){
            if(grep('\\(',unlist(strsplit(x,split="")))[1]<grep(":",unlist(strsplit(x,split="")))[1]  && length(grep('\\(',unlist(strsplit(x,split=""))))==1){
              
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
      
      # On determine le nombre de categorie pour chaque var categorielle
      if(length(vec.factorT) > 0){
        vect.factT <- attr(X_T,"dimnames")[[2]]
        vect.factT <- vect.factT[grep(paste(vec.factorT,collapse="|"),vect.factT)]
        
        occurT <- rep(0,length(vec.factorT))
  
        
        interactionT<-as.vector(apply(matrix(vect.factT,nrow=length(vect.factT)),MARGIN=1,FUN=function(x){length(grep(":",unlist(strsplit(x,split=""))))}))
        which.interactionT <- which(interactionT==1)
        
        for(i in 1:length(vec.factorT)){
          
          if(length(grep(":",unlist(strsplit(vec.factorT[i],split=""))))>0){
            
            
            pos <- grep(":",unlist(strsplit(vec.factorT[i],split="")))
            length.grep <- 0
            for(j in 1:length(vect.factT)){
              if(j%in%which.interactionT){
                if(length(grep(substr(vec.factorT[i],start=1,stop=pos-1),vect.factT[j]))>0 && length(grep(substr(vec.factorT[i],start=pos+1,stop=length(unlist(strsplit(vec.factorT[i],split="")))),vect.factT[j]))>0){
                  length.grep <- length.grep + 1
                  which <- i}
              }}
            occurT[i] <- length.grep
            
          }else{
            
            
            if(length(vect.factT[-which.interactionT])>0){occurT[i] <- length(grep(vec.factorT[i],vect.factT[-which.interactionT]))
            }else{occurT[i] <- length(grep(vec.factorT[i],vect.factT))}
          }
        }
      }
      #=========================================================>
      assign <- lapply(attrassign(X_T, newTermsT)[-1], function(x) x - 1)
      XlevelsT <- .getXlevels(newTermsT, m2)
      contr.saveT <- attr(X_T, 'contrasts')
      #========================================>
      if(length(vec.factorT) > 0){
        positionT <- unlist(assign,use.names=F)
      }
      #========================================>
      if (ncol(X_T) == 1)
      {
        X_T<-X_T-1
        noVarT <- 1
      }else{
        X_T <- X_T[, -1, drop = FALSE]
        noVarT <- 0
      }

      nvarT <- ncol(X_T)

      vartimedepT <- attr(TermsT, "specials")$timedep #nbre de var en fonction de timedep()

      # verifier qu'il y ait du timedep dans la deuxieme formule
      if (!is.null(vartimedepT)) timedepT <- 1
      else timedepT <- 0
      if (timedepT==1) stop("The option 'timedep' is not allowed in this model.")
      
      #varnotdepT <- colnames(X_T)[-grep("timedep",colnames(X_T))]
      #vardepT <- colnames(X_T)[grep("timedep",colnames(X_T))]
      #vardepT <- apply(matrix(vardepT,ncol=1,nrow=length(vardepT)),1,timedep.names)

      #if (length(intersect(varnotdepT,vardepT)) != 0) {
      #  stop("A variable is both used as a constant and time-varying effect covariate in the formula of terminal event")
      #}

      #nvartimedepT <- length(vardepT)

      #filtretpsT <- rep(0,nvarT)
      #filtretpsT[grep("timedep",colnames(X_T))] <- 1

      varT.temp<-matrix(c(X_T),nrow=nrow(X_T),ncol=nvarT)


      if(is.null(nrow(m2)))
      {
        if (length(m2) != nrow(m)){
          stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }
      }else{

        if (nrow(m2) != nrow(m)){
          stop(" There are missing values in the covariates modelling the terminal event. \n Prepare data only with complete cases")
        }

      }


      if (!is.null(ncol(varT.temp))){
        varT<-aggregate(varT.temp[,1],by=list(cluster), FUN=function(x) x[length(x)])[,2]

        if (ncol(varT.temp)>1){

          for (i in 2:ncol(varT.temp)){
            varT.i<-aggregate(varT.temp[,i],by=list(cluster), FUN=function(x) x[length(x)])[,2]
            varT<-cbind(varT,varT.i)
          }
        }
      }else{
        varT<-aggregate(varT.temp,by=list(cluster), FUN=function(x) x[length(x)])[,2]
      }
    }else{
      noVarT <- 1
      varT<-0
    }

   #    cluster <- aggregate(cluster,by=list(cluster),FUN=function(x) x[1])[,2]
    nvarR<-nvar

    if (!missing(formula.terminalEvent)){
      #=======================================>
      #======= Construction du vecteur des indicatrice

      if(length(vec.factorT) > 0){
        k <- 0
        for(i in 1:length(vec.factorT)){
          ind.placeT[i] <- ind.placeT[i]+k
          k <- k + occurT[i]-1
        }
      }

      #==================================

        if(is.null(nrow(varT))){
          nvarEnd<-1
        }else{
          nvarEnd<-ncol(varT)
        }
    }else{
      nvarEnd<-0
    }

    if (sum(as.double(var))==0) nvarR <- 0
    if ( sum(as.double(varT))==0) nvarEnd <- 0

   ng<-length(uni.cluster)

   # ... end preparing data

    #
    # Begin JOINT MODEL
    #


    nvar = nvarR + nvarY + nvarT

    if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0

    if (sum(as.double(varT))==0) nvarT <- 0
    if (sum(as.double(varY))==0) nvarY <- 0

   #if (timedep==0){
  #   npbetatps1 <- 0
   #  npbetatps2 <- 0
  #   npbetatps3 <- 0
   #}else{
   ##  npbetatps1 <- (betaknots+betaorder-1)*nvartimedep
   #  npbetatps2 <- (betaknots+betaorder-1)*nvartimedepT
   #  npbetatps3 <- (betaknots+betaorder-1)*nvartimedepY
   #}
   #npbetatps <- npbetatps1 + npbetatps2 + npbetatps3

   effet <- 1
   indic.alpha <- 1

    np <- switch(as.character(typeof),
                 "0"=(2*(as.integer(n.knots) + 2) + as.integer(nvar) + 1 + ne_re + netadc  + netar +   indic.alpha + effet),
                 "2"=(2*2 + nvar + 1  + ne_re + netadc  + netar + indic.alpha + effet))

   if (all(all.equal(as.numeric(cens),terminal)==T)){
     stop("'Recurrent event' variable and 'Terminal event' variable need to be different")
   }
    # traitement de l'initialisation du Beta rentre par l'utilisateur

    Beta <- rep(0.5,np)
    if (!missing(init.B)) {
      if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
     # if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
      Beta <- c(rep(0.5,np-nvar),init.B)
    }

    if (!missing(init.Random)) {
      if (length(init.Random)!=ne_re+1) stop("init.Random must be of length that corresponds to the number of elements to estimate of the B1 matrix")
      Beta[(np -nvar-ne_re+1):(np-nvar)] <- init.Random[-length(init.Random)]
      Beta[np -nvar-1-ne_re-netadc - netar - indic.alpha] <- init.Random[length(init.Random)]
    }
    if (!missing(init.Eta)) {

      if (length(init.Eta)!=netadc+netar) stop("init.Eta must be of length that corresponds to the dimension of the link function")
      Beta[(np -nvar-1-ne_re-netadc - netar +1):(np-nvar-1-ne_re)] <- init.Eta
    }

   if (!missing(init.Alpha)) {
    if (!is.numeric(init.Alpha)) stop("init.Alpha must be numeric")
     Beta[np -nvar-1-ne_re-netadc - netar ] <- init.Alpha
   }

    xSuT <- matrix(0,nrow=100,ncol=1)
   xSuR <- matrix(0,nrow=100,ncol=1)
    if (typeof==0){
      mt1 <- size1
      mt2 <- size2
    }else{
      mt1 <- 100
      mt2 <- 100
    }

   initialize <- 1

   npinit <- switch(as.character(typeof),
                    "0"=((n.knots + 2) + nvarR + effet),
                    "2"=(2 + nvarR + effet))

    flush.console()
    if (print.times){
      ptm<-proc.time()
      cat("\n")
      cat("Be patient. The program is computing ... \n")
    }


    ans <- .Fortran("joint_longi",
                    as.integer(nsujet),
                    as.integer(nsujety),
                    as.integer(ng),
                    as.integer(n.knots),
                    k0=as.double(kappa), # joint avec generalisation de strate
                    as.double(tt0),
                    as.double(tt1),
                    as.integer(cens),
                    as.integer(cluster),
                    as.double(tt0T),
                    as.double(tt1T),
                    as.integer(terminalEvent),
                    link0 = as.integer(c(link0,link0)),
                    yy0 = as.double(Y),
                    groupey0 = as.integer(clusterY),
                    nb0 = as.integer(nRE),
                    matzy0 =as.double(matzy),
                    cag0 = as.double(cag),
                    as.integer(nvarR),
                    as.double(var),
                    as.integer(nvarT),
                    as.double(varT),
                    nva30 = as.integer(nvarY),
                    vaxy0 = as.double(varY),
                    noVar = as.integer(c(noVarR,noVarT,noVarY)),
                    ag0 = as.integer(AG),
                    as.integer(maxit),
                    np=as.integer(np),
                    neta0 = as.integer(c(netadc,netar)),
                    b=as.double(Beta),
                    H=as.double(matrix(0,nrow=np,ncol=np)),
                    HIH=as.double(matrix(0,nrow=np,ncol=np)),

                    loglik=as.double(0),
                    LCV=as.double(rep(0,2)),
                    xR=as.double(matrix(0,nrow=size1,ncol=1)),
                    lamR=as.double(matrix(0,nrow=size1,ncol=3)),
                    xSuR=as.double(xSuR),
                    survR=as.double(matrix(0,nrow=mt1,ncol=3)),
                    xD=as.double(rep(0,size2)),
                    lamD=as.double(matrix(0,nrow=size2,ncol=3)),
                    xSuD=as.double(xSuT),
                    survD=as.double(matrix(0,nrow=mt2,ncol=3)),
                    as.integer(typeof),
                    as.integer(equidistant),
                    as.integer(c(size1,size2,mt1,mt2)),###
                    counts=as.integer(c(0,0,0)),
                    ier_istop=as.integer(c(0,0)),
                    paraweib=as.double(rep(0,4)),
                    MartinGale=as.double(matrix(0,nrow=ng,ncol=5)),###
                    ResLongi = as.double(matrix(0,nrow=nsujety,ncol=4)),
                    Pred_y  = as.double(matrix(0,nrow=nsujety,ncol=2)),

                    linear.pred=as.double(rep(0,nsujet)),
                    lineardc.pred=as.double(rep(0,as.integer(ng))),
                    zi=as.double(rep(0,(n.knots+6))),

                    paratps=as.integer(c(0,0,0)),# for future developments
                    as.integer(c(0,0,0)),# for future developments
                    BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*0)),# for future developments
                    BetaTpsMatDc=as.double(matrix(0,nrow=101,ncol=1+4*0)),# for future developments
                    BetaTpsMatY = as.double(matrix(0,nrow=101,ncol=1+4*0)),# for future developments
                    EPS=as.double(c(LIMparam,LIMlogl,LIMderiv)),
                    GH = c(as.integer(GH),as.integer(n.nodes)),
                    paGH = data.matrix(cbind(b_lme,invBi_cholDet,as.data.frame(invBi_chol))),
                    PACKAGE = "frailtypack")



    MartinGale <- matrix(ans$MartinGale,nrow=ng,ncol=5)
Residuals <- matrix(ans$ResLongi,nrow=nsujety,ncol=4)

    if (ans$ier_istop[2] == 4){
      warning("Problem in the loglikelihood computation. The program stopped abnormally. Please verify your dataset. \n")
    }

    if (ans$ier_istop[2] == 2){
      warning("Model did not converge.")
    }
    if (ans$ier_istop[2] == 3){
      warning("Matrix non-positive definite.")
    }


    #AD:
    if (noVarT==1 & noVarY==1 & noVarR==1) nvar<-0
    #AD:

    np <- ans$np
    fit <- NULL
     fit$b <- ans$b

    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$n <- nsujet
    fit$groups <- ng

   fit$n.events <- ans$counts[2]
    fit$n.deaths <- ans$counts[3]
    fit$n.measurements <- nsujety

    if(as.character(typeof)=="0"){
      fit$logLikPenal <- ans$loglik
    }else{
      fit$logLik <- ans$loglik
    }
    #AD:
    fit$LCV <- ans$LCV[1]
    fit$AIC <- ans$LCV[2]

   fit$sigma2 <- ans$b[np-nvar-1-ne_re-netadc - netar-indic.alpha]^2
   fit$alpha <- ans$b[np -nvar-1-ne_re-netadc - netar ]
   #random effects
    fit$B1 <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    Ut <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    Utt <- matrix(rep(0,nRE^2),ncol=nRE,nrow=nRE)
    for(j in 1:nRE){
      for(k in 1:j){
        Ut[j,k]=ans$b[np - nvar - ne_re + k + j*(j-1)/2]
        Utt[k,j]=ans$b[np -nvar - ne_re + k + j*(j-1)/2]
      }}
    fit$B1 <- Ut%*%Utt

    fit$ResidualSE <- ans$b[(np  - nvar - 1)]^2
    fit$etaR <- ans$b[(np  - nvar - 1 - ne_re - netadc - netar + 1):(np  - nvar - 1 -ne_re - netadc)]
   fit$etaT <- ans$b[(np - nvar - 1 - ne_re - netadc + 1):(np - nvar - 1 -ne_re)]


    fit$npar <- np

    #AD:
    if ((noVarT==1 & noVarY==1)) {
      fit$coef <- NULL
    }
    else
    {
      fit$coef <- ans$b[(np - nvar + 1):np]
      noms <- c(factor.names(colnames(X)),factor.names(colnames(X_T)),factor.names(colnames(X_L)))
    #  if (timedep == 1){
    #    while (length(grep("timedep",noms))!=0){
    #      pos <- grep("timedep",noms)[1]
    #     noms <- noms[-pos]
    #      fit$coef <- fit$coef[-(pos:(pos+betaorder+betaknots-1))]
    #    }
    #  }
      names(fit$coef) <- noms


    }

    fit$names.re <- names.matzy


    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)

   varH.etaR <- temp1[(np  - nvar - 1 - ne_re - netadc - netar +1):(np  - nvar - 1 - ne_re - netar),
                        (np  - nvar - 1 - ne_re - netadc - netar +1):(np  - nvar - 1 - ne_re - netar )]
   if(netar>1)fit$se.etaR <- sqrt(diag(varH.etaR))
   if(netar==1)fit$se.etaR <- sqrt(varH.etaR)
   varH.etaT<- temp1[(np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re),
                        (np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re )]
   if(netadc>1)fit$se.etaT <- sqrt(diag(varH.etaT))
   if(netadc==1)fit$se.etaT <- sqrt(varH.etaT)
    fit$se.ResidualSE <- sqrt(temp1[(np  - nvar - ne_re- 1 ),(np  - nvar - ne_re- 1)])
    fit$varHtotal <- temp1
    fit$varHIHtotal <- temp2

    fit$varH <- temp1[(np  - nvar +1):np, (np  - nvar +1 ):np]
    fit$varHIH <- temp2[(np  - nvar +1):np, (np  - nvar +1):np]
    noms <- c("alpha","Eta","MeasurementError","B1",factor.names(colnames(X)),factor.names(colnames(X_T)),factor.names(colnames(X_L)))

    #if (timedep == 1){ # on enleve les variances des parametres des B-splines
    #  while (length(grep("timedep",noms))!=0){
    #    pos <- grep("timedep",noms)[1]
    #    noms <- noms[-pos]
    #    fit$varH <- fit$varH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
    #    fit$varHIH <- fit$varHIH[-(pos:(pos+betaorder+betaknots-1)),-(pos:(pos+betaorder+betaknots-1))]
    #  }
    #}
    fit$nvar<-c(nvarR,nvarT,nvarY)
    #fit$nvarnotdep<-c(nvarR-nvartimedep,nvarT-nvartimedepT,nvarY-nvartimedepY)
    fit$formula <- formula(Terms)
    fit$formula.LongitudinalData <- formula(TermsY)

   fit$xR <- matrix(ans$xR, nrow = size1, ncol = 1)
   fit$lamR <- array(ans$lamR, dim = c(size1,3,1))
   fit$xSuR <- matrix(ans$xSuR, nrow = 100, ncol =1)
   fit$survR <- array(ans$survR, dim = c(mt1,3,1))


   fit$xD <-matrix(ans$xD, nrow = size2, ncol = 1)
   fit$lamD <- matrix(ans$lamD, nrow = size2, ncol = 3)
   fit$xSuD <- matrix(ans$xSuD, nrow = 100, ncol =1)
   fit$survD <- matrix(ans$survD, nrow = mt2, ncol = 3)

    fit$link <- link
    fit$type <- type
    fit$n.strat <- 1
    fit$n.iter <- ans$counts[1]
    fit$typeof <- typeof
    if (typeof == 0){
      fit$n.knots<-n.knots
      fit$kappa <- ans$k0#[2]
      fit$n.knots.temp <- n.knots.temp
      fit$zi <- ans$zi
    }

    #AD:
   fit$noVarRec <- noVarR
    fit$noVarEnd <- noVarT
    fit$noVarY <- noVarY

   fit$nvarRec <- nvarR
    fit$nvarEnd <- nvarT
    fit$nvarY <- nvarY
    fit$istop <- ans$ier_istop[2]

   fit$shape.weib <- ans$paraweib[1:2]#ans$shape.weib
   fit$scale.weib <- ans$paraweib[3:4]#ans$scale.weib

    if(ans$cag0[1]==1)fit$leftCensoring <- TRUE
    if(ans$cag0[1]==0)fit$leftCensoring <- FALSE

    if(fit$leftCensoring){fit$leftCensoring.threshold <-ans$cag0[2]
                          fit$prop.censored <- prop.censored}
    #AD:

    # verif que les martingales ont ete bien calculees
    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"




   #            fit$frailty.var <- MartinGale[,4]#ans$frailty.var
  fit$martingale.res <- MartinGale[,1]#ans$martingale.res
  fit$martingaledeath.res <- MartinGale[,2]#ans$martingaledc.res
  fit$conditional.res <- Residuals[,1]
  fit$marginal.res <- Residuals[,3]
  fit$marginal_chol.res <- Residuals[,4]

  fit$conditional_st.res <- Residuals[,2]
  fit$marginal_st.res <- Residuals[,3]/fit$ResidualSE


  fit$random.effects.pred <- MartinGale[,3:(3+nRE-1)]#ans$frailty.pred

  fit$pred.y.marg <- matrix(ans$Pred_y,ncol=2)[,2]
  fit$pred.y.cond <- matrix(ans$Pred_y,ncol=2)[,1]
  fit$linear.pred <- ans$linear.pred
  fit$lineardeath.pred <- ans$lineardc.pred

    fit$frailty.pred <- MartinGale[,3]#ans$frailty.pred



   fit$AG <- recurrentAG

 # fit$BetaTpsMatTY <- matrix(ans$BetaTpsMatY,nrow=101,ncol=1+4*nvartimedepY)
  #  fit$BetaTpsMatR <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedep)
  #  fit$BetaTpsMatT <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedepT)
  #  fit$nvartimedep <- c(nvartimedep,nvartimedepT,nvartimedepY)

  #  fit$Names.vardepR <- vardep
  #  fit$Names.vardepT <- vardepT
  #  fit$Names.vardepY <- vardepY

    fit$EPS <- ans$EPS

   fit$netar<-netar
   fit$netadc<-netadc
   fit$ne_re <- nRE


   #================================> For the reccurrent
   #========================= Test de Wald

   if ((length(vec.factor) > 0) ){
     Beta <- ans$b[(np - nvar + 1):np]
     VarBeta <- diag(diag(fit$varH))#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])

     nfactor <- length(vec.factor)
     p.wald <- rep(0,nfactor)
     ntot <- nvarR + nvarT + nvarY
     fit$global_chisqR <- waldtest(N=nvarR,nfact=nfactor,place=ind.place,modality=occur,b=Beta,Varb=VarBeta,Llast=nvarT+nvarY,Ntot=ntot)
  
     
     fit$dof_chisqR <- occur
     fit$global_chisq.testR <- 1
     # Calcul de pvalue globale
     for(i in 1:length(vec.factor)){
       p.wald[i] <- signif(1 - pchisq(fit$global_chisqR[i], occur[i]), 3)
     }
     fit$p.global_chisqR <- p.wald
     fit$names.factorR <- vec.factor
   }else{
     fit$global_chisq.testR <- 0
   }

    #================================> For the longitudinal
    #========================= Test de Wald

    if ((length(vec.factorY) > 0) ){

      Beta <- ans$b[(np - nvar + 1):np]

      VarBeta <- diag(diag(fit$varH))#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])


      nfactor <- length(vec.factorY)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarY
      fit$global_chisqY <- waldtest(N=nvarY,nfact=nfactor,place=ind.placeY,modality=occurY,b=Beta,Varb=VarBeta,Lfirts=nvarT+nvarR,Ntot=ntot)
      fit$dof_chisqY <- occurY
      fit$global_chisq.testY <- 1
      # Calcul de pvalue globale
      for(i in 1:length(vec.factorY)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqY[i], occurY[i]), 3)
      }
      fit$p.global_chisqY <- p.wald
      fit$names.factorY <- vec.factorY
    }else{
      fit$global_chisq.testY <- 0

    }

    #================================> For the death
    #========================= Test de Wald
 
    if ((length(vec.factorT) > 0) ){
      Beta <- ans$b[(np - nvar - nvarY + 1):np]

      VarBeta <- diag(diag(fit$varH))#[-c(1:(1+indic.alpha+netar+netadc+1+ne_re))])

      nfactor <- length(vec.factorT)
      p.wald <- rep(0,nfactor)
      ntot <- nvarR + nvarT + nvarY
   # print("deces")
      fit$global_chisqT <- waldtest(N=nvarT,nfact=nfactor,place=ind.placeT,modality=occurT,b=Beta,Varb=VarBeta,Llast=nvarY,Lfirts=nvarR,Ntot=ntot)
      fit$dof_chisqT <- occurT
      fit$global_chisq.testT <- 1
      # Calcul de pvalue globale
     
      for(i in 1:length(vec.factorT)){
        p.wald[i] <- signif(1 - pchisq(fit$global_chisqT[i], occurT[i]), 3)
      }
     
      fit$p.global_chisqT <- p.wald
      fit$names.factorT<- vec.factorT
    }else{
      fit$global_chisq.testT <- 0
    }

   fit$max_rep <- max_rep
   fit$joint.clust <- 1
if(intercept)fit$intercept <- TRUE
else fit$intercept <- FALSE
   fit$Frailty <- FALSE
   fit$methodGH <- method.GH
  fit$n.nodes <- n.nodes
    class(fit) <- "trivPenal"

    if (print.times){
      cost<-proc.time()-ptm
      cat("The program took", round(cost[3],2), "seconds \n")
    }
    fit

  }

