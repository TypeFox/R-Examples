
"longiPenal" <-
function (formula, formula.LongitudinalData, data,  data.Longi, random, id, intercept = TRUE, link="Random-effects",left.censoring=FALSE, n.knots, kappa,
             maxit=350, hazard="Splines",   init.B,
             init.Random, init.Eta, method.GH = "Standard", n.nodes, LIMparam=1e-3, LIMlogl=1e-3, LIMderiv=1e-3, print.times=TRUE)
{

m2 <- match.call()
m2$formula <-  m2$data <- m2$random <- m2$id <- m2$link <- m2$n.knots <- m2$kappa <- m2$maxit <- m2$hazard  <- m2$init.B <- m2$LIMparam <- m2$LIMlogl <- m2$LIMderiv <- m2$print.times <- m2$left.censoring <- m2$init.Random <- m2$init.Eta <- m2$method.GH <- m2$intercept <- m2$n.nodes <- m2$... <- NULL
Names.data.Longi <- m2$data.Longi

#### Frailty distribution specification ####
if (!(all(random %in% c(1,names(data.Longi))))) { stop("Random effects can be only related to variables from the longitudinal data or the intercept (1)") }
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
                        nbintervR <- 0
                        nbintervDC <- 0
                }
                ### Weibull
                if (typeof == 2){
                       
                        size1 <- 100
                        size2 <- 100
                        equidistant <- 2
                        nbintervR <- 0
                        nbintervDC <- 0
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
                        equidistant <- 0
                        nbintervR <- 0
                        nbintervDC <- 0
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


        if (missing(formula))stop("The argument formula must be specified in every model")
  if (missing(formula.LongitudinalData))stop("The argument formula.LongitudinalData must be specified in every model") #AK

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

        m <- match.call(expand.dots = FALSE) # recupere l'instruction de l'utilisateur

        m$formula.LongitudinalData <- m$data.Longi <- m$n.knots <- m$random <- m$link  <- m$id <- m$kappa <- m$maxit <- m$hazard  <- m$init.B <- m$LIMparam <- m$LIMlogl <- m$LIMderiv <- m$left.censoring <- m$print.times <- m$init.Random <- m$init.Eta <- m$method.GH <- m$intercept <- m$n.nodes <- m$... <- NULL


        special <- c("strata", "cluster", "subcluster", "terminal","num.id","timedep")

#========= Longitudinal Data preparation =========================

        TermsY <- if (missing(data.Longi)){
                terms(formula.LongitudinalData, special)
        }else{
                terms(formula.LongitudinalData, special, data = data.Longi)
        }

        ord <- attr(TermsY, "order") # longueur de ord=nbre de var.expli

#       if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
#si pas vide tous si il ya au moins un qui vaut 1 on arrete

        m2$formula <- TermsY


        m2[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait


#       m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument

        clusterY <- attr(TermsY, "specials")$cluster # (indice) nbre de var qui sont en fonction de cluster()

#Al : tri du jeu de donnees par cluster croissant
        if (length(clusterY)){
                tempc <- untangle.specials(TermsY, "cluster", 1:10)
                ord <- attr(TermsY, "order")[tempc$terms]
                if (any(ord > 1))stop("Cluster can not be used in an interaction")
                m2 <- m2[order(m2[,tempc$vars]),] # soit que des nombres, soit des caracteres
                ordre <- as.integer(row.names(m2)) # recupere l'ordre du data set
                clusterY <- strata(m2[, tempc$vars], shortlabel = TRUE)
          uni.clusterY <- unique(clusterY)
        }



        if (NROW(m2) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

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
        strats <- attr(TermsY, "specials")$strata #nbre de var qui sont en fonction de strata()
        cluster <- attr(TermsY, "specials")$cluster #nbre de var qui sont en fonction de cluster()
        num.id <- attr(TermsY, "specials")$num.id #nbre de var qui sont en fonction de patkey()
        vartimedep <- attr(TermsY, "specials")$timedep #nbre de var en fonction de timedep()

        #booleen pour savoir si au moins une var depend du tps
        if (is.null(vartimedep)) timedepY <- 0
        else timedepY <- 1


        if (timedepY==1) stop("The option 'timedep' is not allowed in this model.")

        if(is.null(num.id)){
                joint.clust <- 1
        }else{
                joint.clust <- 0
                }

        subcluster <- attr(TermsY, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

  if (length(subcluster))stop("'subcluster' is not an allowed option")
        if (length(cluster))stop("Only the argument 'id' can represent the clusters")

                Names.cluster <- id # nom du cluster

  if (length(num.id))stop("'num.id' is not an allowed option")

        if (length(strats))stop("Stratified analysis is not an allowed option yet")


    

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
		
		  	if(factor.spot<ncol(X_L))  X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1],X_L[(factor.spot+1):ncol(X_L)])
     else X_L <- cbind(X_L[1:(factor.spot-1)],model.matrix(as.formula("~"%+%0%+%"+"%+%paste(vec.factorY[i], collapse= "+")), data.Longi)[,-1])
     
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
        uni.cluster<-as.factor(unique(clusterY))


        if(is.null(id)) stop("grouping variable is needed")

        if(is.null(random))     stop("variable for random effects is needed")


        if(length(uni.cluster)==1){
                stop("grouping variable must have more than 1 level")
        }


        if (length(subcluster))stop("'Subcluster' is not allowed")



                if (typeof==0 && missing(kappa)) stop("smoothing parameter (kappa) is required")



        if ((typeof==0) & (length(kappa)!=1)) stop("wrong length of argument 'kappa'")



#newTerm vaut Terms - les variables dont les position sont dans drop






#========================================>

        nvarY<-ncol(X_L) #nvar==1 correspond a 2 situations:

# au cas ou on a aucune var explicative dans la partie rec, mais X=0
# cas ou on a 1seul var explicative, ici X est en general different de 0

      #  varnotdepY <- colnames(X_L)[-grep("timedep",colnames(X_L))]
      #  vardepY <- colnames(X_L)[grep("timedep",colnames(X_L))]
      #  vardepY <- apply(matrix(vardepY,ncol=1,nrow=length(vardepY)),1,timedep.names)

    #    if (length(intersect(varnotdepY,vardepY)) != 0) {
     #           stop("A variable is both used as a constant and time-varying effect covariate")
     #   }

     #   nvartimedepY <- length(vardepY)

    #    filtretpsY <- rep(0,nvarY)
    #    filtretpsY[grep("timedep",colnames(X_L))] <- 1

  varY <- as.matrix(sapply(X_L, as.numeric))

        nsujety<-nrow(X_L)



#=======================================>
#======= Construction du vecteur des indicatrice
        if(length(vec.factorY) > 0){
#               ind.place <- ind.place -1
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
  if(1%in%random){
    names.matzy<-c("Intercept",random[-which(random==1)])
  }else{
    names.matzy<-random
  }

  matzy <- data.matrix(X_L[,which(names(X_L)%in%names.matzy)])
  if(!intercept && 1%in%random) matzy <- as.matrix(cbind(rep(1,nsujety),matzy))

  if(link0==1)netadc <- ncol(matzy)
  if(link0==2)netadc <- 1


#== Left-censoring ==
  cag <- c(0,0)

  if(!is.null(left.censoring) && is.numeric(left.censoring)){

    if(left.censoring<min(Y))stop("The threshold for the left censoring cannot be smaller than the minimal value of the longitudinal outcome")
    cag[1] <- 1
    cag[2] <- left.censoring
    n.censored <- length(which(Y<=left.censoring))
    prop.censored <- n.censored/nsujety
  }
#============= pseudo-adaptive Gauss Hermite ==============
#m <- lme(measuret ~ time+interact+treatment, data = data, random = ~ 1| idd)
if(method.GH=="Pseudo-adaptive"){
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
B_lme <-pdMatrix(m_lme$modelStruct$reStruct)[[1]]

Bi_tmp  <- vector("list", length(uni.cluster) )
invBi_chol <- matrix(rep(0,ne_re*length(uni.cluster)),nrow=length(uni.cluster) )

invB_lme  <- solve(B_lme )
invB_chol<- chol(B_lme)
invB_cholDet <- det(invB_chol)
for (i in 1:length(uni.cluster )) {
  Zi_lme <- Z_lme[clusterY  == i, , drop = FALSE]

  Bi_tmp[[i]] <- chol(solve(crossprod(Zi_lme) / m_lme$sigma^2 + invB_lme))
  for(j in 1:nRE){
    for(k in 1:nRE){
      if (k>=j)invBi_chol[i,j+k*(k-1)/2] <-  Bi_tmp[[i]][j,k]
      #  else     inv.chol.VC(j,k)=matv[k+j*(j-1)/2]
    }}
 }
#Vs[[i]] <- solve(crossprod(Z.i) / m2$sigma^2 + inv.VC)
#}

#invBi_chol <- lapply(Bi_lme, function (x) solve(chol(solve(x))))
invBi_cholDet <- sapply(Bi_tmp,  det)
}else{

  b_lme <- matrix(rep(0,length(uni.cluster)*nRE),ncol=nRE,nrow=length(uni.cluster))
  invBi_cholDet <-  matrix(rep(0,length(uni.cluster)),ncol=1,nrow=length(uni.cluster))
  invBi_chol <- matrix(rep(0,length(uni.cluster)*ne_re),ncol=ne_re,nrow=length(uni.cluster))

}

#================== Survival Data ===================

  Terms <- if (missing(data)){
    terms(formula, special)
  }else{
    terms(formula, special, data = data)
  }

  ord <- attr(Terms, "order") # longueur de ord=nbre de var.expli

 # if (length(ord) & any(ord != 1))stop("Interaction terms are not valid for this function")
#si pas vide tous si il ya au moins un qui vaut 1 on arrete

  m$formula <- Terms

  m[[1]] <- as.name("model.frame") # m[[1]]=frailtypenal, il le remplace par model.frame en fait

  #model.frame(formula = Surv(time, event) ~ cluster(id) + as.factor(dukes) +
  #as.factor(charlson) + sex + chemo + terminal(death), data = readmission)

  m <- eval(m, sys.parent()) #ici la classe de m est un data.frame donc il recupere ce qu'on lui donne en argument


  cluster <- id # (indice) nbre de var qui sont en fonction de cluster()



  subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

  # booleen pour voir si l'objet Y est reponse avant tri des donnees Surv ou SurvIC
  classofT <- attr(model.extract(m, "response"),"class")
# attention le package pec rajoute un element dans l'attribut "class" des objets de survie
  if (length(classofT)>1) classofT <- classofT[2]

  typeofT <- attr(model.extract(m, "response"),"type") # type de reponse : interval etc..





# verification de la sutructure nested si besoin
  if (length(subcluster))stop("subcluster can not be used in the model")



  # tri par ordre croissant de subcluster a l'interieur des clusters
  ordre <- as.integer(row.names(m))



  if (NROW(m) == 0)stop("No (non-missing) observations") #nombre ligne different de 0

  T <- model.extract(m, "response") # objet de type Surv =Time


  if (classofT != "Surv") stop("Response must be a survival object")


  llT <- attr(Terms, "term.labels")#liste des variables explicatives

#cluster(id) as.factor(dukes) as.factor(charlson) sex chemo terminal(death)

#=========================================================>

  mt <- attr(m, "terms") #m devient de class "formula" et "terms"

  X_T <- if (!is.empty.model(mt))model.matrix(mt, m, contrasts) #idem que mt sauf que ici les factor sont divise en plusieurs variables

  ind.placeT <- unique(attr(X_T,"assign")[duplicated(attr(X_T,"assign"))]) ### unique : changement au 25/09/2014


  vec.factorT <- NULL
  vec.factorT <- c(vec.factorT,llT[ind.placeT])

#=========================================================>
# On determine le nombre de categorie pour chaque var categorielle
  strats <- attr(Terms, "specials")$strata #nbre de var qui sont en fonction de strata()
  cluster <- attr(Terms, "specials")$cluster #nbre de var qui sont en fonction de cluster()
  num.id <- attr(Terms, "specials")$num.id #nbre de var qui sont en fonction de patkey()
  vartimedep <- attr(Terms, "specials")$timedep #nbre de var en fonction de timedep()

#booleen pour savoir si au moins une var depend du tps
  if (is.null(vartimedep)) timedepT <- 0
  else timedepT <- 1
 
if (timedepT==1) stop("The option 'timedep' is not allowed in this model.")

  if(length(num.id)) stop("num.id function is not allowed")

  subcluster <- attr(Terms, "specials")$subcluster #nbre de var qui sont en fonction de subcluster()

  if (length(subcluster))stop("'subcluster' is not allowed")
  if (length(cluster))stop("Clusters are defined by the argument 'id'")

  if (length(strats))stop("Stratified analysis is not allowed")


  mat.factorT <- matrix(vec.factorT,ncol=1,nrow=length(vec.factorT))

# Fonction servant a prendre les termes entre "as.factor"
vec.factorT <-apply(mat.factorT,MARGIN=1,FUN=function(x){
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

# On determine le nombre de categorie pour chaque var categorielle
if(length(vec.factorT) > 0){
  vect.factT <- attr(X_T,"dimnames")[[2]]
  vect.factT <- vect.factT[grep(paste(vec.factorT,collapse="|"),vect.factT)]
  
  occurT <- rep(0,length(vec.factorT))
  #    for(i in 1:length(vec.factordc)){
  #      occurdc[i] <- length(grep(vec.factordc[i],vect.factdc))
  #    }
  #  }
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

dropx <- NULL


  clusterT <- 1:nrow(m) #nrow(data) # valeurs inutiles pour un modele de Cox
  uni.clusterT <- 1:nrow(m) #nrow(data)


  if(length(uni.cluster)==1)stop("grouping variable must have more than 1 level")


  if (length(subcluster))stop("subcluster can not be used in the model")

 if (length(strats))stop("Stratified analysis not allowed")


#type <- attr(Y, "type")
  type <- typeofT

if (type != "right" && type != "counting" && type != "interval" && type != "intervaltronc") { # Cox supporte desormais la censure par intervalle
  stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = ""))
}

#       if ((type == "interval" || type == "interval2" || type == "intervaltronc") && intcens == FALSE) { # rajout
#               stop("You are trying to do interval censoring without intcens = TRUE")
#       }



#newTerm vaut Terms - les variables dont les position sont dans drop

  X_T <- model.matrix(Terms, m)

  assign <- lapply(attrassign(X_T, Terms)[-1], function(x) x - 1)
  Xlevels <- .getXlevels(Terms, m)
  contr.save <- attr(X_T, 'contrasts')


# assigne donne la position pour chaque variables
#ncol(X) : nombre de variable sans sans les fonction speciaux comme terminal()...+id
  if(length(vec.factorT) > 0){
  #========================================>
    position <- unlist(assign,use.names=F)
  }

#========================================>

  if (ncol(X_T) == 1){
    X_T<-X_T-1
    noVarT <- 1
  }else{
   X_T <- X_T[, -1, drop = FALSE]
   noVarT <- 0
  }
# on enleve ensuite la premiere colonne correspondant a id


  nvarT<-ncol(X_T) #nvar==1 correspond a 2 situations:


# au cas ou on a aucune var explicative dans la partie rec, mais X=0
# cas ou on a 1seul var explicative, ici X est en general different de 0

#  varnotdepT <- colnames(X_T)[-grep("timedep",colnames(X_T))]
#  vardepT <- colnames(X_T)[grep("timedep",colnames(X_T))]
#  vardepT <- apply(matrix(vardepT,ncol=1,nrow=length(vardepT)),1,timedep.names)

 # if (length(intersect(varnotdepT,vardepT)) != 0) {
#    stop("A variable is both used as a constant and time-varying effect covariate")
 # }

#  nvartimedepT <- length(vardepT)

#  filtretpsT <- rep(0,nvarT)
#  filtretpsT[grep("timedep",colnames(X_T))] <- 1

  varT<-matrix(c(X_T),nrow=nrow(X_T),ncol=nvarT) #matrix sans id et sans partie ex terminal(death)

  ng<-nrow(X_T)


#add Alexandre 04/06/2012
#lire les donnees differemment si censure par intervalle

  if (type=="right"){
    tt0dc <- rep(0,ng)
    tt1dc <- T[,1]
    cens <- T[,2]
     } else {
    tt0dc <- T[,1]
    tt1dc <- T[,2]
    cens <- T[,3]
  }                   # attention ne pas mettre de 0 sinon en cas de left trunc probleme dans la logV


  if (min(cens)==0) cens.data<-1
  if (min(cens)==1 && max(cens)==1) cens.data<-0

  AG<-0



  #=======================================>
  #======= Construction du vecteur des indicatrice
  if(length(vec.factorT) > 0){
  #             ind.place <- ind.place -1
  k <- 0
  for(i in 1:length(vec.factorT)){
    ind.placeT[i] <- ind.placeT[i]+k
    k <- k + occurT[i]-1
    }
  }

#
# Begin JOINT MODEL
#

# Preparing data ...
  nvar = nvarY + nvarT

  if ((typeof == 0) | (typeof == 2)) indic.nb.int <- 0

  if (sum(as.double(varT))==0) nvarT <- 0
  if (sum(as.double(varY))==0) nvarY <- 0

 

  np <- switch(as.character(typeof),
             "0"=((as.integer(n.knots) + 2) + as.integer(nvar) + 1 + ne_re + netadc  ),
             "2"=(2 + nvar + 1  + ne_re + netadc ))

# traitement de l'initialisation du Beta rentre par l'utilisateur

  Beta <- rep(0.5,np)
  if (!missing(init.B)) {
   if (length(init.B) != nvar) stop("Wrong number of regression coefficients in init.B")
 #  if (timedep) stop("You can hardly know initial regression coefficient while time-varying effect")
   Beta <- c(rep(0.5,np-nvar),init.B)
  }

  if (!missing(init.Random)) {
   if (length(init.Random)!=ne_re) stop("init.Random must be of length that corresponds to the number of elements to estimate of the B1 matrix")
      Beta[(np-nvar-ne_re+1):(np-nvar)] <- init.Random
  }
  if (!missing(init.Eta)) {
    if (length(init.Eta)!=netadc) stop("init.Eta must be of length that corresponds to the dimension of the link function")
   Beta[(np -nvar-1-ne_re-netadc+1):(np-nvar-1-ne_re)] <- init.Eta
  }

  xSuT <- matrix(0,nrow=100,ncol=1)
  if (typeof==0){
   mt1 <- size1
  }else{
    mt1 <- 100
  }
  size2 <- mt1

  flush.console()
  if (print.times){
   ptm<-proc.time()
    cat("\n")
   cat("Be patient. The program is computing ... \n")
  }

        ans <- .Fortran("joint_longi",
      as.integer(1),
                        as.integer(nsujety),
                        as.integer(ng),
                as.integer(n.knots),
                k0=as.double(c(0,kappa)), # joint avec generalisation de strate
                        as.double(0),
                        as.double(0),
      as.integer(0),
      as.integer(0),
                as.double(tt0dc),
                        as.double(tt1dc),
                as.integer(cens),
      link0 = as.integer(c(link0,0)),
      yy0 = as.double(Y),
      groupey0 = as.integer(clusterY),
      nb0 = as.integer(nRE),
      matzy0 =as.double(matzy),
      cag0 = as.double(cag),
      as.integer(1),
      matrix(as.double(0),nrow=1,ncol=1),
                as.integer(nvarT),
                as.double(varT),
      nva30 = as.integer(nvarY),
      vaxy0 = as.double(varY),
      noVar = as.integer(c(0,noVarT,noVarY)),
      ag0 = as.integer(1),
                as.integer(maxit),
                np=as.integer(np),
      neta0 = as.integer(c(netadc,0)),
                b=as.double(Beta),
                H=as.double(matrix(0,nrow=np,ncol=np)),
                        HIH=as.double(matrix(0,nrow=np,ncol=np)),

                        loglik=as.double(0),
                        LCV=as.double(rep(0,2)),
                        xR=as.double(matrix(0,nrow=1,ncol=1)),
                        lamR=as.double(matrix(0,nrow=1,ncol=3)),
                        xSuR=as.double(array(0,dim=100)),
                        survR=as.double(array(0,dim=1)),
                        xD=as.double(rep(0,100)),
                        lamD=as.double(matrix(0,nrow=size1,ncol=3)),
                        xSuD=as.double(xSuT),
                        survD=as.double(matrix(0,nrow=size2,ncol=3)),
                as.integer(typeof),
                        as.integer(equidistant),
                  as.integer(c(1,size1,1,mt1)),###
                        counts=as.integer(c(0,0,0)),
                        ier_istop=as.integer(c(0,0)),
                        paraweib=as.double(rep(0,4)),
                        MartinGale=as.double(matrix(0,nrow=ng,ncol=5)),###
      ResLongi = as.double(matrix(0,nrow=nsujety,ncol=4)),
      Pred_y  = as.double(matrix(0,nrow=nsujety,ncol=2)),

                        linear.pred=as.double(rep(0,ng)),
                        lineardc.pred=as.double(rep(0,as.integer(ng))),
                        zi=as.double(rep(0,(n.knots+6))),

      paratps=as.integer(c(0,0,0)),#for future developments
                        as.integer(c(0,0,0)),#for future developments
                        BetaTpsMat=as.double(matrix(0,nrow=101,ncol=1+4*0)), #for future developments
                        BetaTpsMatDc=as.double(matrix(0,nrow=101,ncol=1+4*0)),#for future developments
      BetaTpsMatY = as.double(matrix(0,nrow=101,ncol=1+4*0)),#for future developments
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
    if (noVarT==1 & noVarY==1) nvar<-0
#AD:

    np <- ans$np
    fit <- NULL
    fit$b <- ans$b
    fit$na.action <- attr(m, "na.action")
    fit$call <- call
    fit$groups <- ng



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

    fit$ResidualSE <- sqrt(ans$b[(np  - nvar - ne_re )]^2)
    fit$eta <- ans$b[(np  - nvar - 1 - ne_re - netadc + 1):(np  - nvar - 1 -ne_re)]

    fit$npar <- np

#AD:
    if ((noVarT==1 & noVarY==1)) {
      fit$coef <- NULL
    }
    else
     {
       fit$coef <- ans$b[(np - nvar + 1):np]
        noms <- c(factor.names(colnames(X_T)),factor.names(colnames(X_L)))
      #  if (timedep == 1){
      #          while (length(grep("timedep",noms))!=0){
      #                  pos <- grep("timedep",noms)[1]
      #                  noms <- noms[-pos]
      #                  fit$coef <- fit$coef[-(pos:(pos-1))]
      #          }
       # }
        names(fit$coef) <- noms


     }

  fit$names.re <- names.matzy


    temp1 <- matrix(ans$H, nrow = np, ncol = np)
    temp2 <- matrix(ans$HIH, nrow = np, ncol = np)

    varH.eta <- temp1[(np  - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re),
                      (np - nvar - 1 - ne_re - netadc +1):(np  - nvar - 1 - ne_re )]

    if(netadc>1)fit$se.eta <- sqrt(diag(varH.eta))
    if(netadc==1)fit$se.eta <- sqrt(varH.eta)

    fit$se.ResidualSE <- sqrt(temp1[(np  - nvar - ne_re ),(np  - nvar - ne_re)])
    fit$varHtotal <- temp1
    fit$varHIHtotal <- temp2

    fit$varH <- temp1[(np  - nvar +1):np, (np  - nvar +1 ):np]
    fit$varHIH <- temp2[(np - nvar +1):np, (np - nvar +1):np]

  noms <- c("MeasurementError","B1",factor.names(colnames(X_T)),factor.names(colnames(X_L)))

   #     if (timedep == 1){ # on enleve les variances des parametres des B-splines
  #              while (length(grep("timedep",noms))!=0){
  #                      pos <- grep("timedep",noms)[1]
  #                      noms <- noms[-pos]
  #                      fit$varH <- fit$varH[-(pos:(pos-1)),-(pos:(pos-1))]
  #                      fit$varHIH <- fit$varHIH[-(pos:(pos-1)),-(pos:(pos-1))]
  #              }
  #      }
    fit$nvar<-c(nvarT,nvarY)
     fit$formula <- formula(Terms)
    fit$formula.LongitudinalData <- formula(TermsY)

    fit$xD <- matrix(ans$xD, nrow = size1, ncol = 1)

    fit$lamD <- array(ans$lamD, dim = c(size1,3,1))
    fit$xSuD <- matrix(ans$xSuD, nrow = 100, ncol = 1)
    fit$survD <- array(ans$survD, dim = c(size2,3,1))

    fit$link <- link
    fit$type <- type
    fit$n.strat <- 1
    fit$n.iter <- ans$counts[1]
    fit$typeof <- typeof
    if (typeof == 0){
        fit$n.knots<-n.knots
        fit$kappa <- ans$k0[2]
        fit$n.knots.temp <- n.knots.temp
        fit$zi <- ans$zi
    }

#AD:

    fit$noVarEnd <- noVarT
    fit$noVarY <- noVarY


    fit$nvarEnd <- nvarT
    fit$nvarY <- nvarY
    fit$istop <- ans$ier_istop[2]

    fit$shape.weib <- ans$paraweib[2]#ans$shape.weib
    fit$scale.weib <- ans$paraweib[4]#ans$scale.weib

    if(ans$cag0[1]==1)fit$leftCensoring <- TRUE
    if(ans$cag0[1]==0)fit$leftCensoring <- FALSE

    if(fit$leftCensoring){fit$leftCensoring.threshold <-ans$cag0[2]
                          fit$prop.censored <- prop.censored}
#AD:

    # verif que les martingales ont ete bien calculees
    msg <- "Problem in the estimation of the random effects (perhaps high number of events in some clusters)"

#       if (any(MartinGale[,2]==0)){

#               fit$martingaledeath.res <- msg


#               fit$linear.pred <- msg

#       }else{

#fit$martingale.res <- MartinGale[,1]#ans$martingale.res
fit$martingaledeath.res <- MartinGale[,2]#ans$martingaledc.res
fit$conditional.res <- Residuals[,1]
fit$marginal.res <- Residuals[,3]
fit$marginal_chol.res <- Residuals[,4]

    fit$conditional_st.res <- Residuals[,2]
    fit$marginal_st.res <- Residuals[,3]/fit$ResidualSE


    fit$random.effects.pred <- MartinGale[,3:(3+nRE-1)]#ans$frailty.pred
          fit$pred.y.marg <- matrix(ans$Pred_y,ncol=2)[,2]
fit$pred.y.cond <- matrix(ans$Pred_y,ncol=2)[,1]
         fit$lineardeath.pred <- ans$lineardc.pred
#       }



#    if (joint.clust==0){
#        fit$kendall <- matrix(ans$kendall,nrow=4,ncol=2)
#    }

  #  fit$BetaTpsMat <- matrix(ans$BetaTpsMat,nrow=101,ncol=1+4*nvartimedepT)
  #  fit$BetaTpsMatDc <- matrix(ans$BetaTpsMatDc,nrow=101,ncol=1+4*nvartimedepY)
  #  fit$nvartimedep <- c(nvartimedepT,nvartimedepY)

  #  fit$Names.vardepT <- vardepT
  #  fit$Names.vardepY <- vardepY

    fit$EPS <- ans$EPS

    fit$ne_re <- nRE
    fit$netadc<-netadc




#================================> For the longitudinal
#========================= Test de Wald


        if ((length(vec.factorY) > 0)){

                Beta <- ans$b[(np-nvar + 1):np]

                VarBeta <- diag(diag(fit$varH))

                nfactor <- length(vec.factorY)
                p.wald <- rep(0,nfactor)
                fit$global_chisq <- waldtest(N=nvarY,nfact=nfactor,place=ind.placeY,modality=occurY,b=Beta,Varb=VarBeta,Lfirts=nvarT,Ntot=nvar)

                fit$dof_chisq <- occurY
                fit$global_chisq.test <- 1
# Calcul de pvalue globale
                for(i in 1:length(vec.factorY)){
                        p.wald[i] <- signif(1 - pchisq(fit$global_chisq[i], occurY[i]), 3)
                }
                fit$p.global_chisq <- p.wald
                fit$names.factor <- vec.factorY
        }else{
                fit$global_chisq.test <- 0

        }

#================================> For the death
#========================= Test de Wald


                if ((length(vec.factorT) > 0) ){
                        Beta <- ans$b[(np - nvar + 1):(np)]

                       # if(npbetatps>0){VarBeta <- diag(diag(fit$varH)[-c((np-npbetatps):np)])
                       # }else{
                          VarBeta <- diag(diag(fit$varH))
                        #}
                        nfactor <- length(vec.factorT)
                        p.waldT <- rep(0,nfactor)
                       fit$global_chisq_d <- waldtest(N=nvarT,nfact=nfactor,place=ind.placeT,modality=occurT,b=Beta,Varb=VarBeta,Llast=nvarY,Ntot=nvar)
                        fit$dof_chisq_d <- occurT
                        fit$global_chisq.test_d <- 1
        # Calcul de pvalue globale
                        for(i in 1:length(vec.factorT)){
                                p.waldT[i] <- signif(1 - pchisq(fit$global_chisq_d[i], occurT[i]), 3)
                        }
                        fit$p.global_chisq_d <- p.waldT
                        fit$names.factordc <- vec.factorT
                }else{
                        fit$global_chisq.test_d <- 0
                }
if(intercept)fit$intercept <- TRUE
else fit$intercept <- FALSE
        fit$Frailty <- FALSE
  fit$max_rep <- max_rep
  fit$joint.clust <- 1
  fit$methodGH <- method.GH
  fit$n.nodes <- n.nodes
  class(fit) <- "longiPenal"

        if (print.times){
                cost<-proc.time()-ptm
                cat("The program took", round(cost[3],2), "seconds \n")
        }
 fit

}

