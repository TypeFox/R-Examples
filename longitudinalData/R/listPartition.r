cat("####################################################################
######################## Class ListPartition ######################
############################## Creation ############################
####################################################################\n")

cat("### Definition ###\n")

ListPartition_validity <- function(object){
#    cat("**************** ListPartition ****************")
    return(TRUE)
}

# nbCluster : cluster number
# clusterIndex :
setClass(
   Class="ListPartition",
   representation=representation(
      criterionActif="character",
      initializationMethod="character",
      sorted="logical",
      c1="list",
      c2="list",
      c3="list",
      c4="list",
      c5="list",
      c6="list",
      c7="list",
      c8="list",
      c9="list",
      c10="list",
      c11="list",
      c12="list",
      c13="list",
      c14="list",
      c15="list",
      c16="list",
      c17="list",
      c18="list",
      c19="list",
      c20="list",
      c21="list",
      c22="list",
      c23="list",
      c24="list",
      c25="list",
      c26="list"
   ),
   prototype=prototype(
      criterionActif="Calinski.Harabatz",
      initializationMethod=character(),
      sorted=logical(),
      c1=list(),
      c2=list(),
      c3=list(),
      c5=list(),
      c6=list(),
      c7=list(),
      c8=list(),
      c9=list(),
      c10=list(),
      c11=list(),
      c12=list(),
      c13=list(),
      c14=list(),
      c15=list(),
      c16=list(),
      c17=list(),
      c18=list(),
      c19=list(),
      c20=list(),
      c21=list(),
      c22=list(),
      c23=list(),
      c24=list(),
      c25=list(),
      c26=list()
   )
)


cat("####################################################################
######################## Class ListPartition ######################
############################ Constructeur ##########################
####################################################################\n")

### Le constructeur ne construit que des listes vides
listPartition <- function(){#criterionActif=""){
    return(new("ListPartition"))#,criterionActif=""))
}




cat("\n####################################################################
######################## Test  ListPartition ######################
############################# Accesseurs ###########################
####################################################################\n")


# Si on veut rendre [<- utilisable pour partition, il faut modifier ICI
cat("### Setteur ###\n")
ListPartition_set <- function(x,i,j,value){
    switch(EXPR=i,
           "initializationMethod"={x@initializationMethod <- value},
           "criterionActif"={
               if(value%in%CRITERION_NAMES){
                   if(x@criterionActif!=value){
                       x@criterionActif<-value
                       x@sorted<-FALSE
                   }else{}
               }else{
                   stop("[ListPartition:set] 'criterionActif' should be one of ",CRITERION_NAMES)
               }
           },
           "add"={
               if(class(value)!="Partition"){
                   stop("[ListPartition:setteur]: a ListPartition object shall contain only Partition object.")
               }else{}
               eval(parse(text=paste('x@c',value['nbClusters'],' <- c(x@c',value['nbClusters'],',list(value))',sep="")))
               x@sorted <- FALSE
           },
#           "clear"={eval(parse(text=paste('x@',value,' <-  list()',sep="")))},
#                if(value=="all"){
 #                   x <- listPartition()#new("ListPartition",criterionActif=x@criterionActif,initializationMethod=x@initializationMethod,sorted=x@sorted)
  #              }else{
#                    eval(parse(text=paste('x@',value,' <-  list()',sep="")))
   #             }
#            },
#            if(i %in% CLUSTER_NAMES){
 #               eval(parse(text=paste('x@',i,' <- ',value,sep="")))
  #          }else{
           if(i %in% CLUSTER_NAMES){
               if(identical(value,"clear")){
                   eval(parse(text=paste('x@',i,' <-  list()',sep="")))
                                        #                   eval(parse(text=paste('x@',value,' <-  list()',sep="")))
               }else{
                   stop("[ListPartition:setteur]: Direct acces to ",i," is not permited.")
               }
           }else{
               stop("[ListPartition:setteur] ",i," is not a ListPartition slot.")
           }
    )
    validObject(x)
    return(x)
}

setReplaceMethod("[","ListPartition",ListPartition_set)


cat("### Getteur ###\n")
ListPartition_get <- function(x,i,j,drop){
##    if(is.numeric(i) & (i<1|i>26)){
  ##      stop("[ListPartition:getteur]: i should be in [1:26]")
    ##}else{}
    if(is.numeric(i)){
        stop("[ListPartition:getteur]: to get a clusters list, use ['ci']")
    }else{}
    if(i%in%CRITERION_NAMES){
        return(x['criterionValuesAsMatrix',i])
    }else{}
    switch(EXPR=i,
           "initializationMethod"={return(x@initializationMethod)},
           "criterionActif"={return(x@criterionActif)},
           "c1"={return(x@c1)},
           "c2"={return(x@c2)},
           "c3"={return(x@c3)},
           "c4"={return(x@c4)},
           "c5"={return(x@c5)},
           "c6"={return(x@c6)},
           "c7"={return(x@c7)},
           "c8"={return(x@c8)},
           "c9"={return(x@c9)},
           "c10"={return(x@c10)},
           "c11"={return(x@c11)},
           "c12"={return(x@c12)},
           "c13"={return(x@c13)},
           "c14"={return(x@c14)},
           "c15"={return(x@c15)},
           "c16"={return(x@c16)},
           "c17"={return(x@c17)},
           "c18"={return(x@c18)},
           "c19"={return(x@c19)},
           "c20"={return(x@c20)},
           "c21"={return(x@c21)},
           "c22"={return(x@c22)},
           "c23"={return(x@c23)},
           "c24"={return(x@c24)},
           "c25"={return(x@c25)},
           "c26"={return(x@c26)},
           "criterionValues"={
               if(missing(j)){j <- x@criterionActif}else{}
               listI <- NULL
               result <- list()
               for(i in CLUSTER_NAMES){
                   eval(parse(text=paste("listI <- lapply(x@",i,",function(x){x['criterionValues']['",j,"']})",sep="")))
                   if(length(listI)!=0){
                       eval(parse(text=paste("result <- c(result,",i,"=list(listI))",sep="")))
                   }else{}
               }
               return(result)
           },
           "criterionValuesAsMatrix"={
               if(missing(j)){j <- x@criterionActif}else{}
               result <- list()
               for(i in CLUSTER_NAMES){
                   eval(parse(text=paste("listI <- lapply(x@",i,",function(x){x['criterionValues']['",j,"']})",sep="")))
                   if(length(listI)!=0){
                       eval(parse(text=paste("result <- c(result,",i,"=list(listI))",sep="")))
                   }else{}
               }
               lengthList <- max(sapply(result , length))
               return(t(sapply(result , function(x) c(x,rep(NA,lengthList-length(x))))))
           },
           "sorted"={return(x@sorted)},
           stop("[ListPartition:getteur]: ",i," is not a ListPartition slot")
    )
}

setMethod("[","ListPartition",ListPartition_get)


cat("####################################################################
######################## Class ListPartition ######################
############################## Affichage ###########################
####################################################################\n")



cat("### Method : 'show' for yPartition ###\n") # Si on ajouter un titre a traj, on pourra afficher 'associate traj ='
ListPartition_show <- function(object){
    cat("\n ~ criterionActif          = ",object@criterionActif)
    cat("\n ~ initializationMethod    = ",object@initializationMethod)
    cat("\n ~ sorted                  = ",object@sorted)
    cat("\n ~ criterion values (",object@criterionActif,"):",sep="")
    allCrit <- object['criterionValues']
    if(length(allCrit)==0){
        cat("\n    <no Partition>\n")
    }else{
        for(i in 1:length(allCrit)){
            cat("\n    - ",names(allCrit)[i]," : ")
            catShort(unlist(allCrit[[i]]))
        }
        cat("\n")
    }
    return(invisible(object))
}

setMethod(f="show",signature="ListPartition",
    definition=function(object){
        cat("   ~~~ Class: ListPartition ~~~ ")
        ListPartition_show(object)
    }
)



cat("\n####################################################################
######################## Class ListPartition ######################
############################### Autres #############################
####################################################################\n")

ListPartition_ordered <- function(x,...){
    nameObject<-deparse(substitute(x))
    criterName <- x['criterionActif']
    matPermut <- list()

 #   if(length(x['criterionActif'])==0){stop("ListPartition:ordered]: 'criterionActif' is not define")}else{}
    listCriterion <- x['criterionValues']
    for(i in 1:length(listCriterion)){
        ##listCriterion <- lapply(x[i],function(x){x['criterionValues'][criterName]})
        ##        if(length(listCriterion)!=0){
        orderCi <- order(unlist(listCriterion[[i]]),decreasing=TRUE,na.last=TRUE)
        eval(parse(text=paste("x@",names(listCriterion)[i]," <- x@",names(listCriterion)[i],"[orderCi]",sep="")))
        matPermut <- c(matPermut,list(orderCi))
    }
#    }
    x@sorted <- TRUE
    assign(nameObject,x,envir=parent.frame())
    lengthList <- max(sapply(matPermut , length))
    return(t(sapply(matPermut , function(x) c(x,rep(NA,lengthList-length(x))))))
#    return(matPermut)
}
setMethod("ordered",signature="ListPartition",definition=ListPartition_ordered)



### Attention, si on regroupe en se fixant que sur la partition, on perd les differences liées aux imputations
regroup <- function(object){
    nameObject<-deparse(substitute(object))
    if(!object['sorted']){ordered(object)}else{}
    for (i in 1:26){
        eval(parse(text=paste("listCi <- object['c",i,"']",sep="")))
        j <- length(listCi)
        keep <- rep(TRUE,j)

        while(j>1){
            if(identical(listCi[[j-1]]['clusters'],listCi[[j]]['clusters'])){
                keep[j] <- FALSE
                listCi[[j-1]]['multiplicity'] <- listCi[[j-1]]['multiplicity']+listCi[[j]]['multiplicity']
            }else{}
            j <- j-1
        }
        eval(parse(text=paste("object@c",i," <- listCi[keep]",sep="")))
    }
    assign(nameObject,object,envir=parent.frame())
}



ListPartition_plotCriterion <- function(x, criterion=x['criterionActif'],nbCriterion = 1000){
    ##    minMax <- criterionMinOrMax(calinski=1,test=-1,test2=1)
    if(length(criterion)!=1){stop("[ListPartition:plot] To plot several criterion, use 'plotAllCriterion'")}else{}
    allCrit <- x["criterionValues",criterion]
    if(length(allCrit)!=0){
        lengthList <- max(sapply(allCrit , length))
        allCrit <- sapply(allCrit , function(x) c(x,rep(NA,lengthList-length(x)+1)))
        lengthList <- min(lengthList,nbCriterion)
        mainTitle <- paste(criterion,ifelse(x["sorted"]&x["criterionActif"]==criterion,"\nSorted","\nUnsorted"),sep="")
        matplot(1:lengthList,allCrit[1:lengthList,,drop=FALSE],type="b",lty=1,pch=c(1:9,letters[1:16])[CLUSTER_NAMES %in% dimnames(allCrit)[[2]]],
                xlab="Rerolling",ylab="",main=mainTitle
                )
    }else{
        plot(1,type="n",xlab="Rerolling",ylab="",main=criterion)
    }
    return(invisible())
}

setMethod("plotCriterion",signature=c(x="ListPartition"),ListPartition_plotCriterion)


ListPartition_plotAllCriterion <- function(x, criterion=CRITERION_NAMES[1:5], standardized = TRUE){
    ##    minMax <- criterionMinOrMax(calinski=1,test=-1,test2=1)
    lengthCrit <- length(criterion)

    if(!identical(x['sorted'],logical())){

        ## On plot plusieurs critères, une ligne par critères, le nombre de groupe en abscisse
        if(!x['sorted'] | !(x['criterionActif']%in%criterion)){
            warning("[ListCriterion:plot]: the Partition are unsorted")
            titleSort <- "Unsorted"
        }else{
            titleSort <- paste("Sorted using '",x['criterionActif'],"'",sep="")
        }
        matCrit <- matrix(NA,lengthCrit,26,dimnames=list(criterion,CLUSTER_NAMES))
        for(i in 1:lengthCrit){
            allCrit <- sapply(x["criterionValues",criterion[i]] , function(x){result <- x[[1]];names(result)<-NULL;result})
            matCrit[i,CLUSTER_NAMES%in%names(allCrit)] <- allCrit
        }

        if(standardized){
            for(i in 1:lengthCrit){
                if(!all(is.na(matCrit[i,]))){
                    matCrit[i,] <- matCrit[i,]-min(matCrit[i,],na.rm=TRUE)
                    matCrit[i,] <- matCrit[i,]/max(matCrit[i,],na.rm=TRUE)
                }else{}
            }
            mainTitle <-c("Standardized criterions",titleSort)
        }else{
            mainTitle <- c("Non standardized criterions",titleSort)
        }
        xlab <- paste(1:lengthCrit,":",criterion,sep="",collapse=" ; ")
        rangeVal <- range(which(apply(matCrit,2,function(x){any(!is.na(x))})))
        matplot(rangeVal[1]:rangeVal[2],t(matCrit[,(rangeVal[1]:rangeVal[2]),drop=FALSE]),type="b",main=mainTitle,lty=1,xlab=xlab,ylab="")
    }else{
        plot(1,type="n",xlab="Rerolling",ylab="")
    }
    return(invisible())
}

setMethod("plotAllCriterion",signature=c(x="ListPartition"),ListPartition_plotAllCriterion)



#setMethod("plotCriterion",signature=c(x="ListPartition"),.ListPartition.plotCriterion)


cat("--------------------------------------------------------------------
------------------------ Class ListPartition ----------------------
------------------------------ Creation ----------------------------
--------------------------------------------------------------------\n")
