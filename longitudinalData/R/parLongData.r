cat("\n####################################################################
######################### Class parLongData ########################
############################# Creation #############################
####################################################################\n")

### Pas de trajectoire totalement vide => maxNA<length(time)

ParLongData_validity <- function(object){
#    cat("**** validity ParLongData <empty> ****\n")
    return(TRUE)
}

setClass(
    Class="ParLongData",
    representation=representation(
        type="character",
        col="character",
        pch="character",
        pchPeriod="numeric",
        cex="numeric",
        xlab="character",
        ylab="character"
    ),
    prototype=prototype(
        type=character(),
        col=character(),
        pch=character(),
        pchPeriod=numeric(),
        cex=numeric(),
        xlab=character(),
        ylab=character()
    ),
    validity=ParLongData_validity
)

cat("\n###################################################################
######################## Class parLongData ########################
########################### Constructeur ##########################
###################################################################\n")

parLongData <- function(type,col,pch,pchPeriod,cex,xlab,ylab){
    if(is.numeric(col)){col<-palette()[col]}
    new("ParLongData",type=type,col=col,pch=as.character(pch),pchPeriod=pchPeriod,cex=cex,xlab=xlab,ylab=ylab)
}

parTRAJ <- function(type="l",col="black",pch="1",pchPeriod=0,cex=1,xlab="Time",ylab=""){
    parLongData(type=type,col=col,pch=pch,pchPeriod=pchPeriod,cex=cex,xlab=xlab,ylab=ylab)
}

parMEAN <- function(type="b",col="clusters",pch="letters",pchPeriod=1,cex=1.2,xlab="Time",ylab=""){
    parLongData(type=type,col=col,pch=pch,pchPeriod=pchPeriod,cex=cex,xlab=xlab,ylab=ylab)
}


cat("### Method : 'show' for ParLongData ###\n") # Si on ajoute un titre a traj, on pourra afficher 'associate traj ='
ParLongData_show <- function(object){
    cat("   ~~~ Class: ParLongData ~~~ ")
    cat("\n ~ type       : ",object@type)
    cat("\n ~ col        : [",length(object@col),"] ",sep="");catShort(object@col)
    cat("\n ~ pch        : [",length(object@pch),"] ",sep="");catShort(object@pch)
    cat("\n ~ pchPeriod  : ",object@pchPeriod)
    cat("\n ~ cex        : ",object@cex)
    cat("\n ~ xlab       : ",object@xlab)
    cat("\n ~ ylab       : ",object@ylab,"\n")
    return(invisible(object))
}
setMethod(f="show",signature="ParLongData",definition=ParLongData_show)



cat("### Getteur ###\n")
setMethod("[","ParLongData",
    function(x,i,j,drop){
        switch(EXPR=i,
            "type"={return(x@type)},
            "col"={return(x@col)},
            "pch"={ifelse(grepl("^\\d+$",x@pch,perl=TRUE),as.numeric(x@pch),x@pch)},
            "pchPeriod"={return(x@pchPeriod)},
            "cex"={return(x@cex)},
            "xlab"={return(x@xlab)},
            "ylab"={return(x@ylab)},
            stop("[ParLongData:get]: there is not such a slot in ParLongData")
        )
    }
)


cat("### Setteur ###\n")
setMethod("[<-","ParLongData",
    function(x,i,j,value){
        switch(EXPR=i,
            "type"={x@type<-value},
            "col"={x@col<-value},
            "pch"={x@pch<-as.character(value)},
            "pchPeriod"={x@pchPeriod<-value},
            "cex"={x@cex<-value},
            "xlab"={x@xlab<-value},
            "ylab"={x@ylab<-value},
            stop("[ParLongData:set]: there is not such a slot in ParLongData")
        )
        validObject(x)
        return(x)
    }
)

### Prépare un ParLongData en fonction d'une partition :
###  - Si les champs 'col' et 'pch' sont de taille 1, ils sont remplacés par des vecteurs.
###  - Si col="clusters", crée un vecteur couleur
###  - Si pch="letters" ou "symbol", un vecteur de lettres ou de symboles est créé?
ParLongData_Partition_expand <- function(xParLongData,y){
    col <- xParLongData['col']
    if(identical(col,"clusters")){
        col <- rainbow(y['nbClusters'])[y['clustersAsInteger']]
    }else{}
    col <- rep_len(col,length(y['clusters']))
    xParLongData['col'] <- col

    if(length(xParLongData['pch'])==1){
        xParLongData['pch'] <- switch(as.character(xParLongData['pch']),
                                      "letters"=LETTERS[1:y['nbClusters']][y['clustersAsInteger']],
                                      "symbols"=as.character(1:y['nbClusters'])[y['clustersAsInteger']],
                                      xParLongData['pch'])
    }else{}
    return(xParLongData)
}
setMethod("expandParLongData",signature=c(xParLongData="ParLongData",y="Partition"),def=ParLongData_Partition_expand)


### Prépare un ParLongData en fonction d'un nombre de clusters :
###  - Si col="clusters", crée un vecteur avec une couleur pour chaque cluster
###  - Si pch="letters" ou "symbols", crée un vecteur avec un pch pour chaque cluster
ParLongData_nbClusters_expand <- function(xParLongData,y){
    col <- xParLongData['col']
    if(identical(col,"clusters")){
        col <- rainbow(y)
    }else{
        if(length(col)==1){
            col <- rep(col,y)
        }else{}
    }
    xParLongData['col'] <- col

    if(length(xParLongData['pch'])==1){
        xParLongData['pch'] <- switch(xParLongData['pch'],
                                      "letters"=LETTERS[1:y],
                                      "symbols"=as.character(1:y),
                                      xParLongData['pch'])
    }else{}

    return(xParLongData)
}
setMethod("expandParLongData",signature=c(xParLongData="ParLongData",y="numeric"),def=ParLongData_nbClusters_expand)


cat("\n--------------------------------------------------------------------
----------------------- Fin Test ParLongData -----------------------
--------------------------------------------------------------------\n")
