cat("\n####################################################################
######################### Class parWindows #########################
############################# Creation #############################
####################################################################\n")

### Pas de trajectoire totalement vide => maxNA<length(time)

ParWindows_validity <- function(object){
#    cat("**** validity ParWindows <empty> ****\n")
    return(TRUE)
}

setClass(
    Class="ParWindows",
    representation=representation(
        nbRow="numeric",
        nbCol="numeric",
        addLegend="logical",
        closeScreen="logical",
        screenMatrix="matrix"
    ),
    prototype=prototype(
        nbRow=numeric(),
        nbCol=numeric(),
        addLegend=logical(),
        closeScreen=logical(),
        screenMatrix=matrix(,0,0)
    ),
    validity=ParWindows_validity
)

cat("\n###################################################################
######################### Class parWindows ########################
########################### Constructeur ##########################
###################################################################\n")

parWindows <- function(nbRow,nbCol,addLegend,closeScreen){
    xDivision <- seq.int(0, 1, length.out = nbCol + 1)
    yDivision <- seq.int(ifelse(addLegend,0.9,1), 0,length.out = nbRow + 1)

    screenMatrix <- matrix(c(rep.int(xDivision[-(nbCol + 1)], nbRow),
                     rep.int(xDivision[-1],nbRow),
                     rep.int(yDivision[-1], rep.int(nbCol, nbRow)),
                     rep.int(yDivision[-(nbRow + 1)], rep.int(nbCol, nbRow))
                     ),ncol = 4)
    if(addLegend){screenMatrix <- rbind(screenMatrix,c(0,1,0.5,1))}else{} ### Pourquoi 0.5 et non 0.9 ?
    return(new("ParWindows",nbRow=nbRow,nbCol=nbCol,addLegend=addLegend,closeScreen=closeScreen,screenMatrix=screenMatrix))
}


windowsCut <- function(x,addLegend=TRUE,closeScreen=TRUE){
    if(length(x)==1){
        nbRow <- ceiling(sqrt(x))
        nbCol <- ceiling(x/nbRow)
    }else{
        nbRow <- x[1]
        nbCol <- x[2]
    }
    return(parWindows(nbRow=nbRow,nbCol=nbCol,addLegend=addLegend,closeScreen=closeScreen))
}

cat("### Method : 'show' for ParWindows ###\n") # Si on ajouter un titre a traj, on pourra afficher 'associate traj ='
ParWindows_show <- function(object){
    cat("   ~~~ Class: ParWindows ~~~ ")
    cat("\n ~ nbRow       :",object@nbRow)
    cat("\n ~ nbCol       :",object@nbCol)
    cat("\n ~ addLegend   :",object@addLegend)
    cat("\n ~ closeScreen :",object@closeScreen)
    cat("\n ~ screenMatrix\n")
    print(object@screenMatrix)
    return(invisible(object))
}
setMethod(f="show",signature="ParWindows",definition=ParWindows_show)

cat("### Getteur ###\n")
setMethod("[","ParWindows",
    function(x,i,j,drop){
        switch(EXPR=i,
            "nbRow"={return(x@nbRow)},
            "nbCol"={return(x@nbCol)},
            "addLegend"={return(x@addLegend)},
            "closeScreen"={return(x@closeScreen)},
            "screenMatrix"={return(x@screenMatrix)},
            stop("[ParWindows:get]: there is not such a slot in ParWindows")
        )
    }
)

cat("### Setteur ###\n")
setMethod("[<-","ParWindows",
    function(x,i,j,value){
        switch(EXPR=i,
            "nbRow"={x@nbRow<-value},
            "nbCol"={x@nbCol<-value},
            "addLegend"={x@addLegend<-value},
            "closeScreen"={x@closeScreen<-value},
            "screenMatrix"={x@screenMatrix<-value},
            stop("[ParWindows:set]: there is not such a slot in ParWindows")
        )
        validObject(x)
        return(x)
    }
)


#cat("### Setteur ###\n")
#setMethod("[<-","ParWindows",
#    function(x,i,j,value){
#        switch(EXPR=i,
#            "type"={x@type<-value},
#            "col"={x@col<-value},
#            "pch"={x@pch<-value},
#            "pchPeriod"={x@pchPeriod<-value},
#            "cex"={x@cex<-value},
#            "xlab"={x@xlab<-value},
#            "ylab"={x@ylab<-value},
#            stop("[ParWindows:set]: there is not such a slot in ParWindows")
#        )
#        validObject(x)
#        return(x)
#    }
#)

cat("\n--------------------------------------------------------------------
-------------------------- Fin ParWindows --------------------------
--------------------------------------------------------------------\n")
