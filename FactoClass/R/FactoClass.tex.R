###########################################################################################
## Crea un archivo tipo tex con los datos y guarda una tabla con el redondeo deseado     ##
##                                                                                       ##
##                                                                                       ##
##                                                                                       ##
## Elaborado por: Campo elias Pardo,                                                     ##
## modificación Pedro Cesar del Campo Neira                                              ##
##                                                                                       ##
## Universidad Nacional de Colombia                                                      ##
##                                                                                       ##
## requiere:xtable      library(xtable)                                                  ##
##                                                                                       ##
## latexDF  ( Obj       := objeto tipo 'data frame'                                      ##
##            job       := nombre del objeto u otro para  las salidas                    ##
##            tit       := Titulo de la tabla                                            ##
##            lab       := (label) etiqueta de la tabla                                  ##
##            append    := Archivo adicionado o no                                       ##
##            dec       := Numero de decimales                                           ##
##            dir       := extension de la carpeta donde desea guardar el archivo tex    ##
##            to.print  := Imprime en consola o no                                       ##
##          )                                                                            ##
##                                                                                       ##
###########################################################################################
# funcion para una tabla
latexDF <- function( obj , job="latex" ,tit="" ,lab="" ,append=TRUE ,dec=1 ,
                     dir = getwd() , to.print = TRUE )
{
        if(!is.data.frame(obj)){return( cat("Only 'data.frame' \n") ) }

    latex <- xtable( obj , digits = rep( dec , ncol(obj) + 1 ) )
    caption(latex) <- paste(job,"-",tit)
    label(latex) <- paste("t:",job,"-",lab,sep="")
    print( latex, type="latex" , file=paste(dir,"/",job,".tex",sep=""),
                  append , caption.placement="top" )

    # impresion en consola en caso de no ser elemento 'FactoClass'

       if(to.print==TRUE){

    cat("\n",tit,"\n")
    print(roundDF(obj,dec))


    if(append) {cat("\n",paste("Table of ",job,":",tit,sep=""),
            "\n was appended in the file:",paste(dir,"/",job,".tex",sep=""),
            "\n")}
    if(!append) {cat("\n",paste("Table of ",job,":",tit,sep=""),
                 "\n was printed in the new file:",paste(dir,"/",job,".tex",sep=""),
                 "\n")}

        }

        invisible( roundDF(obj,dec) )
}
###########################################################################################
##############   FIN DE LA FUNCION            #############################################
###########################################################################################


###########################################################################################
## Función que redondea objetos tipo 'data.frame' sin tener encuenta los factores    ##
##                                                                                   ##
##                                           ##
##                                           ##
## Elaborado por: Pedro Cesar del Campo Neira                        ##
## Universidad Nacional de Colombia                          ##
##                                           ##
##                                                               ##
##                                           ##
## roundDF (  tabla     := objeto 'data frame'                                           ##
##            dec       := numero de decimales                                           ##
##          )                                                                ##
##                                           ##
###########################################################################################

roundDF <- function(tabla,dec=1)
{

           if(!is.data.frame(tabla)){return( cat("Only 'data.frame' \n") ) }

           c.var<- dim(tabla)[2]
           clase <- NULL
           for(i in 1:c.var){ clase <- cbind(clase,class(tabla[,i]))}

           if( all(clase!="numeric") ){return(tabla)}

           numericos <- (1:c.var)[clase=="numeric"]

           SALIDA            <-  tabla
           SALIDA[numericos] <-  round(tabla[numericos],dec)

           SALIDA
}
###########################################################################################
##############FIN DE LA FUNCION            ################################################
###########################################################################################




###########################################################################################
## Crea un archivo tipo tex con los datos y guarda una tabla con el redondeo deseado     ##
## para observalos en consola                                                            ##
##                                                                                   ##
## Elaborado por: Campo Elias Pardo,                             ##
## modificación Pedro Cesar del Campo Neira                          ##
##                                           ##
## Universidad Nacional de Colombia                          ##
##                                           ##
## requiere:ade4      library(ade4)                          ##
##                                           ##
## latexDF  ( FC        := objeto tipo 'FactoClass'                                      ##
##            job       := nombre del objeto u otro para  las salidas            ##
##            append    := Archivo adicionado o no                   ##
##            dec       := Numero de decimales                       ##
##            dir       := extension de la carpeta donde desea guardar el archivo tex    ##
##            p.clust   := (TRUE o FALSE) Opcional para imprimir la clasificación        ##
##          )                                                                ##
##                                           ##
###########################################################################################

FactoClass.tex <- function(FC,job="",append=TRUE, dir = getwd() , p.clust = FALSE ){

 if(class(FC)!="FactoClass"){return(cat("Only object 'FactoClass' \n"))}

###---------------------------------------------------------------------------------



file         <- paste(dir,"/",job,".tex",sep="")

afg          <- FC
afgI         <- inertia.dudi(afg$dudi,row.inertia=TRUE,col.inertia=TRUE)
percent      <- diff(afgI$TOT$ratio,1)*100
percent      <- c(afgI$TOT$ratio[1]*100,percent)
tvalp        <- cbind(afgI$TOT[,1:2]*1000,percent,subset(afgI$TOT,select=3)*100)
names(tvalp) <- c("Eigenvalue","CumInertia","Percent","CumPercent")

carac.cate <- NULL
carac.cont <- NULL
carac.frec <- NULL

if(is.null(afg$carac.cate)==FALSE){carac.cate <- list.to.data(afg$carac.cate)}
if(is.null(afg$carac.cont)==FALSE){carac.cont <- list.to.data(afg$carac.cont)}

cluster    <- NULL

if(p.clust==TRUE){cluster <- data.frame(afg$cluster) }

###---------------------------------------------------------------------------------


e.values    <- latexDF(tvalp,job,"Eigenvalues * 1000","eigenvalues",FALSE,dec=1, dir=dir, to.print=FALSE)
e.vector    <- latexDF(afg$dudi$c1,job,"Eigenvectors","eigenvectors", dir=dir, to.print=FALSE)

co          <- latexDF(afg$dudi$co,job,"Column Coordinates","col-coor",dec=4, dir=dir, to.print=FALSE)
col.abs     <- latexDF(data.frame(afgI$col.abs/100),job,"Column Contributions","col-cont", dir=dir, to.print=FALSE)
col.rel     <- latexDF(data.frame(afgI$col.rel/100),job,"Representation Quality of the Columns","col-qual", dir=dir, to.print=FALSE)
col.cum     <- latexDF(data.frame(afgI$col.cum/100),job,"Cumulated Representation Quality of the Columns","cum-qual-col", dir=dir, to.print=FALSE)

li          <- latexDF(afg$dudi$li,job,"Row Coordinates","row-coor",dec=4, dir=dir, to.print=FALSE)
row.abs     <- latexDF(data.frame(afgI$row.abs/100),job,"Row Contributions","row-cont", dir=dir, to.print=FALSE)
row.rel     <- latexDF(data.frame(afgI$row.rel/100),job,"Representation Quality of the Rows","row-qual", dir=dir, to.print=FALSE)
row.cum     <- latexDF(data.frame(afgI$row.cum/100),job,"Cumulated Representation Quality of the Rows","cum-row-qual", dir=dir, to.print=FALSE)

indices     <- latexDF(afg$indices,job,"Indices for Hierarchical Clustering (WARD)","Indices Ward", dec = 7 , dir=dir, to.print=FALSE)

cor.clus    <- latexDF(afg$cor.clus,job,"cluster coordinates","", dir=dir, to.print=FALSE , dec = 4 )

clus.summ   <- latexDF(afg$clus.summ,job,"cluster description","", dir=dir, to.print=FALSE , dec = 4 )

if(is.null(carac.cate)==FALSE){carac.cate  <- latexDF(carac.cate,job,"Characterization of qualitative variables in the cluster ","Characterization cluster", dir=dir, to.print=FALSE)}
if(is.null(carac.cont)==FALSE){carac.cont  <- latexDF(carac.cont,job,"Characterization of quantitative variables in the cluster","Characterization cluster", dir=dir, to.print=FALSE)}
if(is.null(carac.frec)==FALSE){carac.frec  <- latexDF(carac.frec,job,"Characterization of frequence variables in the cluster","Characterization cluster", dir=dir, to.print=FALSE)}

if(is.null(cluster)==FALSE){cluster <- latexDF(cluster,job,"Table with group of cluster","cluster" , dir=dir, to.print=FALSE)}

###---------------------------------------------------------------------------------

SALIDA <- list( file        =    file        ,
                e.values    =    e.values    ,
        e.vector    =    e.vector    ,
        co          =    co          ,
        col.abs     =    col.abs     ,
        col.rel     =    col.rel     ,
        col.cum     =    col.cum     ,
        li          =    li          ,
        row.abs     =    row.abs     ,
        row.rel     =    row.rel     ,
        row.cum     =    row.cum     ,
        indices     =    indices     ,
                cor.clus    =    cor.clus    ,
                clus.summ   =    clus.summ   ,
        carac.cate  =    carac.cate  ,
        carac.cont  =    carac.cont  ,
        carac.frec  =    carac.frec  ,
        cluster     =    cluster
           )



class(SALIDA)<-"FactoClass.tex"
return(SALIDA)

}
###########################################################################################
############## FIN DE LA FUNCION            ###############################################
###########################################################################################

print.FactoClass.tex <- function(x, ...)
{

 cat("\n")
 cat("FactoClass.tex  (Latex impression) \n\n ")
 cat("The file was printed in: ",x$file ,"\n\n")
 cat("The content of the object and the file is: \n\n")

 sumry <- array("", c(15, 2), list(1:15, c("Table", "Description")))
 sumry[ 1, ] <- c("$e.values  ","Eigenvalues * 1000"                         )
 sumry[ 2, ] <- c("$e.vector  ","Eigenvectors"                           )
 sumry[ 3, ] <- c("$co        ","Column Coordinates"                         )
 sumry[ 4, ] <- c("$col.abs   ","Column Contributions"                       )
 sumry[ 5, ] <- c("$col.rel   ","Representation Quality of the Columns"              )
 sumry[ 6, ] <- c("$col.cum   ","Cumulated Representation Quality of the Columns"        )
 sumry[ 7, ] <- c("$li        ","Row Coordinates"                        )
 sumry[ 8, ] <- c("$row.abs   ","Row Contributions"                      )
 sumry[ 9, ] <- c("$row.rel   ","Representation Quality of the Rows"                 )
 sumry[10, ] <- c("$row.cum   ","Cumulated Representation Quality of the Rows"           )
 sumry[11, ] <- c("$indices   ","Indices for Hierarchical Clustering (WARD)"             )
 sumry[12, ] <- c("$cor.clus  ","Cluster coordinates"            )
 sumry[13, ] <- c("$carac.cate","Characterization of qualitative variables in the cluster "  )
 sumry[14, ] <- c("$carac.cont","Characterization of quantitative variables in the cluster"  )
 sumry[15, ] <- c("$carac.frec","Characterization of frequence variables in the cluster"     )
 class(sumry) <- "table"
 print(sumry)

 if(is.null(x$cluster)==FALSE)
 {

 cat("\n16 $cluster     the cluster to which each row is allocated
")

 }
 cat("\n\n")

}
###########################################################################################
##############FIN DE LA FUNCION            ################################################
###########################################################################################
