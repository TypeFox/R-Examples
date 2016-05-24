# Marzo 4 2007 se nombra dudi.tex
#--------------------------------------------------------------------
# afg.tex - función para crear tablas en Latex de análisis factoriales
# C.E. Pardo Julio/04
# requiere paquete xtable
# entra afc de tipo "dudi",
# job = título del trabajo para las tablas
# append
# se imprimen los resultados en consola
# se escribe en el archivo job.tex
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# funcion para una tabla
latex <- function(obj,job="latex",tit="",lab="",append=TRUE,dec=1)
{
    latex <- xtable(obj,digits=rep(dec,(ncol(obj)+1)))
    caption(latex) <- paste(job,"-",tit)
    label(latex) <- paste("t:",job,"-",lab,sep="")
    print(latex, type="latex", file=paste(job,".tex",sep=""),
            append,caption.placement="top")
    # impresion en consola
    cat("\n",tit,"\n")
    print(round(obj,dec))
    if(append) {cat("\n",paste("Table of ",job,":",tit,sep=""),
            "\n was appended in the file:",paste(job,".tex",sep=""),
            "in the work directory \n")}
    if(!append) {cat("\n",paste("Table of ",job,":",tit,sep=""),
        "\n was printed in the new file: ",paste(job,".tex",sep=""),
        "in the work directory \n")}
}
#---------------------------------------------------------------------
dudi.tex <- function(dudi,job="",aidsC=TRUE,aidsR=TRUE,append=TRUE){
afg <- dudi
afgI <- inertia.dudi(afg,row.inertia=TRUE,col.inertia=TRUE)
percent <- diff(afgI$TOT$ratio,1)*100
percent <- c(afgI$TOT$ratio[1]*100,percent)
tvalp <-cbind(afgI$TOT[,1:2]*1000,percent,subset(afgI$TOT,select=3)*100)
names(tvalp) <- c("Eigenvalue","CumInertia","Percent","CumPercent")
latex(tvalp,job,"Eigenvalues * 1000","eigenvalues",FALSE,dec=1)
latex(afg$c1,job,"Eigenvectors","eigenvectors",dec=4)
if(aidsC)
{
    latex(afg$co,job,"Column Coordinates","col-coor",dec=4)
    latex(afgI$col.abs/100,job,"Column Contributions","col-cont")
    latex(afgI$col.rel/100,job,"Representation Quality of the Columns","col-qual")
    latex(afgI$col.cum/100,job,"Cumulated Representation Quality of the Columns","cum-qual-col")
}
if(aidsR)
{
    latex(afg$li,job,"Row Coordinates","row-coor",dec=4)
    latex(afgI$row.abs/100,job,"Row Contributions","row-cont")
    latex(afgI$row.rel/100,job,"Representation Quality of the Rows","row-qual")
    latex(afgI$row.cum/100,job,"Cumulated Representation Quality of the Rows","cum-row-qual")
}
}
#------------------------------------------------------------------------------
