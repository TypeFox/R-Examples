#Barplot simple

"ENbarplot"=function(dataset,X,spl=FALSE,numr=NULL,numc=NULL,report=FALSE){      #X=vecteur des variables choisies, numr=nb de lignes de la fenêtre graphique, numc=nb de colonnes de la fenête graphique, col=couleurs
#le plus simple est de mettre dans X le numéro des colonnes que l'on veut représenter
#mettre des tests si la variable n'est pas qualitative
if(report==FALSE){
barplot_function(dataset,X,spl,numr,numc)
#s rajouter l'argument colour de barplot.function
#env2=environment(ENbarplot)
#assign("barplot.function",barplot.function)
#barplot.function(dataset,X,spl=FALSE,numr=NULL,numc=NULL)
}

if (report==TRUE){
assign("X",X,envir=.GlobalEnv)
assign("dataset",dataset,envir=.GlobalEnv)
a=getwd()
dir.create(paste(a,"/EnQuireR/",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/fancyvrb.sty",sep=""),paste(a,"/EnQuireR/fancyvrb.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/Sweave.sty",sep=""),paste(a,"/EnQuireR/Sweave.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/upquote.sty",sep=""),paste(a,"/EnQuireR/upquote.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/algorithmic.sty",sep=""),paste(a,"/EnQuireR/algorithmic.sty",sep=""))
setwd(paste(a,"/EnQuireR",sep=""))
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/barplot/Univariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
tools::texi2dvi(paste(a,"/EnQuireR/Univariate_report.tex",sep=""), pdf=TRUE)
setwd(a)
}
}
