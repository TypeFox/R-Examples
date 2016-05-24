#Dataset <- read.table("N:/Nouveau dossier/saumon2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

"chisq.sweave" <- function(data,num.var,num.desc,language){

c=names(data)[num.var]
d=names(data)[num.desc]

assign("language",language,envir=.GlobalEnv)

set.X <- function(dataset,var,desc){
#a=vector()
#b=vector()
#a=names(dataset)[c(var)]
#b=names(dataset)[c(desc)]
assign("varY",var,envir=.GlobalEnv)
assign("varX",desc,envir=.GlobalEnv)
assign("dataset",dataset,envir=.GlobalEnv)
c=names(data)[var]
d=names(data)[desc]
 assign("X",c,envir=.GlobalEnv)
assign("Y",d,envir=.GlobalEnv)
}
env=environment(set.X)

Sweave_f<- function(){
get("varY",envir=.GlobalEnv)
get("varX",envir=.GlobalEnv)
get("dataset",envir=.GlobalEnv)
a=getwd()
dir.create(paste(a,"/EnQuireR/",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/fancyvrb.sty",sep=""),paste(a,"/EnQuireR/fancyvrb.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/Sweave.sty",sep=""),paste(a,"/EnQuireR/Sweave.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/upquote.sty",sep=""),paste(a,"/EnQuireR/upquote.sty",sep=""))
file.copy(paste(.libPaths()[1],"/EnQuireR/Sweave/sty/algorithmic.sty",sep=""),paste(a,"/EnQuireR/algorithmic.sty",sep=""))
setwd(paste(a,"/EnQuireR",sep=""))
if(language=="english"){
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/chi_square/En/Bivariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
}
if(language=="french"){
Sweave(paste(.libPaths()[1],"/EnQuireR/Sweave/chi_square/Fr/Bivariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
}
tools::texi2dvi(paste(a,"/EnQuireR/Bivariate_report.tex",sep=""), pdf=TRUE)
setwd(a)
}
set.X(data,num.var,num.desc)
Sweave_f()
}

#reporting.chisq(Dataset,c(34:38),c(39:43))