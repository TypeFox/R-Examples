#Dataset <- read.table("N:/Nouveau dossier/saumon2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

"barplot_sweave" <- function(dataset,var){
set.X <- function(data,num.var){
a=vector()
a=c(names(data)[c(num.var)])
assign("X",a,envir=env)
return(X)
}
env=environment(set.X)
Sweave_f <- function(){
get("X",envir=env)
a=getwd()
dir.create(paste(a,"/EnQuireR/",sep=""))
file.copy(paste(.libPaths(),"/EnQuireR/Sweave/sty/fancyvrb.sty",sep=""),paste(a,"/EnQuireR/fancyvrb.sty",sep=""))
file.copy(paste(.libPaths(),"/EnQuireR/Sweave/sty/Sweave.sty",sep=""),paste(a,"/EnQuireR/Sweave.sty",sep=""))
file.copy(paste(.libPaths(),"/EnQuireR/Sweave/sty/upquote.sty",sep=""),paste(a,"/EnQuireR/upquote.sty",sep=""))
file.copy(paste(.libPaths(),"/EnQuireR/Sweave/sty/algorithmic.sty",sep=""),paste(a,"/EnQuireR/algorithmic.sty",sep=""))
setwd(paste(a,"/barplot",sep=""))
Sweave(paste(.libPaths(),"/EnQuireR/Sweave/barplot/Univariate_report.Rnw",sep=""), driver = RweaveLatex(),syntax = getOption("SweaveSyntax"))
tools::texi2dvi(paste(a,"/EnQuireR/Univariate_report.tex",sep=""), pdf=TRUE)
setwd(a)
}
set.X(dataset,var)
Sweave_f()
}

