# Some Rcmdr dialogs for the orloca package for graphics

Rcmdr.plot.loca.p <- function(){
   # To ensure that menu name is included in pot file
   gettext("Graphics", domain="R-RcmdrPlugin.orloca")
   gettext("Demands Points Plot", domain="R-RcmdrPlugin.orloca")
   command <- paste("plot(as(", ActiveDataSet(), ", \"loca.p\"), main= \"", gettext("Plot of demand points in", domain="R-RcmdrPlugin.orloca"), " ", ActiveDataSet(), "\")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

Rcmdr.contour.zsum <- function(){
   # To ensure that menu name is included in pot file
   gettext("Contour Plot of sum", domain="R-RcmdrPlugin.orloca")
   command <- paste("contour(as(", ActiveDataSet(), ", \"loca.p\"), main=\"", gettext("Contour Level plot of min-sum objective for", domain="R-RcmdrPlugin.orloca"), " ", ActiveDataSet(), sep="")
   norma <- .RcmdrPlugin.orloca.get.norma(sep="")
   if (norma!="")
      command <- paste(command, "\\n(", gettext("Norm", domain="R-RcmdrPlugin.orloca"), norma, ")\", ", norma, ")", sep="")
   else command <- paste(command, "\")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}


Rcmdr.persp.zsum <- function(){
   # To ensure that menu name is included in pot file
   gettext("3D Plot of sum", domain="R-RcmdrPlugin.orloca")
   command <- paste("persp(as(", ActiveDataSet(), ", \"loca.p\"), main=\"", gettext("3D plot of min-sum objective for", domain="R-RcmdrPlugin.orloca"), " ", ActiveDataSet(), sep="")
   norma <- .RcmdrPlugin.orloca.get.norma(sep="")
   if (norma!="")
      command <- paste(command, "\\n(", gettext("Norm", domain="R-RcmdrPlugin.orloca"), norma, ")\", ", norma, ")", sep="")
   else command <- paste(command, "\")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

Rcmdr.plot.contour.loca.p <- function(){
   # To ensure that menu name is included in pot file
   gettext("Demand & Contour Plot", domain="R-RcmdrPlugin.orloca")
   command <- paste("plot(as(", ActiveDataSet(), ", \"loca.p\"), main= \"", gettext("Plot of demand points and contour plot of", domain="R-RcmdrPlugin.orloca"), " ", ActiveDataSet(), sep="")
   norma <- .RcmdrPlugin.orloca.get.norma(sep="")
   if (norma!="")
      command <- paste(command, "\\n(", gettext("Norm", domain="R-RcmdrPlugin.orloca"), norma, ")\" )", sep="")
   else command <- paste(command, "\")", sep="")
   doItAndPrint(command)
   command <- paste("contour(as(", ActiveDataSet(), ",\"loca.p\"), add=T", sep="")
   command <- paste(command, .RcmdrPlugin.orloca.get.norma(), ")", sep="")
   doItAndPrint(command)
   invisible(NULL)
}

