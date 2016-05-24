# Builds observed over imputed plots in an imputed data frame.
# Arguments:
#   x is of class yai or yai.impute If the object is yai, impute is called
#     with observed=TRUE.
#   vars is a list of variables for the plots, if NULL, all those with imputed and
#     observed values are plotted.
#   pointColor is the color value (or vector) for the xy plots
#   lineColor is the color of the 1:1 plot.
#   spineColor is the color vector of the spineplot 
#   residual, plot residuals when TRUE, imputed over observed when FALSE.
#   ... passed to the impute funciton when it is called and passed to the plot function.

plot.yai = function (x,vars=NULL,pointColor=1,lineColor=2,spineColor=NULL,residual=FALSE,...)
{
   if (missing(x)) stop ("x required.")
   sub=deparse(substitute(x))
   if (class(x)[1] == "yai") x = impute.yai(x,vars=vars,observed=TRUE,...)
   if (is.null(x)) stop ("no imputations found using this x")
   if (is.null(vars)) vars=names(x)
   vi=paste(unique(strsplit(vars,".o",fixed=TRUE)))
   vi=intersect(vi,names(x))
   notFound=setdiff(vars,names(x))
   if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
   if (length(vi) == 0) stop("nothing to plot")
   vo=paste(vi,"o",sep=".")
   notFound=setdiff(vo,names(x))
   if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
   vo=intersect(vo,names(x))
   both=intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)
   if (length(both) == 0) stop("nothing to plot")
   n=length(both)
   rows=floor(sqrt(n))
   if (rows==0) rows=1
   if (rows*rows == n) cols=rows
   else cols=rows+1
   if (rows*cols < n) rows=rows+1
   oldpar=par(mfcol=c(rows,cols))
   on.exit(par(oldpar))

   if (is.null(pointColor)) pointColor=1
   for (imp in both)
   {
      obs=paste(imp,"o",sep=".")
      if ((is.factor(x[,imp])|is.factor(x[,obs])))
      {
         p=try(spineplot(x=x[,imp],y=x[,obs],xlab="Imputed",ylab="Observed",
                   col=spineColor))
         if (class(p)=="try-error") warning ("no plot could be created for ",imp)
      }
      else
      {
         if (residual)
         {   #TODO: figure out how to avoid using suppressWarnings and still pass ... when it contains non-graphical arguments
            suppressWarnings(plot(x=x[,imp],y=x[,obs]-x[,imp],ylab="Residual",xlab="Imputed",col=pointColor,...))
            abline(0,0,col=lineColor)
         }
         else
         {
            suppressWarnings(plot(x=x[,imp],y=x[,obs],xlab="Imputed",ylab="Observed",col=pointColor,...))
            abline(0,1,col=lineColor)
         }
      }
      mtext(imp,font=par("font.main"),line=.7)
   }
   mtext(sub,outer = TRUE,line=-1.6, cex=1.5)
}

plot.impute.yai <- plot.yai

# function to plot the results of running compare
plot.compare.yai=function(x,pointColor=1,lineColor=2,...)
{
   if (!("compare.yai" %in% class(x))) stop("class must include compare.yai")
   addpoints=function(x,y,...) {points(x,y,col=pointColor,...);abline(0,1,col=lineColor)}

   pairs(x,lower.panel=addpoints,upper.panel=addpoints,
         xlim=c(0,max(max(x,na.rm=TRUE),1)),ylim=c(0,max(max(x,na.rm=TRUE),1)),...)
}


