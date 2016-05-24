"summary.pairwiseTest"<-

function( object, digits=4, p.adjust.method="none", ... )
{

args<-list(...)

 byout <- object$byout
 bynames <- object$bynames
 method <- object$method
 by <- object$by


if(is.null(by))
 {
 # no by-factor:

 p.adjust.args<-list()
 p.adjust.args$method <- p.adjust.method
 p.adjust.args$p <-  byout[[1]][,"p.value", drop=TRUE]


 if(!is.null(args$n))
  {p.adjust.args$n<-args$n}

 p.val.adj <- do.call("p.adjust", p.adjust.args )

 
 out<-data.frame(cbind(p.val.adj, byout[[1]]))

 names(out)<-c("p.val.adj", "p.val.raw", "comparison",
 "groupx", "groupy")
}
else
{
# with by-factor(i.e. byout has length>1)

 byv <- factor()

  byoutdat<-data.frame()

  for (i in 1:length(byout))
   {
    byoutdat <- rbind(byoutdat, byout[[i]])
    byv=c( byv, rep( bynames[i], times=length(byout[[i]][,"p.value"]) ))
   }

 p.adjust.args<-list()
 p.adjust.args$method <- p.adjust.method
 p.adjust.args$p <- byoutdat[,"p.value"]
 if(!is.null(args$n))
  {p.adjust.args$n<-args$n}

 p.val.adj <- do.call("p.adjust", p.adjust.args )

byvarout<-paste(by, collapse=".")

 out<-data.frame(cbind(p.val.adj, byoutdat, byv))

names(out)<-c("p.val.adj", "p.val.raw", "comparison",
 "groupx", "groupy", byvarout)
}

attr(x=out, which="testmethod")<-object$method
attr(x=out, which="p.adjust.method")<-p.adjust.method

out[,"p.val.adj"]<-round(out[,"p.val.adj"], digits=digits)
out[,"p.val.raw"]<-round(out[,"p.val.raw"], digits=digits)

class(out)<-c("summary.pairwiseTest", "data.frame")

return(out)
}



"print.summary.pairwiseTest" <- 

function(x, ...)
{
cat("P-values calculated using", attr(x,"testmethod"),",\n")
cat("Adjustment for multiplicity:", attr(x,"p.adjust.method")  , "\n")
cat("\n")
class(x)<-"data.frame"
print(x, ...)
invisible(x)
}


