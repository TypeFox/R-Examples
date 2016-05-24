print.yai = function(x,...)
{
   if (missing(x)) stop ("x required.")
   if (class(x)[1] != "yai") stop("arg class must be yai")
   cat ("\nCall:\n")
   print (x$call)
   if (length(x$obsDropped)== 0) cat ("0 observations dropped\n")
   else cat (length(x$obsDropped),"observations dropped: ",
             x$obsDropped[1:min(15,length(x$obsDropped))],"...\n")
   cat ("method used: ",x$method,"\n")
   if (is.null(x$cancor))
   {
      cat ("Cancor not run\n")
   }
   else
   {
      cat ("Cancor report:\n")
      print(format(data.frame(cor=x$cancor$cor, F=x$ftest$F,
            Pr.F=x$ftest$pgF, Sig=c(" ","NS")[(x$ftest$pgF>x$pVal)+1]),digits=4))
      cat (x$nVec,"vectors used, pVal=",x$pVal,"\n")
      cat ("cancor$xcoef:\n")
      print(x$cancor$xcoef[,1:x$nVec])
   }
   if (!is.null(x$projector))
   {
      cat ("Projector:\n")
      print (x$projector)
   }
   if (is.null(x$ccaVegan))
   {
      cat ("CCA not run\n")
   }
   else
   {
      cat ("CCA analysis:\n")
      if (!requireNamespace ("vegan")) stop("install vegan and try again")
      print (x$ccaVegan)
   }
   if (is.null(x$ranForest))
   {
      cat ("randomForest not run\n")
   }
   else
   {
      cat ("randomForest analysis:\n")
      if (!requireNamespace ("randomForest")) 
        stop("install randomForest and try again")
      print(yaiRFsummary(x))
   }
   cat (sum(x$yDrop),"y variables dropped ",
        paste(names(x$yDrop[x$yDrop]),collapse=","),"\n")
   cat (sum(x$xDrop),"x variables dropped ",
           paste(names(x$xDrop[x$xDrop]),collapse=","),"\n")
   if (x$ann & x$method!="randomForest")  cat ("Note: ann used\n") else cat ("ann not used\n")
   if (length(x$neiDstTrgs)==0) cat ("No target neighbors computed.\n")
   else
   {
      nPr=min(10,nrow(x$neiDstTrgs))
      part=data.frame(x$neiDstTrgs[1:nPr,,drop=FALSE],x$neiIdsTrgs[1:nPr,,drop=FALSE],
                      stringsAsFactors = FALSE)
      names(part)=c(colnames(x$neiDstTrgs),colnames(x$neiIdsTrgs))
      cat ("First",nPr,"targets:\n")
      print (part)
   }
   if (length(x$neiDstRefs)==0) cat ("No reference neighbors computed.\n")
   else
   {
      nPr=min(10,nrow(x$neiDstRefs))
      part=data.frame(x$neiDstRefs[1:nPr,,drop=FALSE],x$neiIdsRefs[1:nPr,,drop=FALSE],
                      stringsAsFactors = FALSE)
      names(part)=c(colnames(x$neiDstRefs),colnames(x$neiIdsRefs))
      cat ("First",nPr,"references:\n")
      print (part)
   }
   if (!is.null(x$biasParameters))
   {
      cat ("Bias correction parameters:\n")
      cat ("trgVal CI =",x$biasParameters$trgValCI,
           " curVal =",x$biasParameters$curVal,
           "\nNumber of passes used =",x$biasParameters$npasses,
           " of ",x$biasParameters$oldk-1,"possible\n")
   }

}

summary.yai = function (object,...) print.yai(object,...)
