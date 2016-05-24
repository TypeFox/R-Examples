# Creates a summary of the variable importance scores over several randomForests (one
# for each variable) when method randomForest is used. These values are then scaled.

yaiVarImp = function(object, nTop=20, plot=TRUE, ...)
{
   if (class(object) != "yai") stop ("arg must be of class yai")
   if (object$method != "randomForest") stop ("method must be randomForest")
   if (!requireNamespace ("randomForest")) 
   {
     stop("install randomForest and try again")
     # the purpose of this line of code is to suppress CRAN check notes
     importance <- function (...) NULL
   } else importance <- randomForest::importance
   
   scaledImportance = matrix(NA, nrow = length(names(object$ranForest)), 
      ncol=length(xvars(object)))
   colnames(scaledImportance) = xvars(object)
   rownames(scaledImportance) = names(object$ranForest)

   i = 0
   for (Rf in object$ranForest)
   {
     i = i+1
     one = importance(Rf)
     scale = FALSE
     attr = if (Rf$type == "regression") "%IncMSE" else "MeanDecreaseAccuracy"
     if (nrow(one)>1) scale = sd(one[,attr])>0
     imports = scale(one[,attr],center=TRUE,scale=scale)
     scaledImportance[i,rownames(imports)] = imports
   }

   if (is.na(nTop) | nTop == 0) nTop=ncol(scaledImportance)
   scaledImportance = data.frame(scaledImportance)
   nTop = min(ncol(scaledImportance), nTop)
   best = sort(apply(scaledImportance, 2, median, na.rm=TRUE), 
               decreasing = TRUE, index.return = TRUE)$ix[1:nTop]
   if (plot)
   {
      plt = par()$plt
      oldplt = plt
      plt[1] = .2
      boxplot(as.data.frame(scaledImportance[,best,drop=FALSE]), 
              horizontal=TRUE, par(plt=plt), las=1,
              main=deparse(substitute(object)), xlab="Scaled Importance",...)
      par(plt=oldplt)
      invisible(scaledImportance[,best,FALSE])
   }
   else return(scaledImportance[,best,FALSE])
}

