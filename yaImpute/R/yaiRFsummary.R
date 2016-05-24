# Creates a summary various aspects of several randomForests (one for each variable)
# when method randomForest is used.


yaiRFsummary = function(object, nTop=0)
{
   if (class(object) != "yai") stop ("arg must be of class yai")
   if (object$method != "randomForest") stop ("method must be randomForest")
   if (!requireNamespace ("randomForest")) stop("install randomForest and try again")
   scaledImportance = yaiVarImp(object, nTop, plot=FALSE)

   error  = vector(mode="numeric",length=length(names(object$ranForest)))
   errtag = vector(mode="character",length=length(names(object$ranForest)))
   levels = vector(mode="integer",length=length(names(object$ranForest)))
   ntree  = vector(mode="integer",length=length(names(object$ranForest)))
   type   = vector(mode="character",length=length(names(object$ranForest)))

   i = 0 
   for (Rf in object$ranForest)
   { 
     i = i+1
     type[i] = Rf$type    
     if(Rf$type == "regression") 
     {
        error [i] = round(100*Rf$rsq[length(Rf$rsq)], digits=2) 
        errtag[i] = "%var explained"
        levels[i] = NA
     }
     else if(Rf$type == "classification") 
     {
       error [i] = Rf$err.rate[Rf$ntree,"OOB"]
       errtag[i] = "OOB error rate"                    
       levels[i] = nrow(Rf$confusion)
     }
     else 
     {
       error [i] = NA
       errtag[i] = "N/A"                    
       levels[i] = NA
     }
     ntree [i] = Rf$ntree
   }
   forestAttributes=data.frame(ntree,error,errtag,levels,type)
   rownames(forestAttributes)=names(object$ranForest)
   list(forestAttributes=forestAttributes,scaledImportance=scaledImportance)
}
