# Computes the RMSD (root mean square difference) between observed
# and imputed observations
# RMSD is like RMSE
# Arguments:
#   object is of class yai or yai.impute If the object is yai, impute is called
#     with observed=TRUE.
#   vars is a list of variables, if NULL, all those with imputed values are processed.
#   ... passed to the impute funciton when it is called
#   scale if true, scale the rmsd by the std dev of the observed.
# Value:
#   A data frame with the rownames as vars and the column as RMSD

rmsd.yai <- function (object,vars=NULL,scale=FALSE,...)
{
   if (missing(object)) stop ("object required.")
   if (class(object)[1] == "yai") object = impute.yai(object,vars=vars,observed=TRUE,...)
   if (is.null(object)) stop ("no imputations found using this object")
   nuke = unlist(lapply(object,function (x) all(is.na(x))))
   nuke=nuke[nuke]
   if (length(nuke) > 0) object = object[,-match(names(nuke),names(object)),drop=FALSE]
   object = na.omit(object)
   if (is.null(vars)) vars=names(object)
   vi=paste(unique(strsplit(vars,".o",fixed=TRUE)))
   vi=intersect(vi,names(object))
   notFound=setdiff(vars,names(object))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",paste(notFound,collapse=", "))
   if (length(vi) == 0) stop("nothing to compute")
   vo=paste(vi,"o",sep=".")
   notFound=setdiff(vo,names(object))
   if (length(notFound)>0) warning ("variables not found or had missing values: ",paste(notFound,collapse=", "))
   vo=intersect(vo,names(object))
   both=intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)
   if (length(both) == 0) stop("nothing to compute")
   vo=paste(both,"o",sep=".")
   rmsd=data.frame(rep(NA,length(vo)),row.names=both)
   names(rmsd)=if (scale || length(scale)>1) "rmsdS" else "rmsd"
   usedScale = list()
   for (i in 1:length(both))
   {
      if (!is.factor(object[,both[i]]))
      {
         rmsd[i,1]=sqrt(mean(((object[,both[i]]-object[,vo[i]])^2)))
         if (scale || length(scale)>1) 
         {
            div = NULL
            if (length(scale) > 1) div = scale[[both[i]]]
            if (is.null(div) || is.na(div)) div=attr(object,"scale")[both[i],"scale"] # in data
            if (is.null(div) || is.na(div)) div = sd(object[,vo[i]]) #use observed when needed.
            usedScale[[both[i]]] = div
            rmsd[i,1] = if (!is.na(div) && div > 0.01) rmsd[i,1]/div else NA
         }
      }
   }
   if (length(usedScale) > 0) attr(rmsd,"scale") = unlist(usedScale)
   rmsd
}

rmsd <- rmsd.yai


