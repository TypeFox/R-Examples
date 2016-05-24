# Computes the correlation between observed and imputed observations
# Arguments:
#   object is of class yai or yai.impute If the object is yai, impute is called
#     with observed=TRUE.
#   vars is a list of variables, if NULL, all those with imputed values are processed.
#   ... passed to the impute funciton when it is called

cor.yai = function (object,vars=NULL,...)
{
   if (missing(object)) stop ("object required.")
   if (class(object)[1] == "yai") object = impute.yai(object,vars=vars,observed=TRUE,...)
   if (is.null(object)) stop ("no imputations found using this object")
   object=na.omit(object)
   if (is.null(vars)) vars=names(object)
   vi=paste(unique(strsplit(vars,".o",fixed=TRUE)))
   vi=intersect(vi,names(object))
   notFound=setdiff(vars,names(object))
   if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
   if (length(vi) == 0) stop("nothing to compute")
   vo=paste(vi,"o",sep=".")
   notFound=setdiff(vo,names(object))
   if (length(notFound)>0) warning ("variables not found: ",paste(notFound,collapse=", "))
   vo=intersect(vo,names(object))
   both=intersect(paste(unique(strsplit(vo,".o",fixed=TRUE))),vi)
   if (length(both) == 0) stop("nothing to compute")
   vo=paste(both,"o",sep=".")
   cors=data.frame(rep(NA,length(vo)),row.names=both)
   names(cors)="r"
   for (i in 1:length(both)) if (!is.factor(object[,both[i]]))
      cors[i,1]=cor(object[,both[i]],object[,vo[i]])
   cors
}
