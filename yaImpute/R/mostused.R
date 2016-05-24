# Takes a "yai" or "foruse.yaImpute" object and returns an array
# of the n reference observations that are most often used as
# data sources.

mostused = function (object,n=20,kth=NULL)
{
   if (is.null(object)) stop ("object required.")
   if (class(object) == "yai") object = foruse(object,kth=kth)
   if (is.null(object)) stop ("no neighbors found using this object")
   if (class(object)[2] != "foruse.yaImpute") stop("class must be yai or foruse.yaImpute")
   tab=table(object[,1])
   class(tab)="vector"
   tab=data.frame(tab,row.names=names(tab))
   n=min(n,nrow(tab))
   sort(tab[,1],decreasing = TRUE)[1:n]
}
