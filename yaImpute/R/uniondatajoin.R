unionDataJoin=function(...,warn=TRUE)
{
#  creates a data frame that has the rows defined by a union of all rownames in the
#  arguments and columns defined by a union of all colnames in the arguments.
#  a single argument can be a list of data frames or matrices
#
#  when warn is TRUE, columns that occur in more than one source are listed as a warning.

   args=list(...)
   if (length(args)==1 && class(args[[1]]) == "list") args=args[[1]]
   for (d in args)
   {
      if (!is.data.frame(d) && !is.matrix(d)) stop ("arguments or list members must be matrices or data frames")
      if (is.matrix(d))
      {
         if (is.null(colnames(d))) stop ("column names are requried within all input matrices")
         if (is.null(rownames(d))) stop ("row names are requried within all input matrices")
         if (length(unique(colnames(d))) != length(colnames(d))) stop("column names must be unique within all input matrices")
      }
   }
   rows=NULL
   cols=NULL
   haveCol=NULL
   for (d in args)
   {
      rows=union(rows,rownames(d))
      haveCol=union(intersect(cols,colnames(d)),haveCol)
      cols=union(cols,colnames(d))
   }
   if (warn & length(haveCol)>0)
      warning ("Columns: \"",paste(haveCol,collapse=", "),
               "\" were defined more than once")
   all=matrix(data=NA,nrow=length(rows),ncol=length(cols))
   all=data.frame(all)
   rownames(all)=rows
   colnames(all)=cols
   factors=rep(FALSE,length(cols))
   names(factors)=cols
   for (d in args)
   {
      theCols=colnames(d)
      if (is.data.frame(d))
      {
         for (var in theCols)
         {
            if (is.factor(d[,var]))
            {
               factors[var] = TRUE
               all[rownames(d),var]=levels(d[,var])[d[,var]]
            }
            else all[rownames(d),var]=d[,var,drop=FALSE]
         }
      }
      else all[rownames(d),theCols]=d
   }
   for (var in colnames(all)) if (factors[var]) all[,var]=as.factor(all[,var])
   all
}
