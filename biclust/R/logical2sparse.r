#Logical2sparse (a bicluster)
logical2sparse=function(bicluster)
  {
  ret=vector("list",1)
  ret=list(rows=row(as.matrix(bicluster[[1]]))[bicluster[[1]]==TRUE],
            cols=row(as.matrix(bicluster[[2]]))[bicluster[[2]]==TRUE]   )
  ret
  }
  
#Logical2sparse (all biclusters)
logical2sparse=function(bicluster)
  {
  ret=vector("list",length(bicluster))
  for(i in 1:length(bicluster))
    {
    ret[[i]]=list(rows=row(as.matrix(bicluster[[1]]))[bicluster[[1]]==TRUE],
            cols=row(as.matrix(bicluster[[2]]))[bicluster[[2]]==TRUE]   )

    }
  ret
  }
  
#Logical2sparse (list of rows or columns)
logical2sparse=function(input)
  {
  row(as.matrix(input))[input==TRUE]
  }
  
  #Works for ccbiclust()
sebastian2rodrigo=function(biclustering)
  {
  ret=c()
  ret$call=biclustering[[1]]
  ret$number=biclustering[[4]]
  ret$warnings=biclustering[[5]]
  #NOTE: warnings() only prints the output, but returns null, last.warning should
  #be used to get warnings
  ret$rows=vector("list", ret$number)
  ret$cols=vector("list", ret$number)
  n=dim(biclustering[[3]])[1]
  for(i in 1:n)
    {
    ret$rows[[i]]=logical2sparse(biclustering[[2]][i,])
    ret$cols[[i]]=logical2sparse(biclustering[[3]][i,])
    }
  ret
  }
  
rodrigo2sebastian=function(biclustering)
  {
  }