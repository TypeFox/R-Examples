gen.partsList3<-function(nobj,cList,ENV)
{     g.part<-NULL
      # generate the elements (parts) for each Cov group
      g.part<-c(g.part,lapply(cList, covpart3, nobj,ENV))
      g.part
}
