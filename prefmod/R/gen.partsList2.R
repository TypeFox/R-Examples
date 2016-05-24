gen.partsList2<-function(nobj,cList,ENV)
{     g.part<-NULL
      # generate the elements (parts) for each Cov group
      g.part<-c(g.part,lapply(cList, covpart2, nobj,ENV))
      g.part
}
