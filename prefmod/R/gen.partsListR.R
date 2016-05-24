gen.partsListR<-function(nobj,cList,ENV)
{     g.part<-NULL
      # generate the elements (parts) for each Cov group
      g.part<-c(g.part,lapply(cList, covpartR, nobj,ENV))
      g.part
}
