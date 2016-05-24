# main function to set up parts for the single cov groups
gen.partsList<-function(nobj,cList,ENV)
{
      g.part<-NULL

      # generate the elements (parts) for each Cov group
      g.part<-c(g.part,lapply(cList, covpart, nobj, ENV))

      g.part
}
