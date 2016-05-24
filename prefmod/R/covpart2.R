covpart2<-function(obj,nobj,ENV)
{
     # splits the data in actual covgroup into NA groups
     blList<-splitData2(obj,nobj,ENV)

     # adds covariate information to actual cov group
     blList<-c(blList,obj[2])
     blList
}
