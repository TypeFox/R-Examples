
################################################################
#________________    Function  MdsDimChoice     ________________
#                  Choice of the dimension for Mds
# Part is an object of class SortingPartition
# dimen is a vector indicating the min and the max of dimensions
# Returns a table of stress as a function of dimension of Mds
################################################################

MdsDimChoice<-function(Part,dimen=c(2,4),metric=FALSE,ties="primary",itmax=5000,eps=1e-06){
  
  if (!class(Part)=="SortingPartition"){
    return("The argument is not an object of class SortingPartition")
  }  
  else
  {
    dimin<-dimen[1]
    dimax<-dimen[2]
    
    TabStress<-matrix(0,dimax-dimin+1,2)
    
    ListDissimil<-Dissimil(Part)
    MatDissim<-apply(simplify2array(ListDissimil),c(1,2),'sum')
    
    #library(smacof)
    
    for (dimen in dimin:dimax){
      
      res<-MdsDiss(MatDissim,ndim=dimen,metric=metric,itmax=itmax)
      TabStress[dimen-dimin+1,1]<-dimen
      TabStress[dimen-dimin+1,2]<-res$Stress
      
    }
    colnames(TabStress)<-c("Dimension","Stress")
    return(TabStress)
  }
}
