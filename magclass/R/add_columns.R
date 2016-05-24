
add_columns<-function(x,addnm=c("new"),dim=3.1){
  dim=old_dim_convention(dim)
  if (dim==1) {
    new_columns<-x[rep(1,length(addnm)),,]
    new_columns<-setCells(new_columns,paste(substr(addnm,1,3),".",(dim(x)[dim]+1):(dim(x)[dim]+length(addnm)),sep=""))
    new_columns[,,]<-NA
  } else if (dim==2) {
    new_columns<-x[,rep(1,length(addnm)),]
    new_columns<-setYears(new_columns,addnm)
    new_columns[,,]<-NA
  } else if (dim>2) {
    new_columns<-x[,,fulldim(x)[[2]][[dim]][[1]]]
    getNames(new_columns,dim=dim-2)<-addnm[1]
    if(length(addnm)>1){
      single_column_x<-new_columns
      for (i in 2:length(addnm)){
        getNames(single_column_x,dim=dim-2)<-addnm[i]
        new_columns<-mbind(new_columns,as.magpie(single_column_x))
      }}
    new_columns[,,]<-NA
  }
  return(mbind(x,new_columns))
}