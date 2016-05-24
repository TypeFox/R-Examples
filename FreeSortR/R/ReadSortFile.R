ReadSortFile<-function(filename,terms=FALSE,septerms=",",sep=";",dec="."){
  #Read a csv file of free sorting data
  
  mat<-read.csv2(filename,sep=sep,dec=dec,row.names=1,header=TRUE)

  nprod<-nrow(mat)
  nsubjects<-ncol(mat)
  
  # Labels of the terms used by the subjects
  lev<-unique(unlist(strsplit(as.character(unlist(mat)),septerms)))
  
  MatTermSubject<-vector("list",length=nsubjects)

  MatSort<-matrix(0,nprod,nsubjects)
  colnames(MatSort)<-colnames(mat)
  rownames(MatSort)<-rownames(mat)
  
  
  
  for (h in 1:nsubjects){
    
    #Sorting of subject h
    S<-as.factor(mat[,h])
    levels(S)<-1:nlevels(S)
    MatSort[,h]<-S
    
    #Terms used by subject h
    S<-strsplit(as.character(mat[,h]),septerms)
    
    Mats<-matrix(0,nprod,length(lev))
    colnames(Mats)<-lev
    rownames(Mats)<-rownames(mat)
    
    for (p in 1:nprod){
      for (l in 1:length(S[[p]])){
        Mats[p,which(S[[p]][l]==lev)]<-1
      }
    }
    MatTermSubject[[h]]<-Mats
  }

  MatTerms<-apply(simplify2array(MatTermSubject),c(1,2),FUN=sum)
    
  if (terms==TRUE){
    return(list(MatSort=MatSort,MatTerms=MatTerms,MatTermSubject=MatTermSubject))
  } else {
    return(MatSort)
  }  
}
