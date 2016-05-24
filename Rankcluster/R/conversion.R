#convertir les mu (multi dim, plusieurs groupes) en une matrice
liste3d2mat=function(liste)
{
  g=length(liste)
  d=length(liste[[1]])
  m=rep(0,d)
  
  for(i in 1:d)
    m[i]=length(liste[[1]][[i]])
  
  mu=matrix(0,g,sum(m))
  rownom=rep(0,g)
  colnom=rep("",sum(m))
  
  for(j in 1:g)
  {
    compteur=1
    rownom[j]=paste0("cl",j)
    for(i in 1:d)	
    {
      colnom[compteur]=paste0("dim",i)
      for(k in 1:m[i])
      {
        mu[j,compteur]=liste[[j]][[i]][k]
        compteur=compteur+1	
      }
    }
  }
  colnames(mu)=colnom
  rownames(mu)=rownom
  
  return(mu)
}



liste2d2matD=function(liste)
{
  d=length(liste)
  m=rep(0,d)
  for(i in 1:d)
    m[i]=length(liste[[i]])
  
  
  mu=matrix(0,1,sum(m))
  colnom=rep("",sum(m))
  
  compteur=1
  for(i in 1:d)	
  {
    colnom[compteur]=paste0("dim",i)
    for(k in 1:m[i])
    {
      mu[1,compteur]=liste[[i]][k]
      compteur=compteur+1	
    }
  }
  colnames(mu)=colnom
  
  
  return(mu)
  
}

#convertir les mu (1 dim, plusieurs groupes) en une matrice
liste2d2matG=function(liste)
{
  g=length(liste)
  m=length(liste[[1]])
  
  
  mu=matrix(0,g,m)
  rownom=rep(0,g)
  colnames(mu)=c("dim1",rep("",m-1))
  for(i in 1:g)	
  {
    rownom[i]=paste0("cl",i)
    
    for(k in 1:m)
    {
      mu[i,k]=liste[[i]][k]	
    }
  }
  rownames(mu)=rownom
  
  return(mu)
  
}



#convertir les mu (1 dim, plusieurs groupes) en une matrice
listedistPartiel=function(liste)#,index)
{
  
  n=length(liste)
  d=length(liste[[1]])
  mu=matrix(ncol=d,nrow=n)
  #rownames(mu)=index
  colnames(mu)=paste0("dim",c(1:d))
  for(i in 1:n)
    mu[i,]=liste[[i]]
  
  
  return(mu)
  
}

partition<-function(row)
{
  return(which(row==max(row))-1)
}


#####convertir liste en array
liste3d2array=function(liste)
{
  
  iter=length(liste)
  d=nrow(liste[[1]])
  g=ncol(liste[[1]])
  tab=array(dim=c(iter,g,d))
  for(i in 1:iter)
  {
    tab[i,,]=t(liste[[i]])
  }
  
  return(tab)
  
}


tliste3d2mat=function(liste)
{
  d=length(liste)
  g=length(liste[[1]])
  m=rep(0,d)
  
  for(i in 1:d)
    m[i]=length(liste[[i]][[1]])
  
  mu=matrix(0,g,sum(m))
  rownom=rep(0,g)
  colnom=rep("",sum(m))
  
  for(j in 1:g)
  {
    compteur=1
    rownom[j]=paste0("cl",j)
    for(i in 1:d)	
    {
      colnom[compteur]=paste0("dim",i)
      for(k in 1:m[i])
      {
        mu[j,compteur]=liste[[i]][[j]][k]
        compteur=compteur+1	
      }
    }
  }
  colnames(mu)=colnom
  rownames(mu)=rownom
  
  return(mu)
}

liste2d2matgd=function(liste)
{
  d=length(liste)
  g=length(liste[[1]])
  
  mat=matrix(ncol=d,nrow=g)
  colnames(mat)=paste0("dim",1:d)
  rownames(mat)=paste0("cl",1:g)
  for(i in 1:d)
    mat[,i]=liste[[i]]
  
  return(mat)
}

liste3d2listematgd=function(liste)
{
  d=length(liste[[1]])
  g=length(liste[[1]][[1]])
  coln=paste0("dim",1:d)
  rown=paste0("cl",1:g)
  outliste=lapply(liste,FUN=function(x)
  {
    mat=matrix(ncol=d,nrow=g)
    colnames(mat)=coln
    rownames(mat)=rown
    for(i in 1:d)
      mat[,i]=x[[i]]
    
    return(mat)
  })
  
  return(outliste)
  
}

