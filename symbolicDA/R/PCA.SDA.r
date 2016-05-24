.PCA.classic<-function(t){
  df<-array(0,c(dim(t)[[1]],dim(t)[[2]]))
  for(i in 1:nrow(t)){
    df[i,]<-apply(t[i,,],1,mean)
  }
  row.w = rep(1, nrow(df))/nrow(df)
  col.w = rep(1, ncol(df))
  center = TRUE 
  scale = TRUE 
  scannf = TRUE 
  nf = 2 
  {
       df <- as.data.frame(df)
      nc <- ncol(df)
      if (any(is.na(df))) 
          stop("na entries in table")
      f1 <- function(v) sum(v * row.w)/sum(row.w)
      f2 <- function(v) sqrt(sum(v * v * row.w)/sum(row.w))
      if (is.logical(center)) {
          if (center) {
              center <- apply(df, 2, f1)
              df <- sweep(df, 2, center)
          }
          else center <- rep(0, nc)
      }
      else if (is.numeric(center) && (length(center) == nc)) 
          df <- sweep(df, 2, center)
      else stop("Non convenient selection for center")
      if (scale) {
          norm <- apply(df, 2, f2)
          norm[norm < 1e-08] <- 1
          df <- sweep(df, 2, norm, "/")
      }
      else norm <- rep(1, nc)
      X <- as.dudi(df, col.w, row.w, scannf = FALSE, nf = 2, 
          call = NULL, type = "pca")
      X$cent <- center
      X$norm <- norm
      X
  }
}

PCA.vertices.SDA<-function(t,pc.number=2){
  vertices<-array(0,c(dim(t)[[1]]*2^dim(t)[[2]],dim(t)[[2]]))
  for(i in 1:dim(t)[[1]]){
    for(j in 1:2^dim(t)[[2]]){
      jj<-(j-1)
      for(k in 1:dim(t)[[2]]){
    if(jj%%2==0){
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,1]
    }
    else{
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,2]
    }
    jj<-jj%/%2
      }
    }
  }
#  res<-as.matrix(vertices)%*%eigen(cor(vertices))$vectors[,1:pc.number]
  res<-as.matrix(vertices)%*%eigen(cor(vertices))$vectors[,1:pc.number]
  pc<-t[,1:pc.number,]
  for(i in 1:dim(t)[[1]]){
    pc[i,,1]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,min)
    pc[i,,2]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,max)
  }
  resul<-list(vertices=vertices,pc=pc)
}


PCA.centers.SDA<-function(t,pc.number=2){
  vertices<-array(0,c(dim(t)[[1]]*2^dim(t)[[2]],dim(t)[[2]]))
  for(i in 1:dim(t)[[1]]){
    for(j in 1:2^dim(t)[[2]]){
      jj<-(j-1)
      for(k in 1:dim(t)[[2]]){
    if(jj%%2==0){
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,1]
    }
    else{
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,2]
    }
    jj<-jj%/%2
      }
    }
  }
  centers<-as.matrix(t[,,1]+t[,,2])/2
  #res<-as.matrix(vertices)%*%eigen(cov(centers))$vectors[,1:pc.number]
  res<-as.matrix(vertices)%*%eigen(cor(centers))$vectors[,1:pc.number]
  #res<-as.matrix(vertices)%*%princomp(centers)$loadings[,1:pc.number]
  pc<-t[,1:pc.number,]
  for(i in 1:dim(t)[[1]]){
    pc[i,,1]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,min)
    pc[i,,2]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,max)
  }
  resul<-list(centers=centers,pc=pc)
}

PCA.spca.SDA<-function(t,pc.number=2){
  vertices<-array(0,c(dim(t)[[1]]*2^dim(t)[[2]],dim(t)[[2]]))
  lvertices<-array(0,c(dim(t)[[1]]*2^dim(t)[[2]],dim(t)[[1]]))
  for(i in 1:dim(t)[[1]]){
    for(j in 1:2^dim(t)[[2]]){
      jj<-(j-1)
      for(k in 1:dim(t)[[2]]){
    if(jj%%2==0){
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,1]
    }
    else{
        vertices[j+(i-1)*2^dim(t)[[2]],k]<-t[i,k,2]
    }
    jj<-jj%/%2
      }
      lvertices[j+(i-1)*2^dim(t)[[2]],i]<-1
    }
  }
  centers<-as.matrix(t[,,1]+t[,,2])/2
  #res<-as.matrix(vertices)%*%eigen(cov(centers))$vectors[,1:pc.number]
  mt<-t(as.matrix(vertices))%*%(as.matrix(lvertices)%*%solve(t(as.matrix(lvertices))%*%as.matrix(lvertices))%*%t(as.matrix(lvertices))%*%as.matrix(vertices))
  res<-as.matrix(vertices)%*%eigen(mt)$vectors[,1:pc.number]
#  #res<<-(as.matrix(lvertices)%*%solve(t(as.matrix(lvertices))%*%as.matrix(lvertices))%*%t(as.matrix(lvertices))%*%as.matrix(vertices))%*%eigen(cor(centers))$vectors[,1:pc.number]
  print(dim(res))
  #res<-as.matrix(vertices)%*%princomp(centers)$loadings[,1:pc.number]
  pc<-t[,1:pc.number,]
  for(i in 1:dim(t)[[1]]){
    pc[i,,1]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,min)
    pc[i,,2]<-apply(res[(1+(i-1)*2^dim(t)[[2]]):(i*2^dim(t)[[2]]),],2,max)
  }
  resul<-list(centers=centers,pc=pc)
}


PCA.mrpca.SDA<-function(t,pc.number=2){

  centers<-(t[,,1]+t[,,2])/2
  radii<-(t[,,2]-t[,,1])/2
  cp<-prcomp(centers)$x[,1:pc.number]
  rp<-prcomp(radii)$x[,1:pc.number]
  rp<-procOPA(cp,rp)$Bhat
  pc<-t[,1:pc.number,]
  for(i in 1:dim(t)[1])
	for(j in 1:pc.number){
  pc[i,j,1]<-cp[i,j]-0.5*rp[i,j]
  pc[i,j,2]<-cp[i,j]+0.5*rp[i,j]
	}
  resul<-list(pc=pc)
}

PCA.spaghetti.SDA<-function(t,pc.number=2){
  s<-t[,,1]
  e<-t[,,2]
  n<-dim(t)[1]
  m<-dim(t)[2]
  cor<-array(0,c(m,m))
  pc<-t[,1:pc.number,]
  for(j in 1:(m-1)){
  for(k in (j+1):m){
   #print(paste(j,k))
    covar<-0
    varj<-0
    vark<-0
    for(i in 1:n){
      varj<-varj+(s[i,j]^2+s[i,j]*e[i,j]+e[i,j]^2)/3
      vark<-vark+(s[i,k]^2+s[i,k]*e[i,k]+e[i,k]^2)/3
      covar<-covar+(2*s[i,j]*s[i,k]+s[i,j]*e[i,k]*2+e[i,j]*e[i,k]+e[i,j]*s[i,k])/6
    }
    varj<-varj/n-mean(t[,j,])^2
    vark<-vark/n-mean(t[,k,])^2
    print(paste("varj=",varj))
    print(paste("vark=",vark))
    covar<-covar/n-mean(t[,j,])*mean(t[,k,])
    cor[j,k]<-cor[j,k]<-covar/sqrt(varj*vark)
  }
  }
  pc[,,1]<-as.matrix(t[,,1])%*%eigen(cor)$vectors[,1:pc.number]
  pc[,,2]<-as.matrix(t[,,2])%*%eigen(cor)$vectors[,1:pc.number]
  resul<-list(cor=cor,pc=pc)
  resul
}
