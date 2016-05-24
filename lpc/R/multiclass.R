multiclass.func <- function(x,y,s0=0){

    ##assumes y is coded 1,2...

    nn <- table(y)
      m <- matrix(0,nrow=nrow(x),ncol=length(nn))
      v <- m
      for(j in 1:length(nn)){
            m[,j] <- rowMeans(x[,y==j])
                v[,j] <- (nn[j]-1)*varr(x[,y==j], meanx=m[,j])
          }
      mbar <- rowMeans(x)
      mm <- m-matrix(mbar,nrow=length(mbar),ncol=length(nn))
      fac <- (sum(nn)/prod(nn))
      scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=TRUE)*mm*mm,1,sum)))

      sd <- sqrt(rowSums(v)*(1/sum(nn-1))*sum(1/nn))
      tt <- scor/(sd+s0)
      mm.stand=t(scale(t(mm),center=FALSE,scale=sd))
      return(list(tt=tt, numer=scor, sd=sd,stand.contrasts=mm.stand))

  }
