tdt.cnt<-function(sample){

  #Mother with genotype 1
  idx.100=length(grep(TRUE, (sample[,2]==1&sample[,3]==0&sample[,4]==0)))
  idx.101=length(grep(TRUE, (sample[,2]==1&sample[,3]==0&sample[,4]==1)))
  idx.102=length(grep(TRUE, (sample[,2]==1&sample[,3]==0&sample[,4]==2)))

  idx.110=length(grep(TRUE, (sample[,2]==1&sample[,3]==1&sample[,4]==0)))
  idx.112=length(grep(TRUE, (sample[,2]==1&sample[,3]==1&sample[,4]==2)))

  idx.120=length(grep(TRUE, (sample[,2]==1&sample[,3]==2&sample[,4]==0)))
  idx.121=length(grep(TRUE, (sample[,2]==1&sample[,3]==2&sample[,4]==1)))
  idx.122=length(grep(TRUE, (sample[,2]==1&sample[,3]==2&sample[,4]==2)))

  #Father with genotype 1
  idx.010=length(grep(TRUE, (sample[,3]==1&sample[,2]==0&sample[,4]==0)))
  idx.011=length(grep(TRUE, (sample[,3]==1&sample[,2]==0&sample[,4]==1)))
  idx.012=length(grep(TRUE, (sample[,3]==1&sample[,2]==0&sample[,4]==2)))

  idx.210=length(grep(TRUE, (sample[,3]==1&sample[,2]==2&sample[,4]==0)))
  idx.211=length(grep(TRUE, (sample[,3]==1&sample[,2]==2&sample[,4]==1)))
  idx.212=length(grep(TRUE, (sample[,3]==1&sample[,2]==2&sample[,4]==2)))


  if (dim(sample)[2]==5){

    idx=sample[,5]

    idx.111M=length(grep(TRUE, (sample[,2]==1&sample[,3]==1&sample[,4]==1&idx==1)))
    idx.111F=length(grep(TRUE, (sample[,2]==1&sample[,3]==1&sample[,4]==1&idx==0)))

    b.m=idx.101+idx.111M+idx.112+idx.122
    c.m=idx.100+idx.110+idx.111F+idx.121

    b.f=idx.111F+idx.112+idx.011+idx.212
    c.f=idx.110+idx.111M+idx.010+idx.211

    b=idx.101+idx.111M+idx.111F+2*idx.112+idx.122+idx.011+idx.212
    c=idx.100+2*idx.110+idx.111M+idx.111F+idx.121+idx.010+idx.211

    cnt=c(b.m,c.m,b.f,c.f,b,c)

  } else if (dim(sample)[2]==4){

    idx.111=length(grep(TRUE, (sample[,2]==1&sample[,3]==1&sample[,4]==1)))
    b=idx.101+idx.111+2*idx.112+idx.122+idx.011+idx.212
    c=idx.100+2*idx.110+idx.111+idx.121+idx.010+idx.211
    cnt=c(NA,NA,NA,NA,b,c)

  }

  names(cnt)=c('b.m','c.m','b.f','c.f','b','c')
  cnt
}

