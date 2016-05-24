strata.cnt<-function(sample){

  cnt.222=length(grep(TRUE,(sample[,2]==2&sample[,3]==2&sample[,4]==2)))
  cnt.212=length(grep(TRUE,(sample[,2]==2&sample[,3]==1&sample[,4]==2)))
  cnt.211=length(grep(TRUE,(sample[,2]==2&sample[,3]==1&sample[,4]==1)))
  cnt.122=length(grep(TRUE,(sample[,2]==1&sample[,3]==2&sample[,4]==2)))
  cnt.121=length(grep(TRUE,(sample[,2]==1&sample[,3]==2&sample[,4]==1)))
  cnt.201=length(grep(TRUE,(sample[,2]==2&sample[,3]==0&sample[,4]==1)))
  cnt.021=length(grep(TRUE,(sample[,2]==0&sample[,3]==2&sample[,4]==1)))
  cnt.112=length(grep(TRUE,(sample[,2]==1&sample[,3]==1&sample[,4]==2)))
  cnt.110=length(grep(TRUE,(sample[,2]==1&sample[,3]==1&sample[,4]==0)))
  cnt.101=length(grep(TRUE,(sample[,2]==1&sample[,3]==0&sample[,4]==1)))
  cnt.100=length(grep(TRUE,(sample[,2]==1&sample[,3]==0&sample[,4]==0)))
  cnt.011=length(grep(TRUE,(sample[,2]==0&sample[,3]==1&sample[,4]==1)))
  cnt.010=length(grep(TRUE,(sample[,2]==0&sample[,3]==1&sample[,4]==0)))
  cnt.000=length(grep(TRUE,(sample[,2]==0&sample[,3]==0&sample[,4]==0)))


  if (dim(sample)[2]==5){
    cat('Parent-of-origin information is available.',fill=T)

  idx=sample[,5]
  cnt.111m=length(grep(TRUE,(sample[,2]==1&sample[,3]==1&sample[,4]==1&idx==1)))
  cnt.111f=length(grep(TRUE,(sample[,2]==1&sample[,3]==1&sample[,4]==1&idx==0)))

  cnt=c(cnt.222,
        cnt.212,cnt.211,cnt.122,cnt.121,
        cnt.201,cnt.021,
        cnt.112,cnt.111m,cnt.111f,cnt.110,
        cnt.101,cnt.100,cnt.011,cnt.010,
        cnt.000)
  names(cnt)=c('cnt.222',
               'cnt.212','cnt.211','cnt.122','cnt.121',
               'cnt.201','cnt.021',
               'cnt.112','cnt.111M','cnt.111F','cnt.110',
               'cnt.101','cnt.100','cnt.011','cnt.010',
               'cnt.000')

  } else if (dim(sample)[2]==4){
    cat('No parent-of-origin information is entered.',fill=T)

  cnt.111=length(grep(TRUE,(sample[,2]==1&sample[,3]==1&sample[,4]==1)))

  cnt=c(cnt.222,
        cnt.212,cnt.211,cnt.122,cnt.121,
        cnt.201,cnt.021,
        cnt.112,cnt.111,cnt.110,
        cnt.101,cnt.100,cnt.011,cnt.010,
        cnt.000)
  names(cnt)=c('cnt.222',
               'cnt.212','cnt.211','cnt.122','cnt.121',
               'cnt.201','cnt.021',
               'cnt.112','cnt.111','cnt.110',
               'cnt.101','cnt.100','cnt.011','cnt.010',
               'cnt.000')

  }
  cnt
}
