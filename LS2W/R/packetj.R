packetj <-
function(imwd,level,o){

m <- NULL

J <- imwd$nlevels

n <- c(2^J, 2^J)

  tmpH<-lt.to.name(level,"DC")
  tmpV<-lt.to.name(level,"CD")
  tmpD<-lt.to.name(level,"DD")
  tmpS<-lt.to.name(level,"CC")

  v <- matrix(imwd[[tmpV]],nrow=n[1],ncol=n[2],byrow=TRUE)
  h <- matrix(imwd[[tmpH]],nrow=n[1],ncol=n[2],byrow=TRUE)
  d <- matrix(imwd[[tmpD]],nrow=n[1],ncol=n[2],byrow=TRUE)
  s <- matrix(imwd[[tmpS]],nrow=n[1],ncol=n[2],byrow=TRUE)

m<-rbind(cbind(s,h),cbind(v,d))

m<-m[o,o]		

m
}

