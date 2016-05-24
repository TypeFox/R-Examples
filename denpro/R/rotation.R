rotation<-function(t,d=2,basis=FALSE)
{

if (d==2){

  rota<-matrix(0,2,2)
  rota[1,1]<-cos(t)
  rota[1,2]<-sin(t)
  rota[2,1]<--sin(t)
  rota[2,2]<-cos(t)

}

if ((d==2) && (basis)){

  rota<-matrix(0,2,2)
  basis1<-c(1,1)
  basis2<-c(-1,1)
  basis1<-basis1/sqrt(sum(basis1^2))
  basis2<-basis2/sqrt(sum(basis2^2))
  rota[,1]<-basis1
  rota[,2]<-basis2

}

if (d==4){

  rotxy<-matrix(0,4,4)
  for (i in 1:4) rotxy[i,i]<-1
  rotxy[1,1]<-cos(t)
  rotxy[1,2]<-sin(t)
  rotxy[2,1]<--sin(t)
  rotxy[2,2]<-cos(t)

  rotyz<-matrix(0,4,4)
  for (i in 1:4) rotyz[i,i]<-1
  rotyz[2,2]<-cos(t)
  rotyz[2,3]<-sin(t) 
  rotyz[3,2]<--sin(t)
  rotyz[3,3]<-cos(t)

  rotzx<-matrix(0,4,4)
  for (i in 1:4) rotzx[i,i]<-1
  rotzx[1,1]<-cos(t)
  rotzx[1,3]<--sin(t)
  rotzx[3,1]<-sin(t)
  rotzx[3,3]<-cos(t)

  rotxw<-matrix(0,4,4)
  for (i in 1:4) rotxw[i,i]<-1
  rotxw[1,1]<-cos(t)
  rotxw[1,4]<-sin(t)
  rotxw[4,1]<--sin(t)
  rotxw[4,4]<-cos(t)

  rotyw<-matrix(0,4,4)
  for (i in 1:4) rotyw[i,i]<-1
  rotyw[2,2]<-cos(t)
  rotyw[2,4]<--sin(t)
  rotyw[4,2]<-sin(t)
  rotyw[4,4]<-cos(t)

  rotzw<-matrix(0,4,4)
  for (i in 1:4) rotzw[i,i]<-1
  rotzw[3,3]<-cos(t)
  rotzw[2,4]<--sin(t)
  rotzw[4,3]<-sin(t)
  rotzw[4,4]<-cos(t)

  rota<-rotxy%*%rotyz%*%rotzx%*%rotxw%*%rotyw%*%rotzw

}

return(rota)
}

