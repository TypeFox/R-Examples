"zimuu" <-
function (level=6,lelev=3) {
  M<-2^level
  N<-M^2
  .Q1<-sample(1:4,N,T)
  .Q1<-matrix(.Q1,nrow=M,ncol=M)
  .Q2<-felkev(.Q1,lelev)
  ketegykep()
  attr(.Q1,"cim")<-"C1 tokveletlen"
  attr(.Q2,"cim")<-"C2 a masik kep 8x8-asain belul osszekeverve"
  imaks(.Q1)
  imaks(.Q2)
  bend(.Q1,.Q2)
  teszt(fnev="test")
  return(cat("\nOK"))
}

