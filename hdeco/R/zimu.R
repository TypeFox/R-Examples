"zimu" <-
function (level=6,lelev=3) {
  M <- 2^level
  N <- M^2
  .Q1 <- sample(1:4,N,T)
  .Q1 <- matrix(.Q1,nrow=M,ncol=M)
  MM <- 2^(level-lelev)
  NN <- MM^2
  .Q2 <- sample(1:4,NN,T)
  .Q2 <- matrix(.Q2,nrow=MM,ncol=MM)

  for( i in 1:lelev) {
    .Q2<-negyit(.Q2)
  }

  attr(.Q1,"cim")<-"C1 tokveletlen"
  attr(.Q2,"cim")<-"C2 tokveletlen 8x8 as homogen darabok"
  ketegykep()
  imaks(.Q1)
  imaks(.Q2)
  bend(.Q1,.Q2)
  teszt(fnev="test")
  return(cat("\nOK"))
}

