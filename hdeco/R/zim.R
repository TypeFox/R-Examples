"zim" <-
function (level=6,lelev=3) {

  MM <- 2^(level-lelev)
  NN <- MM^2
  .Q2 <- sample(1:4,NN,T)
  .Q2 <- matrix(.Q2,nrow=MM,ncol=MM)

  for ( i in 1:lelev){
    .Q2 <- negyit(.Q2)
  }

  attr(.Q2,"cim") <- "C1 tokveletlen 8x8-as homogen darabok"
  egykep()
  imaks(.Q2)
  benc(.Q2)
  teszt(fnev="test")
  return(cat("\nOK"))
}

