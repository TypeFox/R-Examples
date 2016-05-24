cai <- function(seq, w, numcode = 1, zero.threshold = 0.0001, zero.to = 0.01){
  stops <- which("Stp" == aaa(translate(s2c(c2s(words())), numcode = numcode)))
  singulets <- which(sapply(syncodons(words(), numcode = numcode), length) == 1)
  exclude <- c(stops, singulets)
  w <- w[-exclude]
  w[w < zero.threshold] <- zero.to # if value is effectively zero make it 0.01
  nncod <- uco(seq)
  nncod <- nncod[-exclude]
  sigma <- nncod %*% log(w)
  exp(sigma/sum(nncod))
}
