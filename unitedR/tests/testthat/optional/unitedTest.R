require(unitedR)
F <- formation(9,NA,c(8,8), c(1,6,7), c(10,6), hardness = c(6,6,0,0,0))
L <- getLineup(F)

start <- sample(.Machine$integer.max,1)
for(i in start:(start+100)) {
  set.seed(i)
  print(i)
  print(simRedCard(F,L))
}

set.seed(187871191)
obj <- F
lineup <- L