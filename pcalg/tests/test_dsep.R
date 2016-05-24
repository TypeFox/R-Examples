#### Compare our dsep() with  dSep() from package "ggm" :
library(pcalg)

set.seed(22)
p <- 8
nreps <- 10
ok <- rep(FALSE,nreps)
for (i in 1:nreps) {
  myDAG <- randomDAG(p, prob = 0.3)
  amat <- as(myDAG,"matrix")
  amat[amat!=0] <- 1

  x <- sample(1:p,1)
  y <- sample(setdiff(1:p,x),1)
  S <- sample(setdiff(1:p,c(x,y)),sample(1:5,1))

  dsepOld <- ggm::dSep(amat,as.character(x),as.character(y),as.character(S))
  dsepRes <- dsep	   (as.character(x),as.character(y),as.character(S),
			    myDAG)
  ok[i] <- (dsepRes == dsepOld)
}

if (!all(ok)) stop("Test dsep wrong: dsep oracle made a mistake!")

#### Test with a graph that is NOT top. sorted (need not be)
amat <- rbind(c(0,1,1,1,0,1,0),
              c(0,0,0,0,0,0,0),
              c(0,1,0,1,0,0,0),
              c(0,0,0,0,0,0,0),
              c(1,0,1,0,0,1,0),
              c(0,0,0,0,0,0,0),
              c(1,0,0,0,1,0,0))
               
colnames(amat) <- rownames(amat) <- as.character(1:7)
g <- as(amat,"graphNEL")

ok <- rep(FALSE, 2)
ok[1] <- ( dsep("7","2", c("1","5"),g) == TRUE ) ## should be TRUE
ok[2] <- ( dsep("7","2", c("3","5"),g) == FALSE ) ## should be FALSE

if (!all(ok)) stop("Test dsep wrong: dsep oracle made a mistake!")
