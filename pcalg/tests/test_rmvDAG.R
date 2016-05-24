library(pcalg)

set.seed(100)

wmat <- rbind(c(0,1,0,0,0),
              c(0,0,0,1,0),
              c(0,0,0,1,0),
              c(0,0,0,0,1),
              c(0,0,0,0,0))
colnames(wmat) <- rownames(wmat) <- c("1","2","3","4","5")
print.table(wmat, zero.print=".")

g <- as(wmat,"graphNEL")

e.true <- 0
var.true <- 5

dat <- rmvDAG(1000,g)
x5 <- dat[,5]

## test mean
if (t.test(x5,alternative="two.sided")$p.value<0.05) {
  stop("Test of rmvDAG: Mean not correct!")
}

## test variance
if (var.test(x5,rnorm(1000,0,sqrt(5)),ratio=1,
             alternative="two.sided")$p.value<0.05) {
  stop("Test of rmvDAG: Variance not correct!")
}

###----- Check  gmG  generation: ---> ../man/gmG.Rd

## Used to generate "gmG"
set.seed(40)
p <- 8
n <- 5000
## true DAG:
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
gGtrue <- randomDAG(p, prob = 0.3, V = vars)
x <- rmvDAG(n, gGtrue, back.compatible=TRUE)

data(gmG)

## gmG, gmI were produced on 64-bit -> very small difference even in weights:
stopifnot(all.equal(gGtrue, gmG$g,  tol=6e-16),
          all.equal(x,      gmG$ x, tol=1e-15))

###----- Check  gmI  generation: ---> ../man/gmI.Rd

## Used to generate "gmI"
set.seed(123)
p <- 7
n <- 10000
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
datI <- rmvDAG(n, myDAG, back.compatible=TRUE)

data(gmI)
stopifnot(all.equal(myDAG, gmI$ g, tol=6e-16),# for 32-bit
          all.equal(datI,  gmI$ x, tol=1e-15))

