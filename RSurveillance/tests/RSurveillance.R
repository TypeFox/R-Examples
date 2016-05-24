library(RSurveillance)

## Define a tolerance
tol = 1e-5

## Test case - sep.binom
tested<- seq(10,100, by=10)
prev<- 0.05
sens<- 0.9
expected_result <- c(0.368993670137906, 0.60183101167397, 0.748752848011488, 0.841461456735425,
                     0.899961175672938, 0.936874868617662, 0.960167642524368, 0.974865530299547,
                     0.984139990521287, 0.989992233627259)
observed_result <- sep.binom(tested, prev, sens)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.exact
expected_result <- c(0.80000, 0.96000, 0.99200, 0.99840, 0.99968)
observed_result <- sep.exact(d=1:5, se=0.8)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.hypergeo
N<- c(10, 50, 100, 250, 500)
expected_result <- c(0.40000, 0.40000, 0.40000, 0.78400, 0.92224)
observed_result <- sep.hypergeo(se=0.8, N=N, n=c(5, 25, 50, 125, 250), d = ceiling(0.01*N))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - spp
expected_result <- c(0.9043821, 0.8179069, 0.6050061, 0.3660323)
observed_result <- spp(c(10, 20, 50, 100), 0.99)
stopifnot(all(abs(observed_result - expected_result) < tol))


## Test case - sep
N<- seq(30, 100, by = 5)
se<- 0.95
pstar<- 0.1
n<- rep(30, length(N))
expected_result <- c(0.9998750, 0.9988105, 0.9931679, 0.9933724, 0.9852992, 
                     0.9874888, 0.9790610, 0.9823942, 0.9742575, 0.9781660, 
                     0.9705057, 0.9746682, 0.9675142, 0.9717525, 0.9650811)
observed_result <- sep(N, n, pstar, se = se)
stopifnot(all(abs(observed_result - expected_result) < tol))


## Test case - sep.var.se
sens<- c(rep(0.9, 50), rep(0.95, 100))
expected_result <- 0.8065082
observed_result <- sep.var.se(N=500, sens, 0.01)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.sys
H<- 500
N<- rep(1000, 150)
N[5]<- NA
n<- rep(30, 150)
expected_result <- 0.9008415
observed_result <- sep.sys(NA, N, n, 0.02, 0.05, 0.95)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.binom
expected_result <- c(328, 164,  65, 32, 15)
observed_result <- n.binom(sep=0.95, pstar=c(0.01, 0.02, 0.05, 0.1, 0.2), 0.91)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.hypergeo
expected_result <- c(195, 282, 324, 369)
observed_result <- n.hypergeo(sep=0.95, N=c(200, 500, 1000, 10000), d=ceiling(0.01*c(200, 500, 1000, 10000)), se=0.8)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.freedom
expected_result <- c(173, 251, 288, 324, 328, 333, 332)
observed_result <- n.freedom(N=c(200, 500, 1000, 5000, 10000, 100000, NA), sep=0.95, pstar=0.01, se=0.9)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.2stage
expected_result.1 <- 518
expected_result.2 <- c(13, 17, 17, 18, 18, 18, 17)
observed_result.1 <- n.2stage(1000, c(50, 100, 200, 500, 1000, 5000, NA), 0.95, 0.5, 0.01, 0.05, 0.8)[[1]]
observed_result.2 <- n.2stage(1000, c(50, 100, 200, 500, 1000, 5000, NA), 0.95, 0.5, 0.01, 0.05, 0.8)[[2]]
stopifnot(all(abs(observed_result.1 - expected_result.1) < tol))
stopifnot(all(abs(observed_result.2 - expected_result.2) < tol))

## Test case - pstar.calc
expected_result <- c(0.11680795, 0.09683112, 0.08256039, 
                     0.07185595, 0.06352901, 0.05747207, 0.05686632)
observed_result <- pstar.calc(500, n=50, sep=0.95, se=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - pfree.equ
expected_result <- c(0.9848485, 0.9876543, 0.9898990, 0.9917355, 0.9932660,
                     0.9945610, 0.9956710, 0.9966330, 0.9974747, 0.9982175,
                     0.9988777, 0.9994684)
observed_result <- pfree.equ(seq(0.4, 0.95, by = 0.05), 0.01)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))
## prior.pfree.equ
expected_result <- c(0.9750000, 0.9777778, 0.9800000, 0.9818182, 0.9833333,
                     0.9846154, 0.9857143, 0.9866667, 0.9875000, 0.9882353,
                     0.9888889, 0.9894737)
observed_result <- pfree.equ(seq(0.4, 0.95, by = 0.05), 0.01)[[2]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - disc.prior
expected_result <- c(0.4950, 0.5940, 0.6930, 0.7920, 0.8910, 0.9405)
observed_result <- disc.prior(c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), 0.01)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - pfree.1
expected_result <- c(0.6364104, 0.6667116, 0.7000424, 
                     0.6203008, 0.7101865, 0.6490526)
observed_result <- pfree.1(c(0.44, 0.51, 0.58, 0.4, 0.6, 0.47), 0.01, 0.5)[,4]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - pfree.calc
expected_result <- c(0.6203008, 0.7261937, 0.8099976, 0.8709093,
                     0.9124970, 0.9396922, 0.9569776, 0.9677666,
                     0.9744246, 0.9785045)
observed_result <- pfree.calc(rep(0.4,10), 0.01, 0.5)[,5]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.pfree
expected_result <- c(0.7500000, 0.8888889, 0.9473684, 0.9898990)
observed_result <- sep.pfree(0.5, c(0.8, 0.9, 0.95, 0.99))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.prior
expected_result <- c(0.1, 0.2, 0.5, 1.0)
observed_result <- sep.prior(c(0.9, 0.95, 0.98, 0.99), 0.01)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.pfree
expected_result <- c(247, 319, 403, 461)
observed_result <- n.pfree(pfree = c(0.9, 0.95, 0.98, 0.99), prior = 0.5, 0.01, 0.01, 0.8, 1000)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - se.series
expected_result <- 0.7524
observed_result <- se.series(c(0.99, 0.95, 0.8))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - se.parallel
expected_result <- 0.9999
observed_result <- se.parallel(c(0.99, 0.95, 0.8))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sp.series
expected_result <- 0.9999
observed_result <- sp.series(c(0.99, 0.95, 0.8))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sp.parallel
expected_result <- 0.7524
observed_result <- sp.parallel(c(0.99, 0.95, 0.8))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.pooled
expected_result <- c(0.3946345, 0.6335326, 0.7781533, 0.8657016,
                     0.9187004, 0.9507840, 0.9702064, 0.9819640,
                     0.9890816, 0.9933904)
observed_result <- sep.pooled(1:10*5, 5, 0.02, 0.9, 0.99)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))
## spp.pooled
expected_result <- c(0.9509900, 0.9043821, 0.8600584, 0.8179069,
                     0.7778214, 0.7397004, 0.7034477, 0.6689718,
                     0.6361855, 0.6050061)
observed_result <- sep.pooled(1:10*5, 5, 0.02, 0.9, 0.99)[[2]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.pooled
expected_result <- c(15, 6, 3, 2)
observed_result <- n.pooled(0.95, c(2, 5, 10, 20), 0.1, c(0.99, 0.98, 0.97, 0.95), 1)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.freecalc
expected_result <- 0.870566
observed_result <- sep.freecalc(150, 30, 2, 0.9, 0.98, 0.1)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.hp
expected_result <- 0.878252
observed_result <- sep.hp(150, 30, 2, 0.9, 0.98, 15)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.binom.imperfect
expected_result <- c(0.1012041, 0.3144826, 0.5181777, 0.6775625,
                     0.7911787, 0.8678935, 0.9178840, 0.9496525,
                     0.9694677, 0.9816500)
observed_result <- sep.binom.imperfect(1:10*5, 2, 0.95, 0.98, 0.1)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sph.hp
expected_result <- c(0.5487535, 0.6594135, 0.7747132, 0.8839634,
                     0.9688938, 1.0000000)
observed_result <- sph.hp(500, 30, 2, 95:100/100)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sph.binom
expected_result <- c(0.99615761, 0.98382236, 0.96466169, 0.94010102,
                     0.87945431, 0.73577139, 0.40327171, 0.08937548)
observed_result <- sph.binom(c(5, 10, 15, 20, 30, 50, 100, 200), 2, 0.98)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.hp
expected_result <- c(10, 0.467678, 0.9043821, 65, 1, 0.05)
observed_result <- n.hp(65,0.95,c=1,se=0.95,sp=0.99,pstar=0.05, minSpH=0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.c.hp
expected_result <- c(52, 0.9544449, 0.9044236, 65, 2, 0.05)
observed_result <- n.c.hp(65,0.95,c=5,se=0.95,sp=0.99,pstar=0.05, minSpH=0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.freecalc
expected_result <- c(10, 0.4411522, 0.9043821, 65, 1, 0.05)
observed_result <- n.freecalc(65,0.95,c=1,se=0.95,sp=0.99,pstar=0.05, minSpH=0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.c.freecalc
expected_result <- c(27, 0.9555394, 0.9999933, 120, 1, 0.1)
observed_result <- n.c.freecalc(120,0.95,c=5,se=0.9,sp=0.99,pstar=0.1, minSpH=0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - adj.risk
expected_result <- c(3.125, 1.875, 0.625)
observed_result <- adj.risk(c(5, 3, 1), c(0.1, 0.1, 0.8))
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - epi.calc
expected_result <- c(0.0625, 0.0375, 0.0125)
observed_result <- epi.calc(0.02, c(5, 3, 1), c(0.1, 0.1, 0.8))[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rp.bin
expected_result <- 0.9430056
observed_result <- sep.rb.bin(0.1, c(5, 3, 1), c(0.1, 0.1, 0.8), c(5, 5, 5), 0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rp.hypergeo
expected_result <- 0.9623211
observed_result <- sep.rb.hypergeo(0.1, c(5, 3, 1), c(10, 10, 80), c(5, 5, 5), 0.9)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rb.bin.varse
rg<- c(1, 1, 2, 2)
se<- c(0.92, 0.85, 0.92, 0.85)
n<- c(80, 30, 20, 30)
df<- data.frame(rg, se, n)
expected_result <- 0.9800017
observed_result <- sep.rb.bin.varse(0.01, c(5, 1), c(0.1, 0.9), df)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rb.hypergeo.varse
rg<- c(1, 1, 2, 2)
se<- c(0.92, 0.85, 0.92, 0.85)
n<- c(80, 30, 20, 30)
df<- data.frame(rg, se, n)
expected_result <- 0.9938065
observed_result <- sep.rb.hypergeo.varse(0.01, c(5, 1), c(200, 1800), df)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rb2.bin
pstar<- 0.01
rr1<- c(3, 1)
ppr1<- c(0.2, 0.8)
rr2<- rbind(c(4,1), c(4,1))
ppr2<- rbind(c(0.1, 0.9), c(0.3, 0.7))
se<- 0.8
n<- rbind(c(50, 20), c(20, 10))
expected_result <- 0.9611147
observed_result <- sep.rb2.binom(pstar, rr1, ppr1, rr2, ppr2, n, se)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rb2.hypergeo
pstar<- 0.01
rr1<- c(3, 1)
rr2<- rbind(c(4,1), c(4,1))
N<- rbind(c(100, 500), c(300, 1000))
n<- rbind(c(50, 20), c(20, 10))
se<- 0.8
expected_result <- 0.9487643
observed_result <- sep.rb2.hypergeo(pstar, rr1, rr2, N, n, se)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sep.rb2stage
pstar.c<- 0.02
pstar.u<- 0.1
rr.c<- c(5, 1)
ppr.c<- c(0.1, 0.9)
rr.u<- c(3, 1)
se<- 0.9
n<- cbind(rep(10, 50), rep(5, 50))    
rg<- c(rep(1, 30), rep(2, 20))
ppr.u<- cbind(rep(0.2, 50), rep(0.8, 50))
N<- cbind(rep(30, 50), rep(120, 50))
C<- 500        
expected_result <- 0.8990356
observed_result <- sse.rb.2stage(C=NA, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N=NA, n, rg, se)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))
expected_result <- 0.9581237
observed_result <- sse.rb.2stage(C, pstar.c, pstar.u, rr.c, ppr.c, rr.u, ppr.u, N, n, rg, se)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sse.combined
C<- c(300, 1200)
pstar<- 0.01
rr<- c(3,1)
ppr<- c(0.2, 0.8)
comp1<- data.frame(id=1:100, rg=c(rep(1,50), rep(2,50)), cse=rep(0.5,100)) 
comp2<- data.frame(id=seq(2, 120, by=2), rg=c(rep(1,25), rep(2,35)), cse=rep(0.8,60))
comp3<- data.frame(id=seq(5, 120, by=5), rg=c(rep(1,10), rep(2,14)), cse=rep(0.9,24))
sep<- list(comp1, comp2, comp3)
expected_result <- matrix(c(0.6858347, 0.6954111, 0.8051234, 0.8121286), ncol=2, byrow=T)
observed_result <- sse.combined(C, pstar, rr, sep = sep)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.rb
expected_result <- c(7,5,2)
observed_result <- n.rb(0.1, c(5, 3, 1), c(0.1, 0.10, 0.80), 
                        c(0.5, 0.3, 0.2), 0.9, 0.95)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.rb.varse
m<- rbind(c(0.8, 0.2), c(0.5, 0.5), c(0.7, 0.3))
expected_result <- matrix(c(52, 12, 64, 32, 32, 64, 22, 9, 31, 106, 53, 159),
                          nrow = 4, byrow = T)
observed_result <- n.rb.varse(0.01, c(5, 3, 1), c(0.1, 0.1, 0.8), 
                              c(0.4, 0.4, 0.2), c(0.92, 0.8), m, 0.95)[[1]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.ap
expected_result <- c(139, 246, 323, 369, 385)
observed_result <- n.ap(seq(0.1, 0.5, by = 0.1), 0.05)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - binom.agresti
expected_result <- c(25, 200, 0.132067, 0.08558957, 0.1785444, 0.95)
observed_result <- binom.agresti(25, 200)[1,1:6]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - binom.jeffreys
expected_result <- c(25, 200, 0.125, 0.08462743, 0.176126, 0.95)
observed_result <- binom.jeffreys(25, 200)[1,1:6]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - binom.cp
expected_result <- c(25, 200, 0.125, 0.08255234, 0.1789738, 0.95)
observed_result <- binom.cp(25, 200)[1,1:6]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - ap
expected_result <- c(25, 200, 0.125, 0.08611974, 0.1780143, 0.95)
observed_result <- ap(25, 200)[1,1:6]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - n.tp
expected_result <- 174
observed_result <- n.tp(0.1, 0.9, 0.99, 0.05)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - sd.tp
expected_result <- c(0.05475727, 0.07537308, 0.08971191, 0.10049744,
                     0.10879167, 0.11513428, 0.11983556, 0.12308372,
                     0.12499211, 0.12562180)
observed_result <- sd.tp(1:10, 20, 0.9, 0.99)
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - tp.normal
expected_result <- c(0.2228464, 0.2144773, 0.2312155)
observed_result <- tp.normal(25, 120, 0.9, 0.99)[[2]]
stopifnot(all(abs(observed_result - expected_result) < tol))

## Test case - tp
expected_result <- c(0.02941176, 0, 0.09291029)
observed_result <- tp(25, 200, 0.95, 0.9, "c-p")[[2]]
stopifnot(all(abs(observed_result - expected_result) < tol))


