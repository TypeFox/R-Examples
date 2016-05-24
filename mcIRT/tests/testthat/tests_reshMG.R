

################## CREATE DATA ################################


set.seed(2046)
x1 <- rnorm(10*5,0,1.7)
zeta <- tapply(x1,rep(1:10,each=5), function(x) round(x - mean(x),3))

# lambda parameters
x2 <- rnorm(10*5,0,1.1)
lam  <- tapply(x2,rep(1:10,each=5), function(x)
{
  sort(round(x - mean(x),3),decreasing=FALSE)
})

# create a parlist as a first step
ParList <- mapply(function(one,two)
  {
    x1 <- c(one,two)
    names(x1) <- paste0(rep(c("zeta","lam"),each=length(one)),1:length(one))
    return(x1)
  },one=zeta,two=lam,SIMPLIFY=FALSE)

nper <- 3000
names(ParList) <- paste0("item",1:length(ParList))

erglist     <- vector(mode="list",length=100)
erglistTIME <- vector(mode="list",length=100)

perp1 <- rnorm(nper,0,1)
perp2 <- rnorm(nper,0.3,1)

simdat1 <- NRM.sim(ParList,perp1)
simdat2 <- NRM.sim(ParList,perp2)
simdatall <- rbind(simdat1,simdat2)

simdatallg <- data.frame(GROUP=factor(rep(c("A","B"),each=nper),levels=c("B","A")),simdatall)
simdatallg1 <- simdatallg

# one category no obs in one group
simdatallg1[simdatallg1$GROUP == "A","item1"] <- ifelse(simdatallg1[simdatallg1$GROUP == "A","item1"] == 2,1,simdatallg1[simdatallg1$GROUP == "A","item1"])

# one category skipped (no "1" cat)
simdatallg2 <- simdatallg
simdatallg2$item2 <- ifelse(simdatallg2$item1 == 1,5,simdatallg2$item1)

### designs
mydesign <- designTemp(ngru=2,nit=10,TYPE="NRM")
mydesign2 <- designTemp(ngru=2,nit=9,TYPE="NRM")
mydesign3 <- designTemp(ngru=2,nit=10,TYPE="NRM")
mydesign4 <- designTemp(ngru=2,nit=10,TYPE="NRM")
mydesign5 <- designTemp(ngru=2,nit=10,TYPE="NRM")

mydesign[[1]][2,1] <- 3 
mydesign3[[1]][2,1] <- 2


mydesign4[[2]][1,2:5] <- 2
mydesign4[[1]][2,7] <- 2

mydesign5[[1]][2,7] <- 2

test <- reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign4)

################## CREATE DATA - fin ################################



test_that("reshMG throws an error/warning when necessary", {
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(2,11),echo=FALSE), throws_error()) 
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(5,10),echo=FALSE), throws_error()) 
  expect_that(reshMG(simdatallg1,items=2:11,groups=1,correct=rep(3,10),echo=FALSE), throws_error())
  expect_that(reshMG(simdatallg2,items=2:11,groups=1,correct=rep(3,10),echo=FALSE), throws_error())
  expect_that(reshMG(simdatallg2,items=2:12,groups=1,correct=rep(3,10),echo=FALSE), throws_error())
  expect_that(reshMG(simdatallg2,items=2:11,groups=1:2,correct=rep(3,10),echo=FALSE), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE), gives_warning("small number"))
  expect_that(reshMG(simdatallg2,items=2:11,groups=12,correct=rep(3,10),echo=FALSE), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign), throws_error("The numbers inside the matrix should refer to the groups!"))
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign2), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign[[1]]), throws_error("design"))
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign3,TYPE="AAA"), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign3,TYPE="AAA"), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign3,TYPE="BOCK"), throws_error())
  expect_that(reshMG(simdatallg,items=2:11,groups=1,correct=rep(3,10),echo=FALSE,design=mydesign4)$design, is_identical_to(mydesign5))
})






