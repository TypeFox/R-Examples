

f1 <- as.factor(LETTERS[floor(runif(100,1,3))])
f1b <- as.factor(LETTERS[floor(runif(100,1,3))])
f2 <- as.factor(LETTERS[floor(runif(100,1,5))])
f3 <- as.factor(LETTERS[floor(runif(100,1,10))])

l1 <- f1=="A"
l2 <- f1b=="A"

o1 <- as.ordered(as.factor(LETTERS[floor(runif(100,1,3))]))
o2 <- as.ordered(as.factor(LETTERS[floor(runif(100,1,5))]))
o3 <- as.ordered(as.factor(LETTERS[floor(runif(100,1,10))]))

i1 <- d1 <- as.numeric(floor(runif(100,1,3)))
i2 <- d2 <- as.numeric(floor(runif(100,1,5)))
i3 <- d3 <- as.numeric(floor(runif(100,1,10)))
class(d1) <- c("discrete","numeric")
class(d2) <- c("discrete","numeric")
class(d3) <- c("discrete","numeric")
class(i1) <- "integer"
class(i2) <- "integer"
class(i3) <- "integer"

n1 <- c1 <- rnorm(100)
n2 <- c2 <- rnorm(100,1:100)
n3 <- c3 <- rnorm(100,1:100)
class(c1) <- c("continuous","numeric")
class(c2) <- c("continuous","numeric")
class(c3) <- c("continuous","numeric")

df <- data.frame(f1,i1,f2,o2,n2,n1)


malade <- rep(c("Vivant","Mort"),c(49,51))
traitement <- rep(c("A","B","A","B"),c(22,31,27,20))

concours <- rep(c("Reussite","Echec"),c(40,160))
sexe <- rep(c("H","F","H","F"),c(10,30,20,140))
bacType <- rep(c("A","B","C","D","A","B","C","D"),c(0,1,35,4,20,40,60,40))
bacType2 <- rep(c("A","B","C","D","H","F","St","E","M","Sttk"),20)
cheveux <- c("noir","blond","chatain","roux")[floor(runif(200,1,5))]

bacMension <- ordered(rep(c("P","AB","B","TB","P","AB","B","TB"),c(5,5,15,13,95,45,15,7),
                          levels=c("P","AB","B","TB")))
bacMension2 <- ordered(rep(c("P","P+","AB","AB+","B","B+","TB","TB+","P","P+","AB","AB+","B","B+","TB","TB+"),c(2,3,3,4,7,8,9,4,40,45,20,25,18,5,3,4),
                           levels=c("P","P+","AB","AB+","B","B+","TB","TB+")))

nbRedoublement <- rpois(200,0.5)
nbRedoublement2 <- rpois(200,3)
nbRedoublement2[nbRedoublement2>8]<-8
class(nbRedoublement) <- c("discrete","numeric")
class(nbRedoublement2) <- c("discrete","numeric")

taille <- rnorm(200,170,5)
poids <- rnorm(200,65,10)
