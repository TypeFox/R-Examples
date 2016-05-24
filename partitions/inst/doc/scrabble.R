### R code from vignette source 'scrabble.Rnw'

###################################################
### code chunk number 1: set_cutdownproblem
###################################################



###################################################
### code chunk number 2: scrabble.Rnw:78-79
###################################################
do_from_scratch <- FALSE


###################################################
### code chunk number 3: load_library
###################################################
library(partitions)
library(polynom)
options(width=85)


###################################################
### code chunk number 4: chess_problem
###################################################
blockparts(c(Bishops=2,Knights=2,Rooks=2,Queens=1),3)


###################################################
### code chunk number 5: scrabble_tiles
###################################################
scrabble <-
c(a=9,b=2,c=2,d=4,e=12,f=2,g=3,h=2,i=9,j=1,k=1,l=4,m=2,n=6,o=8,p=2,q=1,r=6,s=4,t=6,u=4,v=2,w=2,x=1,y=2,z=1," "=2)


###################################################
### code chunk number 6: scrabble.Rnw:197-198
###################################################
scrabble


###################################################
### code chunk number 7: call_S
###################################################
S(scrabble,7)


###################################################
### code chunk number 8: call_S_noblanks
###################################################
S(scrabble[-27],7)/S(scrabble,7)


###################################################
### code chunk number 9: probOfRacks
###################################################
f <- function(a){prod(choose(scrabble,a))/choose(sum(scrabble),7)}
if(do_from_scratch){
  racks <- blockparts(scrabble,7)
  probs <- apply(racks,2,f)
  otarine_ans <- rep(names(scrabble), racks[, which.max(probs)])
  max_probs <- max(probs)
  mp <- min(probs)
  a <- floor(log10(mp))
  how_many_racks <- sum(probs==min(probs))
} else {
  load("answers.Rdata")
}


###################################################
### code chunk number 10: print_otarine_letters
###################################################
otarine_ans


###################################################
### code chunk number 11: print_max_probs
###################################################
max_probs


