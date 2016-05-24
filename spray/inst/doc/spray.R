### R code from vignette source 'spray.Rnw'

###################################################
### code chunk number 1: setup
###################################################
ignore <- require(spray)


###################################################
### code chunk number 2: define_S
###################################################
library("spray")
M <- matrix(c(0,0,0,1,0,0,1,1,1,2,0,3),ncol=3)
M
S1 <- spray(M, 1:4)
S1


###################################################
### code chunk number 3: use_set_method_S
###################################################
S1[diag(3)] <- -3
S1


###################################################
### code chunk number 4: define_S1
###################################################
S2 <- spray(matrix(c(6,0,1,7,0,1,8,2,3),nrow=3), c(17,11,-4))
S2
S1+S2


###################################################
### code chunk number 5: demonstrate_addrepeats
###################################################
spray(matrix(0:5,8,3),addrepeats=TRUE)


###################################################
### code chunk number 6: show_a_multivariate_polynomial
###################################################
options(polyform=TRUE)
S1


###################################################
### code chunk number 7: exhibit.multiplication
###################################################
S1*S2


###################################################
### code chunk number 8: symbolic
###################################################

x <- lone(1,2)
y <- lone(2,2)
(x+y)^3
(1+x+y)^3


###################################################
### code chunk number 9: exhibit.negative.indices
###################################################
S1[0,-1,-2] <- 1
S1


###################################################
### code chunk number 10: exhibit.function.coercion
###################################################
f <- as.function(S1)
f(matrix(1:9,3,3))


###################################################
### code chunk number 11: exhibit_asum
###################################################
asum(S1,2)


###################################################
### code chunk number 12: substitute_show
###################################################
options(polyform=TRUE)
homog(3,3)
subs(homog(3,3),dims=2,1.5)


###################################################
### code chunk number 13: knight_generating_function
###################################################
chess_knight <- 
  spray(matrix(
      c(1,2,1,-2,-1,2,-1,-2,2,1,2,-1,-2,1,-2,-1),
      byrow=TRUE,ncol=2))
options(polyform=FALSE)
chess_knight


###################################################
### code chunk number 14: knight_six_moves
###################################################
constant(chess_knight^6,drop=TRUE)


###################################################
### code chunk number 15: define.d.dimensional.knight
###################################################
knight <- function(d){
  n <- d*(d-1)
  out <- matrix(0,n,d)
  out[cbind(rep(seq_len(n),each=2),c(t(which(diag(d)==0,arr.ind=TRUE))))] <- seq_len(2)
  spray(rbind(out,-out,`[<-`(out,out==1,-1),`[<-`(out,out==2,-2)))
}


###################################################
### code chunk number 16: dnightmoves
###################################################
constant(knight(4)^6, drop=TRUE)


###################################################
### code chunk number 17: dnightmoves_can_wait
###################################################
constant((1+knight(4))^6, drop=TRUE)


