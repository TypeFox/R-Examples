### R code from vignette source 'hyperdirichlet.Rnw'

###################################################
### code chunk number 1: set_seed_chunk
###################################################
set.seed(0)


###################################################
### code chunk number 2: time_saver
###################################################
calc_from_scratch <- FALSE


###################################################
### code chunk number 3: hyperdirichlet.Rnw:105-106
###################################################
ignore <- require(hyperdirichlet,quietly=TRUE)


###################################################
### code chunk number 4: calculatePriors
###################################################
g <- function(p){
  out <- uniform(3)
  out[2] <- p
  out[3] <- p
  out[5] <- p
  out[4] <- -p
  out[6] <- -p
  out[7] <- -p
  out
}

h <- function(p){
  out <- uniform(3)
  out[5] <- p
  out[3] <- p
  out[7] <- -2*p
  out
}

plot_quality <- 233

if(calc_from_scratch){
  z.prior <- triplot(g(0.5),l=plot_quality)
  z.posterior <- triplot(g(0.5) + mult_restricted_obs(3,1:2,c(7,3)),l=plot_quality)
  z.prior2 <- triplot(h(0.5),l=plot_quality)
  z.posterior2 <- triplot(h(0.5) + mult_restricted_obs(3,1:2,c(7,3)),l=plot_quality)
} else {
  load("vig.Rdata") 
}
  


###################################################
### code chunk number 5: plotPriors_openfile
###################################################
png("plotPriors.png",width=800,height=800)


###################################################
### code chunk number 6: plotPriors_drawfigure
###################################################
  layout(matrix(1:4,2,2))
  par(xpd=NA)

f <- function(jj,title, levels = -(0:10)*exp(1)){
  par(mai=c(0,0,0,0)+0.2)
  jj <- jj-max(jj,na.rm=TRUE)
  image(jj, asp=sqrt(3)/2, axes=FALSE,main=title)
  contour(jj, drawlabels=FALSE, add=TRUE, levels = levels)
  segments(0   , 0 , 0.5 , 1 , lwd=2)
  segments(0.5 , 1 , 1   , 0 , lwd=2)
  segments(1   , 0 , 0   , 0 , lwd=2)
  text(0.95,-0.05,"p1")
  text(0.05,-0.05,"p3")
  text(0.4, 0.97 , "p2")
}

f(pmin(z.prior,5), "(a)", levels= -c(2.5 , 3, 10))
f(z.posterior, "(b)")
f(pmin(z.prior2,7), "(c)")
f(z.posterior2, "(d)")


###################################################
### code chunk number 7: plotPriors_close
###################################################
null <- dev.off()


###################################################
### code chunk number 8: loadChess
###################################################
data("chess")
(w <- as.hyperdirichlet(chess))


###################################################
### code chunk number 9: hyperdirichlet.Rnw:568-574
###################################################
if(calc_from_scratch){
  w <- as.hyperdirichlet(w , calculate_NC = TRUE)
  z.chess <- triplot(w,l=plot_quality)
} else {
    ignore <- 0
}


###################################################
### code chunk number 10: plotChess_dostuff
###################################################
png("plotChess.png",width=800,height=800)
  par(xpd=NA)
  image(z.chess, asp=sqrt(3)/2, axes=FALSE)
  contour(z.chess, drawlabels=FALSE, add=TRUE)
  segments(0   , 0 , 0.5 , 1 , lwd=4)
  segments(0.5 , 1 , 1   , 0 , lwd=4)
  segments(1   , 0 , 0   , 0 , lwd=4)
  text(0.95,-0.05,"Topalov")
  text(0.05,-0.05,"Karpov")
  text(0.4, 0.97 , "Anand")
dev.off()


###################################################
### code chunk number 11: printChessNC
###################################################
w


###################################################
### code chunk number 12: calculate_chess_test
###################################################
jj.qwerty <- dhyperdirichlet(rep(1/3,3), w)
f <- function(p){dhyperdirichlet(p, w) > jj.qwerty}
if(calc_from_scratch){
  jj.test.ch <- probability(w,disallowed=f)
  }


###################################################
### code chunk number 13: print_chess_test
###################################################
jj.test.ch


###################################################
### code chunk number 14: cheat.TgtA
###################################################
T.lt.A <- function(p){p[1] < p[2]}
if(calc_from_scratch){
  jjT <- probability(w , disallowed = T.lt.A)
}


###################################################
### code chunk number 15: printTgtA
###################################################
jjT


###################################################
### code chunk number 16: dataIcons
###################################################
data("icons")


###################################################
### code chunk number 17: CalculateIconStats
###################################################
ic <- as.hyperdirichlet(icons)
if(calc_from_scratch){
  m.free <- maximum_likelihood(ic)
}
m.free$MLE <- round(m.free$MLE,5)


###################################################
### code chunk number 18: printIconStats
###################################################
m.free


###################################################
### code chunk number 19: CalculateIconStatsH1
###################################################
f1 <- function(p){p[1] > 1/6}
if(calc_from_scratch){
  m.f1 <- maximum_likelihood(ic , disallowed=f1)
}
m.f1$MLE <- round(m.f1$MLE,5)


###################################################
### code chunk number 20: printIconStatsH1
###################################################
m.f1


###################################################
### code chunk number 21: H1_printer
###################################################
m.free$support - m.f1$support


###################################################
### code chunk number 22: CalculateIconStatsH2
###################################################
f2 <- function(p){p[1] > max(p[-1])}
if(calc_from_scratch){
  m.f2 <- maximum_likelihood(ic , disallowed=f2)
}
m.f2$MLE <- round(m.f2$MLE,5)


###################################################
### code chunk number 23: printIconStatsH1
###################################################
m.f2


###################################################
### code chunk number 24: H2_printer
###################################################
m.free$support - m.f2$support


###################################################
### code chunk number 25: CalculateIconStatsH3
###################################################
f3 <- function(p){sum(p[5:6]) < 1/3}
if(calc_from_scratch){
  m.f3 <- maximum_likelihood(ic , start_p = c(0,0,0,0,1,1), disallowed=f3)
}
m.f3$MLE <- round(m.f3$MLE, 5)


###################################################
### code chunk number 26: printIconStatsH3
###################################################
m.free$support - m.f3$support


###################################################
### code chunk number 27: CalculateIconStatsH4
###################################################
f4 <- function(p){min(p[1:4])>max(p[5:6])}
if(calc_from_scratch){
m.f4 <- maximum_likelihood(ic , start_p = c(0,0,0,0,1,1), disallowed=f4)
}
m.f4$MLE <- round(m.f4$MLE,5)


###################################################
### code chunk number 28: loadVB
###################################################
data("volleyball")


###################################################
### code chunk number 29: cheat_Volleyball_H0
###################################################
if(calc_from_scratch){
v.HF <- maximum_likelihood(vb_synthetic , start_p = 
                 c(0.30627, 0.17491, 0.10036,
                   0.09195, 0.07228, 0.08449,
                   0.05898, 0.04303, 0.06773)
                 )
}


###################################################
### code chunk number 30: cheat_Volleyball_HA
###################################################
o <- function(p){all(order(p[1:4])==4:1)}
v.HA <- maximum_likelihood(vb_synthetic, disallowed=o , start_p =
                 c(0.14522, 0.14522, 0.14522,
                   0.14522, 0.14523, 0.08623, 
                   0.06744, 0.05365, 0.06657)
                 )


###################################################
### code chunk number 31: showVB
###################################################
round(v.HF$MLE , 4)


###################################################
### code chunk number 32: calcVBerror
###################################################
zipf <- function(n){jj <- 1/seq_len(n); return(jj/sum(jj))}
vb.error <- round(zipf(9) - v.HF$MLE , 4)


###################################################
### code chunk number 33: printVBerror
###################################################
vb.error


###################################################
### code chunk number 34: print_VB_HA
###################################################
v.HF$support - v.HA$support


###################################################
### code chunk number 35: print_data_doubles
###################################################
data("doubles")


###################################################
### code chunk number 36: CalculateDoublesStats
###################################################
if(calc_from_scratch){
  doubles_result <- maximum_likelihood(doubles)$support - maximum_likelihood(doubles,zero=5)$support
}


###################################################
### code chunk number 37: printIconStats
###################################################
doubles_result


