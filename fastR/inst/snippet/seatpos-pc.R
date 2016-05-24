pc=with(seatpos,princomp(cbind(HtShoes,Ht,Seated,Arm,Thigh,Leg),
    scores=T))
summary(pc, loadings=T)
seatpos.lmpc <-lm(hipcenter ~ Age + Weight + pc$scores[,1], seatpos )
###hop:3-9
summary(seatpos.lmpc)
vif(seatpos.lmpc)
