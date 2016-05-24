seatpos.lm2=lm(hipcenter ~ Age + Weight + Ht, seatpos )
###hop:3-9
summary(seatpos.lm2)
vif(seatpos.lm2)
