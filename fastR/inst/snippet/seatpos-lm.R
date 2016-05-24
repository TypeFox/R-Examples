require(faraway)
data(seatpos, package="faraway")
seatpos.lm1=lm(hipcenter ~ ., data=seatpos)
###hop:3-9
summary(seatpos.lm1)
