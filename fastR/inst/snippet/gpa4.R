###hop:3-9
gpa.lm5<- lm(gpa~act+satv,gpa); summary(gpa.lm5)
###hop:3-9
gpa.lm6<- lm(satv~act,gpa); summary(gpa.lm6)
