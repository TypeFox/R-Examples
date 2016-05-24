###hop:3-9
gpa.lm2<- lm(satm~satv+act,gpa); summary(gpa.lm2)
###hop:3-9
gpa.lm3<- lm(satm~satv,gpa); summary(gpa.lm3)
###hop:3-9
gpa.lm4<- lm(satm~act,gpa); summary(gpa.lm4)
