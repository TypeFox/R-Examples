prop.test(39,215,.15)
binom.test(39,215,.15)
lewitt.machin.success <- c(9,4)
lewitt.machin.total <- c(12,13)
prop.test(lewitt.machin.success,lewitt.machin.total)
matrix(c(9,4,3,9),2) 
lewitt.machin <- matrix(c(9,4,3,9),2)
fisher.test(lewitt.machin)
chisq.test(lewitt.machin)
caesar.shoe
caesar.shoe.yes <- caesar.shoe["Yes",]
caesar.shoe.total <- margin.table(caesar.shoe,2)    
caesar.shoe.yes                     
caesar.shoe.total                          
prop.test(caesar.shoe.yes,caesar.shoe.total)
prop.trend.test(caesar.shoe.yes,caesar.shoe.total)
caff.marital <- matrix(c(652,1537,598,242,36,46,38,21,218
,327,106,67),
nrow=3,byrow=T)
colnames(caff.marital) <- c("0","1-150","151-300",">300")
rownames(caff.marital) <- c("Married","Prev.married","Single")
caff.marital
chisq.test(caff.marital)
chisq.test(caff.marital)$expected
chisq.test(caff.marital)$observed
E <- chisq.test(caff.marital)$expected   
O <- chisq.test(caff.marital)$observed   
(O-E)^2/E
attach(juul)
chisq.test(tanner,sex)            
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
