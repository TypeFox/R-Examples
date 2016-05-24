# Bees data:
data(Bees)
# contrasts to resolve trtime into treat, time, treat:time

cont<-matrix(scan(zz<-textConnection("
   10 -1 -1 -1 -1 -1  -1 -1 -1 -1 -1      # '0-treat'  
    0  1  1  1  1  1  -1 -1 -1 -1 -1      # 'treat'
    0 -2 -1  0  1  2  -2 -1  0  1  2      # 'time.1'   
    0  2 -1 -2 -1  2   2 -1 -2 -1  2      # 'time.2'   
    0 -1  2  0 -2  1  -1  2  0 -2  1      # 'time.3'   
    0  1 -4  6 -4  1   1 -4  6 -4  1      # 'time.4'   
"), comment.char="#"),11)
close(zz)

cnames <- c( '0-treat', 'treat', 'time.1', 'time.2', 'time.3', 'time.4')  

 
# generate interaction contrasts for treat:time
for (t in 3:6) {
	 cont <- cbind(cont, cont[,2]*cont[,t])
	 cnames <-c(cnames, paste("treat", cnames[t],sep=":"))
}
colnames(cont)<-cnames
rownames(cont)<- levels(Bees$trtime)
cont

contrasts(Bees$trtime) <-cont

# Tests of linear hypotheses using these contrasts
bees.mod1 <- lm(cbind(Iz,Iy) ~ caste*trtime, data=Bees)
# get coefficient names for linearHypothesis
coefs <- rownames(coef(bees.mod1))
linearHypothesis(bees.mod1,"trtimetreat" , title="Treat")

print(linearHypothesis(bees.mod1, coefs[grep("^trtimetime", coefs)], title="Time"),SSP=FALSE)
print(linearHypothesis(bees.mod1, coefs[grep("^trtimetreat:time", coefs)], title="Treat:Time"),SSP=FALSE)

heplot(bees.mod1, xlab="Iz: Ovarian development", ylab="Iz: Ovarian reabsorption",
		size="effect", main="Bees: ~caste*trtime",
		hypotheses=list(Time=coefs[grep("^trtimetime", coefs)])
)
