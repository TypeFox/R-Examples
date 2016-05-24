load("/data/SEER/00/pops.RData") # this loads in pops
pyf=pym=vector(3,mode="list"); 
for (i in 0:18) { for (r in 1:2) {
		pym[[r]][i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(poprace==r)]))
		pyf[[r]][i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(poprace==r)]))}
	pym[[3]][i+1]=with(pops,sum(population[(popsex==1)&(popage==i)&(poprace>2)]))
	pyf[[3]][i+1]=with(pops,sum(population[(popsex==2)&(popage==i)&(poprace>2)])) }

load("/data/SEER/00/lymyleuk.RData") # this loads in DF
casesf=casesm=vector(3,mode="list"); 
for (i in 1:3) {
	if (i==3) d=DF[(DF$histo2==9863)&(DF$numprims==1)&(DF$race>2)&(DF$race<98),] else 
		d=DF[(DF$histo2==9863)&(DF$numprims==1)&(DF$race==i),] 
	casesm[[i]]=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
	casesf[[i]]=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts}	

(DF=data.frame(age=rep(c(0.5,3,seq(7.5,87.5,5)),6),race=rep(c("white","black","asian"),each=19,times=2),
	sex=rep(c("male","female"),each=57),cases=c(unlist(casesm),unlist(casesf)),py=c(unlist(pym),unlist(pyf))))
lapply(DF,class)
DF=subset(DF,age>20)
DF=transform(DF,agec=age-55)  # center so that intercept is incidence at age 55
head(DF)

summary(lmf<-glm(cases~agec*race*sex+offset(log(py)),family=poisson,data=DF))
# as expected, sex has no significant interactions with anything so move it to additive
summary(lmf<-glm(cases~agec*race+sex+offset(log(py)),family=poisson,data=DF))
# shows that Asians (reference group) do not differ significantly from blacks in slope, but do in intercept
# also shows that Asians differ from whites in both slope and intercept	
exp(0.383664) # confirms average M/F across races of ~1.5
exp(coef(lmf)[5]) # same thing: 5th coeff is male offset

summary(lmf<-glm(cases~agec+sex+offset(log(py)),family=poisson,data=subset(DF,race=="asian") ))
# confirms Asian slope of k= 0.025 in Figure 4 
exp(coef(lmf)[3]) # confirms Asian M/F = 1.7 in Figure 4
