population <-
function(max.repeat, max.time, Na, sex.ratio=FALSE, N.loci=10, N.allele=20, mut.rate=FALSE, one.step=0.9,
n.recruit.fem, skew.recruit.fem=0, skew.recruit.mal=0, delay=1,
surv.juv=FALSE, surv.ad.fem, surv.ad.mal=surv.ad.fem, const=TRUE, quiet=TRUE) {

## check generation length calculation
## include maturity delay

# missing parameters
if (missing(max.repeat) | missing(max.time) | missing(Na) | missing(n.recruit.fem) | missing(surv.ad.fem)) {
stop(cat(message <- paste("important parameters are missing; \n",
"obligatory parameters are: population(max.repeat, max.time, Na, N.loci, N.allele, n.recruit.fem, surv.ad.fem, surv.ad.mal) \n",
"see ?population for further information \n")))
}

# check number of repeats
if (max.repeat<2 & quiet == TRUE) {
stop("max.repeat must be greater than 1; for results without repetitions set quiet=FALSE")
}

# check sex ratio I
if (sex.ratio>1) {
stop("sex.ratio is the proportion of females, please insert a value between 0 and 1")
}

# check sex ratio II
if (const==TRUE & sex.ratio!=FALSE & sex.ratio*Na*n.recruit.fem+sex.ratio*Na*surv.ad.fem+(1-sex.ratio)*Na*surv.ad.mal < Na) {
stop("WARNING: too few females in the population to keep population size constant")
}

# check sex ratio III
vec <- vector(length=max.time)

if (sex.ratio==FALSE){
	vec[1] <- Na*0.5
	for(i in 1:length(vec)){
		vec[i+1]<-trunc(vec[i]*(surv.ad.fem+n.recruit.fem*0.5*surv.juv))
	}
}

else{
	vec[1] <- Na*sex.ratio
	for(i in 1:length(vec)){
		vec[i+1]<-trunc(vec[i]*(surv.ad.fem+n.recruit.fem*sex.ratio*surv.juv))
	}
}

#if (const==FALSE & length(vec[vec==0])>0) {
#stop(cat(message <- paste("WARNING: after ",length(vec)-length(vec[vec==0])+1," years no females are in the population any more \n")))
#}


# check survival
if (surv.juv>1 | surv.ad.fem>1 | surv.ad.mal>1) {
stop("survival rates express the proportion of individuals that survive from one year to the next, please insert values between 0 and 1")
}

# check constant population I
if (surv.juv==FALSE & const == FALSE) {
stop("If the population size is not assumed to be constant, please insert a value for surv.juv")
}

# check constant population II
if (surv.juv>0 & const == TRUE) {
stop("If surv.juv is included, the population size cannot be assumed to be constant. Set const=FALSE.")
}


## expected generation length at start point (approx.)
s <- mean(c(surv.ad.fem,surv.ad.mal))
glength <- sum(surv.ad.fem^c(0:100))

whole.rep <- as.data.frame(matrix(ncol=max.repeat*3,nrow=max.time))

## starting the repeat loop
for (r in 1:max.repeat) {


## initializing the population

# adult matrix

adults <- as.data.frame(matrix(NA, ncol=2*N.loci+2, nrow=Na))
seq.loc <- rep(1:N.loci, each = 2)
seq.all <- rep(c(1,2),N.loci)
print.loc <- paste("l",seq.loc,"all",seq.all, sep=".")
colnames (adults) <- c("sex", "age", print.loc)
rownames (adults) <- NULL

if(sex.ratio==FALSE)
{
	if (Na%%2==0){
		adults[,1] <- rep(c("M","F"),each=Na/2)
	}
	else {
		adults[,1] <- c(rep(c("M","F"),each=Na/2),sample(c("M","F"),1))
	}
	adults[,2] <- round(rtnorm(Na,glength,glength,1,Inf),0)# adults are at least at age 1, mean and sd = 1 generation
	adults[,c(3:dim(adults)[2])] <- round(rtnorm((dim(adults)[2]-2)*Na,sum(1:N.allele)/N.allele,N.allele/5,1,N.allele),0)*10
}

else
{
	adults[,1] <- c(rep("M",round(Na*(1-sex.ratio))),rep("F",Na-round(Na*(1-sex.ratio))))
	adults[,2] <- round(rtnorm(Na,glength,glength,1,Inf),0)# adults are at least at age 1, sd = 1 generation
	adults[,c(3:dim(adults)[2])] <- round(rtnorm((dim(adults)[2]-2)*Na,sum(1:N.allele)/N.allele,N.allele/5,1,N.allele),0)*10
}


# define single parameters and matices

n.adults <-  dim(adults)[1]
genotypes.adults <- adults[,3:dim(adults)[2]]
Het <- 1-sum(genotypes.adults[,seq(1,dim(genotypes.adults)[2],2)] == genotypes.adults[,seq(2,dim(genotypes.adults)[2],2)])/ (n.adults*N.loci)# whole Het
wt1 <- round(rtnorm(sum(adults$sex=="F"& adults$age>=delay),n.recruit.fem,sqrt(skew.recruit.fem),1,Inf))
wt <- wt1/sum(wt1)
Gen.length <- round(weighted.mean(adults$age[adults$sex=="F" & adults$age>=delay],wt))# generation length as mean age of reproducing females (weighted according to repoductive success)
Nind <- n.adults# adult population size


### running the population

for (i in 2:max.time)

{

# sexes of reproducing parents
females <- which(adults$sex=="F" & adults$age>=delay)
males <- which(adults$sex=="M" & adults$age>=delay)


## genotypes of the offspring

# offspring can have different fathers, random mating is assumed

# n.recruit offspring per female
# round(rtnorm(#juv (=mean*# females), mean, skew (=interindividual variance),1,Inf) for females, the reproductive skew is not that strong, most females get at least one offspring
dist.fem <- round(rtnorm(length(females),n.recruit.fem,sqrt(skew.recruit.fem),1,Inf))
genotypes.fem <- genotypes.adults[c(rep(females,dist.fem)),]

# n.recruit offspring per female per year including reproductive skew
num.juv <- sum(dist.fem)
juv <- data.frame(sex=rep("J",num.juv), age=0) # sex is determined later, age=0 per def.


# distribution of males: left skewed, between equal dist. and few males get most and many males get no offspring (rtnorm[0,Inf])
# using the probability as the reproductive skew

# males: sample males of data set randomly, amount=as much as offspring produced, with replace = true because of random mating

dist.mal <- (rtnorm(length(males),num.juv/length(males),sqrt(skew.recruit.mal),0,Inf))
prob.dist.mal <- dist.mal/sum(dist.mal)
genotypes.mal <- genotypes.adults[males[sample(1:length(males),num.juv,replace=TRUE,prob=prob.dist.mal)],]


# select alleles from females and build the genotypes of the offspring using bits

sel.alleles.fem <- matrix(sample(c(TRUE,FALSE),num.juv*N.loci, replace=T), ncol=N.loci)# TRUE/FALSE matrix Njuv ro x #loci col
sel.alleles.fem <- cbind(sel.alleles.fem,!sel.alleles.fem)# sel.allels and the opposite
sel.alleles.fem <- sel.alleles.fem[,rep(1:N.loci, each=2)+c(0,N.loci)]# sort cols by sel.alleles1,!sel.alleles1,sel.alleles2,!sel.alleles2...
sel.alleles.fem <- as.vector(t(sel.alleles.fem))# t() to convert the matrix, as.vector is writing columnwise
sel.alleles.fem <- as.bit(sel.alleles.fem)

genotypes.juv <-as.matrix(t(genotypes.fem))
gen.juv.vec <- as.vector(genotypes.juv)


# by now, the offspring have genotypes only from the females. We take that allele from the female, where the T/F Matrix is true. The other one has to be raplaced by an allele from the father. 
# That can be either the first or the second one. Therefore, I create a new T/F matrix for males.

sel.alleles.mal <- matrix(sample(c(TRUE,FALSE),num.juv*N.loci, replace=T), ncol=N.loci)# TRUE/FALSE matrix Njuv ro x #loci col
sel.alleles.mal <- cbind(sel.alleles.mal,!sel.alleles.mal)# sel.allels and the opposite
sel.alleles.mal <- sel.alleles.mal[,rep(1:N.loci, each=2)+c(0,N.loci)]# sort cols by sel.alleles1,!sel.alleles1,sel.alleles2,!sel.alleles2...
sel.alleles.mal <- as.vector(t(sel.alleles.mal))# t() to convert the matrix, as.vector is writing columnwise
sel.alleles.mal <- as.bit(sel.alleles.mal)

gen.mal.vec <- as.vector(as.matrix(t(genotypes.mal)))
gen.mal.sel <- gen.mal.vec[as.which(sel.alleles.mal)]

# replace always one allele by the one from the father
gen.juv.vec[as.which(!sel.alleles.fem)] <- gen.mal.sel

# include mutation
if(mut.rate==FALSE){
	genotypes.juv <- as.data.frame(t(matrix(gen.juv.vec,nrow=2*N.loci)))
	colnames(genotypes.juv) <- print.loc
	juv <- cbind(juv, genotypes.juv)
}

else {
	if(mut.rate==TRUE){
		mut.rate=5e-3	# standard mutation rate
	}
	else{
		mut.rate=mut.rate
	}
	help.vec <- 1:length(gen.juv.vec)
	prob.vec <- sample(c(one.step,1-one.step),length(gen.juv.vec),prob=c(one.step,1-one.step),replace=TRUE)
		
	if(mut.rate*length(gen.juv.vec)*(i-1)==as.integer(mut.rate*length(gen.juv.vec)*(i-1)) & mut.rate*length(gen.juv.vec)<1){
		g <- sample(help.vec,1)
		if (prob.vec[g]==one.step){
			gen.juv.vec[g] <- sample(c(gen.juv.vec[g]-1,gen.juv.vec[g]+1),1)
		}
		else{
			mut <- rpois(100,gen.juv.vec[g])
			mut.prob <- 1-table(mut)/sum(table(mut))
			gen.juv.vec[g] <- sample(sort(unique(mut)),1,prob=mut.prob)
		}
		
	}
	else if (mut.rate*length(gen.juv.vec)*(i-1)==as.integer(mut.rate*length(gen.juv.vec)*(i-1)) & mut.rate*length(gen.juv.vec)>=1){
		g <- sample(help.vec,mut.rate*length(gen.juv.vec),replace=FALSE)
		g.one.step <- help.vec[g][prob.vec[g]==one.step]
		mut.mat <- rbind(gen.juv.vec[g.one.step]+1,gen.juv.vec[g.one.step]-1)
		gen.juv.vec[g.one.step] <- apply(mut.mat,2,sample,1)
		
		g.multi.step <- help.vec[g][prob.vec[g]==1-one.step]
		for(m in g.multi.step){
			mut <- rpois(100,gen.juv.vec[m])
			mut.prob <- 1-table(mut)/sum(table(mut))
			gen.juv.vec[m] <- sample(sort(unique(mut)),1,prob=mut.prob)
		}
				
	}
	else{
		gen.juv.vec <- gen.juv.vec
	}
	genotypes.juv <- as.data.frame(t(matrix(gen.juv.vec,nrow=2*N.loci)))
	colnames(genotypes.juv) <- print.loc
	juv <- cbind(juv, genotypes.juv)
		
}

## survival to the next year and status change

# after reproduction, some juv and adults will survive, others wont, whereby survival will be set as constant over the years but different between age

# new adults (seperatet by sex, because of the assumption of a specific sex ratio)
all.fem <- which(adults$sex=="F")
all.mal <- which(adults$sex=="M")
ad.mal <- adults[all.mal[sample(1:length(all.mal),sample(round(rtnorm(100,surv.ad.mal*length(all.mal),sd=1,1,length(all.mal))),1))],]# surv.ad of the adult males from the previous year
ad.fem <- adults[all.fem[sample(1:length(all.fem),sample(round(rtnorm(100,surv.ad.fem*length(all.fem),sd=1,1,length(all.fem))),1))],]# surv.ad of the adult females from the previous year
adults <- rbind(ad.mal,ad.fem)

# survived juveniles
if(const==TRUE & sex.ratio==FALSE)
{
	n.juv <- Na-dim(adults)[1]
	juv.sample <- juv[sample(1:dim(juv)[1],n.juv),]
	
	if (n.juv%%2==0){
		juv.sample$sex <- c(rep(c("M","F"),each=n.juv/2))
	}
	else {
		juv.sample$sex <- c(rep(c("M","F"),each=n.juv/2),sample(c("M","F"),1))
	}
}

else if(const==TRUE & sex.ratio>0)
{
	n.juv <- Na-dim(adults)[1]
	juv.sample <- juv[sample(1:dim(juv)[1],n.juv),]
	juv.sample$sex <- c(rep("M",round(n.juv*(1-sex.ratio))),rep("F",n.juv-round(n.juv*(1-sex.ratio))))
}

else if(const==FALSE & sex.ratio==FALSE)
{
	n.juv <- sample(round(rtnorm(100,surv.juv*nrow(juv),sd=1,1,nrow(juv))),1)
	juv.sample <- juv[sample(1:nrow(juv),n.juv),]
	if (n.juv%%2==0){
		juv.sample$sex <- c(rep(c("M","F"),each=n.juv/2))
	}
	else {
		juv.sample$sex <- c(rep(c("M","F"),each=n.juv/2),sample(c("M","F"),1))
	}
}

else 
{
	n.juv <- sample(round(rtnorm(100,surv.juv*nrow(juv),sd=1,1,nrow(juv))),1)
	juv.sample <- juv[sample(1:nrow(juv),n.juv),]
	juv.sample$sex <- c(rep("M",round(n.juv*(1-sex.ratio))),rep("F",n.juv-round(n.juv*(1-sex.ratio))))
}

# change status, former juveniles become adults, adults age 1 year
adults <- rbind(adults,juv.sample)# new adults = survived adults + survived juv from last year
n.adults <-  nrow(adults)[1]# new # adults
adults$age <- adults$age+1# change status. age+1
rownames(adults) <- NULL

# remove everything important that is marked
rm(juv)
rm(juv.sample)
rm(n.juv)
rm(num.juv)

## calculating Heterozygosity
genotypes.adults <- adults[,3:dim(adults)[2]]

Het[i] <- 1-sum(genotypes.adults[,seq(1,dim(genotypes.adults)[2],2)] == genotypes.adults[,seq(2,dim(genotypes.adults)[2],2)])/(n.adults*N.loci)# whole Het
wt1 <- round(rtnorm(sum(adults$sex=="F"& adults$age>=delay),n.recruit.fem,sqrt(skew.recruit.fem),1,Inf))
wt <- wt1/sum(wt1)
Gen.length[i] <- round(weighted.mean(adults$age[adults$sex=="F" & adults$age>=delay],wt))# generation length as mean age of reproducing females (weighted according to repoductive success)
Nind[i] <- n.adults# adult population size (control)

if(nrow(adults[adults$sex=="F",])>4){
	num.fem=TRUE
	next
}
else {
num.fem =FALSE

fill.count <- max.time - length(Het)

Het[(max.time-fill.count+1):max.time] <- NA
Gen.length[(max.time-fill.count+1):max.time] <- NA
Nind[(max.time-fill.count+1):max.time] <- NA

break
}

if(length(adults$age[adults$age>delay&adults$sex=="F"])>4){
	num.fem=TRUE
	next
}
else {
num.fem=FALSE

fill.count <- max.time - length(Het)

Het[(max.time-fill.count+1):max.time] <- NA
Gen.length[(max.time-fill.count+1):max.time] <- NA
Nind[(max.time-fill.count+1):max.time] <- NA

break
}


}

# writing the files for each repeat
mat <- cbind(Het,Gen.length,Nind)

name.seq <- c(1:max.repeat)
name.mat <- paste(paste("mat",name.seq,"s",s,"b",n.recruit.fem, sep="_"),"txt",sep=".")
#write.table(mat,file = name.mat[r], col.names=TRUE, row.names=TRUE)

fill.whole <- seq(r,dim(whole.rep)[2],max.repeat)
whole.rep [,fill.whole] <- mat
colnames(whole.rep)[fill.whole] <- paste(colnames(mat),name.seq[r],sep="_")

message <- paste("I've done the repeat number",r,sep=" ","\n")
cat(message)

}

# extract the data to analyze
Het_pop <- whole.rep[,1:max.repeat]
Gen_pop <- whole.rep[,c((max.repeat+1):(2*max.repeat))]
Nind_pop <- whole.rep[,c((2*max.repeat+1):(3*max.repeat))]

## temporal reference of one year
if (max.repeat<2 & quiet==FALSE)
{
	mean.Nind.y <- Nind_pop
	mean.gen.y <- Gen_pop
	mean.het.y <- Het_pop
	years <- c(1:max.time)
	lm.het.y <- lm(log(mean.het.y[mean.het.y>0])~years[mean.het.y>0])
	slope.y <- coefficients(lm.het.y)[2]
	# calculate Ny from the simulation
	Ny_sim <- -1/(2*slope.y)
	# calculate Ny from the formula
	Ny_cal <- (-n.recruit.fem^2*Na)/(1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*s+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*s^2)
	#ranges
	slope.y.low <- NA
	slope.y.high <- NA
	Ny_sim.low <- NA
	Ny_sim.high <- NA
}

if (max.repeat>1)
{
	mean.Nind.y <- rowMeans(Nind_pop,na.rm=TRUE)
	mean.gen.y <- rowMeans(Gen_pop,na.rm=TRUE)
	mean.het.y <- rowMeans(Het_pop,na.rm=TRUE)
	years <- c(1:max.time)
	lm.het.y <- lm(log(mean.het.y[mean.het.y>0])~years[mean.het.y>0])
	slope.y <- coefficients(lm.het.y)[2]
	# calculate Ny from the simulation
	Ny_sim <- -1/(2*slope.y)
	# calculate Ny from the formula
	Ny_cal <- (-n.recruit.fem^2*Na)/(1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*s+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*s^2)
	# ranges 
	slope.y.range <- vector(length=max.repeat)
	for (r in 1:max.repeat) slope.y.range [r] <- coefficients(lm(log(Het_pop[,r][Het_pop[,r]>0])~years[Het_pop[,r]>0]))[2]
	slope.y.low <- slope.y-1.96*sd(slope.y.range)
	slope.y.high <- slope.y+1.96*sd(slope.y.range)
	Ny_sim.range <- -1/(2*slope.y.range)
	Ny_sim.low <- Ny_sim-1.96*sqrt(sum(Ny_sim.range)/(length(Ny_sim.range)-1)^2)
	Ny_sim.high <- Ny_sim+1.96*sqrt(sum(Ny_sim.range)/(length(Ny_sim.range)-1)^2)
}

## temporal reference one generation
if (max.repeat<2 & quiet==FALSE)
{
	mean.gen.g <- Gen_pop
	mean.Nind.g <- mean.Nind.y[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE)))]
	mean.het.g <- mean.het.y[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE)))]
	years.g <- 1:length(mean.het.g)
	gen.length <- mean(mean.gen.g,na.rm=TRUE)
	lm.het.g <- lm(log(mean.het.g[mean.het.g>0])~years.g[mean.het.g>0])
	slope.g <- coefficients(lm.het.g)[2]
	# calculate Ne from the simulation
	Ne_sim <- -1/(2*slope.g)
	# calculate Ne from the formula
	Ne_cal <- (-n.recruit.fem^2*Na)/((1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*s+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*s^2)*(gen.length))
	# ranges
	slope.g.low <- NA
	slope.g.high <- NA
	Ne_sim.low <- NA
	Ne_sim.high <- NA
	gen.length.low <- NA
	gen.length.high <- NA
}

if (max.repeat>1)
{
	mean.gen.g <- rowMeans(Gen_pop,na.rm=TRUE)
	mean.Nind.g <- mean.Nind.y[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE)))]
	mean.het.g <- mean.het.y[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE)))]
	years.g <- 1:length(mean.het.g)
	gen.length <- mean(mean.gen.g,na.rm=TRUE)
	lm.het.g <- lm(log(mean.het.g[mean.het.g>0])~years.g[mean.het.g>0])
	slope.g <- coefficients(lm.het.g)[2]
	# calculate Ne from the simulation
	Ne_sim <- -1/(2*slope.g)
	# calculate Ne from the formula
	Ne_cal <- (-n.recruit.fem^2*Na)/((1-n.recruit.fem^2+(-1.33+0.44*log(n.recruit.fem))*n.recruit.fem*s+(1+n.recruit.fem^2-0.001*n.recruit.fem^3)*s^2)*(gen.length))
	# ranges 
	slope.g.range <- vector(length=max.repeat)
	for (r in 1:max.repeat) slope.g.range [r] <- coefficients(lm(log(Het_pop[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE))),r][Het_pop[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE))),r]>0])~years.g[Het_pop[seq(1,max.time,round(mean(mean.gen.g,na.rm=TRUE))),r]>0]))[2]
	slope.g.low <- slope.g-1.96*sd(slope.g.range)
	slope.g.high <- slope.g+1.96*sd(slope.g.range)
	Ne_sim.range <- -1/(2*slope.g.range)
	Ne_sim.low <- Ne_sim-1.96*sqrt(sum(Ne_sim.range)/(length(Ne_sim.range)-1)^2)
	Ne_sim.high <- Ne_sim+1.96*sqrt(sum(Ne_sim.range)/(length(Ne_sim.range)-1)^2)
	gen.length.low <- gen.length-1.96*sqrt(sum(mean.gen.g)/(length(mean.gen.g)-1)^2)
	gen.length.high <- gen.length+1.96*sqrt(sum(mean.gen.g)/(length(mean.gen.g)-1)^2)
}


par(mfcol=c(3,2),mar = c(5,5,4,2) + 0.1)

# population size
plot(mean.Nind.y,
xlab="",
ylab="mean population size",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)

# original heterozygosity data
plot(mean.het.y,
xlab="",
ylab="mean Heterozygosity",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)

# ln(heterozygosity) and the best fit of the linear model
plot(log(mean.het.y),
xlab="time in years",
ylab="ln(mean Heterozygosity)",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)
lines(coefficients(lm.het.y)[1] + coefficients(lm.het.y)[2]*years, lwd=2, col="red")

# population size
plot(mean.Nind.g,
xlab="",
ylab="",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)

# original heterozygosity per generation
plot(mean.het.g,
xlab="",
ylab="",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)

# ln(heterozygosity) and the best fit of the linear model per generation
plot(log(mean.het.g),
xlab="time in generations",
ylab="",
type="b",pch=1,lwd=2,
cex.main=2,cex.lab=2,cex.axis=2)
lines(coefficients(lm.het.g)[1] + coefficients(lm.het.g)[2]*years.g, lwd=2, col="red")

## results table
res_table <- matrix(ncol=4,nrow=7)
colnames(res_table) <- c("parameter","value","95% lower CI","95% upper CI")
res_table[,1] <- c("slope of heterozygosity loss per year","slope of heterozygosity loss per generation","mean generation length", "Ny[simulation]", "Ny[calc]", "Ne[simulation]","Ne[calc]")
res_table[,2] <- c(round(slope.y,4), round(slope.g,4), round(gen.length,4), round(Ny_sim,4), round(Ny_cal,4), round(Ne_sim,4), round(Ne_cal,4))
res_table[,3] <- c(round(slope.y.low,4), round(slope.g.low,4), round(gen.length.low,4), round(Ny_sim.low,4), NA, round(Ne_sim.low,4), NA)
res_table[,4] <- c(round(slope.y.high,4), round(slope.g.high,4), round(gen.length.high,4), round(Ny_sim.high,4), NA, round(Ne_sim.high,4), NA)
res_table <- as.data.frame(res_table)

## result list
process <- list(parameters=res_table,mean.heterozygosity=mean.het.y,mean.generation.length=mean.gen.y,mean.population.size=mean.Nind.y)

if (max.repeat<2 & quiet==FALSE & num.fem==TRUE)
{
	message3 <- paste ("WARNING: You did not use any repitition. Results might be highly statistically biased.\n",
	"\n",
	"\n",
	"If you want to analyze another population, recall population().\n",
	"\n",
	"\n",
	"NEff Vers. 1.1 2015",
	"\n",
	"written by Annegret Grimm :-)\n",
	"\n",sep="")
}

if (max.repeat<2 & quiet==FALSE & num.fem==FALSE)
{
	message3 <- paste ("WARNING: You did not use any repitition. Results might be highly statistically biased.\n",
	"WARNING: The maximum number of years was cut in some of the repeats since there were no females in the population any more.\n",
	"\n",
	"\n",
	"If you want to analyze another population, recall population().\n",
	"\n",
	"\n",
	"NEff Vers. 1.1 2015",
	"\n",
	"written by Annegret Grimm :-)\n",
	"\n",sep="")
}

if (max.repeat>1 & num.fem==TRUE)
{
	message3 <- paste ("If you want to analyze another population, recall population().\n",
	"\n",
	"\n",
	"NEff Vers. 1.1 2015",
	"\n",
	"written by Annegret Grimm :-)\n",
	"\n",sep="")
}


if (max.repeat>1 & num.fem==FALSE)
{
	message3 <- paste ("WARNING: The maximum number of years was cut in some of the repeats since there were no females in the population any more.\n",
	"\n",
	"\n",	
	"If you want to analyze another population, recall population().\n",
	"\n",
	"\n",
	"NEff Vers. 1.1 2015",
	"\n",
	"written by Annegret Grimm :-)\n",
	"\n",sep="")
}

message2 <- paste ("Please switch to the Graphics window to have a look at the heterozygozity behavior over time.\n")
message1 <- paste ("Summary table of an optimal behaviour of your population:\n")

cat(paste("\n","\n"))
cat (message1)
return(process)
cat(paste("\n","\n"))
cat (message2)
cat(paste("\n","\n"))
cat (message3)
}
