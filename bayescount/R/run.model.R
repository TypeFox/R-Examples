run.model <- function(data=stop("No data supplied"), model=stop("No model specified"), call.jags = FALSE, alt.prior = FALSE, monitor.lambda=FALSE, monitor.deviance=FALSE, ...){
	
###  Currently only the IP model requires gamma monitored for the likelihood bit - others are integrated

if(monitor.lambda==TRUE && model!="IP"){
	monitor.lambda <- FALSE
}

if(call.jags){
	warning('The call.jags argument is deprecated and will be removed from version 1 of the bayescount package')
}

updates <- sample

if(class(data)=="list"){
	
	if(suppressWarnings(class(data$totals) != "NULL" & class(data$repeats) != "NULL")){
		# data is a list of counts and repeats
		totals <- data$totals
		repeats <- data$repeats
		
		if(length(totals)!=length(repeats)) stop("The number of repeats does not match the length of totals")
		
		if(any(is.na(totals)) | any(is.na(repeats))) stop("Missing data is not allowed in repeats or totals")
		
		data <- matrix(NA, ncol=max(repeats), nrow=length(totals))
		
		for(i in 1:length(totals)){
				
			if(round(totals[i])!=totals[i]) stop("All totals must be an integer")
			if(round(repeats[i])!=repeats[i] | repeats[i] < 1) stop("All repeats must be an integer greater than 0")
			
			if(repeats[i] > 1){

				done <- FALSE
	
				while(!done){
					data[i,1:repeats[i]] <- rpois(repeats[i], totals[i]/repeats[i])
					if(sum(data[i,],na.rm=TRUE)==totals[i] & all(na.omit(data[i,]) >= 0)) done <- TRUE
		
				}
			}else{
				data[i,1] <- totals[i]
			}
		}
		
	}else{

		delist <- function(list){
		
		N <- length(list)
		R <- numeric(length=N)
	
		for(i in 1:N){
			R[i] <- length(na.omit(as.numeric(list[[i]])))
		
		}
	
		if(sum(R>0)==0) stop("No real values in the data")
	
		data <- matrix(NA, ncol=max(R), nrow=N)
	
		for(i in 1:N){
			data[i,(min(R[i], 1):R[i])] <- na.omit(as.numeric(list[[i]]))
		}
	
		data <- data[apply(!is.na(data), 1, sum)>0,]

		if(sum(R>0)==1) return(matrix(data, nrow=1)) else return(data)
		}
	
		data <- delist(data)	
	}
}

if(class(data)!="matrix"){
	counts <- matrix(as.numeric(data[!is.na(data)]), ncol=1)
	N <- length(counts)
	R <- numeric(length=N)
	R[] <- 1
}else{
	counts <- data
	counts[] <- as.numeric(data)
	sums <- function(data) return(sum(!is.na(data)))
	R <- apply(counts, 1, sums)
	use <- R!=0
	R <- R[use]
	counts <- matrix(counts[use,], nrow=length(R))
	N <- length(R)
}

model <- toupper(model)

models <- c("SP", "ZISP", "GP", "ZIGP", "WP", "ZIWP", "LP", "ZILP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "independant Poisson")

#models <- c("SP", "ZISP", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP", "IP", paste("G", 1:10, sep=""), paste("L", 1:10, sep=""))
#modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "independant Poisson")

if(sum(model==models) != 1){
	cat("Invalid model selection.  Please choose from ONE of the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}






chains <- 2

inits <- matrix("", ncol=5, nrow=2)
	
if(alt.prior!=FALSE & (any(model==c("IP", "SP", "ZISP")))){
	cat("Warning:  Alternate or custom prior distribution for dispersion for the ", model, " model is not allowed.  The standard prior distribution will be used\n", sep="")
	alt.prior <- FALSE
}

dataformean = dataformean2 <- apply(counts, 1, mean, na.rm=TRUE)

dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))
meanofall <- pmax(0.2, round(mean(dataformean2, na.rm=TRUE), 2), na.rm=TRUE)

initstring <- character(length=chains)

smallestmean <- as.integer(meanofall / 2)
largestmean <- as.integer(meanofnonzeros * 2)

if(any(!is.na(dataformean2))){
	smallestdispersion <- ((var(dataformean2, na.rm=TRUE) - mean(dataformean2, na.rm=TRUE))^0.5 / mean(dataformean2, na.rm=TRUE)) / 2
}else{
	smallestdispersion <- 0.01
}

if(any(!is.na(dataformean))){
	largestdispersion <- ((var(dataformean, na.rm=TRUE) - mean(dataformean, na.rm=TRUE))^0.5 / mean(dataformean, na.rm=TRUE)) * 2 
}else{
	largestdispersion <- 5
}

if(is.na(smallestmean)==TRUE){
	smallestmean <- 1
}
if(is.na(largestmean)==TRUE){
	largestmean <- 10
}

if((smallestmean < 1)==TRUE){
	smallestmean <- 1
}
if((largestmean < 10)==TRUE){
	largestmean <- 10
}
if((smallestmean > 20)==TRUE){
	smallestmean <- 20
}
if((largestmean > 200)==TRUE){
	largestmean <- 200
}

if(is.na(smallestdispersion)) smallestdispersion <- 0.01
if(is.na(largestdispersion)) largestdispersion <- 5

smallestdispersion <- max(smallestdispersion, 0.1)
largestdispersion <- min(largestdispersion, 9)
if(largestdispersion < 0.1) largestdispersion <- 0.1
if(smallestdispersion > 5) smallestdispersion <- 5


# FOR TESTING DIFFERENT MODELS:
if(any(model==c(paste("G", 1:10, sep=""), paste("L", 1:10, sep="")))){

	true.model <- model
	model <- ""

	probpos <- dataformean2
	zeros <- probpos==0
	probpos[] <- 1
	gammas <- dataformean2
	gammas[gammas < 1] <- 1
	
	initstring[1] <- paste(dump.format(list(prob=0.95, probpos=probpos, gamma=gammas)))
	probpos[zeros] <- 0
	initstring[2] <- paste(initstring[2], dump.format(list(prob=min(round(sum(data!=0)/length(data),2), 0.05, na.rm=TRUE), probpos=probpos, gamma=gammas)))
		
}else{
	true.model <- ""
}



if(model=="SP" | model=="ZISP"){
	initstring[1] <- dump.format(list(mean=smallestmean))
	initstring[2] <- dump.format(list(mean=largestmean))
	
}

if(model=="LP" | model=="ZILP"){

	lgammas1 = lgammas2 <- pmax(dataformean, 0.1, na.rm=TRUE)
	
	lninits1 <- lnormal.params(mean=smallestmean, NA, largestdispersion)
	lninits2 <- lnormal.params(mean=largestmean, NA, smallestdispersion)
	#not actually var but i'm not changing code:
	varinitl <- pmax(c(round(lninits1[[2]],3)),0.1, na.rm=TRUE)
	varinitu <- pmax(c(round(lninits2[[2]],3)),0.1,na.rm=TRUE)
	
	initstr1 <-list(lmu=round(lninits1[[1]],3), gamma=lgammas1)
	initstr2 <- list(lmu=round(lninits2[[1]],3), gamma=lgammas2)
	
	if(class(alt.prior)=="character"){
		priorstring <- paste("lsd ~ ", alt.prior, ";\n", sep="")
		inistr1 <- c(initstr1, lsd=varinitl)
		inistr2 <- c(initstr2, lsd=varinitu)
		
	}else{
		if(alt.prior==TRUE){  #CHANGED FROM FALSE TO TRUE
			priorstring <- "llsd ~ dunif(-4.61, 0.765);\nlsd <- exp(llsd);"
			varinitl <- log(varinitl)
			varinitu <- log(varinitu)
			initstr1 <- c(initstr1, llsd=varinitl)
			initstr2 <- c(initstr2, llsd=varinitu)
		}else{
			priorstring <- "lsd ~ dunif(0.01,2.149);\n"
			initstr1 <- c(initstr1, lsd=varinitl)
			initstr2 <- c(initstr2, lsd=varinitu)
		}	
	}
	
	#if(exists("zilp.type")){
	#	if(zilp.type==3 | zilp.type==4){
	#		initstr1$lmu <- smallestmean
	#		initstr2$lmu <- largestmean
	#		names(initstr1)[1] = "mean"
	#		names(initstr2)[1] = "mean"
	#	}
	#}
	initstring[1] <- dump.format(initstr1)
	initstring[2] <- dump.format(initstr2)
	
	
}

if(model=="IP"){

	gammas1 <- dataformean2
	gammas2 <- dataformean2
	gammas1[gammas1 == 0] <- 0.1
	gammas2[gammas2 == 0] <- 0.1

	initstring[1] <- dump.format(list(gamma= gammas1))
	initstring[2] <- dump.format(list(gamma=gammas2))
}

if(model=="GP" | model=="ZIGP"){

	gammas1 = gammas2 <- numeric(length=length(dataformean2))
	gammas1[] <- 1
	gammas2[] <- 1
	
	initstring[1] <- dump.format(list(mean=smallestmean, gamma=gammas1))
	initstring[2] <- dump.format(list(mean=largestmean, gamma=gammas2))			
	
	if(class(alt.prior)=="character"){
		priorstring <- paste("ia ~ ", alt.prior, ";\n", sep="")
		initstring[1] <- paste(initstring[1], dump.format(list(ia=largestdispersion^2)), sep="")
		initstring[2] <- paste(initstring[2], dump.format(list(ia=smallestdispersion^2)), sep="")
	}else{
		if(alt.prior==FALSE){
			priorstring <- "ia <- exp(logia);\nlogia ~ dunif(-9.21,4.6);\n" #-9.21, 4.6
			#initstring[1] <- paste(initstring[1], dump.format("loga", -4.6), sep="")
			#initstring[2] <- paste(initstring[2], dump.format("loga", 4.6), sep="")
			initstring[1] <- paste(initstring[1], dump.format(list(logia=log(largestdispersion^2))), sep="")
			initstring[2] <- paste(initstring[2], dump.format(list(logia=log(smallestdispersion^2))), sep="")
		}else{
			priorstring <- "ia ~ dunif(0.0001,100);\n"
			initstring[1] <- paste(initstring[1], dump.format(list(ia=largestdispersion^2)), sep="")
			initstring[2] <- paste(initstring[2], dump.format(list(ia=smallestdispersion^2)), sep="")
		}
	}
}


#print(initstring[1])
#print(initstring[2])

if(model=="WP" | model=="ZIWP"){

	if(class(alt.prior)=="character"){
		priorstring <- paste("a ~ ", alt.prior, ";\nb ~ ", alt.prior, ";\n", sep="")
	}else{
		if(alt.prior==FALSE){
			priorstring <- "a ~ dgamma(1,0.01);\nb ~ dgamma(1,0.01);\n"
		}else{
			priorstring <- "a ~ dunif(0.001,1000);\nb ~ dunif(0.001,1000);\n"
		}
	}
	
	gammas1 = gammas2 <- numeric(length=length(dataformean2))
	gammas1[] <- 1
	gammas2[] <- 100
	
	initstring[1] <- dump.format(list(a=0.01, gamma=gammas1, b=0.1))
	initstring[2] <- dump.format(list(a=100, gamma=gammas2, b=100))
	
}

if(model=="ZISP" | model=="ZILP" | model=="ZIGP" | model=="ZIWP"){

	probpos <- dataformean2
	zeros <- probpos==0
	probpos[] <- 1
	
	initstring[1] <- paste(initstring[1], dump.format(list(prob=0.95, probpos=probpos)))
	probpos[zeros] <- 0
	initstring[2] <- paste(initstring[2], dump.format(list(prob=min(round(sum(data!=0)/length(data),2), 0.05, na.rm=TRUE), probpos=probpos)))
	
}


##### make model string
# L1-7 and G1-6 still here but can't be accessed 

if(model=="SP"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(mean);
}
}

# Priors
mean ~ dunif(0.001,1000);
}")
monitors <- "mean"
}
if(model=="ZISP"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * mean;
probpos[row] ~ dbern(prob);
}

# Priors
mean ~ dunif(0.001,1000);
prob ~ dbeta(1,1);
}")
monitors <- c("mean", "prob")
}

if(model=="IP"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(gamma[row]);
}
# Priors
gamma[row] ~ dunif(0.001,1000);
}

mean <- mean(gamma[])
sd <- sd(gamma[])

}")
monitors <- c("mean", "sd")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}


if(model=="LP"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(gamma[row]);
}
gamma[row] ~ dlnorm(lmu, lprec);

}

lprec <- 1 / lsd^2;

meanl <- log(0.001) - ((lsd^2)/2);
meanu <- log(1000) - ((lsd^2)/2);

# Priors
lmu ~ dunif(meanl,meanu);
", priorstring, "}", sep="")
monitors <- c("lmu", "lprec")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZILP"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lsd^2;

meanl <- log(0.001) - ((lsd^2)/2);
meanu <- log(1000) - ((lsd^2)/2);

# Priors
lmu ~ dunif(meanl,meanu);
prob ~ dbeta(1,1);
", priorstring, "}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


#if(exists("zilp.type") & FALSE){
#	if(zilp.type==3 | zilp.type==4){
#if(model=="ZILP"){
#modelstring <- paste("model {
#for(row in 1 : N){
#for(repeat in 1:R[N]){
#Count[row,repeat] ~ dpois(lambda[row]);
#}
#lambda[row] <- probpos[row] * mean * gamma[row];
#gamma[row] ~ dlnorm(0, lprec);
#probpos[row] ~ dbern(prob);
#}

#lprec <- 1 / lsd^2;

# Priors
#mean ~ dunif(0.001,1000)
#prob ~ dbeta(1,1);
#", priorstring, "}", sep="")

#monitors <- c("mean", "lprec", "prob")
#if(monitor.lambda==TRUE){
#	monitors <- c(monitors, "gamma")#(, "probpos")
#}
#}}
#if(zilp.type==2 | zilp.type==4){
#	monitors[2] <- "lsd"
#}
#}


if(true.model=="L1"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lvar;

meanl <- log(0.001) - ((lvar)/2);
meanu <- log(1000) - ((lvar)/2);

# Priors
lmu ~ dunif(meanl,meanu);
prob ~ dbeta(1,1);
lvar ~ dunif(0.0001,4.618);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(true.model=="L2"){  #WINNER
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lsd^2;
lvar <- lsd^2;
meanl <- log(0.001) - ((lsd^2)/2);
meanu <- log(1000) - ((lsd^2)/2);

# Priors
lmu ~ dunif(meanl,meanu);
prob ~ dbeta(1,1);
lsd ~ dunif(0.01,2.149);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}



if(true.model=="L3"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lvar <- log((sd/mean)^2 + 1); 
lprec <- 1 / lvar;

lmu <- log(mean) - ((lvar) / 2)

# Priors
mean ~ dunif(0.001,1000);
prob ~ dbeta(1,1);
sd ~ dunif(0.000001,10000);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(true.model=="L4"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lvar;

lmu <- log(mean) - ((lvar) / 2)

# Priors
mean ~ dunif(0.001,1000);
prob ~ dbeta(1,1);
lvar ~ dunif(0.0001,4.618);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(true.model=="L5"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lsd^2;
lvar <- lsd^2

lmu <- log(mean) - ((lsd^2) / 2)

# Priors
mean ~ dunif(0.001,1000);

prob ~ dbeta(1,1);
lsd ~ dunif(0.01,2.149);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(true.model=="L6"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lvar;

meanl <- log(0.001) - ((lvar)/2);
meanu <- log(1000) - ((lvar)/2);

# Priors
lmu ~ dunif(meanl,meanu);
prob ~ dbeta(1,1);
llvar ~ dunif(-9.21, 1.61);
lvar <- exp(llvar);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(true.model=="L7"){
modelstring <- paste("model {
for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1 / lsd^2;
lvar <- lsd^2;
meanl <- log(0.001) - ((lsd^2)/2);
meanu <- log(1000) - ((lsd^2)/2);

# Priors
lmu ~ dunif(meanl,meanu);
prob ~ dbeta(1,1);
llsd ~ dunif(-4.61,0.81);
lsd <- exp(llsd);
}", sep="")

monitors <- c("lmu", "lprec", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}



if(model=="GP"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- mean * gamma[row];
gamma[row] ~ dgamma(a, a)T(10^-200,);

}
a <- 1 / ia;

# Priors
mean ~ dunif(0.001,1000);
", priorstring, "}", sep="")
monitors <- c("mean", "a")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZIGP"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, a)T(10^-200,);

}

a <- 1 / ia;

# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
", priorstring, "}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(true.model=="G1"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, b);

}


# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
a ~ dunif(0,100)
b ~ dunif(0,100)
}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(true.model=="G2"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, a);

}

# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
a <- 1/ia;
ia ~ dunif(0.0001,100);
}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(true.model=="G3"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, b);

}

b <- a / mean;

# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
a <- 1/ia;
ia ~ dunif(0.0001,100);
}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}



if(true.model=="G4"){  #WINNER
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, a);

}

a <- 1/ia;

# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
ia <- exp(logia);
logia ~ dunif(-9.21,4.6);
}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(true.model=="G5"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(a, b);

}

b <- a / mean;
a <- 1/ia;

# Priors
prob ~ dbeta(1,1);
mean ~ dunif(0.001,1000);
ia <- exp(logia);
logia ~ dunif(-9.21,4.6);

}", sep="")
monitors <- c("mean", "a", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}



if(model=="WP"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(gamma[row]);
}
gamma[row] ~ dweib(a,nb);
}

nb <- exp(-(log(b) * a));

# Priors
", priorstring, "}", sep="")
monitors <- c("a", "b")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZIWP"){
modelstring <- paste("model {

for(row in 1 : N){
for(repeat in 1:R[N]){
Count[row,repeat] ~ dpois(lambda[row]);
}
lambda[row] <- probpos[row] * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dweib(a, nb);

}

nb <- exp(-(log(b) * a));

# Priors
prob ~ dbeta(1,1);
", priorstring, "}", sep="")
monitors <- c("a", "b", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(monitor.deviance) monitors <- c(monitors, "deviance")

datastring <- dump.format(list(N=N, R=R, Count=counts))

if(!call.jags){
	return(run.jags(model=modelstring, monitor=monitors, data=datastring, n.chains=2, inits=initstring, burnin=0, adapt=0, sample=0, ...))
}else{
	return(run.jags(model=modelstring, monitor=monitors, data=datastring, n.chains=2, inits=initstring, ...))
}

}

count.model <- run.model
fec.model <- run.model
FEC.model <- run.model
