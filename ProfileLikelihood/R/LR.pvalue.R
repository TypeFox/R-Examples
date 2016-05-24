LR.pvalue <-
function(y1, y2, n1, n2, interval=0.01){
s <- y1 + y2
n <-  n1 + n2
s.failure <- n - s
x1 <- n1-y1
x2 <- n2-y2
Y <- c(rep(1, s), rep(0, s.failure))
treat <- c(rep(1, y1), rep(0, y2), rep(1, x1), rep(0, x2))
fit.glm <- glm( Y ~ treat, family=binomial)
mle <- summary.glm(fit.glm)$coefficients[2,1]
se <- summary.glm(fit.glm)$coefficients[2,2]

if(y1==0 | y2==0){
	cat("Warning message: Likelihood intervals, LRs and the corresonding p-values are not reliable with empty cells in 2-by-2 tables. \n")
lo.psi <- mle - 2*abs(mle) - 2; hi.psi<- mle + 2*abs(mle) + 2
} else{
lo.psi <- mle - 2*se - 2; hi.psi<- mle + 2*se + 2
}

if(round(lo.psi,1) < 0 & round(hi.psi,1) > 0){
psi <- c(seq(from=round(lo.psi,1), to=-interval, by=interval), 0, seq(from=interval, to=round(hi.psi,1), by=interval))
} else if(round(lo.psi,1) >= 0 & round(hi.psi,1) >= 0){
psi <- c(0, seq(from=interval, to=round(hi.psi,1), by=interval))
} else{
psi <- c(seq(from=round(lo.psi,1), to=-interval, by=interval), 0)
} 

L <- function(lambda){
y1*psi[i] + (y1 + y2)*lambda - n1*log(1 + exp(psi[i] + lambda)) - n2*log(1 + exp(lambda))
}

Epi2 <- (y2+0.5)/n2
Odds2 <- Epi2/(1-Epi2)
int <- log(Odds2)

log.lik <- rep(NA, length(psi))
par <- matrix(NA, ncol = length(int), nrow = length(psi))

for(i in 1:length(psi)){
	
	if(int == -Inf){
		low.int <- - 10000
		up.int <- 1000
		ff <- optimize(L, interval=c(low.int, up.int), maximum=TRUE)
		} else if(int == Inf){
		low.int <- - 1000
		up.int <- 10000
		ff <- optimize(L, interval=c(low.int, up.int), maximum=TRUE)
		}	else{
		low.int <- - 1000
		up.int <- 1000
		ff <- optimize(L, interval=c(low.int, up.int), maximum=TRUE)
		}
	par[i,] <- ff$maximum
	log.lik[i] <- ff$objective
}
lik <-  exp(log.lik)
mm <- max(lik, na.rm=TRUE)
profile.lik.norm <- lik/mm

norm <- 0.146 
psi1.x.p.norm <- psi[profile.lik.norm >=norm]
psi.p.ci.norm <- range(psi1.x.p.norm)
if(y1==0 | y2==0){
profile.LI.norm <- NULL
} else{
	profile.LI.norm <- psi.p.ci.norm
	}

H0.profile <- unique(profile.lik.norm[psi == 0])
LR.profile <- H0.profile/1
Pvalue.profile.LR <- 1- pchisq(-2*log(LR.profile), df=1)

cond.func <- function(y1, y2, n1, n2, OR){
den.func <- function(y1, y2, n1, n2){
s <- y1 + y2
N <- n1 + n2
L <- max(0, s+n1-N)
H <- min(n1, s)
dd <- 0
for(u in L:H){
dd <- dd + choose(n1, u)*choose(n2, (s-u))*OR^u
}
return(dd)
}
dd <- den.func(y1=y1, y2=y2, n1=n1, n2=n2)
choose(n1, y1)*choose(n2, y2)*OR^y1 / dd
}

cond.lik <- rep(NA, length(psi))

for(i in 1:length(cond.lik)){
OR <- exp(psi[i])
cond.lik[i] <- cond.func(y1=y1, y2=y2, n1=n1, n2=n2, OR=OR)
}

H0 <- cond.func(y1=y1, y2=y2, n1=n1, n2=n2, OR=1)

mm.cond <- max(cond.lik, na.rm=TRUE)
cond.lik.norm <- cond.lik/mm.cond

if(y1==0 | y2==0){
mle.cond.lor <- NULL
mle <- NULL
} else{
	mle.cond.lor <- max(psi[cond.lik.norm==max(cond.lik.norm, na.rm=TRUE)], na.rm=TRUE)
	}
	
psi1.x.cond.norm <- psi[cond.lik.norm >= norm]
psi.cond.ci.norm <- range(psi1.x.cond.norm)
cond.LI.norm <- psi.cond.ci.norm

if(y1==0 | y2==0){
cond.LI.norm <- NULL
} else{
	cond.LI.norm <- psi.cond.ci.norm
	}

LR.cond <- H0/mm.cond
Pvalue.cond.LR <- 1- pchisq(-2*log(LR.cond), df=1)

tt <- matrix(c(y1, y2, (n1-y1), (n2-y2)), nr = 2, dimnames = list(c("treat", "control"), c("Success", "No success")))
Pvalue.fisher.test <- fisher.test(tt)$p.value  
Pvalue.chisq.cont.correction <- chisq.test(tt)$p.value 
Pvalue.chisq.test <- chisq.test(tt, correct=FALSE)$p.value  

return(list(mle.lor.uncond=mle, mle.lor.cond=mle.cond.lor, LI.norm.profile=profile.LI.norm, LI.norm.cond=cond.LI.norm, LR.profile=1/LR.profile, LR.cond=1/LR.cond, Pvalue.LR.profile=Pvalue.profile.LR, Pvalue.LR.cond=Pvalue.cond.LR, Pvalue.chisq.test=Pvalue.chisq.test, Pvalue.fisher.test=Pvalue.fisher.test, Pvalue.chisq.cont.correction=Pvalue.chisq.cont.correction))
}

