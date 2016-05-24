simfam <-
function(N.fam=500, design="pop", variation="none", depend=1,
			base.dist="Weibull", frailty.dist="gamma",
           base.parms=c(0.016,3), vbeta=c(-1.13, 2.35, 0.5), allelefreq=c(0.02, 0.2), 
           dominant.m=TRUE, dominant.s=TRUE, mrate=0, hr=0, age1=c(65,2.5), age2=c(45,2.5), 
           agemin=20)
  {
	


if(!is.element(frailty.dist, c("gamma", "lognormal"))) stop("Unrecognized frailty distribution") 

if(!is.element(variation, c("none", "frailty", "secondgene"))) stop("variation should be one of none, frailty, or secondgene")
 
if(!is.element(design, c("pop","pop+","cli","cli+","twostage"))) stop("Unrecognized design; should be one of pop, pop+, cli, cli+ or twostage")  

if(design=="pop") {affectnum=1; m.carrier=0}
if(design=="pop+"){affectnum=1; m.carrier=1}
if(design=="cli") {affectnum=3; m.carrier=0}
if(design=="cli+"){affectnum=3; m.carrier=1}
if(design=="twostage"){
	affectnum = 2
	m.carrier = 0 #proband is not necessary to be a carrier.
	if(hr==0) stop("Please specify the sampling rate of high risk families (0<hr<1)")
	}

dat <- data.frame(familyDesign(n=N.fam, affectnum=affectnum, m.carrier=m.carrier, 
        base.dist=base.dist, frailty.dist=frailty.dist,
        dominant.m=dominant.m, dominant.s=dominant.s, depend=depend, vbeta=vbeta, 
        parms=base.parms, variation=variation, allelefreq=allelefreq, 
        mrate=mrate, age1=age1, age2=age2, agemin=agemin))

if(hr==0){ # One stage sampling 
dat$weight <- 1
        }
else { # Two stage sampling 
	if(design!="twostage") stop ("hr can be used only with twostage design" )
  en.hr <- round(hr*N.fam)  # expected number of high risk families 
  en.lr <- N.fam - en.hr	# expected number of low risk families
  datp <- dat[dat$proband==1,]
  n.hr <- sum(datp$naff>1) # number of high risk families in the simulated families
	fam.id <- datp$famID 
	lfam.id <- fam.id[datp$naff==1] # family IDs of low risk families
	hfam.id <- fam.id[datp$naff>1] # family IDs of high risk families 

	if(length(lfam.id) < en.lr)	id <- lfam.id 
	else id <- sample(lfam.id, en.lr, replace=FALSE) # sampled low risk families

	n.more <- ifelse(en.hr < n.hr, 0, en.hr-n.hr)

	# dat2: additinal high risk families to generate

	dat1 <- dat[is.element(dat$famID, c(id, hfam.id)), ]
	dat2 <- data.frame(familyDesign(n=n.more, affectnum=affectnum, m.carrier=m.carrier,
			base.dist=base.dist, frailty.dist=frailty.dist,
			dominant.m=dominant.m, dominant.s=dominant.s, depend=depend, vbeta=vbeta,
          	parms=base.parms, variation=variation, allelefreq=allelefreq, 
          	mrate=mrate,age1=age1, age2=age2, agemin=agemin))

	dat <- rbind(dat1, dat2)
	samrate <- en.lr/(N.fam-n.hr)/(en.hr/n.hr) 
	# samrate=sampling rate for low risk families when we fix sampling rate for high risk as 1
	# Eg: Suppose we want to include 30% high risk families i.e. HR: LR = 30:70
	# From random sampling, HR:LR ratio is 15:85
	# As we fix sampling rate for HR families as 1, the sampling rate for LR families should be (70/30) / (85/15) 

	dat$weight <- ifelse(dat$naff>1, 1, samrate) 
	# weight=1 for HR families and Weight=1/samrate for LR families

} 

class(dat)<-c("simfam","data.frame")
attr(dat, "design") <- design
attr(dat, "variation") <- variation
attr(dat, "frailty.dist") <- frailty.dist
attr(dat, "base.dist") <- base.dist
attr(dat, "base.parms") <- base.parms
attr(dat, "vbeta") <- vbeta
attr(dat, "agemin") <- agemin
return(dat)
}
