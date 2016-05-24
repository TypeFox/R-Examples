lira <- function(
x, y,
delta.x, delta.y, covariance.xy,
y.threshold, delta.y.threshold, x.threshold, delta.x.threshold,
z, z.ref=0.01, time.factor="Ez", distance=c("luminosity","angular"),
Omega.M0=0.3, Omega.L0=0.7,
alpha.YIZ="dunif", beta.YIZ="dt", gamma.YIZ="dt", delta.YIZ=0.0,
sigma.YIZ.0="prec.dgamma", gamma.sigma.YIZ.Fz=0.0, gamma.sigma.YIZ.D=0.0,
alpha.XIZ=0.0, beta.XIZ=1.0, gamma.XIZ=0.0, delta.XIZ=0.0,
sigma.XIZ.0=0.0, gamma.sigma.XIZ.Fz=0.0, gamma.sigma.XIZ.D=0.0,
rho.XYIZ.0=0.0,  gamma.rho.XYIZ.Fz=0.0, gamma.rho.XYIZ.D=0.0,
n.mixture=1, pi="ddirch",
mu.Z.0="dunif", gamma.mu.Z.Fz="dt", gamma.mu.Z.D="dt",
sigma.Z.0="prec.dgamma", gamma.sigma.Z.Fz=0.0, gamma.sigma.Z.D=0.0,
mu.Z.0.mixture="dunif", sigma.Z.0.mixture="prec.dgamma",
Z.knee="dunif", l.knee=1e-04,
beta.YIZ.knee="beta.YIZ", delta.YIZ.knee="delta.YIZ", sigma.YIZ.0.knee = "sigma.YIZ.0",
mu.Z.min.0="-n.large", gamma.mu.Z.min.Fz="dt", gamma.mu.Z.min.D="dt",
sigma.Z.min.0=0.0, gamma.sigma.Z.min.Fz=0.0, gamma.sigma.Z.min.D=0.0,
Z.max= "n.large",
n.large=1e04,
inits, n.chains=1, n.adapt=1e03, quiet=TRUE,
n.iter=1e04, thin=1,
print.summary=TRUE, print.diagnostic=TRUE,
export=FALSE, export.mcmc="", export.jags="",
X.monitored=FALSE, export.X="",
XZ.monitored=FALSE, export.XZ="",
Z.monitored=FALSE, export.Z="",
Y.monitored=FALSE, export.Y="",
YZ.monitored=FALSE, export.YZ=""){
	
	# ====================== consistency checks ========================
	#===================================================================
	n.data <- length(x)
	delta.min <- 0.001
	arg.match.call <- match.call(expand.dots = FALSE)
	names.arg <- names(arg.match.call)
	
	#================== x and y =======================================
	if(n.data==0)stop("no x data: computation suppressed")
	if(n.data!=length(y))stop("x and y of inequal lengths: computation suppressed")
	
	#===== z ==========================================================
	if(is.na(match("z",names.arg))){z.logical <- FALSE}
	else{
		z.logical <- TRUE
		if(n.data!=length(z)){stop("x and z of inequal lengths: computation suppressed")}
		if(length(which(z< 0))>0){stop("negative values of z: computation suppressed")}
		if(length(unique(z))==1){warning("all objects at the same z: no time evolution");z.logical <- FALSE}
	}
	
	# =================== delta.x =========================
	if(is.na(match("delta.x",names.arg))  || length(which(delta.x==0))==n.data ){
		delta.x.logical <- FALSE
	} 
	else if (length(delta.x)!=n.data && length(delta.x)>0){stop("x and delta.x of inequal lengths: computation suppressed")}
	else if (n.data==length(delta.x) && length(which(delta.x==0))>0) {
		ifelse(delta.x==0, min(delta.min,delta.x[which(delta.x!=0)]), delta.x)
		warning("null uncertainties on x set to a minimum error") 
		delta.x.logical <- TRUE
		}
	else{delta.x.logical <- TRUE}
	
	if(!delta.x.logical && X.monitored && !anyNA(x)){
		warning("x without uncertainties: x and X identified")
		X.monitored=FALSE
		}
	
	# =================== delta.y =========================
	if(is.na(match("delta.y",names.arg))  || length(which(delta.y==0))==n.data){
		delta.y.logical <- FALSE
	}
	else if(n.data!=length(delta.y))(stop("y and delta.y of inequal lengths: computation suppressed"))
	else if(n.data==length(delta.y) && length(which(delta.y==0))>0){
		ifelse(delta.y==0, min(delta.min,delta.y[which(delta.y!=0)]), delta.y)
		warning("null uncertainties on y set to a minimum error") 
		delta.y.logical <- TRUE
		}
	else{delta.y.logical <- TRUE}
	
	if(!delta.y.logical && Y.monitored && !anyNA(y)){
		warning("y without uncertainties: y and Y identified")
		Y.monitored=FALSE
		}
	
	# ===================== rho.xy =============================
#	rho.max<-0.999
	if(delta.x.logical && delta.y.logical){
		if(is.na(match("covariance.xy",names.arg))){rho.xy <- rep(0.0,n.data)}
		else{
			if (length(covariance.xy)!=n.data) {stop("x and covariance.xy of inequal lengths: computation suppressed")}
			rho.xy<-covariance.xy/delta.x/delta.y
			if (max(rho.xy)>=1 || min(rho.xy)<=-1) {stop("correlation factors out of range -1<=rho.xy<=1: computation suppressed")}
		}
	}
	else{
		if(!is.na(match("covariance.xy",names.arg)) && !delta.x.logical){warning("delta.x are not provided; covariance.xy set to 0 ")}
		if(!is.na(match("covariance.xy",names.arg)) && !delta.y.logical){warning("delta.y are not provided; covariance.xy set to 0 ")}
	}
	
	#===== y.threshold ==========
	if(!is.na(match("y.threshold",names.arg)) ) {
		if(length(y.threshold)<n.data){stop("y and y.threshold of inequal lengths: computation suppressed")}
		y.threshold.logical <- TRUE
	}
	else {y.threshold.logical <- FALSE}
	
	#===== x.threshold ==========
	if(!is.na(match("x.threshold",names.arg)) ) {
		if(length(x.threshold)<n.data){stop("x and x.threshold of inequal lengths: computation suppressed")}
		x.threshold.logical <- TRUE
	}
	else {x.threshold.logical <- FALSE}
	
	#===== delta.y.threshold ==========
	if(is.na(match("delta.y.threshold",names.arg))){
		delta.y.threshold.logical <- FALSE
	}
	else{
		delta.y.threshold.logical <- TRUE
		if(length(delta.y.threshold)!=n.data){
			stop("y and delta.y.threshold of inequal lengths: computation suppressed")
			}
		else{ifelse(delta.y.threshold==0, min(delta.min,delta.y.threshold[which(delta.y.threshold!=0)]), delta.y.threshold)}
		if(!y.threshold.logical){delta.y.threshold.logical <- FALSE; warning("delta.y.threshold cannot be used since y.threshold's are missing")}
		if(length(which(delta.y.threshold==0))){delta.y.threshold.logical <- FALSE}
	}
	
	#===== delta.x.threshold ==========
	if(is.na(match("delta.x.threshold",names.arg))){
		delta.x.threshold.logical <- FALSE
	}
	else{
		delta.x.threshold.logical <- TRUE
		if(length(delta.x.threshold)!=n.data){
			stop("x and delta.x.threshold of inequal lengths: computation suppressed")
			}
		else{ifelse(delta.x.threshold==0, min(delta.min,delta.x.threshold[which(delta.x.threshold!=0)]), delta.x.threshold)}
		if(!x.threshold.logical){delta.x.threshold.logical <- FALSE; warning("delta.x.threshold cannot be used since x.threshold's are missing")}
		if(length(which(delta.x.threshold==0))){delta.x.threshold.logical <- FALSE}
	}



	#=========JAGS master script ======================================================
	ScriptJAGS.input <- "#= JAGS script
#===========================================================================================
#=
var
# x[n.data], 
y[n.data]
# , prec.delta.x[n.data], prec.delta.y[n.data], rho.xy[n.data]
# , XY[n.data,2], covariance.mat.XY[n.data,2,2]

# data	{
#  for	(i in 1:n.data)	{...}
# }

model {
#============== P(x,y|X,Y) ========================
for	(i in 1:n.data)	{

x[i]                 ~  dnorm(X[i], prec.delta.x[i])                                      # default_if
mean.yIx[i]          <- Y[i]+rho.xy[i]*(x[i]-X[i])*sqrt(prec.delta.x[i]/prec.delta.y[i])  # default_if
prec.delta.yIx[i]    <- min( prec.delta.y[i]/(1.0-rho.xy[i]^2), n.large)                  # default_if
y[i]                 ~  dnorm(Y[i], prec.delta.y[i])                                      # default_if
# (if_delta.x!=0) #                               x[i]              ~  dnorm(X[i], prec.delta.x[i])
# (if_delta.x!=0) #   # (if_delta.y!=0) #         mean.yIx[i]       <- Y[i]+rho.xy[i]*(x[i]-X[i])*sqrt(prec.delta.x[i]/prec.delta.y[i])
# (if_delta.x!=0) #   # (if_delta.y!=0) #         prec.delta.yIx[i] <- prec.delta.y[i]/(1.0-rho.xy[i]^2)
# # # # (if_delta.x!=0) #   # (if_delta.y!=0) #   prec.delta.yIx[i] <- min( prec.delta.y[i]/(1.0-rho.xy[i]^2), n.large)
# (if_delta.x!=0) #   # (if_delta.y!=0) #         y[i]              ~  dnorm(mean.yIx[i],prec.delta.yIx[i])
# (if_delta.x==0) #   # (if_delta.y!=0) #         y[i]              ~  dnorm(Y[i],prec.delta.y[i])

#=
#=============== Malmquist bias======================

# MB as a step function -- y[i]~dnorm(mean.yIx[i],prec.delta.yIx[i])T(y.min[i],) -- , or MB as 0.5(1-erf(y)) if threshold has uncertainty
# (if_delta.y!=0) #   # (if_delta.y.threshold=={}) #	y.min[i]  <-  y.threshold[i]	
# (if_delta.y.threshold!={}) #	                        y.min[i]  ~   dnorm(y.threshold[i],delta.y.threshold[i]^(-2))	

# (if_delta.y!=0) #   # (if_y.threshold!={}) #          Y.min[i]  ~  dnorm(y.min[i],prec.delta.y[i])
# (if_delta.y==0) #   # (if_y.threshold!={}) #          Y.min[i]  <- y.threshold[i]

# (if_delta.x!=0) #   # (if_delta.x.threshold=={}) #	x.min[i]  <-  x.threshold[i]	
# (if_delta.x.threshold!={}) #	                        x.min[i]  ~   dnorm(x.threshold[i],delta.x.threshold[i]^(-2))	

# (if_delta.x!=0) #  # (if_sigma.XIZ.0!=0) # # (if_x.threshold!={}) #          X.min[i]  ~  dnorm(x.min[i],prec.delta.x[i])
# (if_delta.x==0) #  # (if_sigma.XIZ.0!=0) # # (if_x.threshold!={}) #          X.min[i]  <- x.threshold[i]


#=
#=============== X-Z scaling  ===================

fXZ[i]              <- alpha.XIZ+beta.XIZ*Z[i]+gamma.XIZ*T[i]+delta.XIZ*T[i]*Z[i]                        # default_if
X[i]	            ~  dnorm(fXZ[i],prec.sigma.XIZ.z[i])                                                 # default_if
prec.sigma.XIZ.z[i]	<- prec.sigma.XIZ.0*(Fz[i]^-(2.*gamma.sigma.XIZ.Fz))*(D[i]^(-2.*gamma.sigma.XIZ.D))  # default_if

# (if_sigma.XIZ.0!=0) #                         fXZ[i]              <-  alpha.XIZ+beta.XIZ*Z[i] *flag.z* +gamma.XIZ*T[i] *flag.delta.XIZ* +delta.XIZ*T[i]*Z[i]  
# (if_sigma.XIZ.0!=0) #                         X[i]	            ~	dnorm(fXZ[i],prec.sigma.XIZ.z[i])
# (if_sigma.XIZ.0!=0) #                         prec.sigma.XIZ.z[i]	<-  prec.sigma.XIZ.0 *flag.z* *(Fz[i]^-(2.*gamma.sigma.XIZ.Fz))*(D[i]^(-2.*gamma.sigma.XIZ.D))
# (if_delta.x!=0) #   # (if_sigma.XIZ.0==0) #	X[i]                <-  alpha.XIZ+beta.XIZ*Z[i] *flag.z* +gamma.XIZ*T[i] *flag.delta.XIZ* +delta.XIZ*T[i]*Z[i] 

#=
#=============== Y-Z scaling  ===================
fYZ[i]   <-  alpha.YIZ+beta.YIZ*Z[i]+gamma.YIZ*T[i]  +delta.YIZ*Z[i]*T[i]    # default_if 

# (if_beta.YIZ.knee!=beta.YIZ) # fStep[i]  <-  1.0/(1.0+exp((Z[i]-Z.knee)/l.knee)) 

# (if_beta.YIZ.knee==beta.YIZ) # fYZ[i]    <-  alpha.YIZ+beta.YIZ*Z[i] *flag.z* +gamma.YIZ*T[i] *flag.delta.YIZ* +delta.YIZ*Z[i]*T[i] 
# (if_beta.YIZ.knee!=beta.YIZ) # fYZ[i]    <-  alpha.YIZ+beta.YIZ*Z[i] +(beta.YIZ-beta.YIZ.knee)*(Z.knee-Z[i])*fStep[i] *flag.z* +gamma.YIZ*T[i] *flag.delta.YIZ* +delta.YIZ*Z[i]*T[i]+(delta.YIZ-delta.YIZ.knee)*(Z.knee-Z[i])*T[i]*fStep[i]

# (if_sigma.YIZ.0.knee==sigma.YIZ.0) #   fprec.sigma.YIZ.0[i] <-  prec.sigma.YIZ.0
# (if_sigma.YIZ.0.knee!=sigma.YIZ.0) #   fprec.sigma.YIZ.0[i] <-  prec.sigma.YIZ.0+(prec.sigma.YIZ.0.knee-prec.sigma.YIZ.0)*fStep[i]

Y[i]  ~ dnorm(fYZ[i], prec.sigma.YIZ.z[i])         # default_if
# (if_rho.XYIZ==0) #  Y[i]      ~   dnorm(fYZ[i], prec.sigma.YIZ.z[i])

# (if_beta.YIZ.knee==beta.YIZ) # prec.sigma.YIZ.z[i]	<-	prec.sigma.YIZ.0 *flag.z* *(Fz[i]^(-2.*gamma.sigma.YIZ.Fz))*(D[i]^(-2.*gamma.sigma.YIZ.D))
# (if_beta.YIZ.knee!=beta.YIZ) # prec.sigma.YIZ.z[i]	<-	fprec.sigma.YIZ.0[i] *flag.z* *(Fz[i]^(-2.*gamma.sigma.YIZ.Fz))*(D[i]^(-2.*gamma.sigma.YIZ.D))

#=
#=============== rho.XYIZ =================================

# (if_rho.XYIZ!=0) #  mean.YIX[i]          <- fYZ[i]+rho.XYIZ.z[i]*(X[i]-fXZ[i])*sqrt(prec.sigma.XIZ.z[i]/prec.sigma.YIZ.z[i])
# (if_rho.XYIZ!=0) #  prec.sigma.YIX[i]    <- prec.sigma.YIZ.z[i]/(1.0-rho.XYIZ.z[i]^2.0)
# (if_rho.XYIZ!=0) #  rho.XYIZ.z[i]        <- rho.XYIZ.0 *flag.z* *(Fz[i]^gamma.rho.XYIZ.Fz)*(D[i]^gamma.rho.XYIZ.D)
# (if_rho.XYIZ!=0) #  Y[i]                 ~  dnorm(mean.YIX[i],prec.sigma.YIX[i])

#=
#=================== p(Z) =================================
Z[i]              ~  dnorm(mu.Z.z[i],prec.sigma.Z.z[i])		 # default_if 
# Z[i]            ~  dnorm(mu.Z.z[i],prec.sigma.Z.z[i]) T(Z.min.z[i],Z.max)   #   Truncated normal distribution  
mu.Z.z[i]	      <- mu.Z.0+gamma.mu.Z.D*logD[i]+gamma.mu.Z.Fz*T[i]		# default_if 
prec.sigma.Z.z[i] <- prec.sigma.Z.0*(Fz[i]^(-2.0*gamma.sigma.Z.Fz))*(D[i]^(-2.0*gamma.sigma.Z.D)) # default_if 

# (if_n.mixture==1) # Z[i]                ~	  dnorm(mu.Z.z[i],prec.sigma.Z.z[i])
# (if_n.mixture==1) # mu.Z.z[i]           <-  mu.Z.0 *flag.z* +gamma.mu.Z.D*logD[i]+gamma.mu.Z.Fz*T[i]		
# (if_n.mixture==1) # prec.sigma.Z.z[i]   <-  prec.sigma.Z.0 *flag.z* *(Fz[i]^(-2.0*gamma.sigma.Z.Fz))*(D[i]^(-2.0*gamma.sigma.Z.D))

# (if_n.mixture>1) # Z[i]				  ~	  dnormmix(mu.Z.mixture.z[i,], prec.sigma.Z.mixture.z[i,], pi)	
# (if_n.mixture>1) # for (j in 1:n.mixture) {
# (if_n.mixture>1) #     mu.Z.mixture.z[i,j]          <- mu.Z.0.mixture[j] *flag.z* +gamma.mu.Z.D*logD[i]+gamma.mu.Z.Fz*T[i]
# (if_n.mixture>1) #     prec.sigma.Z.mixture.z[i,j]  <- prec.sigma.Z.0.mixture[j] *flag.z* *(Fz[i]^(-2.0*gamma.sigma.Z.Fz))*(D[i]^(-2.0*gamma.sigma.Z.D))
# (if_n.mixture>1) #     }

# (if_mu.Z.min.0!=-n.large) # # (if_sigma.Z.min.0!=0) #  Z.min.z[i]              ~  dnorm(mu.Z.min.z[i],prec.sigma.Z.min.z[i])
# (if_mu.Z.min.0!=-n.large) # # (if_sigma.Z.min.0!=0) #  mu.Z.min.z[i]           <- mu.Z.min.0 *flag.z* +gamma.mu.Z.min.D*logD[i]+gamma.mu.Z.min.Fz*T[i]
# (if_mu.Z.min.0!=-n.large) # # (if_sigma.Z.min.0!=0) #  prec.sigma.Z.min.z[i]   <- prec.sigma.Z.min.0 *flag.z* *(Fz[i]^(-2.0*gamma.sigma.Z.min.Fz))*(D[i]^(-2.0*gamma.sigma.Z.min.D))

# (if_mu.Z.min.0!=-n.large) # # (if_sigma.Z.min.0==0) #  Z.min.z[i]              <- mu.Z.min.0 *flag.z* +gamma.mu.Z.min.D*logD[i]+gamma.mu.Z.min.Fz*T[i]


}

#=
#================= priors ============================
#=====================================================
#=
#====== Y-Z scaling =================
alpha.YIZ  ~  dunif(-n.large,n.large)
beta.YIZ   ~  dt(0,1,1)
*flag.z* gamma.YIZ  ~  dt(0,1,1)
*flag.delta.YIZ* delta.YIZ  ~  dt(0,1,1)
#=
#====== Y-Z scatter =================
prec.sigma.YIZ.0	~  dgamma(1/n.large,1/n.large) 
sigma.YIZ.0	        <- 1.0/sqrt(prec.sigma.YIZ.0)              
*flag.z* gamma.sigma.YIZ.Fz	~	dt(0,1,1)
*flag.z* gamma.sigma.YIZ.D	<-	0.0
#=
#====== X-Z scaling =================
alpha.XIZ  <-  0.0
beta.XIZ   <-  1.0
*flag.z* gamma.XIZ  <-  0.0
*flag.delta.XIZ* delta.XIZ  <-  0.0
#=
#====== X-Z scatter =================
prec.sigma.XIZ.0	~	dgamma(1/n.large,1/n.large) # default_if
sigma.XIZ.0         <-	1.0/sqrt(prec.sigma.XIZ.0)      # default_if
# (if_sigma.XIZ.0!=0) #  prec.sigma.XIZ.0 ~	 dgamma(1/n.large,1/n.large) 
# (if_sigma.XIZ.0!=0) #  sigma.XIZ.0      <- 1.0/sqrt(prec.sigma.XIZ.0)
# (if_sigma.XIZ.0==0) #  sigma.XIZ.0      <- 0.0        
*flag.z* gamma.sigma.XIZ.Fz ~  dt(0,1,1)                 
*flag.z* gamma.sigma.XIZ.D	<- 0.0                   
#=
##====== rho.XYIZ ====================
rho.XYIZ.0        <- 0.0  # default_if
*flag.z* gamma.rho.XYIZ.Fz <- 0.0  # default_if
*flag.z* gamma.rho.XYIZ.D  <- 0.0  # default_if
# (if_rho.XYIZ!=0) # rho.XYIZ.0  ~  dunif(-1.,1.)
# (if_rho.XYIZ!=0) # *flag.z* gamma.rho.XYIZ.Fz ~  dt(0,1,1)
# (if_rho.XYIZ!=0) # *flag.z* gamma.rho.XYIZ.D  <- 0.0
#=
#====== knee ===========================
# (if_beta.YIZ.knee!=beta.YIZ) #  beta.YIZ.knee         ~  dt(0,1,1)
# (if_beta.YIZ.knee!=beta.YIZ) #  Z.knee                ~  dunif(-n.large,n.large)
# (if_beta.YIZ.knee!=beta.YIZ) #  l.knee                <-  1e-04
# (if_sigma.YIZ.0.knee!=sigma.YIZ.0) # prec.sigma.YIZ.0.knee ~  dgamma(1/n.large,1/n.large)
# (if_sigma.YIZ.0.knee!=sigma.YIZ.0) # sigma.YIZ.0.knee	     <- 1.0/sqrt(prec.sigma.YIZ.0.knee)
# (if_delta.YIZ.knee!=delta.YIZ) # *flag.delta.YIZ* delta.YIZ.knee    ~  dt(0,1,1)
# (if_delta.YIZ.knee==delta.YIZ) # *flag.delta.YIZ* delta.YIZ.knee    <- delta.YIZ
#=
#====== p(Z) ============================
# (if_n.mixture==1) # pi[1] <- 1.0

# (if_n.mixture>1) # pi ~ ddirch(alpha.dirch)

# (if_n.mixture>1) # mu.Z.0.mixture[1]         <- mu.Z.0
# (if_n.mixture>1) # prec.sigma.Z.0.mixture[1] <- prec.sigma.Z.0
# (if_n.mixture>1) # for(j in 2:n.mixture) {
# (if_n.mixture>1) #     mu.Z.0.mixture[j]         ~  dunif(-n.large,n.large)
# (if_n.mixture>1) #     prec.sigma.Z.0.mixture[j] ~  dgamma(1.0/n.large,1.0/n.large)
# (if_n.mixture>1) #     sigma.Z.0.mixture[j]      <- 1.0/sqrt(prec.sigma.Z.0.mixture[j])
# (if_n.mixture>1) #     }

mu.Z.0	~	dunif(-n.large,n.large)
*flag.z* gamma.mu.Z.Fz ~  dt(0,1,1)
*flag.z* gamma.mu.Z.D  ~  dt(0,1,1)

prec.sigma.Z.0  ~  dgamma(1.0/n.large,1.0/n.large) 
sigma.Z.0       <- 1.0/sqrt(prec.sigma.Z.0)
*flag.z* gamma.sigma.Z.Fz	~	dt(0,1,1)
*flag.z* gamma.sigma.Z.D	<-  0.0
#=
# (if_mu.Z.min.0!=-n.large) #    mu.Z.min.0 ~ dunif(-n.large,n.large)
# (if_mu.Z.min.0!=-n.large) #    *flag.z* gamma.mu.Z.min.Fz ~  dt(0,1,1)
# (if_mu.Z.min.0!=-n.large) #    *flag.z* gamma.mu.Z.min.D  ~  dt(0,1,1)
#=
# (if_sigma.Z.min.0!=0) #       prec.sigma.Z.min.0 ~  dgamma(1.0/n.large,1.0/n.large)
# (if_sigma.Z.min.0!=0) #       sigma.Z.min.0      <- 1.0/sqrt(prec.sigma.Z.min.0)
# (if_sigma.Z.min.0!=0) #       *flag.z* gamma.sigma.Z.min.Fz	~	dt(0,1,1)
# (if_sigma.Z.min.0!=0) #       *flag.z* gamma.sigma.Z.min.D	<-  0.0
# (if_sigma.Z.min.0==0) #       sigma.Z.min.0      <- 0.0

}"


    #============ from string to matrix =================================
    ScriptJAGSTmp1 <- unlist(strsplit(ScriptJAGS.input,split="\n"))
    ScriptJAGSTmp2 <- NULL
    for (Ind in 1:length(ScriptJAGSTmp1)){ScriptJAGSTmp2[[Ind]] <- unlist(strsplit(ScriptJAGSTmp1[[Ind]], "[[:space:]]+"))}
    n.col.max <- max(sapply(ScriptJAGSTmp2,length))
    ScriptJAGS <- matrix(" ",length(ScriptJAGSTmp2),n.col.max)
    for (i in 1:length(ScriptJAGSTmp2)) { ScriptJAGS[i,] <- c(ScriptJAGSTmp2[[i]],rep("",n.col.max-length(ScriptJAGSTmp2[[i]] ) )) }

	#=========================== JAGS data ==============================
	data.jags<-list('x' = x, 'y' = y, 'n.data'=n.data, 'n.large'=n.large)
	
	# =========== model parameters ===============
	ParStringAllList <- c("alpha.YIZ", "beta.YIZ", "gamma.YIZ", "delta.YIZ", "sigma.YIZ.0", "gamma.sigma.YIZ.Fz","gamma.sigma.YIZ.D",  
	"alpha.XIZ", "beta.XIZ", "gamma.XIZ", "delta.XIZ", "sigma.XIZ.0", "gamma.sigma.XIZ.Fz","gamma.sigma.XIZ.D", 	
	"rho.XYIZ.0", "gamma.rho.XYIZ.Fz", "gamma.rho.XYIZ.D",
	"Z.knee", "l.knee", "beta.YIZ.knee", "delta.YIZ.knee", "sigma.YIZ.0.knee",
	"pi[1]","mu.Z.0",  "gamma.mu.Z.Fz", "gamma.mu.Z.D",
	"sigma.Z.0", "gamma.sigma.Z.Fz", "gamma.sigma.Z.D",
	"mu.Z.min.0","gamma.mu.Z.min.Fz", "gamma.mu.Z.min.D",
	"sigma.Z.min.0", "gamma.sigma.Z.min.Fz", "gamma.sigma.Z.min.D")
	ParStringJAGSList <- c("alpha.YIZ", "beta.YIZ", "sigma.YIZ.0","pi[1]","mu.Z.0", "sigma.Z.0")
	
	if(n.mixture>1){
		for(Ind in 2:n.mixture){
		par.to.add <-  c(paste0("pi[", Ind,"]" ) , paste0("mu.Z.0.mixture[", Ind,"]" ), paste0("sigma.Z.0.mixture[", Ind,"]" ) )
		ParStringAllList <- c(ParStringAllList, par.to.add )
		ParStringJAGSList <- c(ParStringJAGSList, par.to.add )
		}
	}
	
	n.mcmc<-n.chains*ceiling(n.iter/thin)
	mcmc.all<-data.frame(matrix(double(), n.mcmc, length(ParStringAllList), dimnames=list(c(), ParStringAllList)), stringsAsFactors=F)
	

	
	# =============== editing functions==============================================================
		
	fIfScriptJAGS <- function(if.string){
		for (j in  which(ScriptJAGS[,2]==if.string) ) {
			ScriptJAGS[j,] <- c(ScriptJAGS[j,4:n.col.max],"","","")
			}
			assign('ScriptJAGS',ScriptJAGS,environment(fPriorScriptJAGS))
			}
			
	fAppendScriptJAGS <- function(var,string.to.append){
		for (iRow in c(which(ScriptJAGS[,1]==var),which(ScriptJAGS[,4]==var))){
			iCol<- max(which(ScriptJAGS[iRow,]!=""))
			ScriptJAGS[iRow,iCol]<-paste0(ScriptJAGS[iRow,iCol],string.to.append)
			}
			assign('ScriptJAGS',ScriptJAGS,environment(fPriorScriptJAGS))
			}
	
	fPriorScriptJAGS <- function(par, prior) {
		
		if(substr(par,1,5)=="sigma" && prior!="prec.dgamma"){
			par.prec<-sub("s","prec.s",par)
			prior.string.prec<- c(par.prec,"<-",paste0("1.0/",par,"^2.0",sep=""))
			for (Ind in which(ScriptJAGS[,1]==par.prec)){ScriptJAGS[Ind,1:3]<-prior.string.prec}
			for (Ind in which(ScriptJAGS[,4]==par.prec)){ScriptJAGS[Ind,4:6]<-prior.string.prec}
		}
		
		if  (prior == "dt") {
				prior.string  <- c(par,"~","dt(0.0,1.0,1.0)")
				ParStringJAGSList <- c(ParStringJAGSList,par)
			}
		else if (prior == "dunif"){
				prior.string  <- c(par,"~","dunif(-n.large,n.large)")
				ParStringJAGSList <- c(ParStringJAGSList,par)
    		}
    	else if (prior == "dgamma"){
    			prior.string  <- c(par,"~","dgamma(1.0/n.large,1.0/n.large)")
    			ParStringJAGSList <- c(ParStringJAGSList,par)
    			}
    	else if ( class(prior)=="numeric"){
    			ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != par]
    			prior.string <- c(par, "<-", prior)
    			mcmc.all[[which(ParStringAllList==par)]]<<-prior   # mcmc.all$par.tmp1<-prior
    			}
    	else if (prior == "ddirch"){
				prior.string  <- c(par,"~","ddirch(alpha.dirch)")
				ParStringJAGSList <- c(ParStringJAGSList,par)
				}
    	else {
    			prior.string  <- c(par,"~",prior)
    			ParStringJAGSList <- c(ParStringJAGSList,par)
    			}
    	    
        if (prior!="prec.dgamma"){
        	for (Ind in which(ScriptJAGS[,1]==par)){ScriptJAGS[Ind,1:3]<-prior.string}
        	for (Ind in which(ScriptJAGS[,4]==par)){ScriptJAGS[Ind,4:6]<-prior.string}
        }
    	
    	assign('mcmc.all',mcmc.all,environment(fPriorScriptJAGS))
    	assign('ScriptJAGS',ScriptJAGS,envir=environment(fPriorScriptJAGS))
    	assign('ParStringJAGSList',ParStringJAGSList,environment(fPriorScriptJAGS))	
    }
    
       
  	# ============ time evolution ============================================== 
	ind.z.flag<-which(ScriptJAGS=="*flag.z*",arr.ind=TRUE)
	if(z.logical){	# evolution
		#ScriptJAGS[ind.z.flag]<-""
		for (Ind in 1:nrow(ind.z.flag) ) { ScriptJAGS[ind.z.flag[Ind,1],ind.z.flag[Ind,2]:n.col.max] <-c(ScriptJAGS[ind.z.flag[Ind,1],(ind.z.flag[Ind,2]+1):n.col.max],"")}
		fPriorScriptJAGS("gamma.YIZ", gamma.YIZ)
		fPriorScriptJAGS("gamma.XIZ", gamma.XIZ)
		fPriorScriptJAGS("gamma.sigma.YIZ.D", gamma.sigma.YIZ.D)
		fPriorScriptJAGS("gamma.sigma.YIZ.Fz", gamma.sigma.YIZ.Fz)
		fPriorScriptJAGS("gamma.sigma.XIZ.D", gamma.sigma.XIZ.D)
		fPriorScriptJAGS("gamma.sigma.XIZ.Fz",gamma.sigma.XIZ.Fz)
		fPriorScriptJAGS("gamma.sigma.Z.D", gamma.sigma.Z.D)
		fPriorScriptJAGS("gamma.sigma.Z.Fz", gamma.sigma.Z.Fz)
		fPriorScriptJAGS("gamma.mu.Z.D", gamma.mu.Z.D)
		fPriorScriptJAGS("gamma.mu.Z.Fz", gamma.mu.Z.Fz) 
		fPriorScriptJAGS("gamma.mu.Z.min.Fz", gamma.mu.Z.min.Fz)
		fPriorScriptJAGS("gamma.mu.Z.min.D", gamma.mu.Z.min.D)
		fPriorScriptJAGS("gamma.sigma.Z.min.Fz",gamma.sigma.Z.min.Fz)
		fPriorScriptJAGS("gamma.sigma.Z.min.D", gamma.sigma.Z.min.D)
    
		switch(match.arg(distance),
		       "angular"={gamma.Dist<-0},
		       "luminosity"={gamma.Dist<-2}
		)
		ind.flag.delta.YIZ<-which(ScriptJAGS=="*flag.delta.YIZ*",arr.ind=TRUE)
		if(delta.YIZ!=0){
			for (Ind in 1:nrow(ind.flag.delta.YIZ) ) { ScriptJAGS[ind.flag.delta.YIZ[Ind,1],ind.flag.delta.YIZ[Ind,2]:n.col.max] <-c(ScriptJAGS[ind.flag.delta.YIZ[Ind,1],(ind.flag.delta.YIZ[Ind,2]+1):n.col.max],"")}
			fPriorScriptJAGS("delta.YIZ", delta.YIZ)
			}
		else{
			for (Ind in 1:nrow(ind.flag.delta.YIZ) ) { ScriptJAGS[ind.flag.delta.YIZ[Ind,1],ind.flag.delta.YIZ[Ind,2]:n.col.max] <-""}
			fPriorScriptJAGS("delta.YIZ", 0.0)
			}
		ind.flag.delta.XIZ<-which(ScriptJAGS=="*flag.delta.XIZ*",arr.ind=TRUE)
		if(delta.XIZ!=0){
			for (Ind in 1:nrow(ind.flag.delta.XIZ) ) { ScriptJAGS[ind.flag.delta.XIZ[Ind,1],ind.flag.delta.XIZ[Ind,2]:n.col.max] <-c(ScriptJAGS[ind.flag.delta.XIZ[Ind,1],(ind.flag.delta.XIZ[Ind,2]+1):n.col.max],"")}
			fPriorScriptJAGS("delta.XIZ", delta.XIZ)
			}
		else{
			for (Ind in 1:nrow(ind.flag.delta.XIZ) ) { ScriptJAGS[ind.flag.delta.XIZ[Ind,1],ind.flag.delta.XIZ[Ind,2]:n.col.max] <-""}
			fPriorScriptJAGS("delta.XIZ", 0.0)
			}
		
		D.ref <- (1.+z.ref)^gamma.Dist * distance.LCDM(0.,z.ref,Omega.M0, Omega.L0)
		D <- sapply(z, function(x)(1.+x)^gamma.Dist*distance.LCDM(0., x, Omega.M0, Omega.L0)/D.ref)
		if (time.factor=="1+z"){Fz <- (1.+z)/(1+z.ref) } 
		else if (time.factor=="Ez"){
			H.ref <- Hubble.LCDM(z.ref, Omega.M0, Omega.L0)
			Fz <- Hubble.LCDM(z, Omega.M0, Omega.L0)/H.ref
			}
		else {
			warning("Time factor should be one of \"1+z\" or \"Ez\". Set to default \"Ez\""); 
			H.ref <- Hubble.LCDM(z.ref, Omega.M0, Omega.L0)
			Fz <- Hubble.LCDM(z, Omega.M0, Omega.L0)/H.ref
			}
		logD<-log10(D)
		logFz<-log10(Fz)
		data.jags<-c(data.jags,list('D'=D, 'Fz'=Fz, 'logD'=logD, 'T'=logFz))
	} else{ 
		for (Ind in 1:nrow(ind.z.flag) ) { ScriptJAGS[ind.z.flag[Ind,1],ind.z.flag[Ind,2]:n.col.max] <-""}
		ind.flag.delta.YIZ<-which(ScriptJAGS=="*flag.delta.YIZ*",arr.ind=TRUE)
		for (Ind in 1:nrow(ind.flag.delta.YIZ) ) { ScriptJAGS[ind.flag.delta.YIZ[Ind,1],ind.flag.delta.YIZ[Ind,2]:n.col.max] <-""}
		ind.flag.delta.XIZ<-which(ScriptJAGS=="*flag.delta.XIZ*",arr.ind=TRUE)
		for (Ind in 1:nrow(ind.flag.delta.XIZ) ) { ScriptJAGS[ind.flag.delta.XIZ[Ind,1],ind.flag.delta.XIZ[Ind,2]:n.col.max] <-""}
		mcmc.all$"gamma.YIZ" <- 0.0; mcmc.all$"delta.YIZ" <- 0.0; mcmc.all$"delta.YIZ.knee" <- 0.0; mcmc.all$"gamma.XIZ" <- 0.0;  mcmc.all$"delta.XIZ" <- 0.0;
		mcmc.all$"gamma.sigma.YIZ.D" <- 0.0; mcmc.all$"gamma.sigma.YIZ.Fz" <- 0.0; mcmc.all$"gamma.sigma.XIZ.D" <- 0.0; mcmc.all$"gamma.sigma.XIZ.Fz" <- 0.0; 
		mcmc.all$"gamma.rho.XYIZ.Fz" <- 0.0; mcmc.all$"gamma.rho.XYIZ.D" <- 0.0;
		mcmc.all$"gamma.sigma.Z.D" <- 0.0; mcmc.all$"gamma.sigma.Z.Fz" <- 0.0; mcmc.all$"gamma.mu.Z.D" <- 0.0; mcmc.all$"gamma.mu.Z.Fz" <- 0.0;
		mcmc.all$"gamma.sigma.Z.min.D" <- 0.0; mcmc.all$"gamma.sigma.Z.min.Fz" <- 0.0; 
		mcmc.all$"gamma.mu.Z.min.D" <- 0.0; mcmc.all$"gamma.mu.Z.min.Fz" <- 0.0;
    }
    
    #================== priors =================================================
    
    fPriorScriptJAGS("alpha.YIZ",alpha.YIZ)
    fPriorScriptJAGS("beta.YIZ", beta.YIZ)
    fPriorScriptJAGS("sigma.YIZ.0", sigma.YIZ.0)
   
    fPriorScriptJAGS("alpha.XIZ", alpha.XIZ)
    fPriorScriptJAGS("beta.XIZ", beta.XIZ)
    fPriorScriptJAGS("sigma.XIZ.0", sigma.XIZ.0)
    
    fPriorScriptJAGS("mu.Z.0", mu.Z.0)
    fPriorScriptJAGS("sigma.Z.0", sigma.Z.0)
      
    #================ delta.x=0 ==============
	if(delta.x.logical){fIfScriptJAGS("(if_delta.x!=0)")
		data.jags<-c(data.jags,list('prec.delta.x'=delta.x^-2))
	}
	else{
		fIfScriptJAGS("(if_delta.x==0)")
		names(data.jags)[1] <- 'X'
	}
	
	if(delta.y.logical){
		fIfScriptJAGS("(if_delta.y!=0)")
		data.jags<-c(data.jags,list('prec.delta.y'=delta.y^-2))
		if (delta.x.logical) {data.jags<-c(data.jags,list('rho.xy'=rho.xy))}
	}
	else {
		fIfScriptJAGS("(if_delta.y==0)")
		names(data.jags)[2] <- 'Y'
	}
      
   #========== knee=======================================================   
   if (beta.YIZ.knee!="beta.YIZ") {
   	   fIfScriptJAGS("(if_beta.YIZ.knee!=beta.YIZ)")
   	   fPriorScriptJAGS("Z.knee", Z.knee)
   	   fPriorScriptJAGS("l.knee", l.knee)
   	   fPriorScriptJAGS("beta.YIZ.knee", beta.YIZ.knee)
   	   if(delta.YIZ.knee!="delta.YIZ" && z.logical){
   	   	if(delta.YIZ!=0){
   	   	   	fIfScriptJAGS("(if_delta.YIZ.knee!=delta.YIZ)")
   	   		fPriorScriptJAGS("delta.YIZ.knee", delta.YIZ.knee)
   	   		} else{
   	   			warning("delta.YIZ.knee set to zero as delta.YIZ")
   	   			fIfScriptJAGS("(if_delta.YIZ.knee==delta.YIZ)")
   	   			delta.YIZ.knee="delta.YIZ"
   	   			}
   	   	} else {
		fIfScriptJAGS("(if_delta.YIZ.knee==delta.YIZ)")
		}
		if (sigma.YIZ.0.knee != "sigma.YIZ.0") {
				fIfScriptJAGS("(if_sigma.YIZ.0.knee!=sigma.YIZ.0)")
				fPriorScriptJAGS("sigma.YIZ.0.knee", sigma.YIZ.0.knee)
				} else {
					fIfScriptJAGS("(if_sigma.YIZ.0.knee==sigma.YIZ.0)")
			}
		} else {
		fIfScriptJAGS("(if_beta.YIZ.knee==beta.YIZ)")
		delta.YIZ.knee="delta.YIZ"
		
		#fIfScriptJAGS("(if_gamma.YIZ.knee==gamma.YIZ)")
		#fIfScriptJAGS("(if_delta.YIZ.knee==delta.YIZ)")
		#fIfScriptJAGS("(if_sigma.YIZ.0.knee==sigma.YIZ.0)")
	}
	 
	 
	# =========X-Z scatter==================================================
	if (sigma.XIZ.0 != 0.0){
		fIfScriptJAGS("(if_sigma.XIZ.0!=0)")
		fPriorScriptJAGS("sigma.XIZ.0", sigma.XIZ.0)
		if (rho.XYIZ.0 == 0.0){
			fIfScriptJAGS("(if_rho.XYIZ==0)")
			mcmc.all$"rho.XYIZ.0" <- 0.0
			mcmc.all$"gamma.rho.XYIZ.Fz" <- 0.0
			mcmc.all$"gamma.rho.XYIZ.D" <- 0.0
		} else if (rho.XYIZ.0 != 0.0){
			fIfScriptJAGS("(if_rho.XYIZ!=0)")
			fPriorScriptJAGS("rho.XYIZ.0", rho.XYIZ.0)
			if(z.logical){
			fPriorScriptJAGS("gamma.rho.XYIZ.Fz", gamma.rho.XYIZ.Fz)
			fPriorScriptJAGS("gamma.rho.XYIZ.D", gamma.rho.XYIZ.D)
			}
		} 
	} else {
		if(!delta.x.logical){
			names(data.jags)[1] <- 'Z';
			if(alpha.XIZ!=0.0 || beta.XIZ!=1.0 || gamma.XIZ!=0.0 || delta.XIZ!=0.0){ 
			warning("x rescaled to x=Z")
			fPriorScriptJAGS("alpha.XIZ", 0.0); fPriorScriptJAGS("beta.XIZ", 1.0); fPriorScriptJAGS("gamma.XIZ", 0.0); fPriorScriptJAGS("delta.XIZ", 0.0)
			}
			if(Z.monitored && !anyNA(x)){
				warning("x without uncertainties and X not scattered: x and Z identified")
				Z.monitored=FALSE
				}
			}
		fIfScriptJAGS("(if_sigma.XIZ.0==0)") 
		fIfScriptJAGS("(if_rho.XYIZ==0)")
		if(rho.XYIZ.0!=0)warning("correlation rho.XYIZ.0 set to 0")
		for (Par in c("gamma.sigma.XIZ.Fz","gamma.sigma.XIZ.D","gamma.rho.XYIZ.Fz","gamma.rho.XYIZ.D")) {fPriorScriptJAGS(Par,0.0)}
		mcmc.all$"rho.XYIZ.0" <- 0.0
		mcmc.all$"gamma.rho.XYIZ.Fz" <- 0.0
		mcmc.all$"gamma.rho.XYIZ.D" <- 0.0
		mcmc.all$"sigma.XIZ.0" <- 0.0
   }
   
   # =========== Mixture =============================================
   if (n.mixture > 1){
   	   data.jags <- c(data.jags,list('n.mixture' = n.mixture,'alpha.dirch' = rep(1.0,n.mixture)))
   	   if (Z.min!="-n.large"){warning("JAGS distribution dnormmix cannot be truncated. Z.min superseded"); Z.min<-"-n.large"}
   	   if (Z.max!="n.large"){warning("JAGS distribution dnormmix cannot be truncated. Z.max superseded"); Z.max<-"n.large"}
   	   fIfScriptJAGS("(if_n.mixture>1)")
   	   fPriorScriptJAGS("mu.Z.0.mixture[j]", mu.Z.0.mixture)
   	   ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "mu.Z.0.mixture[j]"]
   	   fPriorScriptJAGS("sigma.Z.0.mixture[j]", sigma.Z.0.mixture)
   	   ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "sigma.Z.0.mixture[j]"]
   	   fPriorScriptJAGS("pi", pi)
   	   ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "pi"]
	} else {
	fIfScriptJAGS("(if_n.mixture==1)")
	mcmc.all$"pi.1." <- 1.0
	ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "pi[1]"]
	}
		
	# =========== Truncated p(Z) distribution =============================================
	mcmc.all$"mu.Z.min.0" <- -n.large
	mcmc.all$"gamma.mu.Z.min.D" <- 0.0
	mcmc.all$"gamma.mu.Z.min.Fz" <- 0.0
	mcmc.all$"sigma.Z.min.0" <- 0.0
	mcmc.all$"gamma.sigma.Z.min.D" <- 0.0
	mcmc.all$"gamma.sigma.Z.min.Fz" <- 0.0
		
	if (mu.Z.min.0 != "-n.large" || Z.max != "n.large") {
		fAppendScriptJAGS("Z[i]",paste0("T( Z.min.z[i], ", toString(Z.max), ")"))
		}
			
	if (mu.Z.min.0 != "-n.large") {
		fIfScriptJAGS("(if_mu.Z.min.0!=-n.large)")
		fPriorScriptJAGS("mu.Z.min.0", mu.Z.min.0)
		if(sigma.Z.min.0!=0){
			fIfScriptJAGS("(if_sigma.Z.min.0!=0)");
			fPriorScriptJAGS("sigma.Z.min.0", sigma.Z.min.0)
			}
		else {fIfScriptJAGS("(if_sigma.Z.min.0==0)")}
	}
	else {
		ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "gamma.mu.Z.min.D"]
		ParStringJAGSList <- ParStringJAGSList[ParStringJAGSList != "gamma.mu.Z.min.Fz"]
	}
	
		
		#if (is.na(match(mu.Z.0,names.arg)) || mu.Z.0=="dunif") {fPriorScriptJAGS("mu.Z.0", paste0("dunif(",toString(Z.min),",",toString(Z.max),")",sep="") ) } 
		#else if (class(mu.Z.0)!="numeric" && !is.na(match(mu.Z.0,names.arg)) )
		#{fAppendScriptJAGS("mu.Z.0",paste0("T(",toString(Z.min), ",", toString(Z.max), ")"))}
	

	# ============== Malmquist bias =============================================
	if (y.threshold.logical)  {
		fIfScriptJAGS("(if_y.threshold!={})")
		fAppendScriptJAGS("y[i]","T(y.min[i],)")
		fAppendScriptJAGS("Y[i]","T(Y.min[i],)")
		if (!delta.y.threshold.logical){
			data.jags<-c(data.jags,list('y.threshold'=y.threshold))
			fIfScriptJAGS("(if_delta.y.threshold=={})")
		}
		else if (delta.y.threshold.logical){
			data.jags<-c(data.jags,list('y.threshold'=y.threshold, 'delta.y.threshold'=delta.y.threshold))
			fIfScriptJAGS("(if_delta.y.threshold!={})")
		}
	}     
    
    # ============== Malmquist bias in X ===========================================
	if (x.threshold.logical)  {
		fIfScriptJAGS("(if_x.threshold!={})")
		fAppendScriptJAGS("x[i]","T(x.min[i],)")
		if (sigma.XIZ.0 != 0.0){ fAppendScriptJAGS("X[i]","T(X.min[i],)")}
		if (!delta.x.threshold.logical){
			data.jags<-c(data.jags,list('x.threshold'=x.threshold))
			fIfScriptJAGS("(if_delta.x.threshold=={})")
		}
		else if (delta.x.threshold.logical){
			data.jags<-c(data.jags,list('x.threshold'=x.threshold, 'delta.x.threshold'=delta.x.threshold))
			fIfScriptJAGS("(if_delta.x.threshold!={})")
		}
	}     

    
    #==================== export JAGS Script============================
    ScriptJAGS[1,4]<-paste0("created with lira on ",date(),sep="")
    ScriptJAGS[which(ScriptJAGS=="default_if",arr.ind = TRUE)[,1], ]<-rep("",n.col.max)
    ScriptJAGS[which(ScriptJAGS=="#",arr.ind = TRUE)[,1], ]<-rep("",n.col.max)
    
    IndRowList<-NULL
    for (iRow in 1:nrow(ScriptJAGS)){
	if(
	gsub('([[:space:]]+)\\1+', '\\1', paste0(unique(ScriptJAGS[iRow,]),collapse=""))!=" " &&
	gsub('([[:space:]]+)\\1+', '\\1', paste0(unique(ScriptJAGS[iRow,]),collapse=""))!=""  &&
	ScriptJAGS[iRow,1]!="#" ){IndRowList<-c(IndRowList,iRow)}
	}
	
	ScriptJAGSText <- paste( apply(ScriptJAGS[IndRowList,], 1, paste, collapse = " "), collapse="\n")
	
	#======= check =====================================
#	write.table(ScriptJAGS[IndRowList,],file = "/Users/maurosereno/Documents/calcoli/R_JAGS/tmp/lira_check.jags", quote = FALSE,row.names = FALSE,col.names = FALSE)
	
	#================rjags====================================
	#========================================================
	
	ParStringJAGSListScaling <- sort(ParStringJAGSList[which(!duplicated(ParStringJAGSList))])
		
	if(Z.monitored) { ParStringJAGSList <- c(ParStringJAGSList, paste0("Z[", 1:n.data,"]") ) }
	if(X.monitored) { ParStringJAGSList <- c(ParStringJAGSList, paste0("X[", 1:n.data,"]") ) }
	if(XZ.monitored){ ParStringJAGSList <- c(ParStringJAGSList, paste0("fXZ[", 1:n.data,"]") ) }
	if(Y.monitored) { ParStringJAGSList <- c(ParStringJAGSList, paste0("Y[", 1:n.data,"]") ) }
    if(YZ.monitored){ ParStringJAGSList <- c(ParStringJAGSList, paste0("fYZ[", 1:n.data,"]") ) }


	ParStringJAGSList <- sort(ParStringJAGSList[which(!duplicated(ParStringJAGSList))])
	                 
	                 	
	#=========== chains ======================================
	
	load.module("mix")
	set.factory("mix::TemperedMix", 'sampler', FALSE)
	
	jm <- jags.model(textConnection(ScriptJAGSText), data = data.jags, inits=inits, n.adapt=n.adapt, n.chains=n.chains, quiet=quiet)
	update(jm)
	mcmc.samples <- coda.samples(jm, ParStringJAGSList, n.iter=n.iter, thin=thin)
	
	mcmc.jags <- as.data.frame(mcmc.samples[[1]])
	if(n.chains>1){for (Ind in 2:n.chains){mcmc.jags<-rbind(mcmc.jags,as.data.frame(mcmc.samples[[Ind]]))} }
	
	
	mcmc.X <- NULL
	if(X.monitored){mcmc.X <- mcmc.jags[,paste0("X[", 1:n.data,"]")]}
	mcmc.Z <- NULL
	if(Z.monitored){mcmc.Z <- mcmc.jags[,paste0("Z[", 1:n.data,"]")]}
	mcmc.Y <- NULL
	if(Y.monitored){mcmc.Y <- mcmc.jags[,paste0("Y[", 1:n.data,"]")]}
	mcmc.YZ <- NULL
	if(YZ.monitored){mcmc.YZ <- mcmc.jags[,paste0("fYZ[", 1:n.data,"]")]}
	mcmc.XZ <- NULL
	if(XZ.monitored){mcmc.XZ <- mcmc.jags[,paste0("fXZ[", 1:n.data,"]")]}
	
#	write.table(mcmc.X, file = "/Users/maurosereno/Documents/calcoli/JAGS/mcmc_X.jags", quote = FALSE,row.names = FALSE,col.names = FALSE)
	
	for (par in ParStringJAGSListScaling){mcmc.all[[which(ParStringAllList==par)]] <- mcmc.jags[[which(colnames(mcmc.jags)==par)]] }
	
	if (beta.YIZ.knee == "beta.YIZ"){
		mcmc.all$Z.knee <- 0.0
		mcmc.all$l.knee <- l.knee 
		mcmc.all$beta.YIZ.knee <- mcmc.all$beta.YIZ 
		mcmc.all$delta.YIZ.knee <- mcmc.all$delta.YIZ
		mcmc.all$sigma.YIZ.0.knee <- mcmc.all$sigma.YIZ.0
	}
	else {
		if (delta.YIZ.knee == "delta.YIZ"){mcmc.all$delta.YIZ.knee <- mcmc.all$delta.YIZ}
		if (sigma.YIZ.0.knee == "sigma.YIZ.0") {mcmc.all$sigma.YIZ.0.knee <- mcmc.all$sigma.YIZ.0}
	}
	
	if(export){
		if(export.jags==""){
			warning("No file name specified: JAGS script saved as \"lira.jags\" in the working directory")
			export.jags<-paste0(getwd(),"/lira.jags",sep = "")
			}
		if(export.mcmc==""){
			warning("No file name specified for merged chains: regression parameters saved as \"mcmc_all.dat\" in the working directory")
			export.mcmc<-paste0(getwd(),"/mcmc_all.dat")
			}
		if(X.monitored){
			if(export.X==""){
				warning("no file name specified for merged chains: no export of the X variables")
				#warning("X merged chains saved as \"mcmc_X.dat\" in the working directory")
				#export.X<-paste0(getwd(),"/mcmc_X.dat")
				}
			else {write.table(mcmc.X,file=export.X, sep="\t",row.names=FALSE, col.names=TRUE)}
		}
		if(Z.monitored){
			if(export.Z==""){
				warning("no file name specified for merged chains: no export of the Z variables")}
			else {write.table(mcmc.Z,file=export.Z, sep="\t",row.names=FALSE, col.names=TRUE)}
		}
		if(YZ.monitored){
			if(export.YZ==""){
				warning("no file name specified for merged chains: no export of the YZ variables")}
			else {write.table(mcmc.YZ,file=export.YZ, sep="\t",row.names=FALSE, col.names=TRUE)}
		}
		if(XZ.monitored){
			if(export.XZ==""){
				warning("no file name specified for merged chains: no export of the XZ variables")}
			else {write.table(mcmc.XZ,file=export.XZ, sep="\t",row.names=FALSE, col.names=TRUE)}
		}
		if(Y.monitored){
			if(export.Y==""){
				warning("no file name specified for merged chains: no export of the Y variables")}
			else {write.table(mcmc.Y,file=export.Y, sep="\t",row.names=FALSE, col.names=TRUE)}
		}
		write.table(mcmc.all,file=export.mcmc, sep="\t",row.names=FALSE, col.names=TRUE)
		write.table(ScriptJAGS[IndRowList,],file = export.jags, quote = FALSE,row.names = FALSE,col.names = FALSE)
		}
	
	if(print.summary){
		print(summary(mcmc.samples))
		}
		
	if(print.diagnostic){
		if(n.chains>1){
			print("Gelman and Rubin's convergence diagnostic")
			print(gelman.diag(mcmc.samples))}
		}
	
	return(list(list(mcmc.all, mcmc.Z, mcmc.X, mcmc.YZ, mcmc.Y), mcmc.samples))

}
