`s.check` <-
function(x,K,n,EMstart,EMoption,method,Sdist)
{

#sanity checks for method, EMoption and Sdist
if (!any(EMoption==c("maximization","classification","randomization")))
  stop("EMoption must be element of (maximization, classification, randomization).")

if (!any(method==c("separate","int.gp","main.gp","main.g","main.p")))
  stop("Method must be element of (separate, int.gp, main.gp, main.g, main.p).")
  
if ((EMoption=="maximization") && (method!="separate")) {
  warning("Maximization EM is not applicable to this method! Method is set to 'separate'.") 
  method <- "separate"
}

if (!any(Sdist==c("weibull", "exponential", "rayleigh")))
  stop("Sdist must be element of (weibull, exponential, rayleigh).")

# sanity checks for x

x[is.na(x)==TRUE] <- 0                                        #replace NA's by 0

if (any(rowSums(x)==0)) {                                     #eliminate sessions with all 0 dwell times
  x <- x[rowSums(x)!=0,]
  warning("Subjects with no visit were eliminated!")
  }

if (any(colSums(x)==0)) {                                    #eliminate pages not visited
  x <- x[,colSums(x)!=0]
  warning("Variables with complete 0 dwell times were eliminated!")
  }

#starting values for EMstart
if (length(EMstart)==1) {
  EMstart <- sample(1:K, n, replace=TRUE)                     #same starting values for all models
} else if (length(EMstart) != n) {
   stop("Mismatch between EMstart and number of subjects!")
} else if (max(EMstart)!=K) {
   stop("Number of defined groups in EMstart doesn't correspond to K!")
}
  

if (EMoption == "maximization") {
  nullmat <- matrix(rep(0,n*K),ncol=K)  
  nullmat[cbind(1:n,EMstart)] <- 1                            #converting EMstart vector into 0/1 matrix
  EMstart <- abs(jitter(nullmat))                             
}

list(x=x,EMstart=EMstart,method=method)
}

