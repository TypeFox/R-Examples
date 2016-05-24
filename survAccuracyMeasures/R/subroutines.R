#function to estimate the SE using monte carlo methods
logit = function(x) log(x/(1-x))
expit = function(x) exp(x)/(1+exp(x))








EstRTall<-function(data)
{

  ck      = data[,1]   ## aka Y
  nc      = length(ck) 
  CondSck = data[,2]   ## aka Sy
  Fck     = data[,3]      
		
  dFck    = Fck - c(0,Fck[-nc])
  St0.Fck = cumsum(CondSck*dFck)## St0.Fck = P(T> t0,Y<=ck)
  Ft0.Fck = Fck-St0.Fck         ## Ft0.Fck = P(T<=t0,Y<=ck)

  
  St0     = max(St0.Fck)        ## St0     = P(T> t0      )
  Ft0     = 1-St0               ## Ft0     = P(T<=t0      )
  
  
  FPR.c   = (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)
  TPR.c   = (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
  NPV.c   = St0.Fck/Fck           ## P(T> t0|Y<=ck)
  PPV.c   = (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)

  RT.out  = data.frame("cutoff" = ck,"RiskT" = 1-CondSck,"v"= Fck, "FPR" = FPR.c,"TPR" = TPR.c, "PPV" = PPV.c,"NPV" = NPV.c)
 # TG    = sum(abs(1-CondSck-Ft0)*(-Fck+c(Fck[-1],1)))
  AUC   = sum(TPR.c*(FPR.c-c(FPR.c[-1],0)))
#  RiskT = 1-CondSck
			
#  mmm  = length(RiskT)
#  ITPR = sum(TPR.c*(RiskT-c(0,RiskT[-mmm])))
#  IFPR = sum(FPR.c*(RiskT-c(0,RiskT[-mmm])))
#  IDI  = ITPR - IFPR

  list(RT.out = RT.out, AUC = AUC)
}




EstRTvp<-function(RT.out, uu0, typex, typey) {

  ind0.x = match(typex,c("cutoff","RiskT","v","FPR","TPR","PPV","NPV"))
  ind0.y = match(typey,c("cutoff","RiskT","v","FPR","TPR","PPV","NPV"))

  uuk = RT.out[,ind0.x]

  tmpind <- which.min(abs(uuk - uu0))[1]


  RT.out[tmpind, ind0.y]
}


sum.I<-function(yy,FUN,Yi,Vi=NULL){

	if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}

	pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')

	if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos

	if (!is.null(Vi)) {

      	if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)

		Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
          
		return(rbind(0,Vi)[pos+1,])

	} else return(pos)
}



VTM <- function(vc, dm)
{
    matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}







Est.Var.CCH.trueweights = function(N,Wexp,data,stratum, subcohort) {

  cohort.variance = colSums(data$wi*(Wexp/N)^2)
 # robust.variance = colSums((data$weights*Wexp/N)^2)

  if(subcohort){

  ## strata by cases and conrols 
  stra = sort(unique(stratum))
  nstra=length(stra);
  np = dim(Wexp)[2]
  strvar = rep(0,np);
  for (i in 1:nstra) {
    straWt = data$wi[stratum==stra[i]]	
    straWexp = as.matrix(Wexp[stratum==stra[i],])
    ns = length(straWt)	
    tempstratavar = (ns-1)*(straWt[1]-1)*straWt[1]*apply(straWexp/N,2,var)
    strvar = strvar + tempstratavar
  }
     return(cohort.variance+strvar)
  }else{
    return(cohort.variance)
  }
  
}



