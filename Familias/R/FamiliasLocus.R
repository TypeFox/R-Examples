FamiliasLocus <- function (frequencies, allelenames, name, 
			   MutationModel = "Stepwise", 
                     MutationRate  = 0, 
                     MutationRange = 0.5, 
                     MutationRate2 = 0, 
			   MutationMatrix, 
			   Stabilization = "None", 
			   MaxStabilizedMutrate = 1,
			   femaleMutationModel, 
			   femaleMutationRate, 
			   femaleMutationRange, 
			   femaleMutationRate2, 
			   femaleMutationMatrix, 
			   maleMutationModel, 
			   maleMutationRate, 
			   maleMutationRange, 
			   maleMutationRate2, 
			   maleMutationMatrix)
{
    if (missing(frequencies)) 
        stop("The first argument must be a list of frequencies or a FamiliasLocus object.")
    if(class(frequencies)=="FamiliasLocus") {
        if (!missing(name) || !missing(allelenames))
            stop("Only mutation parameters can be edited with the FamiliasLocus function.")
        x <- frequencies
        frequencies <- x$alleles
        name        <- x$locusname
        allelenames <- names(frequencies)
        if (missing(MutationModel) & missing(MutationRate) & 
            missing(MutationRange) & missing(MutationRate2) & 
            missing(MutationMatrix) & missing(femaleMutationModel) & 
            missing(femaleMutationRate) & missing(femaleMutationRange) & 
            missing(femaleMutationRate2) & missing(femaleMutationMatrix) & 
            missing(maleMutationModel) & missing(maleMutationRate) & 
            missing(maleMutationRange) & missing(maleMutationRate2) & 
            missing(maleMutationMatrix)) {
            if (missing(Stabilization)) 
                stop("When a FamiliasLocus object is given, at least one mutation parameter or 'Stabilization' must be specified.")
            MutationModel <- "Custom"
            maleMutationMatrix <- x$maleMutationMatrix
            femaleMutationMatrix <- x$femaleMutationMatrix
        }
    } else {
        if (!is.numeric(frequencies) || any(frequencies <= 0)) 
            stop("frequencies must be a vector of positive numbers.")
        if (round(sum(frequencies), 6) != 1) 
            stop("The frequencies must sum to 1.")
        if (missing(name)) 
            name <- deparse(substitute(frequencies))
        if (missing(allelenames)) {
            if (is.null(names(frequencies)))  
                allelenames <- as.character(1:length(frequencies))
            else 
                allelenames <- names(frequencies)
            }
        else if (length(allelenames) != length(frequencies)) 
            stop("The number of allele names must be the same as the number of frequencies.")
        if (anyDuplicated(allelenames)) 
            stop("There cannot be duplicates among the allele names.")
        if (any(allelenames[-length(allelenames)] == "silent") || 
            any(allelenames[-length(allelenames)] == "Silent")) 
            stop("Only the last allele can be specified as silent.")
        names(frequencies) <- allelenames
    }
    nAlleles <- length(frequencies)
    hasSilentAllele <- (allelenames[nAlleles] == "silent" || allelenames[nAlleles] == "Silent")
    nAll <- nAlleles - hasSilentAllele
    freq <- frequencies[1:nAll]

    if (missing(femaleMutationModel)) femaleMutationModel <- MutationModel
    if (missing(femaleMutationRate)) femaleMutationRate <- MutationRate
    if (missing(femaleMutationRange)) femaleMutationRange <- MutationRange
    if (missing(femaleMutationRate2)) femaleMutationRate2 <- MutationRate2
    if (missing(femaleMutationMatrix) && !missing(MutationMatrix)) femaleMutationMatrix <- MutationMatrix
    if (missing(maleMutationModel)) maleMutationModel <- MutationModel
    if (missing(maleMutationRate)) maleMutationRate <- MutationRate
    if (missing(maleMutationRange)) maleMutationRange <- MutationRange
    if (missing(maleMutationRate2)) maleMutationRate2 <- MutationRate2
    if (missing(maleMutationMatrix) && !missing(MutationMatrix)) maleMutationMatrix <- MutationMatrix

    if ((femaleMutationModel != "Equal"        & femaleMutationModel != "equal" &
         femaleMutationModel != "Proportional" & femaleMutationModel != "proportional" & 
         femaleMutationModel != "Stepwise"     & femaleMutationModel != "stepwise" & 
         femaleMutationModel != "Custom"       & femaleMutationModel != "custom") | 
        (maleMutationModel   != "Equal"        & maleMutationModel   != "equal" &
         maleMutationModel   != "Proportional" & maleMutationModel   != "proportional" &
         maleMutationModel   != "Stepwise"     & maleMutationModel   != "stepwise" & 
         maleMutationModel   != "Custom"       & maleMutationModel   != "custom")) 
        stop("ERROR: Mutation models must be either 'Equal', 'Proportional', 'Stepwise', or 'Custom'.")
    if (femaleMutationRate < 0 | maleMutationRate < 0 | femaleMutationRate > 1 | maleMutationRate > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRate2 < 0 | maleMutationRate2 < 0 | femaleMutationRate2 > 1 | maleMutationRate2 > 1) 
        stop("ERROR: Mutation rates must be between 0 and 1.")
    if (femaleMutationRange <= 0 | maleMutationRange <= 0) 
        stop("ERROR: Mutation ranges must be positive.")

    if (femaleMutationModel == "Custom" | femaleMutationModel == "custom") {
        if (missing(femaleMutationMatrix)) 
            stop("When the female mutation model is 'Custom' the female mutation matrix must be specified.")
        if (!is.matrix(femaleMutationMatrix) | dim(femaleMutationMatrix)[1] != 
            nAlleles | dim(femaleMutationMatrix)[2] != nAlleles) 
            stop("The female mutation matrix must be of a dimension corresponding to the vector of frequencies.")
        if (any(as.vector(femaleMutationMatrix) < 0)) 
            stop("The female mutation matrix cannot have negative entries.")
        if (any(round(apply(femaleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the female mutation matrix must sum to 1.")
        femaleMutationType <- "A 'Custom' specified mutation matrix"
    }
    else 
    {
        femaleMutationMatrix <- matrix(0, nAlleles, nAlleles)
        diag(femaleMutationMatrix) <- 1
        if (femaleMutationRate == 0 ) 
            femaleMutationType <- "No mutations"
        else if (femaleMutationModel == "Equal" | femaleMutationModel == "equal") {
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                femaleMutationMatrix[i, j] <- 1 - femaleMutationRate
            else femaleMutationMatrix[i, j] <- femaleMutationRate/(nAll - 1)
            femaleMutationType <- paste("An 'Equal' mutation model with mutation rate", 
                femaleMutationRate)
        }
        else if (femaleMutationModel == "Proportional" | femaleMutationModel == "proportional") {
            sumfreq <- sum(freq)
            frq     <- freq/sumfreq
            alpha <- femaleMutationRate/sumfreq/sum(frq * (1 - frq))
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                femaleMutationMatrix[i, j] <- 1 - alpha + alpha * frq[j]
            else femaleMutationMatrix[i, j] <- alpha * frq[j]
            femaleMutationType <- paste("A 'Proportional' mutation model with expected mutation rate", 
                femaleMutationRate)
        }
        else if (femaleMutationModel == "Stepwise" | femaleMutationModel == "stepwise") {
            numfreq <- as.numeric(names(freq))
            if (any(is.na(numfreq)))
                stop("The 'Stepwise' mutation model requires all non-silent alleles to have numerical names.")
            if (any(round(numfreq, 1)!=numfreq))
                stop("Microvariants must be named as a decimal number with one decimal.")
            microgroup <- (numfreq - round(numfreq))*10
            for (i in 1:nAll) {
                microcompats <- (microgroup == microgroup[i])
                for (j in 1:nAll) {
                    if (i==j) {
                        if (all(microcompats)) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate
                        else if (sum(microcompats)==1) femaleMutationMatrix[i,j] <- 1 - femaleMutationRate2
                        else femaleMutationMatrix[i,j] <- 1 - femaleMutationRate - femaleMutationRate2
                    } else if (microcompats[j])
                        femaleMutationMatrix[i,j] <- femaleMutationRange^abs(numfreq[i]-numfreq[j])
                    else
                        femaleMutationMatrix[i,j] <- femaleMutationRate2/(nAll-sum(microcompats))
                }
                microcompats[i] <- FALSE
                if (any(microcompats))
                    femaleMutationMatrix[i,microcompats] <- femaleMutationMatrix[i,microcompats]/
						            sum(femaleMutationMatrix[i,microcompats])*
						            femaleMutationRate
            }
            if (all(microgroup==0))
                femaleMutationType <- paste("A 'Stepwise' mutation model with mutation rate", 
                femaleMutationRate, "and mutation range", femaleMutationRange)
            else 
                femaleMutationType <- paste("A 'Stepwise' mutation model with mutation rate ", 
                femaleMutationRate, ", range ", femaleMutationRange, ", and fractional mutation rate ", 
		femaleMutationRate2, sep="")
        }
    }

    if (maleMutationModel == "Custom" | maleMutationModel == "custom") {
        if (missing(maleMutationMatrix)) 
            stop("When the male mutation model is 'Custom' the male mutation matrix must be specified.")
        if (!is.matrix(maleMutationMatrix) | dim(maleMutationMatrix)[1] != 
            nAlleles | dim(maleMutationMatrix)[2] != nAlleles) 
            stop("The male mutation matrix must be of a dimension corresponding to the vector of frequencies.")
        if (any(as.vector(maleMutationMatrix) < 0)) 
            stop("The male mutation matrix cannot have negative entries.")
        if (any(round(apply(maleMutationMatrix, 1, sum), 6) != 1)) 
            stop("The rows in the male mutation matrix must sum to 1.")
        maleMutationType <- "A 'Custom' specified mutation matrix"
    }
    else 
    {
        maleMutationMatrix <- matrix(0, nAlleles, nAlleles)
        diag(maleMutationMatrix) <- 1
        if (maleMutationRate == 0 ) 
            maleMutationType <- "No mutations"
        else if (maleMutationModel == "Equal" | maleMutationModel == "equal") {
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                maleMutationMatrix[i, j] <- 1 - maleMutationRate
            else maleMutationMatrix[i, j] <- maleMutationRate/(nAll - 1)
            maleMutationType <- paste("An 'Equal' mutation model with mutation rate", 
                maleMutationRate)
        }
        else if (maleMutationModel == "Proportional" | maleMutationModel == "proportional") {
            sumfreq <- sum(freq)
            frq     <- freq/sumfreq
            alpha <- maleMutationRate/sumfreq/sum(frq * (1 - frq))
            for (i in 1:nAll) for (j in 1:nAll) if (i == j) 
                maleMutationMatrix[i, j] <- 1 - alpha + alpha * frq[j]
            else maleMutationMatrix[i, j] <- alpha * frq[j]
            maleMutationType <- paste("A 'Proportional' mutation model with expected mutation rate", 
                maleMutationRate)
        }
        else if (maleMutationModel == "Stepwise" | maleMutationModel == "stepwise") {
            numfreq <- as.numeric(names(freq))
            if (any(is.na(numfreq)))
                stop("The 'Stepwise' mutation model requires all non-silent alleles to have numerical names.")
            if (any(round(numfreq, 1)!=numfreq))
                stop("Microvariants must be named as a decimal number with one decimal.")
            microgroup <- (numfreq - round(numfreq))*10
            for (i in 1:nAll) {
                microcompats <- (microgroup == microgroup[i])
                for (j in 1:nAll) {
                    if (i==j) {
                        if (all(microcompats)) maleMutationMatrix[i,j] <- 1 - maleMutationRate
                        else if (sum(microcompats)==1) maleMutationMatrix[i,j] <- 1 - maleMutationRate2
                        else maleMutationMatrix[i,j] <- 1 - maleMutationRate - maleMutationRate2
                    } else if (microcompats[j])
                        maleMutationMatrix[i,j] <- maleMutationRange^abs(numfreq[i]-numfreq[j])
                    else
                        maleMutationMatrix[i,j] <- maleMutationRate2/(nAll-sum(microcompats))
                }
                microcompats[i] <- FALSE
                if (any(microcompats))
                    maleMutationMatrix[i,microcompats] <- maleMutationMatrix[i,microcompats]/
						          sum(maleMutationMatrix[i,microcompats])*
						          maleMutationRate
            }
            if (all(microgroup==0))
                maleMutationType <- paste("A 'Stepwise' mutation model with mutation rate", 
                maleMutationRate, "and mutation range", maleMutationRange)
            else 
                maleMutationType <- paste("A 'Stepwise' mutation model with mutation rate ", 
                maleMutationRate, ", range ", maleMutationRange, ", and fractional mutation rate ", 
		maleMutationRate2, sep="")
        }
    }

    # do the stabilization 
    Stabilization <- toupper(Stabilization)
    if (Stabilization != "NONE" & Stabilization != "DP" & 
        Stabilization != "RM" & Stabilization != "PM") 
        stop("The Stabilization parameter must be 'None', 'DP', 'RM', or 'PM'.")
    if (Stabilization != "NONE") {
        if (hasSilentAllele & all(maleMutationMatrix[1:nAll,nAlleles]==0) & 
                              all(maleMutationMatrix[nAlleles,1:nAll]==0) & 
                              all(femaleMutationMatrix[1:nAll,nAlleles]==0) & 
                              all(femaleMutationMatrix[nAlleles,1:nAll]==0)) {
            res <- stabilize(femaleMutationMatrix[1:nAll,1:nAll], frequencies[1:nAll]/sum(frequencies[1:nAll]), 
                             Stabilization, MaxStabilizedMutrate) 
            if (res$error!="") {
                print(paste("WARNING:", res$error))
                print("WARNING: Female mutation matrix not stabilized.")
            } else {
                print(paste("Female mutation matrix f ratio:", res$fratio))
                print(paste("Female mutation matrix max specific mutation rate:", 1-res$mindiag))
                femaleMutationMatrix[1:nAll,1:nAll] <- res$stabilized
            }
            res <- stabilize(maleMutationMatrix[1:nAll,1:nAll], frequencies[1:nAll]/sum(frequencies[1:nAll]), 
                             Stabilization, MaxStabilizedMutrate) 
            if (res$error!="") {
                print(paste("WARNING:", res$error))
                print("WARNING: Male mutation matrix not stabilized.")
            } else {
                print(paste("Male mutation matrix f ratio:", res$fratio))
                print(paste("Male mutation matrix max specific mutation rate:", 1-res$mindiag))
                maleMutationMatrix[1:nAll,1:nAll] <- res$stabilized
            }
        } else {
            res <- stabilize(femaleMutationMatrix, frequencies, 
                             Stabilization, MaxStabilizedMutrate) 
            if (res$error!="") {
                print(paste("WARNING:", res$error))
                print("WARNING: Female mutation matrix not stabilized.")
            } else {
                print(paste("Female mutation matrix f ratio:", res$fratio))
                print(paste("Female mutation matrix max specific mutation rate:", 1-res$mindiag))
                femaleMutationMatrix <- res$stabilized
            }
            res <- stabilize(maleMutationMatrix, frequencies, 
                             Stabilization, MaxStabilizedMutrate) 
            if (res$error!="") {
                print(paste("WARNING:", res$error))
                print("WARNING: Male mutation matrix not stabilized.")
            } else {
                print(paste("Male mutation matrix f ratio:", res$fratio))
                print(paste("Male mutation matrix max specific mutation rate:", 1-res$mindiag))
                maleMutationMatrix <- res$stabilized
            }
        }
    }
    rownames(femaleMutationMatrix) <- names(frequencies)
    colnames(femaleMutationMatrix) <- names(frequencies)
    rownames(maleMutationMatrix) <- names(frequencies)
    colnames(maleMutationMatrix) <- names(frequencies)

    simpleMutationMatrices <- TRUE
    for (j in 1:nAlleles) {
        v <- femaleMutationMatrix[-j, j]
        if (any(round(v - v[1], 6) != 0)) 
            simpleMutationMatrices <- FALSE
    }
    for (j in 1:nAlleles) {
        v <- maleMutationMatrix[-j, j]
        if (any(round(v - v[1], 6) != 0)) 
            simpleMutationMatrices <- FALSE
    }
    result <- list(locusname = name, alleles = frequencies, femaleMutationType = femaleMutationType, 
        femaleMutationMatrix = femaleMutationMatrix, maleMutationType = maleMutationType, 
        maleMutationMatrix = maleMutationMatrix, simpleMutationMatrices = simpleMutationMatrices, 
	Stabilization = Stabilization)
    class(result) <- "FamiliasLocus"
    result
}

stabilize = function(M,pe,stabilizationMethod="DP",t=1){
  #library('Rsolnp')
  R = 1-sum(diag(M)*pe)
  n = dim(M)[1]
  m = n^2
  xM = as.vector(M)
  tol = 1e-10
  if (all(M==diag(n))) 
     return(list(stabilized = M, fratio=1, mindiag=min(diag(M)), error=""))
  if (any(M==0))
     return(list(stabilized = M, fratio=1, mindiag=min(diag(M)), error="Cannot stabilize non-identity matrices with zero entries."))
  if (stabilizationMethod == "DP"){
    if (2*max(pe*(1-diag(M))) > sum(pe*(1-diag(M)))){
      return(list(stabilized=M,fratio=1,mindiag=min(diag(M)),error="DP stabilization doesn't exist."))
    }
    if (n>30)
      print("NOTE: Stabilization comuptations may take long time for large systems.")
    pnew = pKCompute(pe,diag(M))
    P0 = (diag(n)-diag(diag(M)))%*%diag(1/(1-pnew))%*%(outer(rep(1,n),pnew)-diag(pnew)) + diag(diag(M))
    P = theOpting(M,P0,pe)
    if (max(abs(pe%*%P-pe))>tol){
      return(list(error="The proposed stabilization doesn't have the desired stationary distribution."))
    }
    if (max(abs(P%*%rep(1,n)-rep(1,n)))>tol){
      return(list(error="The proposed stabilization isn't a proper mutation matrix."))
    }
    if (min(P)<0){
      return(list(error="The proposed stabilization has negative elements."))
    } else if (min(P)<tol){
      print("WARNING: The proposed stabilization has very small elements.")
    }
    fratio = max(max(P/M),max(M/P))
    minS = min(diag(P))
    return(list(stabilized=P,fratio=fratio,mindiag=minS,error=""))
  } else if (stabilizationMethod == "RM"){
    if (t<R+tol | t>1) 
      return(list(error="MaxStabilizedMutrate parameter out of bounds."))
    if (n>30)
      print("NOTE: Stabilization comuptations may take long time for large systems.")
    C = matrix(0,2*n,m)
    for (i in 1:n){
      C[seq(1,n),seq(n*(i-1)+1,n*i)] = diag(n)
      C[n+i,seq(n*(i-1)+1,n*i)] = pe
      C[2*n,seq(n*(n-1)+1,m)] = 0
      C[2*n,n*(i-1)+i] = pe[i]
    }
    b = c(rep(1,n),pe[-n],1-R)
    xP0 = solnp(pars=xM,
                fun=function(x) max(sum(x/xM),sum(xM/x)),
                eqfun = function(x) C%*%x,
                eqB=b,
                LB=as.vector((1-t)*diag(n)),
                UB=rep(1,m),
                control=list("trace"=FALSE))
    xP = solnp(pars=xP0$pars,
               fun=function(x) max(abs(log(x)-log(xM))),
               eqfun=function(x) C%*%x,
               eqB=b,
               LB=as.vector((1-t)*diag(n)),
               UB=rep(1,m),
               control=list("trace"=FALSE))
    P0 = matrix(xP$pars,n,n)
    P = theOptingRM(M,P0,pe,t)
    if (max(abs(pe%*%P-pe))>tol){
      warning("The proposed stabilization doesn't have the desired stationary distribution.")
    }
    if (max(abs(P%*%rep(1,n)-rep(1,n)))>tol){
      warning("The proposed stabilization isn't a proper mutation matrix.")
    }
    if (min(P)<0){
      stop("The proposed stabilization has negative elements.")
    } else if (min(P)<tol){
      warning("The proposed stabilization has very small elements.")
    }
    fratio = max(max(P/M),max(M/P))
    minS = min(diag(P))
    return(list(stabilized=P,fratio=fratio,mindiag=minS,error=""))
  } else if (stabilizationMethod == "PM"){
    # No optimization needed here, the stabilization is unique (if it exists).
    X = t(M)-diag(n)
      X[n,] = rep(1,n)
      v = solve(X,c(rep(0,n-1),1))
      d = R*v/((1-sum(diag(M)*v))*pe)
      if (any(d*(1-diag(M))>t)){
        return(list(statbilized=M,fratio=1,mindiag=min(diag(M)),error="PM stabilization doesn't exist."))
      }
      P = diag(d)%*%(M-diag(n)) + diag(n)
    fratio = max(max(P/M),max(M/P))
    minS = min(diag(P))
    return(list(stabilized=P,fratio=fratio,mindiag=minS,error=""))
  } else {
    stop("Stabilization method must be either \"DP\",\"RM\" or \"PM\".")
  }
}

pKCompute = function(paj,d){
  # Input paj is the desired stationary distribution and d is a vector
  # with the desired diagonal elements.
  # Output should be a probability vector p satisfying
  # p*(I-D(p)) = K*pie*(I-D(d)).

  paj = as.vector(paj)
  n = length(paj)
  
  w = paj*(1-d)
  if (2*max(w) > sum(w)){
    stop('Task impossible with this input')
  }
  
  ord = rev(order(w))
  Perm = matrix(0,n,n)
  for (i in 1:n){
    Perm[i,ord[i]] = 1
  }
  w = Perm%*%w
  
  eps = 1e-5
  p = rep(0,n)
  Kmax = 0.25/w[1]
  ffun = function (x) 0.5*sum(1-sqrt(1-4*x*w))-1
  gfun = function (x) ffun(x)+sqrt(1-4*x*w[1])
  if (ffun(Kmax)>0){
    K = uniroot(ffun,c(0,Kmax))$root
    p = 0.5*(1-sqrt(1-4*K*w))
  } else {
    while (gfun(eps)<0){
      eps = eps/2
    }
    K = uniroot(gfun,c(eps,Kmax))$root
    p[1] = 0.5*(1+sqrt(1-4*K*w[1]))
    p[-1] = 0.5*(1-sqrt(1-4*K*w[-1]))
  }
  
  p = t(Perm)%*%p
  p = as.vector(p)
  return(p)
  
}

fratio = function(A,B){
  max(max(A/B),max(B/A))
}

SgivenMopt = function(M,S0,p){
  
  n1 = dim(M)[1]
  n2 = dim(M)[2]
  n3 = length(p)
  if ( (n1!=n2) || (n1!=n3) ){
    stop('matrix and/or vector dimensions don\'t agree')
  }
  
  n = n1
  tolerance = 1e-10
  factorTolerance = 1
  steg = 1e-2
  S = S0
  best = fratio(M,S)
  change = FALSE
  if (max(M/S0) >= max(S0/M)){
    ii = which(M/S0 == best) %% n
    jj = ceiling(which(M/S0 == best)/n)
  } else {
    ii = which(S0/M == best) %% n
    jj = ceiling(which(S0/M == best)/n)
    steg = -steg
  }
  for (ki in 1:length(ii)){
    i = ii[ki]
    if (i==0) i=n
    for (kj in 1:length(jj)){
      j = jj[kj]
      if (j==0) j=n
      if (j==i) next
      for (i1 in 1:n){
        if ( (i1 == i) || (i1 == j) ) next
        for (j1 in 1:n){
          if ( (j1 == i) || (j1 == i1) || (j1 == j) ) next
          while (TRUE){
            Stemp = S
            Stemp[i,j] = (1+steg)*S[i,j]
            Stemp[i,j1] = S[i,j] + S[i,j1] - Stemp[i,j]
            Stemp[i1,j] = (S[i,j]*p[i] + S[i1,j]*p[i1] - Stemp[i,j]*p[i])/p[i1]
            Stemp[i1,j1] = S[i1,j] + S[i1,j1] - Stemp[i1,j]
            f1 = fratio(M[c(i,i1),c(j,j1)],Stemp[c(i,i1),c(j,j1)])
            f2 = fratio(M[c(i,i1),c(j,j1)],S[c(i,i1),c(j,j1)])
            if ( (min(Stemp)>0) && (f1-f2 < 0) ){
              change = TRUE
              S = Stemp
              best = fratio(M,S)
            } else {
              break
            }
          }
        }
      }
    }
  }
  ut = list(newS=S,changed=change)
  return(ut)
}

SgivenMoptRM = function(M,S0,p,t){
  
  n1 = dim(M)[1]
  n2 = dim(M)[2]
  n3 = length(p)
  if ( (n1!=n2) || (n1!=n3) ){
    stop('matrix and/or vector dimensions don\'t agree')
  }
  
  n = n1
  tolerance = 1e-10
  factorTolerance = 1
  steg = 1e-3
  S = S0
  best = fratio(M,S)
  change = FALSE
  if (max(M/S0) >= max(S0/M)){
    ii = which(M/S0 == best) %% n
    jj = ceiling(which(M/S0 == best)/n)
  } else {
    ii = which(S0/M == best) %% n
    jj = ceiling(which(S0/M == best)/n)
    steg = -steg
  }
  for (ki in 1:length(ii)){
    i = ii[ki]
    if (i==0) i=n
    for (kj in 1:length(jj)){
      j = jj[kj]
      if (j==0) j=n
      for (i1 in 1:n){
        if (i1 == i) next
        for (j1 in 1:n){
          if (j1 == j) next
          while (TRUE){
            Stemp = S
            Stemp[i,j] = (1+steg)*S[i,j]
            Stemp[i,j1] = S[i,j] + S[i,j1] - Stemp[i,j]
            Stemp[i1,j] = (S[i,j]*p[i] + S[i1,j]*p[i1] - Stemp[i,j]*p[i])/p[i1]
            Stemp[i1,j1] = S[i1,j] + S[i1,j1] - Stemp[i1,j]
            f1 = fratio(M[c(i,i1),c(j,j1)],Stemp[c(i,i1),c(j,j1)])
            f2 = fratio(M[c(i,i1),c(j,j1)],S[c(i,i1),c(j,j1)])
            if ( (min(diag(Stemp))>1-t) && (min(Stemp)>0) && (f1-f2 < 0) ){
              change = TRUE
              S = Stemp
              best = fratio(M,S)
            } else {
              break
            }
          }
        }
      }
    }
  }
  return(list(newS=S,changed=change))
}

theOpting = function(M,S,p){
  
  go = TRUE
  while(go){
    aRun = SgivenMopt(M,S,p)
    S = aRun$newS
    go = aRun$changed
    #print(fratio(M,S))
  }
  return(S)
}

theOptingRM = function(M,S,p,t){
  
  go = TRUE
  while(go){
    aRun = SgivenMoptRM(M,S,p,t)
    S = aRun$newS
    go = aRun$changed
    #print(fratio(M,S))
  }
  return(S)
}

