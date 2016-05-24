## This code implements the BAC method with a prior on (\alpha^X, \alpha^Y) that 
##requires the main effect term to be included in the outcome model whenever an interaction is included*/

### Some Suggested Input Arguments
## omega = -1.0;     // 0< omega <= infinity. we use omega = -1.0 to represent omega = infinity.
## num_its = 10000; // number of iterations in the MCMC algorithm
## burnM = 10000; // burin-in iterations for sampling (\alpha^X, \alpha^Y)
## burnB = 10000; //burn-in iterations for sampling \beta
## thin = 250;  //thinning parameter for sampling \beta

bac <- function(data, exposure, outcome, confounders, interactors, familyX, familyY, omega=Inf, num_its, burnM, burnB, thin, population=NULL){
    interactions = match(interactors, confounders)
    if (length(interactions)==0) {interactions = NULL}
    nits = num_its + burnM
    predictorsY = get.predictorsY(exposure, confounders, interactions)
    predictorsX = confounders
    n = nrow(data)
    m = length(confounders)
    mm = m+length(interactions)
    if (length(population)==0) {population = rep(T, n)}
    MX = matrix(nrow=nits, ncol=m)
    MY = matrix(nrow=nits, ncol=mm+1)
    M0X = rep(1, m)
    M0Y = rep(1, mm+1)
    MX[1,] = M0X
    MY[1,] = M0Y
    a0X = calclogpost(exposure, predictorsX[M0X==1], familyX, data)
    a0Y = calclogpost(outcome, predictorsY[M0Y==1], familyY, data)

    ptime = proc.time()[3]
    ## MCMC algorithm
    for (k in 2:nits){
        
        M0X = MX[k-1,]
        M1X = M0X
        if (omega == Inf){
            #update exposure model
            candidate = (1:m)[MY[k-1,1:m]==1]
            if (length(candidate)>0){
                if (length(candidate)>1){
  	                j = sample(candidate, 1)
                } else{
                    j = candidate
                }
                M1X[j] = 1 - M1X[j]
                tuning = 1.0
                temp = update_alpha(exposure, predictorsX, M0X, M1X, a0X, tuning, familyX, data)
                M0X = temp$M; a0X = temp$a0
            } else{
	            j = sample(1:m,1)
                M1X[j] = 1 - M1X[j]
                if (MY[k-1,j] == 1){
                    tuning = 1.0
                } else{
	                tuning = omega^(M0X[j] - M1X[j])
                }
                temp = update_alpha(exposure, predictorsX, M0X, M1X, a0X, tuning, familyX, data)
                M0X = temp$M; a0X = temp$a0
            }
            MX[k,] = M0X
            #update outcome model
            M0Y = MY[k-1,]
            M1Y = M0Y
            ind = rep(0, m)
            ind[interactions] = M0Y[-c(1:m, mm+1)]
        }
        if (omega == Inf){
            candidate = c((1:m)[(ind==0) & (MX[k,]==0)], (1:mm)[-(1:m)][M0Y[interactions]==1])
            if (length(candidate)>0){
                if (length(candidate)>1){
  	                j = sample(candidate, 1)
                } else{
                    j = candidate
                }
                M1Y[j] = 1 - M1Y[j]
                tuning = 1.0
                temp = update_alpha(outcome, predictorsY, M0Y, M1Y, a0Y, tuning, familyY, data)
                M0Y = temp$M; a0Y = temp$a0
            } else{
                candidate = c((1:m)[(ind==0)], (1:mm)[-(1:m)])
                if (length(candidate)>0){
                    if (length(candidate)>1){
  	                    j = sample(candidate, 1)
                    } else{
                        j = candidate
                    }
                M1Y[j] = 1 - M1Y[j]
                tuning = 1.0
                if (j<=m){
                    if (MX[k,j] == 1){
	                    tuning = omega^(M1Y[j] - M0Y[j])
                    }
                }
                temp = update_alpha(outcome, predictorsY, M0Y, M1Y, a0Y, tuning, familyY, data)
                M0Y = temp$M; a0Y = temp$a0
                }
            }
            MY[k,] = M0Y
        }
    }
    MX = MX[-(1:burnM),]
    MY = MY[-(1:burnM),]
    #calculate ACE
    tmp = calc.ACE(population, exposure, outcome, predictorsY, interactions, confounders, MY, familyY, thin, burnB, data)
    ACE.boot = tmp$ACE.boot
    #count the time to run the program
    ptime = round(proc.time()[3] - ptime, digits=2)
    print(paste("It took you ", ptime, " seconds to run BAC",sep=""))
    
    #save results
    para = list(familyX=familyX, familyY=familyY, thin=thin, omega=omega, num_its=num_its, burnM=burnM, burnB=burnB, population=population)
    models = find.models(MY)
    result = list(data = data, MX = MX, MY = MY, models=models, exposure=exposure, outcome=outcome, confounders=confounders, interactions=interactions, predictorsY=predictorsY, ACE=ACE.boot, para=para)
    class(result)='bacr'
    return(result)
}
