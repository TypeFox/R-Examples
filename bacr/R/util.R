get.predictorsY <- function(exposure, confounders, interactions){
    if (!is.null(interactions)){
      predictorsY = c(confounders, paste(exposure, confounders[interactions], sep=':'), exposure)
    } else {
      predictorsY = c(confounders, exposure)
    }
    return(predictorsY)
}


calclogpost <- function(Yvar, Xvar, family, data){
    if (length(Xvar)>0){
        formu = paste(Yvar, " ~ ", paste(Xvar, collapse="+"), sep="")
    } else {
        formu = paste(Yvar, " ~ 1", sep="")
    }
    fit = glm(formu, data = data, family=family)
    return((-1)*BIC(fit)/2)
}

update_alpha<- function(Yvar, Xvar, M0, M1, a0, tuning, family, data){
    a1 = calclogpost(Yvar, Xvar[M1==1], family, data)
    BF = tuning * exp(a1 - a0)
    flag = 1
    if (BF<1){
        if (BF < runif(1)){
            flag = 0
        }
    }
    if (flag==1){
        return(list(M = M1, a0 = a1))
    } else{
        return(list(M = M0, a0 = a0))
    } 
}

find.models <- function(MY){
    bb = apply(MY, 1, paste, collapse=" ")
    count = table(bb)
    models = matrix(as.numeric((sapply(strsplit(names(count), split=" "), unlist))), byrow=T, nrow=length(count))
    return(list(models=models, counts=count))
}

GetACE <- function(nsample, model, population, exposure, outcome, predictorsY, interactions, confounders, familyY, thin, burnB, data){
    startvalue = NA
    formu = paste(outcome, " ~ ", paste(predictorsY[model==1], collapse="+"), sep="")
    kk = coef(glm(formu, data=data, family=familyY))
    model[(which(model==1))[which(is.na(kk[(names(kk) != "(Intercept)") & (names(kk) != exposure)]))]] = 0
    formu = paste(outcome, " ~ ", paste(predictorsY[model==1], collapse="+"), sep="")
    b0 = 0; B0 = 0

    if (familyY=="gaussian"){
        para = data.matrix(MCMCregress(formu, data=data, mcmc=nsample*thin, b0=b0, B0=B0, burnin=burnB, thin=thin)[,1:(sum(model==1)+1)])  ### the last col of MCMCreg is sigma2
        if (nsample>1){
            para = t(para)
        }
    }
    if (familyY=="binomial"){
        para = t(data.matrix(MCMClogit(formu, data=data, mcmc=nsample*thin, b0=b0, B0=B0, burnin=burnB, thin=thin, beta.start = startvalue, tune=0.4)))
    }
    if (familyY=="poisson"){
        para = t(data.matrix(MCMCpoisson(formu, data=data, mcmc=nsample*thin, b0=b0, B0=B0, burnin=burnB, thin=thin, beta.start = startvalue, tune=0.5))) ### MCMCpoisson may re-order variables
    }
    data1 = data.matrix(data)
    n = nrow(data)
    m = length(confounders)
    mm = m+length(interactions)
    if (!is.null(interactions)){
        data_c = cbind(rep(1,n), data1[,confounders[model[1:m]==1]], rep(0,n), rep(0,n)*data1[,confounders[interactions[model[(m+1):mm]==1]]])
        data_t = cbind(rep(1,n), data1[,confounders[model[1:m]==1]], rep(1,n), rep(1,n)*data1[,confounders[interactions[model[(m+1):mm]==1]]])
    } else{
        data_c = cbind(rep(1,n), data1[,confounders[model[1:m]==1]], rep(0,n))
        data_t = cbind(rep(1,n), data1[,confounders[model[1:m]==1]], rep(1,n))
    }

    data_t = data_t[population,]; data_c = data_c[population,]
    predict1 = data_t%*%para; predict0 = data_c%*%para
    if (familyY=="binomial"){
        predict1 = 1/(1+exp((-1)*predict1))   
        predict0 = 1/(1+exp((-1)*predict0))  
    }
    if (familyY=="poisson"){
        predict1 = exp(predict1)
        predict0 = exp(predict0)
    }
    weights = t(rdirichlet(nsample, rep((-1)+2,sum(population)))) * sum(population)
    ACE.boot = colMeans(weights * (predict1 - predict0))
    return(list(ACE.boot=ACE.boot, para=para))
}

calc.ACE <- function(population, exposure, outcome, predictorsY, interactions, confounders, MY, familyY, thin, burnB, data){
    temp = find.models(MY)
    ACE.boot = NULL
    for(h in 1:length(temp$count)){
        model = temp$models[h,]
        vv = GetACE(temp$counts[h], model, population, exposure, outcome, predictorsY, interactions, confounders, familyY, thin, burnB, data)
        ACE.boot = c(ACE.boot, vv$ACE.boot)
    }
    return(list(ACE.boot=ACE.boot))
}


