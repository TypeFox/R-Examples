eblup.mse.f.wrap <-
    function(domain.data = NULL
             , lme.obj = NULL
             , debug=F
             , ...)
{
        #get the name of the domain ID. This is the random variable in lme and thus the second column of the fitted object in lme.obj
    domain.ID.name <- colnames(lme.obj$fitted)[2]
    if(debug){print(domain.ID.name)}
        #create new $domain.ID if this is not already the name
    if(domain.ID.name != "domain.ID"){domain.data$domain.ID <- domain.data[,domain.ID.name]}

        #get sample level data from the lme object
    sample.data <- lme.obj$data
        #fit lm.obj
    lm.obj <- lm(formula(lme.obj), data=sample.data)
    if(debug){print(summary(lm.obj))}

        #vector of unique domain IDs (can be different to domain.ID further down)
    domain.IDs <- domain.data$domain.ID

        #get an overview of the domains
            #variables of the model
    variabs <- rownames(attr(terms(lme.obj), which="factors"))
    if(debug){print(variabs)}

        #get response columns name
    response.nam <- variabs[attr(terms(lme.obj), which="response")]

            #mean of the response and predictor variables from the sample. For the response this is the sample mean estimator.
    mean.response.predicors <-####here it needs to stay sample.data[,domain.ID.name]!!
        aggregate(sample.data[,variabs]#[,c("biomass.ha", "mean.canopy.ht")]
                  , by=list(domain.ID=c(sample.data[,domain.ID.name])), mean)
                #add "sample" in col names
    names(mean.response.predicors)[-which(names(mean.response.predicors)==domain.ID.name)] <-
        paste(names(mean.response.predicors)
              [-which(names(mean.response.predicors)==domain.ID.name)], ".sample.mean", sep="")

            #number of samples within the domains
    n.samples <- aggregate(cbind(n.i.sample=sample.data[,variabs[1]])#[,c("biomass.ha", "mean.canopy.ht")]
                                , by=list(domain.ID=c(sample.data[,domain.ID.name])
                                  ), length)

        #mean lm residual -- needed for GREG
    resids.tmp <- aggregate(cbind(mean.resid.lm=resid(lm.obj))
                            , by=list(domain.ID=sample.data[,domain.ID.name])
                            , mean)

        #mean lme residual -- needed for EBLUP.var
    resids.tmp$mean.resid.lme <- aggregate(resid(lme.obj, level=0)#
                                                 , by=list(domain.ID=sample.data[,domain.ID.name])
                                                 , mean)[,2]

        #synthetic estimate
    estimates.tmp <- data.frame(Synth=predict(lm.obj
                                      , newdata= domain.data))

        #EBLUP estimate
    estimates.tmp$EBLUP <- predict(lme.obj
                                         , newdata=domain.data
                                         , level=1)

#    estimates.tmp[,domain.ID.name] <- domain.data[,domain.ID.name]\
    estimates.tmp$domain.ID <- domain.data$domain.ID

        #add ".domain" in col names
    domain.data.tmp <- domain.data
    names(domain.data.tmp)[-which(names(domain.data.tmp)=="domain.ID")] <-
        paste(names(domain.data.tmp)[-which(names(domain.data.tmp)=="domain.ID")], ".domain", sep="")

    if(debug){print(domain.ID.name)}

            #combine population information of the domains
    tmp2 <- merge(domain.data.tmp, mean.response.predicors, by="domain.ID")
    if(debug){print(domain.ID.name)}
    tmp3 <- merge(tmp2, n.samples, by="domain.ID")
    tmp4 <- merge(tmp3, resids.tmp, by="domain.ID")
    overview.domains <- merge(tmp4, estimates.tmp, by="domain.ID")

    if(debug){print(domain.ID.name)}

        #GREG estimate
    overview.domains$GREG <- overview.domains$Synth +  overview.domains$mean.resid.lm

        #gamma
    overview.domains$gamma.i <- eblup.mse.f.gamma.i(lme.obj=lme.obj
                                              , n.i=overview.domains$n.i.sample)

        #variance of the SRS total estimate (devision by n_i for variance of the mean)
    sample.var.tmp <-
        aggregate(cbind(sample.var.tot=sample.data[,response.nam])
              , by=list(domain.ID=sample.data[,domain.ID.name]), var)

        #variance of the GREG total estimate
    GREG.var.tmp <-
        aggregate(cbind(GREG.var.tot=resid(lm.obj))
                  , by=list(domain.ID=sample.data[,domain.ID.name]), var)



            #add variances to the overview df
    overview.domains <- merge(overview.domains, sample.var.tmp, by="domain.ID")
    overview.domains <- merge(overview.domains, GREG.var.tmp, by="domain.ID")

        #variance of the mean estimates
    overview.domains$sample.var.mean <- overview.domains$sample.var.tot/overview.domains$n.i.sample
    overview.domains$GREG.var.mean <- overview.domains$GREG.var.tot/overview.domains$n.i.sample

        #add fake response column for domains to be able to use model.matrix
    domain.data[,response.nam] <- 0

        #variance of the EBLUP
            #compute the A.i matrices for all domains (only needed once)
    if(debug){print("compute the A.i matrices for all domains")}
            #vector of domain IDs from the overview
    domain.ID <- overview.domains$domain.ID
    if(debug){
        print("domain.IDs:")
        print(domain.ID)}
                #initialize the result vector
    a.i.mats <- vector(mode="list", length=length(domain.ID))
    for(i in 1:length(domain.ID)){
        if(debug){print(i)}
        a.i.mats[[i]] <- eblup.mse.f.c2.ai(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                                     , n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                     , lme.obj=lme.obj
                                     , X.i=model.matrix(formula(lme.obj)
                                       , sample.data[sample.data[,domain.ID.name]==domain.ID[i],])
             )
    }
                    #add all the matrices
    sum.A.i.mats <- Reduce("+", a.i.mats)

            #the assymptotic var-cov matrix
    asy.var.cov.mat <- eblup.mse.f.c3.asyvarcovarmat(n.i=overview.domains$n.i
                                               , lme.obj=lme.obj)

    if(debug){print("Calculate eblup variance components")}
            #put together the variance components
    result <- NULL
    for(i in 1:length(domain.ID)){
        if(debug){print(i)}
        #first comp
        mse.c1.tmp <- eblup.mse.f.c1(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                               , n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                               , lme.obj=lme.obj)
        #second comp
        mse.c2.tmp <- eblup.mse.f.c2(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                               , X.i=model.matrix(formula(lme.obj), sample.data[sample.data[,domain.ID.name]==domain.ID[i],])
                               , X.bar.i =
                                     t(model.matrix(formula(lme.obj)
                                                    , domain.data[domain.data$domain.ID==domain.ID[i],]))
                               , sum.A.i = sum.A.i.mats
                               )
        #third comp
        mse.c3.tmp <- eblup.mse.f.c3(n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                               , lme.obj=lme.obj
                               , asympt.var.covar=asy.var.cov.mat)
        #third star comp
        mse.c3.star.tmp <- eblup.mse.f.c3.star( n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                         , lme.obj=lme.obj
                                         , mean.resid.i=overview.domains$mean.resid.lme[overview.domains$domain.ID==domain.ID[i]]
                                         , asympt.var.covar=asy.var.cov.mat)
        #save result
        result <- rbind(result, data.frame(domain.ID=domain.ID[i]
                                           #, n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                           , c1=as.numeric(mse.c1.tmp), c2=as.numeric(mse.c2.tmp)
                                           , c3=as.numeric(mse.c3.tmp), c3star=as.numeric(mse.c3.star.tmp)))
    }
        #add components to the EBLUP variance
    result$EBLUP.var.1 <- result$c1 + result$c2 + 2* result$c3star
    result$EBLUP.var.2 <- result$c1 + result$c2 + result$c3 + result$c3star
    result.eblup.mse.comp <- result

        #add EBLUP variances
    overview.domains <- merge(overview.domains, result.eblup.mse.comp, by="domain.ID")
            #calculate standard errors
    tmp <- sqrt(overview.domains[,c("sample.var.mean", "GREG.var.mean", "EBLUP.var.1", "EBLUP.var.2")])
    names(tmp) <- c("sample.se", "GREG.se", "EBLUP.se.1", "EBLUP.se.2")#paste(names(tmp), ".se", sep="")
        #delete misleading columns
    overview.domains$sample.var.tot <- NULL
    overview.domains$GREG.var.tot <- NULL
        #rename some columns
    overview.domains <- cbind(overview.domains, tmp)
    return(overview.domains)
}
