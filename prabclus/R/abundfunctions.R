# simulation with regpop.sar
# p.nb=NULL: ignore neighbors
# groupvector is needed for teststat="groups", has to be a factor
# sarestimate=NULL: presence-absence
abundtest <- function (prabobj, teststat = "distratio", tuning = 0.25,
                       times = 1000, p.nb = NULL, 
                       prange = c(0, 1), nperp = 4, step = 0.1, step2 = 0.01, 
                       twostep = TRUE, species.fixed=TRUE, prab01=NULL,
                       groupvector=NULL,
                       sarestimate=prab.sarestimate(prabobj),
                       dist = prabobj$distance,
                       n.species = prabobj$n.species)
#    species.fixed = TRUE 
{
    if (is.null(prab01))
      prab01 <- prabinit(prabmatrix=toprab(prabobj),rows.are.species=FALSE,
                      distance="none",neighborhood=prabobj$nb)
    if (is.null(p.nb) & prabobj$spatial){        
        ac <- autoconst(prab01, twostep = twostep, prange = prange, 
            nperp = nperp, step1 = step, step2 = step2,
                        species.fixed = species.fixed)
        p.nb <- ac$pd
    }
    if (is.null(p.nb)) 
        p.nb <- 1    
    statres <- rep(0, times)
    if (teststat=="groups"){
      groupvector <- as.factor(groupvector)
      ng <- length(levels(groupvector))
      lg <- levels(groupvector)
      nsg <- numeric(0)
      for (i in 1:ng) nsg[i] <- sum(groupvector==lg[i])
      pa <- pb <- rep(1,ng)
      groupinfo <- list(lg=lg,ng=ng,nsg=nsg)
      statreslist <- list(overall=numeric(0),mean=numeric(0),
                          gr=matrix(0,nrow=ng,ncol=times))
    }
    else
      groupinfo <- NULL
#    print(groupinfo)
#    rmatrix <- rankmatrix(prabobj)
    for (i in 1:times) {
        cat("Simulation run ", i)
        if (is.null(sarestimate) || teststat == "inclusions")
          mat <- randpop.nb(neighbors=prabobj$nb,
                             p.nb = p.nb, n.species = prabobj$n.species, 
                             vector.species = prab01$regperspec,
                             species.fixed = species.fixed, 
                             pdf.regions =
                             prab01$specperreg/sum(prab01$specperreg),
                             count = FALSE)
         else
           mat <- regpop.sar(prabobj, prab01, sarestimate,
                          p.nb, count = FALSE)
        if (teststat != "inclusions"){
            if (dist == "jaccard") 
                distm <- jaccard(mat)
            if (dist == "kulczynski") 
                distm <- kulczynski(mat)
            if (dist == "qkulczynski") 
                distm <- qkulczynski(mat)
            if (dist == "logkulczynski") 
                distm <- qkulczynski(mat,log.distance=TRUE)
        }
        else statres[i] <- incmatrix(mat)$ninc
        if (teststat == "isovertice") {
            test <- homogen.test(distm, ne = tuning)
            statres[i] <- test$iv
        }
        if (teststat == "lcomponent") 
            statres[i] <- lcomponent(distm, ne = tuning)$lc
        if (teststat == "mean")
            statres[i] <- mean(as.dist(distm))
        if (teststat == "distratio") 
            statres[i] <- distratio(distm, prop = tuning)$dr
        if (teststat == "nn") 
            statres[i] <- nn(distm, ne = tuning)
        if (teststat == "groups"){
            slist <- specgroups(distm, groupvector, groupinfo)
#            print(slist)
            statreslist$overall[i] <- slist$overall
            statreslist$mean[i] <- mean(as.dist(distm))
            statreslist$gr[,i] <-  slist$gr
#            print(statreslist)
            cat(" statistics value=", statreslist$overall[i], "\n")
        }
        else
          cat(" statistics value=", statres[i], "\n")
    }
    regmat <- prabobj$prab
    if (teststat != "inclusions") {
      if (dist==prabobj$distance)
        distm <- prabobj$distmat
      else{
        if (dist == "jaccard") 
            distm <- jaccard(regmat)
        if (dist == "kulczynski") 
            distm <- kulczynski(regmat)
        if (dist == "qkulczynski") 
            distm <- qkulczynski(regmat)
        if (dist == "logkulczynski") 
            distm <- qkulczynski(regmat,log.distance=TRUE)
      }
    }
    else {
        regmat <- prab01$prab
        test <- incmatrix(regmat)$ninc
        p.above <- (1 + sum(statres >= test))/(1 + times)
        p.below <- (1 + sum(statres <= test))/(1 + times)
        datac <- test
        tuning <- NA
    }
    if (teststat == "mean"){
      test <- mean(as.dist(distm))
      p.above <- (1 + sum(statres >= test))/(1 + times)
      p.below <- (1 + sum(statres <= test))/(1 + times)
      datac <- test
      tuning <- NA
    }

    if (teststat == "isovertice") {
        test <- homogen.test(distm, ne = tuning)
        p.above <- (1 + sum(statres >= test$iv))/(1 + times)
        p.below <- (1 + sum(statres <= test$iv))/(1 + times)
        pb <- min(p.above, p.below) * 2
        p.above <- max(p.above, p.below)
        p.below <- pb
        datac <- test$iv
        tuning <- test$ne
    }
    if (teststat == "lcomponent") {
        test <- lcomponent(distm, ne = tuning)
        p.above <- (1 + sum(statres >= test$lc))/(1 + times)
        p.below <- (1 + sum(statres <= test$lc))/(1 + times)
        datac <- test$lc
        tuning <- test$ne
    }
    if (teststat == "nn") {
        test <- nn(distm, ne = tuning)
        p.above <- (1 + sum(statres >= test))/(1 + times)
        p.below <- (1 + sum(statres <= test))/(1 + times)
        datac <- test
    }
    if (teststat == "distratio") {
        test <- distratio(distm, prop = tuning)
        p.above <- (1 + sum(statres >= test$dr))/(1 + times)
        p.below <- (1 + sum(statres <= test$dr))/(1 + times)
        datac <- test$dr
        tuning <- test$prop
    }
    if (teststat=="groups"){
        test <- specgroups(distm, groupvector, groupinfo)
        testm <- mean(as.dist(distm))
        p.above <- (1 + sum(statreslist$overall >= test$overall))/(1 + times)
        p.below <- (1 + sum(statreslist$overall <= test$overall))/(1 + times)
        p.m.above <- (1 + sum(statreslist$mean >= testm))/(1 + times)
        p.m.below <- (1 + sum(statreslist$mean <= testm))/(1 + times)
        for (i in 1:ng){
          pa[i] <- (1 + sum(statreslist$gr[i,] >= test$gr[i]))/(1 + times)
          pb[i] <- (1 + sum(statreslist$gr[i,] <= test$gr[i]))/(1 + times)
        }
        datac <- test
        tuning <- NA
        groupinfo$testm <- testm
        groupinfo$pa <- pa
        groupinfo$pb <- pb
        groupinfo$pma <- p.m.above
        groupinfo$pmb <- p.m.below
        cat("Data value: ", datac$overall, "\n")
    }
    else
      cat("Data value: ", datac, "\n")
    if (!prabobj$spatial || is.null(sarestimate)) sarlambda <- NULL
    else sarlambda <- sarestimate$lambda*sarestimate$nbweight
    if (teststat=="groups")
      results <- statreslist
    else
      results=statres
    out <- list(results = results, p.above = p.above, p.below = p.below, 
        datac = datac, tuning = tuning, distance=dist, times=times,
                teststat=teststat, pd=p.nb,
                abund=!is.null(sarestimate),
                sarlambda=sarlambda, sarestimate=sarestimate,
                groupinfo=groupinfo)
    class(out) <- "prabtest"
    out
}


    
toprab <- function(prabobj)
  prabobj$prab>0

    
build.nblist <- function(prabobj,prab01=NULL,style="C"){
#  require(spdep)
  if (is.null(prab01))
    prab01 <- prabinit(prabmatrix=toprab(prabobj),rows.are.species=FALSE,
                      distance="none")
  nblist <- list()
  q <- 1
  ijsum <- 0
  for (i in 1:prabobj$n.species){
#    print(i)
    iregs <- (1:prabobj$n.regions)[prab01$prab[,i]]
    for (j in 1:prab01$regperspec[i]){
#      print(j)
      nblist[[q]] <-
        (1:length(iregs))[iregs %in% prabobj$nb[[iregs[j]]]]+ijsum
#      if (sum(nblist[[q]]==q)>0) cat(q," ",i," ",j," ",iregs[j]," ",ijsum,"\n")
      q <- q+1
    }
    ijsum <- ijsum+prab01$regperspec[i]
  }
  nblist <- lapply(nblist,as.integer)
  nblist[sapply(nblist, length) == 0L] <- 0L
#  print(nblist)
  class(nblist) <- "nb"
  out <- spdep::nb2listw(nblist,style=style,zero.policy=TRUE)
  invisible(out)
}

# columns are species
# fixed species
# if neighborhood="none", sampling=TRUE is faster
# abmat is a prab object
prab.sarestimate <- function(abmat, prab01=NULL,sarmethod="eigen",
                             weightstyle="C",
                             quiet=TRUE, sar=TRUE,
                             add.lmobject=TRUE){
#  if (sar) require(spdep)
  if (is.null(prab01))
    prab01 <- prabinit(prabmatrix=toprab(abmat),rows.are.species=FALSE,
                      distance="none")
  # single and regions total abundances log transformed
#  print("sarestimate logs")
  logabund <- log(abmat$prab[prab01$prab[,1],1])
  species <- rep(1,sum(prab01$prab[,1]))
  region <- (1:abmat$n.regions)[prab01$prab[,1]]
  for (j in 2:abmat$n.species){
    logabund <- c(logabund,log(abmat$prab[prab01$prab[,j],j]))
    species <- c(species,rep(j,sum(prab01$prab[,j])))
    region <- c(region,(1:abmat$n.regions)[prab01$prab[,j]])
  }
  species <- as.factor(species)
  region <- as.factor(region)
#  nblistw <- build.nblist(abmat,prab01=prab01,style=weightstyle)
  abundreg <- data.frame(logabund,species,region,
                         row.names=sapply(1:length(species),toString))
  if (sar){
      nblistw <- build.nblist(abmat,prab01=prab01,style=weightstyle)
#    print("errorsarlm")
      abundlm <- spdep::errorsarlm(logabund~region+species,data=abundreg,
                        listw=nblistw,quiet=quiet,zero.policy=TRUE,
                        method=sarmethod)
      interc <- coef(abundlm)[2] # was 1
      sigma <- sqrt(summary(abundlm)$s2)
      regeffects <- c(0,coef(abundlm)[3:(abmat$n.regions+1)]) # was one down
      speffects <- c(0,coef(abundlm)[(abmat$n.regions+2):
                                 (abmat$n.regions+abmat$n.species)]) # dito
      lambda <- abundlm$lambda
      nbweight <- mean(c(nblistw[[3]],recursive=TRUE))
      if (!add.lmobject) abundlm <- NULL
#    print("sarestimate end")
      out <- list(sar=sar,intercept=interc,sigma=sigma,regeffects=regeffects,
              speffects=speffects,lambda=lambda,size=length(nblistw[[3]]),
              nbweight=nbweight,lmobject=abundlm)
  }
  else{
#    print("lm")
    abundlm <- lm(logabund~region+species,data=abundreg)
    interc <- coef(abundlm)[1]
    sigma <- summary(abundlm)$sigma
    regeffects <- c(0,coef(abundlm)[2:abmat$n.regions])
    speffects <- c(0,coef(abundlm)[(abmat$n.regions+1):
                                 (abmat$n.regions+abmat$n.species-1)])
    if (!add.lmobject) abundlm <- NULL
#    print("sarestimate end")
    out <- list(sar=sar,intercept=interc,sigma=sigma,regeffects=regeffects,
              speffects=speffects,lmobject=abundlm)
  }
  out
}
 
# columns are species
# fixed species
# if neighborhood="none", sampling=TRUE is faster
# If sarestimate$sar=FALSE, this does also simulation without
# SAR-neighborhood dependent abundances.
# p.nb=NULL > ignore neighbors
regpop.sar <- function(abmat, prab01=NULL,
                       sarestimate=prab.sarestimate(abmat),
                    p.nb=NULL,
                    vector.species=prab01$regperspec,
                    pdf.regions=prab01$specperreg/(sum(prab01$specperreg)),
                   count=FALSE){
#  require(spdep)
#  require(mvtnorm)
  if (is.null(prab01)){
    prab01 <- prabinit(prabmatrix=toprab(abmat),rows.are.species=FALSE,
                      distance="none")
    vector.species=prab01$regperspec
    pdf.regions=prab01$specperreg/(sum(prab01$specperreg))
  }
  proble <- function(v,val) mean(v<=val, na.rm=TRUE)
  # single and regions total abundances log transformed
#  print("Computing the logs")
  logabund <- log(abmat$prab[prab01$prab[,1],1])
  species <- rep(1,sum(prab01$prab[,1]))
  region <- (1:abmat$n.regions)[prab01$prab[,1]]
  for (j in 2:abmat$n.species){
    logabund <- c(logabund,log(abmat$prab[prab01$prab[,j],j]))
    species <- c(species,rep(j,sum(prab01$prab[,j])))
    region <- c(region,(1:abmat$n.regions)[prab01$prab[,j]])
  }
  species <- as.factor(species)
  region <- as.factor(region)
  neighbors <- abmat$nb
  m01 <- matrix(FALSE, ncol = abmat$n.species, nrow = abmat$n.regions)
  out <- matrix(0, ncol = abmat$n.species, nrow = abmat$n.regions)
  cdf.local <- cdf.regions <- c()
  for (i in 1:abmat$n.regions) cdf.regions[i] <- sum(pdf.regions[1:i])
  for (i in 1:abmat$n.species){
    if (count) 
        cat("Species ", i, "\n")
    spec.regind <- spec.neighb <- rep(FALSE, abmat$n.regions)
    nsize <- vector.species[i]
    if (is.null(p.nb)){
      m01[,i] <- rep(FALSE,abmat$n.regions)
      m01[sample(abmat$n.regions,nsize,prob=pdf.regions),i] <- rep(TRUE,nsize)
    }
    else{
#      print(p.nb)
      r1 <- runif(1)
      reg <- 1 + sum(r1 > cdf.regions)
      spec.regind[reg] <- TRUE
      for (k in neighbors[[reg]]) spec.neighb[k] <- TRUE
      m01[reg, i] <- TRUE
      if (nsize > 1) 
        for (j in 2:nsize)
          if (all(!spec.neighb) | all(pdf.regions[spec.neighb] < 
           1e-08) | all(spec.neighb | spec.regind) |
             all(pdf.regions[!(spec.regind | spec.neighb)] < 1e-08)) {
              nreg <- sum(!spec.regind)
              pdf.local <- pdf.regions[!spec.regind]
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg) cdf.local[l] <- sum(pdf.local[1:l])
              r1 <- runif(1)
              zz <- 1 + sum(r1 > cdf.local[1:nreg])
              reg <- (1:abmat$n.regions)[!spec.regind][zz]
              spec.regind[reg] <- TRUE
              spec.neighb[reg] <- FALSE
              for (k in neighbors[[reg]]) spec.neighb[k] <- !(spec.regind[k])
              m01[reg, i] <- TRUE
          }
          else if (runif(1) < p.nb) {
              regs <- !(spec.regind | spec.neighb)
              nreg <- sum(regs)
              pdf.local <- pdf.regions[regs]
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg) cdf.local[l] <- sum(pdf.local[1:l])
              r1 <- runif(1)
              zz <- 1 + sum(r1 > cdf.local[1:nreg])
              reg <- (1:abmat$n.regions)[regs][zz]
              spec.regind[reg] <- TRUE
              for (k in neighbors[[reg]]) spec.neighb[k] <- !(spec.regind[k])
              m01[reg, i] <- TRUE
          }
          else {
              nreg <- sum(spec.neighb)
              pdf.local <- pdf.regions[spec.neighb]
              pdf.local <- pdf.local/sum(pdf.local)
              for (l in 1:nreg) cdf.local[l] <- sum(pdf.local[1:l])
              r1 <- runif(1)
              zz <- 1 + sum(r1 > cdf.local[1:nreg])
              reg <- (1:abmat$n.regions)[spec.neighb][zz]
              spec.regind[reg] <- TRUE
              spec.neighb[reg] <- FALSE
              for (k in neighbors[[reg]]) spec.neighb[k] <- !(spec.regind[k])
              m01[reg, i] <- TRUE
          }
    }
    iregions <- (1:abmat$n.regions)[m01[,i]]
    if (sarestimate$sar){
      inbmatrix <- matrix(0,nrow=nsize,ncol=nsize)
      for (j in 1:nsize)
        inbmatrix[j,(1:nsize)[iregions %in%
                                       abmat$nb[[iregions[j]]]]] <- 1
      inbmatrix <- sarestimate$lambda*sarestimate$nbweight*inbmatrix
#    print(iregions)
#    print(inbmatrix)
      invmatrix <- solve(diag(nsize)-inbmatrix)
      icov <- sarestimate$sigma^2*invmatrix %*% invmatrix
#    print(icov)
      ierror <- mvtnorm::rmvnorm(1,sigma=icov)
#    print(ierror)
    }
    else
      ierror <- rnorm(nsize,sd=sarestimate$sigma)      
    for (j in 1:nsize){
     abundmean <- sarestimate$intercept+sarestimate$speffects[i]+
            sarestimate$regeffects[iregions[j]]
     out[iregions[j],i] <- exp(abundmean+ierror[j])
#      print(abundmean)
    }
#    break
  }
  out
}

specgroups <- function (distmat,groupvector, groupinfo) 
{
    distmat <- as.matrix(distmat)
    nc <- ncol(distmat)
    sgd <- mgd <- numeric(0)
    sni <- 0
    for (i in 1:groupinfo$ng){
      gd <- distmat[groupvector==groupinfo$lg[i],
                    groupvector==groupinfo$lg[i]]
      ni <- groupinfo$nsg[i]
      ni <- ni*(ni-1)/2
      sni <- sni+ni
#      print(i)
#      print(gd)
      sgd[i] <- sum(gd[upper.tri(gd)])
      mgd[i] <- sgd[i]/ni
    }
#    print(gd)
#    print(sgd)
#    print(ni)
#    print(sni)
    overall <- sum(sgd)/sni
    out <- list(overall=overall,gr=mgd)
    out
}

