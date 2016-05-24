# Includes ignore.richness-parameter
# dist="dice": dice distance
prabtest <- function (prabobject, teststat = "distratio", tuning = switch(teststat, 
    distratio = 0.25, lcomponent = floor(3 * ncol(prabobject$distmat)/4), 
    isovertice = ncol(prabobject$distmat), nn = 4, NA), times = 1000, 
    pd = NULL, prange = c(0, 1), nperp = 4, step = 0.1, step2=0.01,
                      twostep = TRUE, 
    sf.sim = FALSE, sf.const = sf.sim, pdfnb = FALSE, ignore.richness=FALSE) 
{
    if (is.null(pd) & prabobject$spatial) 
        ac <- autoconst(prabobject, twostep = twostep, prange = prange, 
            nperp = nperp, step1 = step, step2=step2, species.fixed = sf.const,
                        ignore.richness=ignore.richness)
    else {
        if (is.null(pd)) 
            pd <- 1
        ac <- list(pd = pd, coef = NA)
    }
    psim <- pop.sim(prabobject$prab, prabobject$nb, teststat = teststat, h0c = ac$pd, 
        dist = prabobject$distance, times = times, geodist = prabobject$geodist, gtf = prabobject$gtf,
                    testc = tuning, n.species = prabobject$n.species, 
        specperreg = prabobject$specperreg, regperspec = prabobject$regperspec, 
        species.fixed = sf.sim, pdfnb = pdfnb, ignore.richness=ignore.richness)
    out <- list(results = psim$results, datac = psim$datac, p.value = ifelse(teststat == 
        "inclusions", psim$p.above, psim$p.below), tuning = tuning, 
        pd = ac$pd, reg = ac$coef, teststat = teststat, distance = prabobject$distance, 
        times = times, pdfnb = pdfnb, ignore.richness=ignore.richness)
    class(out) <- "prabtest"
    out
}

"summary.prabtest" <-
function(object, above.p=object$teststat %in% c("groups","inclusions","mean"),
         group.outmean=FALSE,...){
  if (object$teststat!="groups"){
    rrange <- range(object$results)
    rmean <- mean(object$results)
    groupinfo <- NULL
  }
  else{
    rrange <- range(object$results$overall)
    rmean <- mean(object$results$overall)
    rrangem <- range(object$results$mean)
    rmeanm <- mean(object$results$mean)
    rrangeg <- matrix(0,ncol=2,nrow=object$groupinfo$ng)
    rmeang <- numeric(0)
    for (i in 1:object$groupinfo$ng){
      rrangeg[i,] <- range(object$results$gr[i,])
      rmeang[i] <- mean(object$results$gr[i,])
    }    
    groupinfo <- c(object$groupinfo,list(rrangeg=rrangeg,rmeang=rmeang,
                                         rrangem=rrangem,rmeanm=rmeanm))
  }
  # Aendern in print.summary.prabtest: Verhalten bei is.null(pd), pdfnb=FALSE
  # und Ueberschrift!
  p.value <- if (is.null(object$abund)) object$p.value
             else{
               if (above.p) object$p.above
               else object$p.below
             }
  out <- list(rrange=rrange, rmean=rmean, datac=object$datac,
              p.value=p.value,
              pd=object$pd, tuning=object$tuning,
              teststat=object$teststat, distance=object$distance,
              times=object$times,
              pdfnb=object$pdfnb,
              abund=object$abund, sarlambda=object$sarlambda,
              groupinfo=groupinfo, group.outmean=group.outmean,
              ignore.richness=object$ignore.richness, nbcols=object$nbcols,
              allele=object$allele)
  class(out) <- "summary.prabtest"
  out
}


"print.summary.prabtest" <-
function(x, ...){
  if (!is.null(x$abund))
    cat("* Parametric bootstrap test for spatial abundance data *\n\n")
  else{
    if(!is.null(x$allele))
      cat("* Parametric bootstrap test for allele/loci data *\n\n")
    else
      cat("* Parametric bootstrap test for presence-absence data *\n\n")
  }
  cat("Test statistics: ",x$teststat,", Tuning constant=",x$tuning,"\n")
  cat("Distance: ",x$distance,"\n")
  cat("Simulation runs: ",x$times,"\n")
  if(is.null(x$pd))    
    cat("A model without spatial autocorrelation was used.\n")
  else{
    cat("Disjunction parameter for presence-absence pattern: ",x$pd,"\n")
    if (!is.null(x$abund))
      cat("Neighborhood parameter lambda for SAR-model: ",x$sarlambda,"\n")
  }
  if (!is.null(x$ignore.richness))
    if (x$ignore.richness)
      cat("Rows (regions) richness has been ignored by the null model.\n")
  if (!is.null(x$nbcol))
    if (x$nbcol)
      cat("Neighborhood belongs to columns (species) instead of rows.\n")
  if (!is.null(x$pdfnb))
    if (x$pdfnb)
      cat("Neighbor-based correction of region probabilities was used.\n")
  if(x$teststat=="groups"){
    cat("Mean within group distances for original data: ",x$datac$overall,"\n")
    cat("Mean of mean within group distances for null data: ",x$rmean,
        ", range: ",x$rrange,"\n")
    cat("p= ",x$p.value,"\n")
    if (x$group.outmean){
      cat("Overall mean for original data: ",x$groupinfo$testm,"\n")
      cat("Mean of overall mean for null data: ",x$groupinfo$rmeanm,
          ", range: ",x$groupinfo$rrangem,"\n")
      cat("p= ",x$groupinfo$pma,"\n")
    }
    for (i in 1:x$groupinfo$ng)
      if (x$groupinfo$nsg[i]>1){
        cat("  Group ",x$groupinfo$lg[i],
            " statistics value for original data: ",
          x$datac$gr[i],"\n")
        cat("  Group ",x$groupinfo$lg[i],
            " mean for null data: ",x$groupinfo$rmeang[i],
          ", range: ",x$groupinfo$rrangeg[i,],"\n")
        cat("  p= ",x$groupinfo$pa[i],"\n")
      }
      
  }
  else{
    cat("Statistics value for original data: ",x$datac,"\n")
    cat("Mean for null data: ",x$rmean,", range: ",x$rrange,"\n")
    cat("p= ",x$p.value,"\n")
  }
}



autoconst <- function (x, prange = c(0, 1), twostep = TRUE, step1 = 0.1, step2 = 0.01, 
    plot = TRUE, nperp = 4, ejprob = NULL, species.fixed = TRUE, 
    pdfnb = FALSE, ignore.richness=FALSE) 
{
    probs <- prange[1] + step1 * (0:round((prange[2] - prange[1])/step1))
    if (is.null(ejprob)) {
        cat("  Calculating disjunction probability for original data ")
# str(x$prab)
# print(x$nb)        
        positivespecies <- x$regperspec>0
        cn <- con.regmat(x$prab, x$nb)[positivespecies]
        ejumps <- sum(cn - 1)
        den <- sum(x$regperspec[positivespecies] - 1)
        if (den == 0 | ejumps == 0) 
            ejprob <- 0
        else ejprob <- ejumps/den
        cat(ejprob, "\n")
    }
    out <- list()
    if (ejprob > 0) 
        out <- autoreg(x, probs, ejprob, plot, nperp, species.fixed = species.fixed, 
            pdfnb = pdfnb, ignore.richness=ignore.richness)
    if (ejprob == 0) 
        out$pd <- 0
    if (out$pd < 0) 
        out$pd <- 0
    if (out$pd > 1) 
        out$pd <- 1
    if (twostep){
        out1 <- out
        if (out$pd==0)
          prange2 <- c(0,10*step2)
        else 
          prange2 <- c(max(0, out$pd - 5 * step2), min(1, out$pd + 
            5 * step2))
        probs <- prange2[1] + step2 * (0:round((prange2[2] - 
            prange2[1])/step2))
        out <- autoreg(x, probs, ejprob, plot, nperp, species.fixed = species.fixed, 
            pdfnb = pdfnb, ignore.richness=ignore.richness)
        if (out$pd < 0 | out$pd > 1) 
            out <- out1
    }
    out <- c(out, list(ejprob = ejprob))
    out
}

autoreg <- function (x, probs, ejprob, plot = TRUE, nperp = 4,
                     species.fixed = TRUE, 
    pdfnb = FALSE, ignore.richness=FALSE) 
{
    tjumps <- matrix(nrow = nperp * length(probs), ncol = 3)
    for (j in 1:length(probs)) {
        cat("    Estimating disj. parameter: Simulations for p= ", 
            probs[j], "\n")
        for (i in 1:nperp) {
            if (ignore.richness)
              test <- randpop.nb(x$nb, p.nb = probs[j], n.species = x$n.species, 
                vector.species = x$regperspec, species.fixed = species.fixed,  
                count = FALSE, pdfnb = pdfnb)
            else
              test <- randpop.nb(x$nb, p.nb = probs[j], n.species = x$n.species, 
                vector.species = x$regperspec, species.fixed = species.fixed,
                                 pdf.regions = x$specperreg/sum(x$specperreg), 
                count = FALSE, pdfnb = pdfnb)
            nst <- apply(test, 2, sum)
            tcn <- con.regmat(test, x$nb, count = FALSE)
            ind <- (j - 1) * nperp + i
            tjumps[ind, 1] <- probs[j]
            tjumps[ind, 2] <- sum(tcn - 1)
            tjumps[ind, 3] <- tjumps[ind, 2]/sum(nst - 1)
        }
    }
    ejlm <- lm(tjumps[, 3] ~ tjumps[, 1])
    if (plot) {
        plot(tjumps[, 1], tjumps[, 3], xlab = "pdisj", ylab = "qdisj")
        abline(ejlm$coef)
        abline(c(ejprob, 0), lty = 2)
    }
    pd <- (ejprob - ejlm$coef[1])/ejlm$coef[2]
    cat("  Estimated disjunction parameter =", pd, "\n")
    out <- list(pd = pd, coef = ejlm$coef)
    out
}

pop.sim <-
  function (regmat, neighbors, h0c = 1, times = 200, dist = "kulczynski", 
    teststat = "isovertice", testc = NULL, geodist = NULL, gtf = 0.1,
            n.species = ncol(regmat), 
    specperreg = NULL, regperspec = NULL, species.fixed = FALSE, 
    pdfnb = FALSE, ignore.richness=FALSE) 
{
    statres <- rep(0, times)
    if (is.null(specperreg)) 
        nregions <- apply(regmat, 1, sum)
    else nregions <- specperreg
    if (is.null(regperspec)) 
        nspecies <- apply(regmat, 2, sum)
    else nspecies <- regperspec
    for (i in 1:times) {
        cat("Simulation run ", i)
        if (ignore.richness)
          mat <- randpop.nb(neighbors, p.nb = h0c, n.species = n.species, 
            vector.species = nspecies, species.fixed = species.fixed, 
            count = FALSE, 
            pdfnb = pdfnb)
        else
          mat <- randpop.nb(neighbors, p.nb = h0c, n.species = n.species, 
            vector.species = nspecies, species.fixed = species.fixed, 
            pdf.regions = nregions/sum(nregions), count = FALSE, 
            pdfnb = pdfnb)
#        str(geodist)
#        print(max(geodist))
        if (teststat != "inclusions") {
            if (dist == "jaccard") 
                distm <- jaccard(mat)
            if (dist == "kulczynski") 
                distm <- kulczynski(mat)
            if (dist=="dice")
                distm <- dicedist(mat)
            if (dist == "geco") 
                distm <- geco(mat, geodist, tf = gtf)
#            print("distm")
#            str(distm)
        }
        else statres[i] <- incmatrix(mat)$ninc
        if (teststat == "isovertice") {
            test <- homogen.test(distm, ne = testc)
            statres[i] <- test$iv
        }
        if (teststat == "lcomponent") 
            statres[i] <- lcomponent(distm, ne = testc)$lc
        if (teststat == "distratio") 
            statres[i] <- distratio(distm, prop = testc)$dr
        if (teststat == "nn") 
            statres[i] <- nn(distm, ne = testc)
        cat(" statistics value=", statres[i], "\n")
    }
    if (teststat != "inclusions") {
        if (dist == "jaccard") 
            distm <- jaccard(regmat)
        if (dist == "kulczynski") 
            distm <- kulczynski(regmat)
        if (dist=="dice")
                distm <- dicedist(regmat)
        if (dist == "geco") 
            distm <- geco(regmat, geodist, tf = gtf)
    }
    else {
        test <- incmatrix(regmat)$ninc
        p.above <- (1 + sum(statres >= test))/(1 + times)
        p.below <- (1 + sum(statres <= test))/(1 + times)
        datac <- test
    }
    if (teststat == "isovertice") {
        test <- homogen.test(distm, ne = testc)
        p.above <- (1 + sum(statres >= test$iv))/(1 + times)
        p.below <- (1 + sum(statres <= test$iv))/(1 + times)
        pb <- min(p.above, p.below) * 2
        p.above <- max(p.above, p.below)
        p.below <- pb
        datac <- test$iv
        testc <- test$ne
    }
    if (teststat == "lcomponent") {
        test <- lcomponent(distm, ne = testc)
        p.above <- (1 + sum(statres >= test$lc))/(1 + times)
        p.below <- (1 + sum(statres <= test$lc))/(1 + times)
        datac <- test$lc
        testc <- test$ne
    }
    if (teststat == "nn") {
        test <- nn(distm, ne = testc)
        p.above <- (1 + sum(statres >= test))/(1 + times)
        p.below <- (1 + sum(statres <= test))/(1 + times)
        datac <- test
    }
    if (teststat == "distratio") {
        test <- distratio(distm, prop = testc)
        p.above <- (1 + sum(statres >= test$dr))/(1 + times)
        p.below <- (1 + sum(statres <= test$dr))/(1 + times)
        datac <- test$dr
        testc <- test$prop
    }
    cat("Data value: ", datac, "\n")
    out <- list(results = statres, p.above = p.above, p.below = p.below, 
        datac = datac, testc = testc)
    out
}


# "pop.sim" <-
# function(regmat, neighbors, h0c=1, times=200, dist="kulczynski",
#                     teststat="isovertice",testc=NULL,geodist=NULL,gtf=0.1,
#                     n.species=ncol(regmat), specperreg=NULL,
#                     regperspec=NULL, species.fixed=FALSE, pdfnb=FALSE){
#   statres <- rep(0,times)
#   if (is.null(specperreg))
#       nregions <- apply(regmat,1,sum)
#   else
#       nregions <- specperreg
#   if (is.null(regperspec))
#       nspecies <- apply(regmat,2,sum)
#   else
#       nspecies <- regperspec
#   for (i in 1:times)
#   {
#     cat("Simulation run ",i)
#     mat <- randpop.nb(neighbors,p.nb=h0c,n.species=n.species,
#                         vector.species=nspecies, species.fixed=species.fixed,
#                         pdf.regions=nregions/sum(nregions),count=FALSE,
#                       pdfnb=pdfnb)
#     if (teststat!="inclusions"){
#       if (dist=="jaccard")
#         distm <- jaccard(mat)
#       if (dist=="kulczynski")
#         distm <- kulczynski(mat)
#       if (dist=="geco")
#         distm <- geco(mat,geodist,tf=gtf)
#     }
#     else
#       statres[i] <- incmatrix(mat)$ninc
#     if (teststat=="isovertice")
#     {
#       test <- homogen.test(distm,ne=testc)
#       statres[i] <- test$iv
#     }
#     if (teststat=="lcomponent")
#       statres[i] <- lcomponent(distm,ne=testc)$lc
#     if (teststat=="distratio")
#       statres[i] <- distratio(distm,prop=testc)$dr
#     if (teststat=="nn")
#       statres[i] <- nn(distm,ne=testc)
#     cat(" statistics value=",statres[i],"\n")
#   }
#   if (teststat!="inclusions"){
#     if (dist=="jaccard")
#       distm <- jaccard(regmat)
#     if (dist=="kulczynski")
#       distm <- kulczynski(regmat)
#     if (dist=="geco")
#       distm <- geco(regmat,geodist,tf=gtf)
#   }
#   else{
#     test <- incmatrix(regmat)$ninc
#     p.above <- (1+sum(statres>=test))/(1+times)
#     p.below <- (1+sum(statres<=test))/(1+times)
#     datac <- test
#   }
#   if (teststat=="isovertice")
#   {
#     test <- homogen.test(distm,ne=testc)
#     p.above <- (1+sum(statres>=test$iv))/(1+times)
#     p.below <- (1+sum(statres<=test$iv))/(1+times)
#     pb <- min(p.above,p.below)*2
#     p.above <- max(p.above,p.below)
#     p.below <- pb
#     datac <- test$iv
#     testc <- test$ne
#   }
#   if (teststat=="lcomponent")
#   {
#     test <- lcomponent(distm,ne=testc)
#     p.above <- (1+sum(statres>=test$lc))/(1+times)
#     p.below <- (1+sum(statres<=test$lc))/(1+times)
#     datac <- test$lc
#     testc <- test$ne
#   }
#   if (teststat=="nn")
#   {
#     test <- nn(distm,ne=testc)
#     p.above <- (1+sum(statres>=test))/(1+times)
#     p.below <- (1+sum(statres<=test))/(1+times)
#     datac <- test
#   }
#   if (teststat=="distratio")
#   {
#     test <- distratio(distm,prop=testc)
#     p.above <- (1+sum(statres>=test$dr))/(1+times)
#     p.below <- (1+sum(statres<=test$dr))/(1+times)
#     datac <- test$dr
#     testc <- test$prop
#   }
#   cat("Data value: ",datac,"\n")
#   out <- list(results=statres,p.above=p.above,p.below=p.below,datac=datac,
#               testc=testc)
#   out
# }
