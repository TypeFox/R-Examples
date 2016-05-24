marker_h2 <-
function(data.vector,geno.vector,covariates=NULL,K,alpha=0.05,eps=1e-06,max.iter=100,fix.h2=FALSE,h2=0.5)   {

    her.frame            <- data.frame(dat=data.vector,genotype=geno.vector)
    her.frame$genotype   <- as.character(her.frame$genotype)

    if (!is.null(covariates)) {
      covariates            <- as.data.frame(covariates)
      n.cov                 <- ncol(covariates)
      names(covariates)     <- paste('c',1:n.cov,sep='')
      row.names(covariates) <- row.names(her.frame)
      her.frame             <- cbind(her.frame,covariates)
    }

    her.frame     <- her.frame[!is.na(her.frame$dat),]
    n.rep.vector  <- as.integer(table(her.frame$genotype))

    # re-scale the kinship matrix, for those genotypes for which phenotypic data are available
    K <- K[unique(her.frame$genotype),unique(her.frame$genotype)]
    K <- K / KinshipTransform(K)

    if (is.null(covariates)) {
      VALUE <- as.data.frame(her.frame$dat)
    } else {
      cov.part <- paste(paste('c',1:n.cov,sep=''),collapse='+')
      VALUE <- cbind(her.frame$dat,her.frame[,(ncol(her.frame)-n.cov+1):(ncol(her.frame))])
    }
    names(VALUE)[1] <- 'value'

    reml.obj <- AI_algorithm(Y=data.frame(genotype=her.frame$genotype,VALUE),K=K,eps=eps,max.iter=max.iter,fix.h2=fix.h2,h2=h2)

    reml.obj$gammas        <- c(reml.obj$varcomp[1],1,reml.obj$varcomp[2])

    names(reml.obj$gammas) <- c("giv(genotype, var = TRUE).giv","R!variance","R!units.var")
    reml.obj$gammas.type   <- c(2,1,2)
    names(reml.obj$gammas.type) <- names(reml.obj$gammas)
    reml.obj$ai            <- reml.obj$inv.ai[lower.tri(reml.obj$inv.ai,diag=TRUE)] #reml.obj$inv.ai


    reml.obj$ai            <- c(reml.obj$ai[1],0,0,reml.obj$ai[2],0,reml.obj$ai[3])
    var.comp <- reml.obj$varcomp

    inv.ai <- reml.obj$inv.ai.new


    h2.estimate <- var.comp[1] / (var.comp[1] + var.comp[2])

    st.error1 <- as.numeric(pin(reml.obj,h2 ~ V1 / (V1 + V3))[2])

    conf.int1 <- c(h2.estimate - qnorm(1-alpha/2)*st.error1,h2.estimate + qnorm(1-alpha/2)*st.error1)

    st.error2 <- as.numeric(pin(reml.obj,h2 ~ log(V3 / V1))[2])

    conf.int2 <- c(1 / (1 + exp(log(var.comp[2] /var.comp[1]) + qnorm(1-alpha/2)*st.error2)),1 / (1 + exp(log(var.comp[2]/var.comp[1]) - qnorm(1-alpha/2)*st.error2)))

return(list(va=var.comp[1],ve=var.comp[2],h2=h2.estimate,conf.int1=conf.int1,conf.int2=conf.int2,
            inv.ai=reml.obj$inv.ai,loglik=reml.obj$loglik))

}
