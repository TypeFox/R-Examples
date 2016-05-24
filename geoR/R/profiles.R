"proflik" <- 
  function(obj.likfit, geodata, coords = geodata$coords,
           data = geodata$data,
           sill.values,
           range.values, 
           nugget.values,
           nugget.rel.values,
           lambda.values,
           sillrange.values = TRUE,
           sillnugget.values = TRUE,
           rangenugget.values = TRUE, 
           sillnugget.rel.values = FALSE,
           rangenugget.rel.values = FALSE, 
           silllambda.values = FALSE,
           rangelambda.values = TRUE, 
           nuggetlambda.values = FALSE,
           nugget.rellambda.values = FALSE,
           uni.only = TRUE,
           bi.only = FALSE,
           messages,
           ...)
{
  ##
  ## 1. setting arguments
  ##
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  if(missing(geodata))
    geodata <- list(coords = coords, data = data)
#  if(! "package:stats" %in% search()) require(mva)
  call.fc <- match.call()
  n.cov.pars <- obj.likfit$npars - length(obj.likfit$beta)
  if(obj.likfit$transform.info$fix.lambda == FALSE)
    n.cov.pars <- n.cov.pars - 1
  if(missing(sill.values))
    sill.values <- sillrange.values <- sillnugget.values <- FALSE
  if(missing(range.values))
    range.values <- sillrange.values <- rangenugget.values <- rangelambda.values <- FALSE 
  if(!is.null(obj.likfit$call$fix.nugget))
    if(obj.likfit$call$fix.nugget == TRUE)
      nugget.values <-  nugget.rel.values <- FALSE
  if(missing(nugget.values))
    nugget.values <- sillnugget.values <- rangenugget.values <- nuggetlambda.values <- FALSE
  if(missing(nugget.rel.values))
    nugget.rel.values <- sillnugget.rel.values <- rangenugget.rel.values <- nugget.rellambda.values <- FALSE
  if(missing(lambda.values) | obj.likfit$transform.info$fix.lambda) lambda.values <- FALSE  
  if(uni.only){
    sillrange.values <- sillnugget.values <- rangenugget.values <-
      sillnugget.rel.values <- rangenugget.rel.values <- 
        silllambda.values <- rangelambda.values <- 
          nugget.rellambda.values <- nuggetlambda.values <- FALSE
  }
  else{
    if(all(sillrange.values == TRUE)){
      if(all(sill.values == FALSE) | all(range.values == FALSE)){
        sillrange.values <- FALSE
        stop("if argument sillrange.values = TRUE sill.values and range.values must be provided. Alternatively a matrix can be provided in sillrange.values  or set this to FALSE")
      }
      else
        sillrange.values <- as.matrix(expand.grid(sill.values, range.values))
    }
    if(n.cov.pars == 2){
        sillnugget.values <- rangenugget.values <-
          sillnugget.rel.values <- rangenugget.rel.values <- 
            nugget.rellambda.values <- nuggetlambda.values <- FALSE
    }
    else{
      if(all(sillnugget.values == TRUE)){
        if(all(sill.values == FALSE) | all(nugget.values == FALSE)){
          sillnugget.values <- FALSE
          stop("if argument sillnugget.values = TRUE sill.values and nugget.values must be provided. Alternatively a matrix can be provided in sillnugget.values or set this to FALSE")
        }
        else
          sillnugget.values <- as.matrix(expand.grid(sill.values, nugget.values))
      }
      if(all(rangenugget.values == TRUE)){
        if(all(range.values == FALSE) | all(nugget.values == FALSE)){
          rangenugget.values <- FALSE
          stop("if argument rangenugget.values = TRUE range.values and nugget.values must be provided. Alternatively a matrix can be provided in rangenugget.values or set this to FALSE")
        }
        else
          rangenugget.values <- as.matrix(expand.grid(range.values, nugget.values))
      }
      if(all(sillnugget.rel.values == TRUE)){
        if(all(sill.values == FALSE) | all(nugget.rel.values == FALSE)){
          sillnugget.rel.values <- FALSE
          stop("if argument sillnugget.rel.values = TRUE sill.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in sillnugget.rel.values or set this to FALSE")
        }
        else
          sillnugget.rel.values <- as.matrix(expand.grid(sill.values, nugget.rel.values))
      }
      if(all(rangenugget.rel.values == TRUE)){
        if(all(range.values == FALSE) | all(nugget.rel.values == FALSE)){
          rangenugget.rel.values <- FALSE
          stop("if argument rangenugget.rel.values = TRUE range.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in rangenugget.rel.values or set this to FALSE")
        }
        else
          rangenugget.rel.values <- as.matrix(expand.grid(range.values, nugget.rel.values))
      }
      if(obj.likfit$transform.info$fix.lambda == TRUE){
        if(all(nuggetlambda.values == TRUE)){
          if(all(lambda.values == FALSE) | all(nugget.values == FALSE)){
            nuggetlambda.values <- FALSE
            stop("if argument nuggetlambda.values = TRUE lambda.values and nugget.values must be provided. Alternatively a matrix can be provided in nuggetlambda.values or set this to FALSE")
          }
          else
            nuggetlambda.values <- as.matrix(expand.grid(lambda.values, nugget.values))
        }
        if(all(nugget.rellambda.values == TRUE)){
          if(all(lambda.values == FALSE) | all(nugget.rel.values == FALSE)){
            nugget.rellambda.values <- FALSE
            stop("if argument nugget.rellambda.values = TRUE lambda.values and nugget.rel.values must be provided. Alternatively a matrix can be provided in nugget.rellambda.values or set this to FALSE")
          }
          else
            nugget.rellambda.values <- as.matrix(expand.grid(lambda.values, nugget.rel.values))
        }
      }
    }
    if(obj.likfit$transform.info$fix.lambda == TRUE)
      silllambda.values <- rangelambda.values <- FALSE
    else{
      if(all(silllambda.values == TRUE)){
        if(all(sill.values == FALSE) | all(lambda.values == FALSE)){
          silllambda.values <- FALSE
          stop("if argument silllambda.values = TRUE sill.values and lambda.values must be provided. Alternatively a matrix can be provided in silllambda.values or set this to FALSE")
        }
        else
          silllambda.values <- as.matrix(expand.grid(sill.values, lambda.values))
      }
      if(all(rangelambda.values == TRUE)){
        if(all(range.values == FALSE) | all(lambda.values == FALSE)){
          rangelambda.values <- FALSE
          stop("if argument rangelambda.values = TRUE range.values and lambda.values must be provided. Alternatively a matrix can be provided in rangelambda.values or set this to FALSE")
        }
        else
          rangelambda.values <- as.matrix(expand.grid(range.values, lambda.values))
      }      
    }
  }
  ##
  ## 2. data preparation
  ##
  trend <- unclass(trend.spatial(trend=obj.likfit$trend, geodata = geodata))
  if (nrow(trend) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  data <- as.vector(data)
  dimnames(trend) <- list(NULL, NULL)
  if(obj.likfit$transform.info$fix.lambda) {
    if(obj.likfit$lambda != 1) {
      if(any(data <= 0))
        stop("Data transformation not allowed when there are zeros or negative data"
             )
      if(obj.likfit$lambda == 0)
        data <- log(data)
      else data <- ((data^obj.likfit$lambda) - 1)/obj.likfit$lambda
    }
  }
  n <- length(data)
  dists.vec <- as.vector(dist(coords))
  d <- range(dists.vec)
  min.dist <- d[1]
  max.dist <- d[2]
  tausq <- obj.likfit$nugget
  sigmasq <- obj.likfit$cov.pars[1]
  tausq.rel <- tausq/sigmasq
  phi <- obj.likfit$cov.pars[2]
  lambda <- obj.likfit$lambda
  loglik <- obj.likfit$loglik
  sill.total <- sigmasq + tausq
  n.uni <- 0
  n.bi <- 0
  lower.phi <- 0.01 * (min.dist/max.dist) 
  upper.phi <- 1000 * max.dist
  lower.sigmasq <- 0.01 * sill.total
  result <- list()
  assign(".temp.list", list(n = n,
                            z = data,
                            beta.size = dim(trend)[2],
                            kappa = obj.likfit$kappa,
                            xmat = trend,
                            ## txmat = t(trend),
                            method.lik = obj.likfit$method.lik,
                            dists.lowertri = dists.vec,
                            cov.model = obj.likfit$cov.model,
                            fix.lambda = obj.likfit$transform.info$fix.lambda,
                            lambda = obj.likfit$lambda,
                            lower.phi = lower.phi,
                            upper.phi = upper.phi,
                            lower.sigmasq = lower.sigmasq, 
                            phi.est = phi,
                            tausq.rel.est = tausq.rel,
                            tausq.est = tausq,
                            sigmasq.est = sigmasq), pos=1)
  if(obj.likfit$transform.info$fix.lambda == TRUE)
    eval(substitute(.temp.list$log.jacobian <- xxx, list(xxx=obj.likfit$transform.info$log.jacobian)), envir=.GlobalEnv)
  ##
  ## 3. One-dimentional profile likelihoods
  ##
  ##  
  ## 3.1 Profile for sigmasq
  ##  
  if(bi.only == FALSE) {
    if(any(sill.values != FALSE)) {
      n.uni <- n.uni + 1
      if(messages.screen) cat("proflik: computing profile likelihood for the sill\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          if(obj.likfit$transform.info$fix.lambda == FALSE) {
            ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                  max(range.values),
                                                            l = 5),
                                                  seq(-1,1,l = 5)))
          }
          else {
            ini.grid <- as.matrix(seq(min(range.values),
                                                max(range.values),
                                                l = 10))
          }
          dimnames(ini.grid) <- list(NULL, NULL)
          eval(substitute(.temp.list$ini.grid <- xxx, list(xxx= ini.grid)), envir=.GlobalEnv)
          pl.sigmasq <- apply(matrix(sill.values,
                                     ncol = 1), 1, .proflik.aux2, ...)
          eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
        }
        else {
          stop("not yet implemented for fixed nugget != 0")
        }
      }
      if(n.cov.pars == 3) {
        if(any(lambda.values != FALSE)) {
          ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                max(range.values), l = 6),
                                            seq(0, 2 * tausq.rel, l = 4),
                                            seq(-1, 1, l = 5)))
        }
        else {
          ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                                max(range.values), l = 10),
                                            seq(0, 2 * tausq.rel, l = 4)))
        }
        dimnames(ini.grid) <- list(NULL, NULL)
        eval(substitute(.temp.list$ini.grid <- xxx, list(xxx= ini.grid)), envir=.GlobalEnv)
        pl.sigmasq <- apply(matrix(sill.values, ncol = 
                                   1), 1, .proflik.aux9, ...)
        eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
      }
      v.ord <- order(c(sigmasq, sill.values))
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.sigmasq <- pl.sigmasq + obj.likfit$transform.info$log.jacobian
      result$sill <- list(sill = c(sigmasq, sill.values)[
                            v.ord], proflik.sill = c(loglik, pl.sigmasq)[
                                      v.ord], est.sill = c(sigmasq, loglik))
    }
    ##  
    ## 3.2 Profile for phi
    ##
    if(any(range.values != FALSE)) {
      n.uni <- n.uni + 1
      if(messages.screen) cat("proflik: computing profile likelihood for the range\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          eval(expression(.temp.list$nugget <-  0), envir=.GlobalEnv)
          pl.phi <- apply(matrix(range.values,
                                 ncol = 1), 1, .proflik.aux0, ...)
          eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        pl.phi <- apply(matrix(range.values, ncol = 1),
                        1, .proflik.aux7, ...)
      }
      v.ord <- order(c(phi, range.values))
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.phi <- pl.phi + obj.likfit$transform.info$log.jacobian
      result$range <- list(range = c(phi, range.values)[
                             v.ord], proflik.range = c(loglik, pl.phi)[
                                       v.ord], est.range = c(phi, loglik))
    }
    ##  
    ## 3.3 Profile for \tau^2
    ##  
    if(n.cov.pars == 3) {
      if(any(nugget.values != FALSE)) {
        n.uni <- n.uni + 1
        if(messages.screen) cat("proflik: computing profile likelihood for the nugget\n"
              )
        pl.tausq <- apply(matrix(nugget.values, ncol = 
                                 1), 1, .proflik.aux11, ...)
        v.ord <- order(c(tausq, nugget.values))
        if(obj.likfit$transform.info$fix.lambda == TRUE)
          pl.tausq <- pl.tausq + obj.likfit$transform.info$log.jacobian
        result$nugget <- list(nugget = c(tausq, 
                                nugget.values)[v.ord], proflik.nugget
                              = c(loglik, pl.tausq)[v.ord], 
                              est.nugget = c(tausq, loglik))
      }
      ##  
      ## 3.4 Profile for relative \tau^2
      ##
      if(any(nugget.rel.values != FALSE)) {
        if(messages.screen) cat("proflik: computing profile likelihood for the relative nugget\n"
              )
        n.uni <- n.uni + 1
        pl.tausq.rel <- apply(matrix(nugget.rel.values,
                                     ncol = 1), 1, .proflik.aux5, ...)
        v.ord <- order(c(tausq.rel, nugget.rel.values))
        if(obj.likfit$transform.info$fix.lambda == TRUE)
          pl.tausq.rel <- pl.tausq.rel + obj.likfit$transform.info$log.jacobian
        result$nugget.rel <- list(nugget.rel = c(
                                    tausq.rel, nugget.rel.values)[v.ord],
                                  proflik.nugget = c(loglik, pl.tausq.rel
                                    )[v.ord],
                                  est.nugget.rel = c(tausq.rel,loglik))
      }
    }
    ##  
    ## 3.5 Profile for \lambda
    ##
    if(any(lambda.values != FALSE)) {
      assign(".temp.temp.list", get(".temp.list", pos=1), pos=1)
      eval(substitute(.temp.temp.list$coords <- xxx, list(xxx= coords)), envir=.GlobalEnv)
      n.uni <- n.uni + 1
      if(messages.screen) cat("proflik: computing profile likelihood for lambda\n"
            )
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          eval(expression(.temp.temp.list$fixtau <-  TRUE), envir=.GlobalEnv)
          eval(substitute(.temp.temp.list$ini <- xxx, list(xxx= c(sigmasq,phi))), envir=.GlobalEnv)
          pl.lambda <- apply(as.matrix(lambda.values), 1, .proflik.aux23, ...)
        }
        else {
          stop("not yet implemented for fixed nugget != 0")
        }
      }
      if(n.cov.pars == 3) {
        eval(expression(.temp.temp.list$fixtau <-  FALSE), envir=.GlobalEnv)
        eval(substitute(.temp.temp.list$ini <- xxx, list(xxx= phi)), envir=.GlobalEnv)
        pl.lambda <- apply(matrix(lambda.values,
                                  ncol = 1), 1, .proflik.aux23, ...)
      }
      v.ord <- order(c(lambda, lambda.values))
      result$lambda <- list(lambda = c(lambda, 
                              lambda.values)[v.ord], proflik.lambda
                            = c(loglik, pl.lambda)[v.ord], 
                            est.lambda = c(lambda, loglik))
      remove(".temp.temp.list", inherits=TRUE, pos=1)
    }
  }
  ##
  ## 4. Two-dimentional profile likelihoods
  ##
  ##  
  ## 4.1 Profile for \sigma^2 and \phi
  ##
  if(uni.only == FALSE){
    if(any(sillrange.values != FALSE)) {
      n.bi <- n.bi + 1
      if(messages.screen) cat("proflik: computing 2-D profile likelihood for the sill and range parameters\n")
      if(n.cov.pars == 2) {
        if(tausq == 0) {
          eval(expression(.temp.list$nugget <-  0), envir=.GlobalEnv)
          if(get(".temp.list", pos=1)$fix.lambda == TRUE) {
            pl.sigmasqphi <- apply(cbind(0, sillrange.values, 1), 1, loglik.spatial, ...)
          }
          else {
            pl.sigmasqphi <- apply(sillrange.values,
                                   1, .proflik.aux28, ...)
          }
          eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
        }
        else {
          stop("not yet implemented for fixed nugget != 0"
               )
        }
      }
      if(n.cov.pars == 3) {
        pl.sigmasqphi <- apply(sillrange.values, 1, .proflik.aux13, ...)
      }
      names(pl.sigmasqphi) <- NULL
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.sigmasqphi <- pl.sigmasqphi + obj.likfit$transform.info$log.jacobian
      result$sillrange <- list(sill = as.numeric(levels(as.factor(sillrange.values[,1]))), range = 
                               as.numeric(levels(as.factor(sillrange.values[,2]))), proflik.sillrange = pl.sigmasqphi, 
                               est.sillrange = c(sigmasq, phi, loglik))
    }
    ##  
    ## 4.2 Profile for \sigma^2 and \tau^2
    ##  
    if(any(sillnugget.values != FALSE)) {
      n.bi <- n.bi + 1
      if(messages.screen) cat("proflik: computing 2-D profile likelihood for the sill and nugget\n")
      if(obj.likfit$transform.info$fix.lambda == FALSE)
        ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                              max(range.values), l = 
                                              10), seq(-1, 1, l = 5)))
      else
        ini.grid <- as.matrix(seq(min(range.values),
                                  max(range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      eval(substitute(.temp.list$ini.grid <-  xxx, list(xxx=ini.grid)), envir=.GlobalEnv)
      pl.sigmasqtausq <- apply(sillnugget.values, 1, .proflik.aux15, ...)
      eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
      names(pl.sigmasqtausq) <- NULL
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.sigmasqtausq <- pl.sigmasqtausq + obj.likfit$transform.info$log.jacobian
      result$sillnugget <- list(sill = as.numeric(levels(as.factor(sillnugget.values[,1]))), nugget = as.numeric(levels(as.factor(sillnugget.values[,2]))), proflik.sillnugget = 
                                pl.sigmasqtausq, est.sillrange = c(sigmasq,
                                                   tausq, loglik))
    }
    ##  
    ## 4.3 Profile for \phi and \tau^2
    ##
    if(any(rangenugget.values != FALSE)) {
      n.bi <- n.bi + 1
      if(messages.screen) cat("proflik: computing 2-D profile likelihood for the range and nugget\n"
            )
      eval(substitute(.temp.list$ini.grid <- xxx, list(xxx= as.matrix(seq(sigmasq/4, 5 * 
                                          sigmasq, l = 15)))), envir=.GlobalEnv)
      pl.phitausq <- apply(rangenugget.values, 1, .proflik.aux17, ...)
      eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
      names(pl.phitausq) <- NULL
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.phitausq <- pl.phitausq + obj.likfit$transform.info$log.jacobian
      result$rangenugget <- list(range = as.numeric(levels(as.factor(rangenugget.values[,1]))), nugget
                               = as.numeric(levels(as.factor(rangenugget.values[,1]))), proflik.rangenugget = 
                               pl.phitausq, est.rangenugget = c(phi, tausq,
                                              loglik))
    }
    ##  
    ## 4.4 Profile for \sigma^2 and \tau^2_{rel}
    ##
    if(any(sillnugget.rel.values != FALSE)) {
      n.bi <- n.bi + 1
      if(messages.screen) cat("proflik: computing 2-D profile likelihood for the sill and relative nugget parameters\n"
            )
      if(get(".temp.list", pos=1)$fix.lambda == FALSE)
        ini.grid <- as.matrix(expand.grid(seq(min(range.values), max(range.values), l = 
                                              10), seq(-1, 1, l = 5)))
      else
        ini.grid <- as.matrix(seq(min(range.values),
                                  max(range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      eval(substitute(.temp.list$ini.grid <- xxx, list(xxx= ini.grid)), envir=.GlobalEnv)
      pl.sigmasqtausq.rel <- apply(sillnugget.rel.values, 1, 
                                   .proflik.aux19, ...)
      eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
      names(pl.sigmasqtausq.rel) <- NULL
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.sigmasqtausq.rel <- pl.sigmasqtausq.rel + obj.likfit$transform.info$log.jacobian
      result$sillnugget.rel <- list(sill = as.numeric(levels(as.factor(sillnugget.rel.values[,1]))), 
                                    nugget.rel = as.numeric(levels(as.factor(sillnugget.rel.values[,2]))), 
                                    proflik.sillnugget.rel = pl.sigmasqtausq.rel,
                                    est.sillrange.rel = c(sigmasq, tausq.rel, 
                                      loglik))
    }
    ##  
    ## 4.5 Profile for \phi and \tau^2_{rel}
    ##
    if(any(rangenugget.rel.values != FALSE)) {
      n.bi <- n.bi + 1
      if(messages.screen) cat("proflik: computing 2-D profile likelihood for the range and relative nugget parameters\n"
            )
      pl.phitausq.rel <- apply(rangenugget.rel.values, 1, .proflik.aux30, ...)
      names(pl.phitausq.rel) <- NULL
      if(obj.likfit$transform.info$fix.lambda == TRUE)
        pl.phitausq.rel <- pl.phitausq.rel + obj.likfit$transform.info$log.jacobian
      result$rangenugget.rel <- list(range = as.numeric(levels(as.factor(rangenugget.rel.values[,1]))),
                                   nugget.rel = as.numeric(levels(as.factor(rangenugget.rel.values[,2]))), 
                                   proflik.rangenugget.rel = pl.phitausq.rel,
                                   "est.rangenugget.rel" = c(phi, tausq.rel,
                                     loglik))
    }
  }
  ##  
  ## 4.6 Profile for \sigma^2 and \lambda
  ##
  if(any(silllambda.values != FALSE)) {
    n.bi <- n.bi + 1
    if(messages.screen) cat("proflik: computing 2-D profile likelihood for the sill and transformation parameters\n"
          )
    if(n.cov.pars == 2) {
      ini.grid <- as.matrix(seq(min(range.values), max(
                                                       range.values), l = 10))
      dimnames(ini.grid) <- list(NULL, NULL)
      eval(substitute(.temp.list$ini.grid <- xxx, list(xxx= ini.grid)), envir=.GlobalEnv)
      if(tausq == 0) {
        eval(expression(.temp.list$nugget <-  0), envir=.GlobalEnv)
        pl.sigmasqlambda <- apply(silllambda.values, 1, .proflik.aux24, ...)
        eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
        eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
      }
      else {
        stop("not yet implemented for fixed nugget != 0"
             )
      }
    }
    if(n.cov.pars == 3) {
      ini.grid <- as.matrix(expand.grid(seq(min(range.values),
                                            max(range.values), l = 10), seq(0, 1, l = 5)))
      dimnames(ini.grid) <- list(NULL, NULL)
      eval(substitute(.temp.list$ini.grid <-  xxx, list(xxx=ini.grid)), envir=.GlobalEnv)
      pl.sigmasqlambda <- apply(silllambda.values, 1, .proflik.aux27, ...)
      eval(expression(.temp.list$ini.grid <-  NULL), envir=.GlobalEnv)
    }
    names(pl.sigmasqlambda) <- NULL
    result$silllambda <- list(sill = as.numeric(levels(as.factor(silllambda.values[,1]))), lambda = as.numeric(levels(as.factor(silllambda.values[,2]))), proflik.silllambda = pl.sigmasqlambda,
                              est.silllambda = c(sigmasq, lambda, loglik))
  }
  ##  
  ## 4.7 Profile for \phi and \lambda
  ##
  if(any(rangelambda.values != FALSE)) {
    eval(substitute(.temp.list$data <-  xxx, list(xxx=.temp.list$z)), envir=.GlobalEnv)
    n.bi <- n.bi + 1
    cat("proflik: computing 2-D profile likelihood for the range and transformation parameters\n"
              )
    if(n.cov.pars == 2) {
      if(tausq == 0) {
        eval(expression(.temp.list$nugget <-  0), envir=.GlobalEnv)
        pl.philambda <- apply(rangelambda.values, 1, 
                              .proflik.aux1, ...)
        eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
      }
      else {
        stop("not yet implemented for fixed nugget != 0"
             )
      }
    }
    if(n.cov.pars == 3) {
      pl.philambda <- apply(rangelambda.values, 1, .proflik.aux31, ...)
    }
    names(pl.philambda) <- NULL
    result$rangelambda <- list(range = as.numeric(levels(as.factor(rangelambda.values[,1]))), lambda = as.numeric(levels(as.factor(rangelambda.values[,2]))), proflik.rangelambda = pl.philambda,
                               est.rangelambda = c(phi, lambda, loglik))
  }
  ##  
  ## 4.8 Profile for \tau^2 and \lambda
  ##                                        
  if(any(nuggetlambda.values != FALSE)) {
    n.bi <- n.bi + 1
    cat("proflik: computing 2-D profile likelihood for the nugget and transformation parameters\n"
          )
    pl.nuggetlambda <- apply(nuggetlambda.values, 1, .proflik.aux32, ...)
      names(pl.nuggetlambda) <- NULL
    result$nuggetlambda <- list(nugget = as.numeric(levels(as.factor(nuggetlambda.values[,1]))), lambda = as.numeric(levels(as.factor(nuggetlambda.values[,2])))
                                , proflik.nuggetlambda = pl.nuggetlambda,
                                est.nuggetlambda = c(tausq, lambda, loglik))
  }
  ##  
  ## 4.9 2-D Profile for \tau^2_{rel} and \lambda
  ##
  if(any(nugget.rellambda.values != FALSE)) {
    n.bi <- n.bi + 1
    pl.nugget.rellambda <- apply(nugget.rellambda.values, 1, .proflik.aux33, ...)
    names(pl.nugget.rellambda) <- NULL
    result$nugget.rellambda <- list(nugget.rel = as.numeric(levels(as.factor(nugget.rellambda.values[,1]))),
                                    lambda = as.numeric(levels(as.factor(nugget.rellambda.values[,2]))), proflik.nugget.rellambda = 
                                    pl.nugget.rellambda, est.nugget.rellambda = c(tausq.rel,
                                                           lambda, loglik))
  }
  result$n.uni <- n.uni
  result$n.bi <- n.bi
  result$method.lik <- obj.likfit$method.lik
  result$call <- call.fc
  oldClass(result) <- "proflik"
  return(result)
}

".proflik.aux0" <-
  function(phi, ...)
{
  ## This function computes the value of the profile likelihood for the correlation parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi for each value of \lambda, if this transformation parameter is included in the model
  ## This is an auxiliary function called by likfit.proflik
  ##
  if(get(".temp.list", pos=1)$fix.lambda == TRUE)
    proflik <- .proflik.aux1(philambda = phi)
  else {
    eval(substitute(.temp.list$phi <-  xxx, list(xxx=phi)), envir=.GlobalEnv)
#    proflik <-  - (optim(.temp.list$lambda, .proflik.aux1.1, method="L-BFGS-B", lower
#                         = -2, upper = 2, ...)$value)
    proflik <-  - (optimise(.proflik.aux1.1, lower = -5, upper = 5, ...)$objective)
    return(proflik)
  }
}
".proflik.aux1" <-
  function(philambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi when nugget effect = 0
  .temp.list <- get(".temp.list", pos=1)
  if(length(philambda) == 2) lambda <- philambda[2]
  else lambda <- 1
  n <- .temp.list$n
  main <- .proflik.main(tausq=.temp.list$nugget, sigmasq=1, phi=philambda[1], lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    proflik <-  - (n/2) * log(2 * pi) - main$log.det.to.half -
      (n/2) * log(main$ssresmat/n) - (n/2) + main$
    log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    proflik <-  - ((n - .temp.list$beta.size)/2) * log(2 * pi) - main$
    log.det.to.half - ((n - .temp.list$beta.size)/2) * log(main$ssresmat/
                                                n) - (n/2) + 0.5 * sum(log(eigentrem$values)) + 
                                                  main$log.jacobian
  }
  return(proflik)
}

".proflik.aux10" <-
  function(phitausq.rel.lambda, ...)
{
  .temp.list <- get(".temp.list", pos=1)
  if(length(phitausq.rel.lambda) == 3)
    lambda <- phitausq.rel.lambda[3]
  else lambda <- 1
  phitausq.rel.lambda <- as.vector(phitausq.rel.lambda)
  n <- .temp.list$n
  phi <- phitausq.rel.lambda[1]
  tausq <- phitausq.rel.lambda[2]
  sigmasq <- .temp.list$sigmasq
  main <- .proflik.main(tausq=tausq, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half +
      (n/2) * log(sigmasq) + (0.5/sigmasq) * main$ssresmat -
        main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + ((n - .temp.list$beta.size)/2) * log(sigmasq) +
      (0.5/sigmasq) * main$ssresmat - 0.5 * sum(log(eigentrem$
                                               values)) - main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux11" <-
  function(tausq, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \tau^2.
  ## It requires the minimisation of the function wrt \sigma^2, \phi and \lambda (if the case)  for each value of \tau^2.
  ## This is an auxiliary function called by proflik.
  eval(substitute(.temp.list$nugget <-  xxx, list(xxx=as.vector(tausq))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE) {
    sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$phi.est),
                            .proflik.aux12, method="L-BFGS-B",
                            lower = c(.temp.list$lower.sigmasq,
                              .temp.list$lower.phi),
                            upper=c(+Inf, .temp.list$upper.phi), ...)$value
  }
  else {
    sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$
                              phi.est, .temp.list$lambda), .proflik.aux12,method="L-BFGS-B",  lower = c(.temp.list$lower.sigmasq, .temp.list$lower.phi, -2),
                            upper = c( + Inf, .temp.list$upper.phi, 2), ...)$value
  }
  eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
  return( - sigmasqphi.res)    
}

".proflik.aux1.1" <-
  function(lambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi when nugget effect = 0
  .temp.list <- get(".temp.list", pos=1)
  phi <- .temp.list$phi
  n <- .temp.list$n
  main <- .proflik.main(tausq=.temp.list$nugget, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half +
      (n/2) * log(main$ssresmat/n) + (n/2) - main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + ((n - .temp.list$beta.size)/2) * log(main$ssresmat/n) +
      (n/2) - 0.5 * sum(log(eigentrem$values)) - 
        main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}

".proflik.aux12" <-
  function(sigmasqphi.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the nugget parameter \tau^2, minimizing the likelihood wrt correlation function scale parameter \phi (range), the random field scale parameter \sigma^2 (sill) and the transformation parameter \lambda. 
  .temp.list <- get(".temp.list", pos=1)
  if(length(sigmasqphi.lambda) == 3) lambda <-  sigmasqphi.lambda[3]
  else lambda <- 1
  sigmasqphi.lambda <- as.vector(sigmasqphi.lambda)
  n <- .temp.list$n
  sigmasq <- sigmasqphi.lambda[1]
  phi <- sigmasqphi.lambda[2]
  main <- .proflik.main(tausq=.temp.list$nugget, sigmasq=sigmasq, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * (main$ssresmat) -
          main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * (main$ssresmat) -
          0.5 * sum(log(eigentrem$values)) -
            main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux13" <-
  function(sigmasqphi, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters \sigma^2 and \phi when the nugget is included.
  ## It requires the minimisation of the function wrt \tau^2 and \lambda (if the case) for each value of (\sigma^2, \phi)
  ## This is an auxiliary function called by likfit.proflik
  eval(substitute(.temp.list$sigmasqphi <-  xxx, list(xxx=as.vector(sigmasqphi))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE) {
##      tausq.res <- optim(.temp.list$tausq.est, .proflik.aux14, method="L-BFGS-B", lower
##			 = 0, ...)$value
      tausq.res <- optimise(.proflik.aux14, lower = 0, upper = .Machine$double.xmax^0.25, ...)$objective
  }
  else {
    tausq.res <- optim(
                       c(.temp.list$tausq.est, .temp.list$lambda), .proflik.aux14, method="L-BFGS-B",lower = c(0, -2
                                                                                                      ), upper = c( +Inf, 2), ...)$value
  }
  eval(expression(.temp.list$sigmasqphi <-  NULL), envir=.GlobalEnv)
  return( - tausq.res)
}

".proflik.aux14" <-
  function(tausq.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \phi), minimizing the likelihood wrt the nugget parameter \tau^2.
  ## This functions is called by the auxiliary function .proflik.aux13
  .temp.list <- get(".temp.list", pos=1)
  if(length(tausq.lambda) == 2) lambda <- tausq.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  tausq <- tausq.lambda[1]
  main <- .proflik.main(tausq=tausq, .temp.list$sigmasqphi[1], phi=.temp.list$sigmasqphi[2], lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) - main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux15" <-
  function(sigmasqtausq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters \sigma^2 and \tau^2
  ## It requires the minimisation of the function wrt \phi and also \lambda (if the case) for each value of (\sigma^2, \tau^2) 
  ## This is an auxiliary function called by likfit.proflik
  eval(substitute(.temp.list$sigmasqtausq <-  xxx, list(xxx=as.vector(sigmasqtausq))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      .proflik.aux16))
  ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
  if(.temp.list$fix.lambda == TRUE) {
#    phi.res <- optim(ini, .proflik.aux16, method="L-BFGS-B", lower = 
#                     .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    phi.res <- optimise(.proflik.aux16, lower = .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$objective
  }
  else {
    phi.res <- optim(ini, .proflik.aux16, method="L-BFGS-B", 
                     lower = c(.temp.list$lower.phi, -2),
                     upper = c(.temp.list$upper.phi, 2), ...)$value
  }
  eval(expression(.temp.list$sigmasqtausq <-  NULL), envir=.GlobalEnv)
  return( - phi.res)
}

".proflik.aux16" <-
  function(phi.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the sill and nugget parameters (\sigma^2,\tau^2), minimising the profile likelihood wrt correlation function scale parameter \phi (and the transformation parameter \lambda
  ## This is an auxiliary function called by likfit.aux15  
  .temp.list <- get(".temp.list", pos=1)
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  phi <- phi.lambda[1]
  main <- .proflik.main(tausq=.temp.list$sigmasqtausq[2], sigmasq=.temp.list$sigmasqtausq[1], phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) -
                                                       main$log.jacobian
  }
    return(as.vector(round(neglik, digits=8)))
}
".proflik.aux17" <-
  function(phitausq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2)
  ## It requires the minimisation of the function wrt \sigma^2 and \lambda (if the case) for each value of (\phi, \tau^2) 
  ## This is an auxiliary function called by likfit.proflik
  eval(substitute(.temp.list$phitausq <-  xxx, list(xxx=as.vector(phitausq))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE) {
##    ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
##                                        .proflik.aux18))
##    ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
##    sigmasq.res <- optim(ini, .proflik.aux18, method="L-BFGS-B", 
##                         lower = .temp.list$lower.sigmasq, ...)$value
    sigmasq.res <- optimise(.proflik.aux18, lower = .temp.list$lower.sigmasq, upper = .Machine$double.xmax^0.25, ...)$objective
  }
  else {
    sigmasq.res <- optim(c(.temp.list$sigmasq.est, .temp.list$lambda
                           ), .proflik.aux18, method="L-BFGS-B", lower = c(.temp.list$lower.sigmasq,
                                                                  -2), upper = c( + Inf, 2), ...)$value
  }
  eval(expression(.temp.list$phitausq <-  NULL), envir=.GlobalEnv)
  return( - sigmasq.res)
}

".proflik.aux18" <-
  function(sigmasq.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the range and nugget parameters (\phi, \tau^2), minimising the likelihood wrt the random field scale parameter \sigma^2 (sill) ant the transformation parameter \lambda. 
  ## This is an auxiliary function called by likfit.aux17.
  .temp.list <- get(".temp.list", pos=1)
  if(length(sigmasq.lambda) == 2) lambda <- sigmasq.lambda[2]
  else lambda <- 1
  n <- .temp.list$n
  sigmasq <- sigmasq.lambda[1]
  main <- .proflik.main(tausq=.temp.list$phitausq[2], sigmasq=sigmasq, phi=.temp.list$phitausq[1], lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + 0.5 *
      main$ssresmat - main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
    log.det.to.half + 0.5 * main$ssresmat - 0.5 * sum(log(
                                                     eigentrem$values)) -
                                                       main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux19" <-
  function(sigmasqtausq.rel, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \tau^2_{rel})
  ## It requires the minimisation of the function wrt \phi and \lambda (if the case) for each value of (\sigma^2, \tau^2_{rel})
  ## This is an auxiliary function called by likfit.proflik
  eval(substitute(.temp.list$sigmasqtausq.rel <-  xxx, list(xxx=as.vector(sigmasqtausq.rel))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE) {
##    phi.res <- optim(.temp.list$phi.est, .proflik.aux20, method="L-BFGS-B", lower = 
##                     .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    phi.res <- optimise(.proflik.aux20, lower = 
                     .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$objective
  }
  else {
    phi.res <- optim(c(.temp.list$phi.est, .temp.list$lambda, ...), .proflik.aux20, method="L-BFGS-B", 
                     lower = c(.temp.list$lower.phi, -2),
                     upper = c(.temp.list$upper.phi, 2), ...)$value
  }
 eval(expression( .temp.list$sigmasqtausq.rel <-  NULL), envir=.GlobalEnv)
  return( - phi.res)
}


".proflik.aux2" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the profile likelihood for the random field scale (variance) parameter \sigma^2 when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi and maybe \lambda for each value of \sigma^2
  ## This is an auxiliary function called by likfit.proflik
  ##
  eval(substitute(.temp.list$sigmasq <-  xxx, list(xxx=as.vector(sigmasq))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      .proflik.aux3))
  ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
  if(.temp.list$fix.lambda == TRUE) {
##    phi.res <- optim(ini , .proflik.aux3, method="L-BFGS-B",
 ##                    lower = .temp.list$lower.phi,
  ##                   upper=.temp.list$upper.phi, ...)$value
    phi.res <- optimise(.proflik.aux3,
                     lower = .temp.list$lower.phi,
                     upper=.temp.list$upper.phi, ...)$objective
  }
  else {
    phi.res <- optim(ini, .proflik.aux3, method="L-BFGS-B",
                     lower = c(.temp.list$lower.phi, -2),
                     upper = c(.temp.list$upper.phi, 2), ...)$value
  }
  eval(expression(.temp.list$sigmasq <-  NULL), envir=.GlobalEnv)
  return( - phi.res)
}

".proflik.aux20" <-
  function(phi.lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the sill and relative nugget parameters (\sigma^2, \tau^2_{rel}), minimising the likelihood wrt the correlation function scale parameter \phi and the transformation parameter \lambda.
  ## This is an auxiliary function called by likfit.aux19.
  .temp.list <- get(".temp.list", pos=1)
  phi.lambda <- as.vector(phi.lambda)
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2]
  else lambda <- 1
  sigmasqtausq.rel <- as.vector(.temp.list$sigmasqtausq.rel)
  sigmasq <- sigmasqtausq.rel[1]
  tausq.rel <- sigmasqtausq.rel[2]
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- .proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq) +
          (0.5/sigmasq) * main$ssresmat -
            main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq) +
          (0.5/sigmasq) * main$ssresmat -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  } 
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux21" <-
function(phitausq.rel, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2_{rel})
  ## This is an auxiliary function called by likfit.proflik
  .temp.list <- get(".temp.list", pos=1)
  phitausq.rel <- as.vector(phitausq.rel)
  phi <- phitausq.rel[1]
  tausq.rel <- phitausq.rel[2]
  n <- .temp.list$n
  main <- .proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = 1)
  sigmasq.hat <- main$ssresmat/n
  if(.temp.list$method.lik == "ML") {
    proflik <-  - (n/2) * log(2 * pi) -
      main$log.det.to.half -
      (n/2) * log(sigmasq.hat) -
        (n/2) -
          main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    proflik <-  - ((n - .temp.list$beta.size)/2) * log(2 * pi) -
      main$log.det.to.half -
        (n/2) * log(sigmasq.hat) -
          (n/2) +
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(proflik)
}

".proflik.aux21.1" <-
  function(lambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\phi, \tau^2_{rel})
  ## This requires minimasation wrt to the transformation parameter \lambda
  ## This is an auxiliary function called by likfit.proflik
  .temp.list <- get(".temp.list", pos=1)
  n <- .temp.list$n
  main <- .proflik.main(tausq = .temp.list$phitausq.rel[2], sigmasq = 1,
                       phi = .temp.list$phitausq.rel[1], lambda = lambda)
  sigmasq.hat <- main$ssresmat/n
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half + (n/2) * log(sigmasq.hat) +
        (n/2) -
          main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(sigmasq.hat) +
          (n/2) -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}

".proflik.aux22" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the range and nugget parameters (\phi, \tau^2), minimising the likelihood wrt the random field scale parameter \sigma^2 (sill) 
  ## This is an auxiliary function called by likfit.aux17
  .temp.list <- get(".temp.list", pos=1)
  n <- .temp.list$n
  main <- .proflik.main(tausq=.temp.list$phitausq[2], sigmasq=sigmasq, phi= .temp.list$phitausq[1], lambda = 1)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * main$ssresmat -
          main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        0.5 * main$ssresmat -
          0.5 * sum(log(eigentrem$values)) -
            main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}

".proflik.aux23" <-
  function(lambda, ...)
{
  ## This function computes the value of the profile likelihood for the transformation parameter \lambda
  ## It requires the minimisation of the function wrt \phi and \tau^2 and sigma^2 for each value of \lambda
  ## This is an auxiliary function called by proflik
  .temp.list <- get(".temp.list", pos=1)
  .temp.temp.list <- get(".temp.temp.list", pos=1)
  lambda <- as.vector(lambda)
  if(.temp.temp.list$fixtau == FALSE) {
    if(lambda == 0)
      data.l <- log(.temp.list$z)
    else data.l <- ((.temp.list$z^lambda) - 1)/lambda
    var.l <- var(data.l)
    ini.cov <- c(var.l, .temp.temp.list$ini)
  }
  else
    ini.cov <- .temp.temp.list$ini
  if(dim(.temp.list$xmat)[2] == 1 & all(.temp.list$xmat == 1))
    trend.mat <- "cte"
  else
    trend.mat <- ~ (.temp.list$xmat[,-1])
  lambda.res <- likfit(coords = .temp.temp.list$coords,
                       data = .temp.list$z,
                       ini.cov.pars = ini.cov, trend = trend.mat,
                       fix.nugget = .temp.temp.list$fixtau,
                       lik.method = .temp.list$method.lik,
                       cov.model = .temp.list$cov.model,
                       kappa = .temp.list$kappa, fix.lambda = TRUE,
                       lambda = lambda,
                       messages = FALSE)$loglik
  assign(".temp.list", .temp.temp.list, pos=1)
  return(lambda.res)
}

".proflik.aux24" <-
  function(sigmasqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the parameters (\sigma^2, \lambda) when there is no nugget effect (\tau^2 = 0, fixed)
  ## It requires the minimisation of the function wrt \phi for each value of (\sigma^2, \lambda)
  ## This is an auxiliary function called by proflik
  sigmasqlambda <- as.vector(sigmasqlambda)
  eval(substitute(.temp.list$sigmasq <-  xxx, list(xxx=sigmasqlambda[1])), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  lambda <- sigmasqlambda[2]
  if(lambda == 1) {
    eval(expression(.temp.list$log.jacobian <-  0), envir=.GlobalEnv)
  }
  else {
    if(any(.temp.list$z <= 0))
      stop("Transformation option not allowed when there are zeros or negative data"
           )
    eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=sum(log(.temp.list$z^(lambda - 1))))), envir=.GlobalEnv)
    if(lambda == 0)
      eval(substitute(.temp.list$z <-  xxx, list(xxx=log(.temp.list$z))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$z <- xxx, list(xxx= ((.temp.list$z^lambda) - 1)/lambda)), envir=.GlobalEnv)
  }
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      .proflik.aux3))
  ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
##  phi.res <- optim(ini, .proflik.aux3, method="L-BFGS-B", lower = .temp.list$
##                   lower.phi, upper = .temp.list$upper.phi, ...)$value
  phi.res <- optimise(.proflik.aux3, lower = .temp.list$lower.phi, upper = .temp.list$upper.phi, ...)$objective
  eval(expression(.temp.list$log.jacobian <-  NULL), envir=.GlobalEnv)
  eval(expression(.temp.list$sigmasq <-  NULL), envir=.GlobalEnv)
  eval(substitute(.temp.list$z <-  xxx, list(xxx=.temp.list$data)), envir=.GlobalEnv)
  return( - phi.res)
}

".proflik.aux27" <-
  function(sigmasqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for sill \sigma^2 and the transformation parameter \lambda
  ## It requires the minimisation of the function wrt \phi and \tau^2 and for each value of (\sigma^2,\lambda)
  ## This is an auxiliary function called by .proflik.
  sigmasqlambda <- as.vector(sigmasqlambda)
  eval(substitute(.temp.list$sigmasq <-  xxx, list(xxx=sigmasqlambda[1])), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  lambda <- sigmasqlambda[2]
  if(lambda == 1) {
    eval(expression(.temp.list$log.jacobian <-  0), envir=.GlobalEnv)
  }
  else {
    eval(expression(.temp.list$fix.lambda <-  TRUE), envir=.GlobalEnv)
    if(any(.temp.list$z^(lambda - 1) <= 0))
      eval(substitute(.temp.list$log.jacobian <- xxx, list(xxx=log(prod(.temp.list$z^(lambda -
                                                         1))))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=sum(log(.temp.list$z^(lambda -
                                                           1))))), envir=.GlobalEnv)
    if(lambda == 0)
      eval(substitute(.temp.list$z <-  xxx, list(xxx=log(.temp.list$z))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$z <-  xxx, list(xxx=((.temp.list$z^lambda) - 1)/lambda)), envir=.GlobalEnv)
  }
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      .proflik.aux10))
  ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
  phitausq.rel.res <- optim(ini, .proflik.aux10, method="L-BFGS-B",
                            lower = c(.temp.list$lower.phi,
                              0), upper=c(.temp.list$upper.phi, 100), ...)$value
  eval(expression(.temp.list$log.jacobian <-  NULL), envir=.GlobalEnv)
  eval(expression(.temp.list$sigmasq <-  NULL), envir=.GlobalEnv)
  eval(substitute(.temp.list$z <-  xxx, list(xxx=.temp.list$data)), envir=.GlobalEnv)
  return( - phitausq.rel.res)
}

".proflik.aux28" <-
  function(sigmasqphi, ...)
{
  ## This function computes the value of the 2-D profile likelihood for the random field scale (variance) parameter \sigma^2  and the correlation function parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \lambda for each value of (\sigma^2, \phi)
  ## This is an auxiliary function called by likfit.proflik
  ##
##  ini.seq <- seq(-1.5, 1.5, l=7)
##  .temp.list$sigmasqphi <<- as.vector(sigmasqphi)
##  lambda.lik <- apply(as.matrix(ini.seq), 1, .proflik.aux4)
##  ini <- ini.seq[lambda.lik == max(lambda.lik)]
##  lambda.res <- optim(ini, .proflik.aux4, method="L-BFGS-B", lower = -2.5, upper = 2.5, ...)$value
  lambda.res <- optimise(.proflik.aux4, lower = -5, upper = 5, ...)$objective
  eval(expression(.temp.list$sigmasqphi <-  NULL), envir=.GlobalEnv)
  return( - lambda.res)
}

".proflik.aux30" <-
  function(phitausq.rel, ...)
{
  ## This function computes the value of the profile likelihood for the correlation parameter \phi when nugget effect is not included in the model.
  ## It requires the minimisation of the function wrt \phi for each value of \lambda, if this transformation parameter is included in the model
  ## This is an auxiliary function called by likfit.proflik
  ##
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE)
    proflik <- .proflik.aux21(phitausq.rel = phitausq.rel)
  else {
    eval(substitute(.temp.list$phitausq.rel <-  xxx, list(xxx=phitausq.rel)), envir=.GlobalEnv)
##    proflik <-  - (optim(.temp.list$lambda, .proflik.aux21.1, method="L-BFGS-B", lower =
##                         -2, upper = 2, ...)$value)
    proflik <-  - (optimise(.proflik.aux21.1, lower = -5, upper = 5, ...)$objective)
    eval(expression(.temp.list$phitausq.rel <-  NULL), envir=.GlobalEnv)
  }
  return(proflik)
}


".proflik.aux3" <-
  function(phi.lambda, ...)
{
  ## This function computer the negative of the likelihood function for the correlation function scale parameter \phi (and maybe the transformation parameter \lambda) only for models with fixed nugget effect (i.e., when it is not a parameter to be estimated) 
  ## This function is used when computing the profile likelihood for \sigma^2
  ## This is an auxiliary function called by .proflik.aux2
  ##  phi <- pmax(phi, .temp.list$lower.phi)
  .temp.list <- get(".temp.list", pos=1)
  if(length(phi.lambda) == 2)
    lambda <- phi.lambda[2]
  else lambda <- 1
  sigmasq <- .temp.list$sigmasq
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- .proflik.main(tausq=0, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
      (n/2) * log(sigmasq) + (0.5/sigmasq) * main$ssresmat - 
        main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        ((n - .temp.list$beta.size)/2) * log(sigmasq) +
          (0.5/sigmasq) * main$ssresmat -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}

".proflik.aux31" <-
  function(philambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for range \phi and the transformation parameter \lambda.
  ## It requires the minimisation of the function wrt \tau^2_{rel} and for each value of (\phi,\lambda).
  ## This is an auxiliary function called by .proflik.
  .temp.list <- get(".temp.list", pos=1)
  philambda <- as.vector(philambda)
  eval(substitute(.temp.list$phi <- xxx, list(xxx= philambda[1])), envir=.GlobalEnv)
  .temp.list$lambda <- philambda[2]
##  tausq.rel.res <- optim(.temp.list$tausq.rel.est, .proflik.aux8, method="L-BFGS-B", lower = 
##                         0, upper=100, ...)$value
  tausq.rel.res <- optimise(.proflik.aux8, lower =  0, upper=1000, ...)$objective
  eval(expression(.temp.list$phi <-  NULL), envir=.GlobalEnv)
  return( - tausq.rel.res)
}

".proflik.aux32" <-
  function(tausqlambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for nugget \tau^2 and the transformation parameter \lambda.
                                        # It requires the minimisation of the function wrt \phi and \sigma^2 and for each value of (\tau^2,\lambda).
                                        # This is an auxiliary function called by .proflik.
  tausqlambda <- as.vector(tausqlambda)
  eval(substitute(.temp.list$nugget <-  xxx, list(xxx=tausqlambda[1])), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  lambda <- tausqlambda[2]
  if(lambda == 1) {
    eval(expression(.temp.list$log.jacobian <-  0), envir=.GlobalEnv)
  }
  else {
    if(any(.temp.list$z^(lambda - 1) <= 0))
      eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=log(prod(.temp.list$z^(lambda -
                                                         1))))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=sum(log(.temp.list$z^(lambda - 1))))), envir=.GlobalEnv)
    if(lambda == 0)
      eval(substitute(.temp.list$z <-  xxx, list(xxx=log(.temp.list$z))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$z <-  xxx, list(xxx=((.temp.list$z^lambda) - 1)/lambda)), envir=.GlobalEnv)
  }
  sigmasqphi.res <- optim(c(.temp.list$sigmasq.est, .temp.list$phi.est), .proflik.aux12, method="L-BFGS-B",
                          lower = c(.temp.list$lower.sigmasq, .temp.list$
                            lower.phi), upper=c(+Inf, .temp.list$upper.phi), ...)$value
  eval(expression(.temp.list$log.jacobian <-  NULL), envir=.GlobalEnv)
  eval(substitute(.temp.list$z <-  xxx, list(xxx=.temp.list$data)), envir=.GlobalEnv)
  eval(expression(.temp.list$nugget <-  NULL), envir=.GlobalEnv)
  return( - sigmasqphi.res)
}

".proflik.aux33" <-
  function(tausq.rellambda, ...)
{
  ## This function computes the value of the 2-D profile likelihood for nugget \tau^2 and the transformation parameter \lambda.
  ## It requires the minimisation of the function wrt \phi for each value of (\tau^2,\lambda).
  ## This is an auxiliary function called by .proflik.
  tausq.rellambda <- as.vector(tausq.rellambda)
  eval(substitute(.temp.list$nugget.rel <-  xxx, list(xxx=tausq.rellambda[1])), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  lambda <- tausq.rellambda[2]
  if(lambda == 1) {
    eval(expression(.temp.list$log.jacobian <-  0), envir=.GlobalEnv)
  }
  else {
    if(any(.temp.list$z^(lambda - 1) <= 0))
      eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=log(prod(.temp.list$z^(lambda -
                                                         1))))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$log.jacobian <-  xxx, list(xxx=sum(log(.temp.list$z^(lambda -
                                                           1))))), envir=.GlobalEnv)
    if(lambda == 0)
      eval(substitute(.temp.list$z <-  xxx, list(xxx=log(.temp.list$z))), envir=.GlobalEnv)
    else eval(substitute(.temp.list$z <-  xxx, list(xxx=((.temp.list$z^lambda) - 1)/lambda)), envir=.GlobalEnv)
  }
##  phi.res <- optim(.temp.list$phi.est, .proflik.aux6, method="L-BFGS-B", lower = .temp.list$
##                    lower.phi, upper=.temp.list$upper.phi, ...)$value
  phi.res <- optimise(.proflik.aux6, lower = .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$objective
  eval(expression(.temp.list$log.jacobian <-  NULL), envir=.GlobalEnv)
  eval(expression(.temp.list$nugget.rel <-  NULL), envir=.GlobalEnv)
  eval(substitute(.temp.list$z <-  xxx, list(xxx=.temp.list$data)), envir=.GlobalEnv)
  return( - phi.res)
}

".proflik.aux4" <-
  function(lambda, ...)
{
  ## This function computer the values of the profile likelihood function for the parameters \phi  and \sigma^2 for models with nugget effect = 0, including the tranformation parameter \lambda
  ## This is an auxiliary function called by .proflik.aux28
  ##
  .temp.list <- get(".temp.list", pos=1)
  sigmasqphi <- as.vector(.temp.list$sigmasqphi)
  sigmasq <- sigmasqphi[1]
  phi <- sigmasqphi[2]
  n <- .temp.list$n
  if(lambda > 0.999 & lambda < 1.001)
    lambda <- 1
  main <- .proflik.main(tausq=.temp.list$nugget, sigmasq = sigmasq, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- ((n/2) * log(2 * pi) +
               main$log.det.to.half +
               0.5 * main$ssresmat - 
               main$log.jacobian)
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE,
                      only.values = TRUE)
    neglik <- (((n - .temp.list$beta.size)/2) * log(2 * pi) -
               0.5 * sum(log(xx.eigen$values)) + main$log.det.to.half +
               (0.5) * main$ssresmat + 0.5 * sum(log(eigentrem$values)) + main$log.jacobian)
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux5" <-
  function(tausq.rel, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \tau^2_{rel}.
  ## It requires the minimisation of the function wrt \phi and \lambda (if the case) for each value of \tau^2_{rel}.
  ## This is an auxiliary function called by .proflik.
  eval(substitute(.temp.list$nugget.rel <-  xxx, list(xxx=as.vector(tausq.rel))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda == TRUE) {
##    phi.res <- optim(.temp.list$phi.est, .proflik.aux6, method="L-BFGS-B", lower = 
##                     .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$value
    phi.res <- optimise(.proflik.aux6, lower = 
                     .temp.list$lower.phi, upper=.temp.list$upper.phi, ...)$objective
  }
  else {
    phi.res <- optim(c(.temp.list$phi.est, .temp.list$lambda), .proflik.aux6, method="L-BFGS-B", 
                       lower = c(.temp.list$lower.phi, -2),
                       upper = c(.temp.list$upper.phi, 2), ...)$value
  }
  eval(expression(.temp.list$nugget.rel <-  NULL), envir=.GlobalEnv)
  return( - phi.res)
}

".proflik.aux6" <-
function(phi.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the relative nugget parameter \tau^2_{rel}, minimizing the likelihood wrt correlation function scale parameter \phi (range) and the transformation parameter \lambda.
  .temp.list <- get(".temp.list", pos=1)
  if(length(phi.lambda) == 2) lambda <- phi.lambda[2] else lambda <- 1
  phi.lambda <- as.vector(phi.lambda)
  phi <- phi.lambda[1]
  n <- .temp.list$n
  main <- .proflik.main(tausq=.temp.list$nugget.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) +
      main$log.det.to.half +
        (n/2) * log(main$ssresmat/n) +
          (n/2) -
            main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- ((n - .temp.list$beta.size)/2) * log(2 * pi) +
      main$log.det.to.half +
        ((n - .temp.list$beta.size)/2) * log(main$ssresmat/n) +
          (n/2) -
            0.5 * sum(log(eigentrem$values)) -
              main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux7" <-
  function(phi, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \phi when the nugget \tau^2 is included in the model
  ## It requires the minimisation of the function wrt relative \tau^2_{rel} for each value of \phi
  ## This is an auxiliary function called by .proflik.
  eval(substitute(.temp.list$phi <-  xxx, list(xxx=as.vector(phi))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  if(.temp.list$fix.lambda) {
    eval(expression(.temp.list$lambda <-  1), envir=.GlobalEnv)
##    tausq.rel.res <- optim(.temp.list$tausq.rel.est, .proflik.aux8, method="L-BFGS-B", 
##                           lower = 0, upper=100, ...)$value
    tausq.rel.res <- optimise(.proflik.aux8, lower = 0, upper=1000, ...)$objective
    eval(expression(.temp.list$lambda <-  NULL), envir=.GlobalEnv)
  }
  else {
    tausq.rel.res <- optim(c(.temp.list$tausq.rel.est, .temp.list$lambda), .proflik.aux8, method="L-BFGS-B", lower = c(0, -2), upper = c(100, 2), ...)$value
  }
  eval(expression(.temp.list$phi <-  NULL), envir=.GlobalEnv)
  return( - tausq.rel.res)
}

".proflik.aux8" <-
  function(tausq.rel.lambda, ...)
{
  ## This function computes the value of the profile likelihood for the correlation function scale parameter \phi (and lambda), minimizing the likelihood wrt relative nugget parameter \tau^2_{rel}
  .temp.list <- get(".temp.list", pos=1)
  if(length(tausq.rel.lambda) == 2)
    lambda <- tausq.rel.lambda[2]
  else lambda <- .temp.list$lambda
  n <- .temp.list$n
  phi <- .temp.list$phi
  tausq.rel <- tausq.rel.lambda[1]
  main <- .proflik.main(tausq=tausq.rel, sigmasq=1, phi=phi, lambda = lambda)
  if(.temp.list$method.lik == "ML") {
    neglik <- (n/2) * log(2 * pi) + main$log.det.to.half + (
                                                              n/2) * log(main$ssresmat/n) + (n/2) - main$log.jacobian
  }
  if(.temp.list$method.lik == "RML") {
    eigentrem <- eigen(main$ixix, symmetric = TRUE, only.values = TRUE)
    neglik <- (((n - .temp.list$beta.size)/2) * log(2 * pi) + main$
               log.det.to.half + ((n - .temp.list$beta.size)/2) * log(main$ssresmat/
                                                           n) + (n/2) - 0.5 * sum(log(eigentrem$values))) - 
                                                             main$log.jacobian
  }
  return(as.vector(round(neglik, digits=8)))
}
".proflik.aux9" <-
  function(sigmasq, ...)
{
  ## This function computes the value of the profile likelihood for the parameter \sigma^2 when \tau^2 is included in the model
  ## It requires the minimisation of the function wrt \phi and \tau^2 for each value of \sigma^2
  ## This is an auxiliary function called by likfit.proflik
  eval(substitute(.temp.list$sigmasq <-  xxx, list(xxx=as.vector(sigmasq))), envir=.GlobalEnv)
  .temp.list <- get(".temp.list", pos=1)
  ini.lik <- round(100000000. * apply(.temp.list$ini.grid, 1,
                                      .proflik.aux10))
  ini <- as.vector(.temp.list$ini.grid[which(ini.lik == min(ini.lik, na.rm = TRUE)),,drop=FALSE][1,])
  if(.temp.list$fix.lambda == TRUE) {
    phitausq.rel.res <- optim(ini, .proflik.aux10, method="L-BFGS-B",
                              lower = c(.temp.list$
                                lower.phi, 0),
                              upper=c(.temp.list$upper.phi, 100), ...)$value
  }
  else {
    if(ini[2] == 0) ini[2] <- 0.01
    phitausq.rel.res <- optim(ini, .proflik.aux10, method="L-BFGS-B", 
                                lower = c(.temp.list$lower.phi, 0,-2),
                                upper = c(.temp.list$upper.phi, 100, 2), ...)$value
  }
  eval(expression(.temp.list$sigmasq <- NULL), envir=.GlobalEnv)
  return( - phitausq.rel.res)
}

"plot.proflik" <-
  function(x, pages = c("user", "one", "two"),
           uni.only, bi.only, type.bi = c("contour", "persp"),
           conf.int = c(0.90,0.95),
           yaxis.lims = c("conf.int", "as.computed"),
           by.col = TRUE, log.scale = FALSE, use.splines = TRUE,
           par.mar.persp = c(0, 0, 0, 0), ask = FALSE, 
           ...)
{
  ##
  ## Saving original par() parameters
  ##
#  if (is.R()) 
#    par.ori <- par(no.readonly = TRUE)
#  else par.ori <- par()
#  on.exit(par(par.ori))
  ##
  parask <- par()$ask
  par(ask = ask)
  on.exit(par(ask = parask))  
  ##
  ## Checking whether to plot 1D and/or 2D profiles
  ##
  if(missing(uni.only)){
    if(x$n.bi == 0) uni.only <- TRUE
    else uni.only <- FALSE
  }
  if(missing(bi.only)){
    if(x$n.uni == 0) bi.only <- TRUE
    else bi.only <- FALSE
  }
  if(!is.logical(uni.only))
    stop("argument uni.only must be logical (TRUE or FALSE)")
  if(!is.logical(bi.only))
    stop("argument bi.only must be logical (TRUE or FALSE)")
  n.uni <- x$n.uni
  n.bi <- x$n.bi
  if((uni.only == FALSE) & (bi.only == FALSE))
    np <- n.uni + n.bi
  if((uni.only == TRUE) & (bi.only == FALSE))
    np <- n.uni
  if((uni.only == FALSE) & (bi.only == TRUE))
    np <- n.bi
  if(n.uni == 0 & np > 0) bi.only <- TRUE
  if(n.bi == 0 & np > 0) uni.only <- TRUE
  ##
  ##
  ##
  if(mode(yaxis.lims) == "character")
    yaxis.lims <- match.arg(yaxis.lims)
  type.bi <- match.arg(type.bi)
  ##
  ## Definig number of pages to place the plots
  ##
  pages <- match.arg(pages)
  if(pages == "one") {
    if(np >= 1 & np < 4)
      par(mfrow = c(np, 1))
    if(np >= 4) {
      if(by.col == TRUE)
        par(mfcol = c(ceiling(np/2), 2))
      else par(mfrow = c(ceiling(np/2), 2))
    }
  }
  if(pages == "two") {
    if(n.uni > 1 & n.uni < 4)
      par(mfrow = c(n.uni, 1))
    if(n.uni >= 4)
      par(mfrow = c(ceiling(n.uni/2), 2))
  }
  ##
  ##
  ##
  if(bi.only == FALSE) {
    for(i in 1:n.uni) {
      if(x$method.lik == "ML")
        ylabm <- "profile log-likelihood"
      else ylabm <- "profile log-(restricted) likelihood"
      if(!is.logical(conf.int) || all(conf.int) != FALSE) {
        if(mode(conf.int) != "numeric"| any(conf.int > 1))
          stop("argument conf.int must be numerical (scalar or vector) with values between 0 and 1")
        conf.int.drop <- x[[i]][[3]][2] - 0.5 * qchisq(conf.int,1)
      }
      if(all(is.character(yaxis.lims))){
        if(yaxis.lims == "conf.int")
          lik.lims <- c(min(conf.int.drop), 
                        x[[i]][[3]][2])
        else lik.lims <- c(min(x[[i]][[2]]),
                           x[[i]][[3]][2])
      }
      else
        lik.lims <- yaxis.lims
      if(log.scale == TRUE) {
        if(use.splines){
          nxpoints <- 5*length(x[[i]][[1]])
          nodups <- which(duplicated(x[[i]][[1]]) == FALSE)
          plot(spline(x = log(x[[i]][[1]][nodups]), 
                      y = x[[i]][[2]][nodups],
                      n = nxpoints,
                      method="natural"), type = "l",
               xlab = paste("log-",
                 .proflik.plot.aux1(names(x[[i]])[1])),
               ylab = ylabm, ylim = lik.lims)
        }
        else{
          plot(log(x[[i]][[1]]), 
               x[[i]][[2]],
               type = "l",
               xlab = paste("log-",
                 .proflik.plot.aux1(names(x[[i]])[1])),
               ylab = ylabm, ylim = lik.lims)
        }
        lines(log(c(x[[i]][[3]][1], x[[i]][[3]][1])),
              c(min(lik.lims), x[[i]][[3]][2]), lty = 2)
      }
      else {
        if(use.splines){
          nxpoints <- 5*length(x[[i]][[1]])
          nodups <- which(duplicated(x[[i]][[1]]) == FALSE)
          plot(spline(x = x[[i]][[1]][nodups],
                      y = x[[i]][[2]][nodups],
                      n = nxpoints,
                      method="natural"),
               type = "l",
               xlab = .proflik.plot.aux1(names(x[[i]])[1]),
               ylab = ylabm, ylim = lik.lims)
        }
        else{
          plot(x[[i]][[1]],
               x[[i]][[2]],
               type = "l", xlab = 
               .proflik.plot.aux1(names(x[[i]])[1]),
               ylab = ylabm, ylim = lik.lims)
        }
        lines(c(x[[i]][[3]][1], 
                x[[i]][[3]][1]),
              c(min(lik.lims), x[[
                                            i]][[3]][2]), lty = 2)
      }
      abline(h = conf.int.drop, lty = 3)
    }
  }
  if(uni.only == FALSE) {
    if(pages == "two") {
      if(n.bi >= 1 & n.bi < 4)
        par(mfrow = c(n.bi, 1))
      if(n.bi >= 4)
        par(mfrow = c(ceiling(n.bi/2), 2))
    }
    for(i in 1:n.bi) {
      if(type.bi == "contour") {
        if(log.scale == TRUE) {
          contour(log(x[[(n.uni + i)]][[1]]),
                  log(x[[(n.uni + i)]][[2]]),
                  matrix(x[[(n.uni + i)]][[3]],
                         ncol = length(x[[(n.uni +i)]][[2]])),
                  xlab = paste("log-", .proflik.plot.aux1(names(x[[(n.uni + i)]][1]))),
                  ylab = paste("log-", .proflik.plot.aux1(names(x[[(n.uni + i)]][2]))),
                    ...)
          points(log(t(x[[(n.uni + i)]][[4]][1:2])))
        }
        else {
          contour(x[[(n.uni + i)]][[1]],
                  x[[(n.uni + i)]][[2]],
                  matrix(x[[(n.uni + i)]][[3]],
                         ncol = length(x[[(n.uni + i)]][[2]])),
                  xlab = .proflik.plot.aux1(names(x[[(n.uni + i)]][1])),
                  ylab = .proflik.plot.aux1(names(x[[(n.uni + i)]][2])),
                  ...)
          points(t(x[[(n.uni + i)]][[4]][1:2]))
        }
      }
      if(type.bi == "persp") {
        cat("For better visualisation arguments for the funtion `persp` can be passed.\nSome relevant argments are: theta, phi, r, d, among others.\n Type help(persp) for a description of the options\n")
        if(x$method.lik == "ML")
          zlabm <- 
            "profile log-likelihood"
        else zlabm <- "profile log-(restricted) likelihood"
        zlimm <- range(x[[(n.uni +
                                     i)]][[3]])
        zlimm[1] <- 1.01 * zlimm[1]
        minlik <- min(x[[(n.uni + i)]][[3]])
        if(log.scale == TRUE) {
          persp(log(x[[(n.uni + i)]][[1]]),
                log(x[[(n.uni + i)]][[2]]),
                matrix(x[[(n.uni + i)]][[3]],
                       ncol = length(x[[(n.uni + i)]][[2]])),
                xlab = .proflik.plot.aux1(paste("log-", names(x[[(n.uni + i)]][1]))),
                ylab = paste("log-", .proflik.plot.aux1(names(x[[(n.uni + i)]][2]))),
                zlab = zlabm, box = TRUE, ...)
                                        #          pp1 <- perspp(x = log(c(x[[(n.uni + i)]][[4]][1],
                                        #                          min( x[[(n.uni + i)]][[1]]))[c(1, 1, 1, 2)]),
                                        #                        y = log(c(x[[(n.uni + i)]][[4]][2],
                                        #                          min(x[[(n.uni + i)]][[2]]))[c(1, 1, 2, 1)]),
                                        #                        z = c(minlik, x[[(n.uni + i)]][[4]][3])[c(1, 2, 1, 1)], pp)
                                        #          segments(log(pp1$x[1]), log(pp1$y[1]), log(pp1$x[2]), log(pp1$y[2]),
                                        #                   lwd = 2)
        }
        else {
          persp(x = x[[(n.uni + i)]][[1]],
                y = x[[(n.uni + i)]][[2]],
                z = matrix(x[[(n.uni + i)]][[3]],
                  ncol = length(x[[(n.uni + i)]][[2]])),
                xlab = .proflik.plot.aux1(names(x[[(n.uni + i)]][1])),
                ylab = .proflik.plot.aux1(names(x[[(n.uni + i)]][2])),
                zlab = zlabm, box = TRUE, ...)
                                        #          pp1 <- perspp(x = c(x[[(n.uni + i)]][[4]][1],
                                        #                          min(x[[(n.uni + i)]][[1]]))[c(1, 1, 1, 2)],
                                        #                        y = c(x[[(n.uni + i)]][[4]][2],
                                        #                          min(x[[(n.uni + i)]][[2]]))[c(1, 1,2, 1)],
                                        #                        z = c(minlik, x[[(n.uni + i)]][[4]][3])[c(1, 2, 1, 1)], pp)
                                        #          segments(pp1$x[1], pp1$y[1], pp1$x[2], pp1$y[2], lwd = 2)
        }
      }
    }
  }
  return(invisible())
}

".proflik.plot.aux1" <-
  function(parameter.name)
{
  switch(parameter.name,
         range = expression(phi),
         sill = expression(sigma^2),
         lambda = expression(lambda),
         nugget = expression(tau^2),
         nugget.rel = expression(tau[rel]^2))
}

".proflik.main" <-
  function(tausq, sigmasq, phi, lambda)
{
  .temp.list <- get(".temp.list", pos=1)
  z <- .temp.list$z
  n <- .temp.list$n
  if(lambda == 1){
    ##    log.jacobian <- .temp.list$log.jacobian
    log.jacobian <- 0
  }
  else {
    if(any(z <= 0))
      stop("Transformation option not allowed when there are zeros or negative data")
    if(any(z^(lambda - 1) <= 0))
      log.jacobian <- log(prod(z^(lambda - 1)))
    else log.jacobian <- sum(log(z^(lambda - 1)))
    if(lambda == 0)
      z <- log(z)
    else z <- ((z^lambda) - 1)/lambda
  }
  kappa <- .temp.list$kappa
  V <- varcov.spatial(dists.lowertri = .temp.list$dists.lowertri,
                      cov.model = .temp.list$cov.model, kappa = kappa,
                      nugget = tausq, cov.pars = c(sigmasq, phi),
                      det = TRUE)  
  xix <- crossprod(.temp.list$xmat,solve(V$varcov,.temp.list$xmat))
  if(length(as.vector(xix)) == 1) choldet <- 0.5 * log(xix)
  else choldet <- sum(log(diag(chol(xix))))
  iz <- solve(V$varcov,z)
  xiy <- crossprod(.temp.list$xmat,iz)
  beta.hat <- drop(solve(xix,xiy))
  yiy <- drop(crossprod(z,iz))
  ssresmat <- as.vector(yiy - crossprod(beta.hat,xiy))
  return(list(log.det.to.half = V$log.det.to.half,
              ssresmat = ssresmat,
              ixix = solve(xix), log.jacobian = log.jacobian))
}


".proflik.lambda" <-
function(lambda)
{
  .temp.list <- get(".temp.list", pos=1)
  .temp.lower.lambda <- get(".temp.lower.lambda", pos=1)
  .temp.upper.lambda <- get(".temp.upper.lambda", pos=1)
  if (any(is.na(lambda)) | any(lambda==Inf) | any(is.nan(lambda)))
    neglik <- 1e+32
  else{
    if(.temp.list$minimisation.function == "nlm"){
      if (exists(".temp.lambda", where=1)) remove(".temp.lambda", pos=1, inherits = TRUE)
      lambda.minimiser <- lambda
      penalty <-  1000 * (.temp.lower.lambda - min(lambda, .temp.lower.lambda))
      lambda <- max(lambda, .temp.lower.lambda)
      penalty <- penalty + 1000 * (.temp.upper.lambda - max(lambda, .temp.upper.lambda))
      lambda <- min(lambda, .temp.upper.lambda)
      if (round(1000 * lambda.minimiser) <= round(1000 * .temp.lower.lambda))
        assign(".temp.lambda", lambda, pos=1)
      if (round(1000 * lambda.minimiser) >= round(1000 * .temp.upper.lambda))
        assign(".temp.lambda", lambda, pos=1)
    }
    z <- .temp.list$z
    n <- .temp.list$n
    if(lambda == 1) {
      .temp.list$log.jacobian <- 0
    }
    else {
      if(any(z < 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      if(any(z^(lambda - 1) <= 0))
        .temp.list$log.jacobian <- log(prod(z^(lambda - 1)))
      else .temp.list$log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
    beta.size <- .temp.list$beta.size
    kappa <- .temp.list$kappa
    xmat <- .temp.list$xmat
    txmat <- .temp.list$txmat
    ixx <- solve(crossprod(xmat))
    tausqhat <- (z %*% (diag(n) - xmat %*% ixx %*% txmat) %*% z)/n
    if(.temp.list$method == "ML")
      neglik <- ((n/2) * log(2 * pi) +
                 (n/2) * log(tausqhat) +
                 (n/2) -
                 .temp.list$log.jacobian
                 )
    if(.temp.list$method == "RML") {
      eigentrem <- eigen(ixx, symmetric = TRUE, only.values = TRUE)
      neglik <- (((n - beta.size)/2) * log(2 * pi) +
                 ((n - beta.size)/2) * log(tausqhat) +
                 (n/2) -
                 0.5 * sum(log(eigentrem$values)) -
                 .temp.list$log.jacobian
                 )
    }
  }
  if(.temp.list$minimisation.function == "nlm")
    return(as.vector(neglik + penalty))
  else
    return(as.vector(neglik))
}




