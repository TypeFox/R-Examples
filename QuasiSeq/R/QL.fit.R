
QL.fit <- function(counts, design.list, test.mat = NULL, log.offset = NULL, Model = "NegBin", print.progress = TRUE, 
                   NBdisp = "trend", bias.fold.tolerance=1.10, ...) {
  if(is.data.frame(counts)) counts<-as.matrix(counts)
  ### Note: First element of design.list should pertain to overall full model.  This is the design used to obtain
  ### dispersion estimates for quasi-likelihood models.
  if(any(rowSums(counts)==0)) stop(paste(sum(rowSums(counts)==0)," genes have 0 counts across all samples.  Please remove genes with zero total counts before analyzing.",sep=""))
  ### Check for errors
  if (any(round(counts) != counts)) 
    stop("Count data contains non-integers.")
  if (any(counts < 0)) 
    stop("Count data contains negative counts.")
  
  if (!Model %in% c("NegBin", "Poisson")) 
    stop("Unidentified Model: Model must be either 'NegBin' or 'Poisson'.")
  if (Model == "NegBin") {
    
    if (!NBdisp[1] %in% c("trend", "common") & length(NBdisp) != nrow(counts)) 
      stop("Unidentified NegBin Dispersion: NBdisp must be set as 'trend' or 'common' to estimate negative binomial dispersion from data using GLM edgeR (McCarthy et al., 2012),\n or it must be a vector providing negative binomial dispersion parameter value to use for each gene.")
    
    if (length(NBdisp) == nrow(counts) & !is.numeric(NBdisp)) 
      stop("NBdisp contains non-numeric values.\n\tAll negative binomial dispersion parameters must be non-negative real numbers.")
    
    if (length(NBdisp) == nrow(counts) & any(NBdisp < 0)) 
      stop("NBdisp contains negative values.\nAll negative binomial dispersion parameters must be non-negative real numbers.")
  }
  
  ### Fit model and evaluate deviance under each design provided in design.list
  deviance.list <- vector("list", length(design.list))
  p <- NULL
  n <- ncol(counts)  # p is used to store the d.f. for each model (it will be a vector) and n is the total number of samples (also the number of observations for each gene)
  
  for (jj in 1:length(design.list)) {
    design <- design.list[[jj]]
    
    if (is.vector(design)) {
      p <- c(p, length(unique(design)))  ## Record the d.f. for current model
      
      ### Check for errors if the current model design is specified as a vector
      if (p[jj] > p[1]) 
        stop(paste("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design."))
      if (length(design) != n) 
        stop(paste("Element", jj, "in 'design.list' has length", length(design), ".\nDesign vectors must have length", 
                   n, "(to match number of columns in data)."))
    }
    
    if (is.matrix(design)) {
      p <- c(p, ncol(design))
      
      ### Check for errors if the current model design is specified as a matrix 
      ### if(prod(design[,1]==1)!=1)
      ### stop(paste('The first column of matrix in element',jj,'of 'design.list' is not a column of 1s for the
      ### intercept. Please include intercept.'))
      if (nrow(design) != n) 
        stop(paste("Element", jj, "in 'design.list' has", nrow(design), "rows.\nDesign matrices must have", 
                   n, "rows (to match number of columns in data)."))
      if (p[jj] > p[1]) 
        stop(paste("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design."))
    }
    
    ### Analyze using quasi-negative binomial model, if chosen
    if (Model == "NegBin") {
      if (is.vector(design)) {
        if (length(unique(design)) > 1) 
          design <- model.matrix(~as.factor(design))
        if (length(unique(design)) == 1) 
          design <- matrix(1, ncol(counts), 1)
      }
      if (jj == 1) {
        ### Obtain negative binomial dispersion estimates from edgeR, if requested
        if(NBdisp %in% c("trend", "common")){
          
          ### Default to using edgeR's suggested normalization routine
          if (is.null(log.offset)) { 
            d <- DGEList(counts = counts, group = design[, 2])
            d <- calcNormFactors(d)
          }
          
          ### If provided, use prespecified library size normalizations instead of edgeR's normalization routine
          if (!is.null(log.offset)) 
            d <- DGEList(counts = counts, group = design[, 2], lib.size = exp(log.offset))
          
          ### If requested, use gene-specific trended dispersion estimates from GLM edgeR (McCarthy et al., 2012).
          if (NBdisp == "trend") 
            nb.disp <- estimateGLMTrendedDisp(d, design)$trended.dispersion
          
          ### If requested, use common dispersion estimate from GLM edgeR (McCarthy et al., 2012).
          if (NBdisp == "common") 
            nb.disp <- rep(estimateGLMCommonDisp(d, design, ...)$common.dispersion, nrow(counts))
        }
        ### If provided, use prespecified dispersion estimates.
        if (length(NBdisp) == nrow(counts)) {
          if (is.numeric(NBdisp) & !any(NBdisp < 0)) 
            nb.disp <- NBdisp
        }
      }
      
      ### Analyze genes with positive dispersion parameters using quasi-negative binomial model
      if (any(nb.disp > 0)) 
        res <- NBDev(counts[nb.disp > 0, ], design, log.offset, nb.disp[nb.disp > 0], print.progress, bias.fold.tolerance=bias.fold.tolerance)
      
      ### If present, analyze genes for which nb.disp==0 using quasi-Poisson model
      if (any(nb.disp == 0)) {
        res2 <- PoisDev(counts[nb.disp == 0, ], design, log.offset, print.progress)
        means <- dev <- rep(NA, nrow(counts))
        parms <- matrix(NA, nrow(counts), ncol(design))
        means[nb.disp == 0] <- res2$means
        dev[nb.disp == 0] <- res2$dev
        parms[nb.disp == 0, ] <- res2$parms
        if (any(nb.disp > 0)) {
          means[nb.disp > 0] <- res$means
          dev[nb.disp > 0] <- res$dev
          parms[nb.disp > 0, ] <- res$parms
        }
        res <- list(dev = dev, means = means, parms = parms)
      }
    }
    
    ### Analyze using quasi-Poisson model, if chosen
    if (Model == "Poisson") res <- PoisDev(counts, design, log.offset, print.progress)
    ### Record means and parameter estimate from full model
    if (jj == 1) {
      means <- res$means
      parms <- res$parms
    }
    deviance.list[[jj]] <- res$dev
  }
  
  LRT <- num.df <- NULL
  
  if (length(design.list) > 1) {
    ### Compute likelihood ratio test statistics. If not otherwise specified, compare each model to the first model
    ### in design.list, which should be the full model
    
    if (is.null(test.mat)) {
      print("Note: 'test.mat' not provided. Comparing each model \nfrom 'design.list' to first model in 'design.list', which must be the full model")
      test.mat <- cbind(1, 2:length(design.list))
      rownames(test.mat) <- paste("Design", 1, " vs Design", 2:length(design.list), sep = "")
    }
    
    for (i in 1:nrow(test.mat)) {
      i1 <- test.mat[i, 1]
      i2 <- test.mat[i, 2]
      num.df <- c(num.df, abs(p[i2] - p[i1]))
      LRT <- cbind(LRT, -(deviance.list[[i2]] - deviance.list[[i1]])/(p[i2] - p[i1]))
    }
    colnames(LRT) <- rownames(test.mat)
  }
  den.df <- (n - p[1])
  
  ### Compute deviance dispersion estimate
  phi.hat.dev <- deviance.list[[1]]/den.df
  
  ### Compute Pearson dispersion estimate
  if (Model == "NegBin") 
    phi.hat.pearson <- (means - counts)^2/(means + means^2 * nb.disp)
  if (Model == "Poisson") 
    phi.hat.pearson <- (means - counts)^2/means
  phi.hat.pearson[means == 0] <- 0
  phi.hat.pearson <- rowSums(phi.hat.pearson)/den.df
  
  if (Model == "Poisson") 
    return(list(LRT = LRT, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts), 
                den.df = den.df, num.df = num.df, Model = Model, fitted.values = means, coefficients = parms))
  if (Model == "NegBin") 
    return(list(LRT = LRT, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts), 
                den.df = den.df, num.df = num.df, Model = Model, NB.disp = nb.disp, fitted.values = means, coefficients = parms))
} 