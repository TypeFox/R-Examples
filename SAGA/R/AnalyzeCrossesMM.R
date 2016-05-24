AnalyzeCrossesMM <- function(data, Cmatrix = "XY", model.sum = .95, 
                             max.models = 300000, even.sex = F, graph=F,
                             cex.axis=1, cex.names=1, cex.main=1, max.pars = NULL){
  # lets store the graphic paratmeters so we leave people 
  # unscathed for their next plots
  old.par <- par()
  old.par <- old.par[c(1:12,14:18,20,24:53,55:72)]
  
  # lets check and make sure that people picked or supplied a cmatrix
  if(is.null(Cmatrix)) stop("Please supply or choose a Cmatrix")
  
  # if they are supplying the Cmatrix lets do a couple of basic checks
  if(!is.null(Cmatrix)){
    if(!is.vector(Cmatrix)){
      # is it even a matrix
      if(!is.matrix(Cmatrix)) stop("Your supplied c-matrix is not a matrix")
      # does it contain the cohorts we need
      if(!sum(data[,1] %in% Cmatrix[,1]) == nrow(data)){
        stop("The cohort IDs in your data don't match those in your c-matrix")
      }
    }
  }
  
  ## load the possible contributions to cohort means
  if(is.vector(Cmatrix)){
    if(Cmatrix == "XY"){
      Cmatrix <- read.csv(file = system.file("cmatrix.xy.csv", package = "SAGA"), row.names=1)[, -1]
      # remove the effects that are not estimable 
      # with mixed sex lines of unequal ratios
      if(even.sex == F){
        Cmatrix <- Cmatrix[, c(4:6,13:19,22:24) * -1]
      }
    } else if(Cmatrix == "XO" | Cmatrix == "X0"){
      Cmatrix <- read.csv(file = system.file("cmatrix.xy.csv", package = "SAGA"), row.names=1)[, -1]
      # remove the Y chromosome effects
      Cmatrix <- Cmatrix[, c(6,17:19,24) * -1]
      # remove the effects that are not estimable 
      # with mixed sex lines of unequal ratios
      if(even.sex == F){
        Cmatrix <- Cmatrix[, c(4:5,12:15,18:19) * -1]
      }
    } else if(Cmatrix == "ZW"){
      Cmatrix <- read.csv(file = system.file("cmatrix.zw.csv", package = "SAGA"), row.names=1)[, -1]
      # remove the effects that are not estimable 
      # with mixed sex lines of unequal ratios
      if(even.sex == F){
        Cmatrix <- Cmatrix[, c(4:6,13:19,22:24) * -1]
      }
    } else if(Cmatrix == "ZO" | Cmatrix == "Z0"){
      Cmatrix <- read.csv(file = system.file("cmatrix.zw.csv", package = "SAGA"), row.names=1)[, -1]
      # remove the W chromosome effects
      Cmatrix <- Cmatrix[, c(6,17:19,24) * -1]
      # remove the effects that are not estimable 
      # with mixed sex lines of unequal ratios
      if(even.sex == F){
        Cmatrix <- Cmatrix[, c(4:5,12:15,18:19) * -1]
      }
    } else if(Cmatrix == "esd"){
      Cmatrix <- read.csv(file = system.file("cmatrix.esd.csv", package = "SAGA"), row.names=1)[, -1]
    } else {
      stop("Your selection for the c-matrix to use does not match any of the options")
    }
  }
  
  # store the identity of the cohorts
  identities <- data[,1]
  # remove them from the dataframe we will be using for analysis
  data <- data[,-1]
  # keep only those lines that correspond to crosses that thee user has data for
  red.Cmatrix <- Cmatrix[identities,]  
  # lets remove variables that have no difference in lines
  red.Cmatrix <- red.Cmatrix[, c(1, which(apply(red.Cmatrix, 2, var) != 0))]      

  #lets look for composite effects that are identical
  drop.counter <- vector()
  for(i in 2:(ncol(red.Cmatrix)-1)){
    for(j in (i+1):ncol(red.Cmatrix)){
      if(sum(red.Cmatrix[,i] == red.Cmatrix[,j]) == nrow(red.Cmatrix)){
        drop.counter <- c(drop.counter, j)
      }
    }
  }
  #drop those composite effects that are equivelant of lower order simpler effects
  if(length(drop.counter)>0){
    leslie <- paste(colnames(red.Cmatrix)[drop.counter], sep=", ", collapse=", ")
    cat(paste("The following composite effects cannot be estimated with the line \n",
              "means available because they estimate identical quantities to \n",
              "lower order effects: \n", leslie, "\n\n", sep=""))
    red.Cmatrix <- red.Cmatrix[,-drop.counter]
  }
  have.data <- paste(colnames(red.Cmatrix)[-1], collapse = ", ")
  cat(paste("The composite genetic effects that will be tested are: \n", 
              have.data, collapse = ", "), "\n\n")
  # calcualte the potential size of model space
  # the final -2 is because we will always be including the mean so we have
  # one less choice to make
  mod.space.size <- sum(choose(ncol(red.Cmatrix) -1 , 
                               1:(nrow(red.Cmatrix) - 2)))
  if(!is.null(max.pars)) mod.space.size <- sum(choose((ncol(red.Cmatrix) -1), 1:max.pars))
  # warn the user if the model space is very large
  if(mod.space.size > 5000){
    cat(paste("Since there are ", mod.space.size, " possible models this may take a bit:\n", sep=""))
  }
  # generate all possible models storing each matrix in a list
  pos.cols <- 2:ncol(red.Cmatrix)             # col that could be used
  eqns <- list()                              # store the eqns
  counter <- 1                                # index for eqns
  max.par <- nrow(red.Cmatrix) - 2            #
  # if a user has very many cohorts model space can become problematically large
  # however I think that actually very few datasets support >>large models with 
  # many important factors.  So one solution is simply to allow users to set a max
  # model size this makes things fairly easy to handle
  if(!is.null(max.pars)) max.par <- max.pars
  
  cat(paste("Generating Models"))
  if(length(pos.cols) < max.par){
    max.par <- length(pos.cols)
  }
  for(i in 1:max.par){                     # different number of par models
    cat(".")
    foo <- combn(pos.cols, i)              # all pos models with i variables
    # this loop just places the models generate with i variables into the list
    # of all possible models.  Models are described by the columns they include
    for(j in 1:ncol(foo)){
      eqns[[counter]] <- as.vector(foo[,j])
      counter <- counter + 1
    }
  }
  # just the setup for a small counter
  if(length(eqns) < 101) x <- 50
  if(length(eqns) > 100) x <- 50
  if(length(eqns) > 1000) x <- 500
  if(length(eqns) > 10000) x <- 5000
  # now we test each model
  mod.results <- list()                  # stores glm results
  num.pars <- dev <- aic <- vector()     # stores various useful values
  # we need a counter because redundant models arrise.  These originate because
  # some components will have high covariance depending on the lines included
  # in the dataset.  The glm function automatically throws these variables
  # resulting in fitting the same model more than once.
  counter <- 1
  for(i in 1:length(eqns)){
    # generate the matrix for the current model
    test.mat <- as.matrix(red.Cmatrix[, c(1, eqns[[i]])])
    # fit the model weight is equal to the inverse of the square of the SE
    temp.mod <- glm(data[, 1] ~ test.mat, weights = data[, 2] ^ - 2)
    # this if statement will bypass a model with a singularity
    # 1 NA will be generated for the line mean any additional are sign of sing.
    if(sum(is.na(temp.mod$coef)) < 2){
      # name model results as eqns
      mod.results[[counter]] <- temp.mod
      names(mod.results)[counter] <- i
      # record the number of parameters in the model
      num.pars[counter] <- length(mod.results[[counter]]$coefficients) - 1
      # record the residual deviances
      dev[counter] <- mod.results[[counter]]$dev
      # record the AIC of the models
      aic[counter] <- mod.results[[counter]]$aic
      counter <- counter + 1
    }
    if(i / x == round(i / x)) cat(paste("\n", i))
  }
  ## need to report the number of models thrown out due to 
  ## high covariance ~ singularity
  if(i > counter){
  cat(paste("\n", i - (counter - 1), 
            " models were removed due to high covariances \n",
              "or linear relationships between predictor variables.  \n", "The remaining ", 
              counter - 1, " models have been evaluated.\n\n", sep = ""))
  }
  # in the unrealistic situation where there was a model that predicted the
  # data perfectly we would get -Inf for the AIC should only be an issue in 
  # simulated data... for these purposes lets just plug in something that is 
  # equal to the lowest AIC value for models in that same parameter range
  for(i in 1:length(aic)){
    if(aic[i] == -Inf){
      aic[i] <- sort(unique(aic))[2] -1
    }
  }
  # calculate aicc and delta aicc
  aicc <- aic + (((2 * num.pars) * (num.pars + 1)) / 
                   (nrow(data) - num.pars))
  daicc <- aicc - min(aicc)
  # this code correctly produces akaike weights 
  waic <- (exp(-.5 * daicc)) / (sum(exp(-.5 * daicc)))
  new.waic <- waic
  new.waic.names <- eqns[as.numeric(names(mod.results))]
  # so now we have a copy of the waics to play with
  new.vars <- matrix(0,(ncol(red.Cmatrix)-1),2)
  new.vars[,1] <- colnames(red.Cmatrix)[2:ncol(red.Cmatrix)]
  for(i in 1:nrow(new.vars)){
    for(j in 1:length(new.waic)){
      if((i+1) %in% new.waic.names[[j]]){
        new.vars[i,2] <- as.numeric(new.vars[i,2]) + new.waic[j]
      }
    }
  }
  # lets calculate the 95% probability set of models
  best.models <- list()
  counter <- i <- 0
  good.model.waics <- vector()
  while(counter < model.sum){
    i <- i + 1
    counter <- counter + waic[order(waic, decreasing = T, na.last = F)][i]
    good.model.waics[i] <- waic[order(waic, decreasing = T, na.last = F)][i]
  }
  best.models.ind <- order(waic, decreasing = T, na.last = F)[1:i]
  best.models <- mod.results[best.models.ind]
  cat(paste("\nAICc weights were used to select the minimum number of models ",
              "whose weights sum \nto greater than ", 
            model.sum * 100, "% this model set includes ", length(best.models), 
              " model(s)\n", sep = ""))
   #lets calculate variable importance
   #which equations are being used
   best.eqns <- eqns[as.numeric(names(best.models))]
   best.eqns.w <- waic[sort(best.models.ind)]
   #lets tell the user the models being included
   #for(i in 1:length(best.eqns)){
   #     cat(paste(colnames(red.Cmatrix)[c(1, best.eqns[[i]])], collapse = ", "), 
   #      "  waic = ",good.model.waics[i],"\n")
   #  }
  # now we need to print the model weighted averages and SE
  # lets make a matrix of the calculated values under each model
  par.est <- matrix(0, length(best.eqns), ncol(red.Cmatrix) + 2)
  colnames(par.est) <- c('eqn', colnames(red.Cmatrix), 'mw')
  par.est[, 1] <- names(best.models)
  par.est[, 2] <- 1
  # now we need a 1 or 0  if the parameter is in the eqn
  for(i in 1:nrow(par.est)){
    bar <- as.numeric(par.est[i, 1])
    par.est[i, eqns[[bar]] + 1] <- 1
  }
  # now replace 1's with the parameter estimate for each variable
  for(i in 1:nrow(par.est)){
    bar <- best.models[[i]]$coefficients[-2]
    counter <- 0
    for(j in 2:ncol(par.est)){
      if(par.est[i, j] == 1){
        counter <- counter + 1
        par.est[i, j] <- bar[counter]
      }
    }
  }
  # add in the aicw
  # best.models has eqn lookup in waic
  names(waic) <- names(mod.results)
  for(i in 1:nrow(par.est)){
    par.est[i, ncol(par.est)] <- waic[names(waic) == par.est[i, 1]]
  }
  # recalculate model waic to sum to 1
  par.est[, 'mw'] <- as.numeric(par.est[, 'mw']) / 
    sum(as.numeric(par.est[, 'mw']))
  # calculate the model weighted parameter estimates
  par.est <- rbind(par.est, rep(0,ncol(par.est)))
  for(i in 2:(ncol(par.est)-1)){
    par.est[nrow(par.est), i] <- sum(as.numeric(par.est[1:nrow(par.est)-1, i]) * 
                             as.numeric(par.est[1:nrow(par.est)-1, 'mw']))
  }
  par.est[nrow(par.est), 1] <- 'mw.avg'
  # now lets calculate the unconditional variances as proposed in burnham and
  # anderson 2002 pg 162
  # so lets duplicate table par.est to use to fill in our values
  var.est <- par.est
  # lets loop through models first with i ... the rows
  for(i in 1:(nrow(var.est) - 1)){
    counter <- 0
    mod.vars <- diag(vcov(best.models[[which(names(best.models) == 
                                               var.est[i, 1])]]))
    # now lets loop through parameters with j ... the columns
    for(j in 2:(ncol(var.est) - 1)){
      if(var.est[i, j] != 0){
        counter <- counter + 1
        var.est[i,j] <- mod.vars[counter] 
      }
    }
  }
  # we now have all of the required variables for the uncond. variance est.
  for(i in 2:(ncol(var.est) - 1)){
    if(as.numeric(par.est[nrow(par.est), i]) != 0){
      foo <- vector()
      for(j in 1:(nrow(var.est) - 1)){
        foo[j] <- as.numeric(par.est[j, 'mw']) * 
                  sqrt(as.numeric(var.est[j, i]) + 
                         ((as.numeric(par.est[j, i]) - 
                             as.numeric(par.est[nrow(par.est), i])) ^ 2))
      }
      var.est[nrow(var.est), i] <- sqrt(sum(foo) ^ 2)
    }
  }
  # now lets make a table with the stuff we want
  results <- matrix(, 2, ncol(red.Cmatrix))
  colnames(results) <- colnames(red.Cmatrix)
  row.names(results) <- c('Model Weighted Average', 
                          'Unconditional Standard Error')
  results[1, ] <- par.est[nrow(par.est), 2:(ncol(par.est) - 1)]
  results[2, ] <- var.est[nrow(var.est), 2:(ncol(var.est) - 1)]
  if(graph == T){
    # make colors for barplot
    # get some extra room
    par(mar=c(2, 2, 2, 6))
    foo.colors <- heat.colors(100)[100:1]
    plot.colors <- vector()
    counter <- 0
    plot.colors <- rep(0, (ncol(results) - 1))
    names(plot.colors) <- colnames(results)[2:ncol(results)]
    for(i in 2:ncol(results)){
      counter <- counter + 1
      if(colnames(results)[i] %in% new.vars[, 1]){
        plot.colors[counter] <- round(as.numeric(new.vars[new.vars[, 1] == colnames(results)[i], 2]) * 100)
      }
      if(plot.colors[counter] == 0){
        plot.colors[counter] <- 1
      }
    }
    maxval <- max(as.numeric(results[1, 2:ncol(results)]) + 
                    as.numeric(results[2, 2:ncol(results)]))
    # this is a little fix for when everything is below 0
    if(maxval < 0){
      maxval <- 0
    }
    minval <- min(as.numeric(results[1, 2:ncol(results)]) - 
                    as.numeric(results[2, 2:ncol(results)]))
    if(minval > 0){
      minval <- 0
    }
    if(length(best.eqns.w) > 1){
      mp <- barplot(as.numeric(results[1,2:ncol(results)]), 
                    names.arg=colnames(results)[2:ncol(results)], 
                    col=foo.colors[as.vector(plot.colors)],
                    ylim = c(minval - .4 * abs(minval), maxval + .4 * abs(maxval)), 
                    main = "Model Weighted Averages and Unc. SE", 
                    cex.axis=cex.axis, cex.names=cex.names, cex.main=cex.main)
      segments(mp, as.numeric(results[1,2:ncol(results)]) - 
                 as.numeric(results[2,2:ncol(results)]), 
               mp, as.numeric(results[1,2:ncol(results)]) + 
                 as.numeric(results[2,2:ncol(results)]), lwd=2)
      # Now plot the horizontal bounds for the error bars
      # 1. The lower bar
      segments(mp - 0.1, as.numeric(results[1,2:ncol(results)]) - 
                 as.numeric(results[2,2:ncol(results)]), mp + 0.1, 
               as.numeric(results[1,2:ncol(results)]) - 
                 as.numeric(results[2,2:ncol(results)]), lwd = 2)
      # 2. The upper bar
      segments(mp - 0.1, as.numeric(results[1,2:ncol(results)]) + 
                 as.numeric(results[2,2:ncol(results)]), mp + 0.1, 
               as.numeric(results[1,2:ncol(results)]) + 
                 as.numeric(results[2,2:ncol(results)]), lwd = 2)
      # add a legend for the variable importance
      locs <- par("usr")  
      color.legend(locs[2] ,  #xl
                   locs[4]- ((locs[4]-locs[3]) * 0.5),  #yb
                   locs[2] + (locs[2]*.05),    #xr
                   locs[4],        #yt
                   legend = c("0.00", "0.25", "0.50", "0.75","1"),
                   rect.col = foo.colors,
                   cex = 0.6,
                   align="rb",
                   gradient = "y")
      par(xpd = TRUE)
      text(x = (((locs[2]) + (locs[2] + (locs[2]*.05)))/2), 
           y = locs[4], 
           labels = "Var. Imp.", cex = 0.7, pos = 3)
    }else{
      result.1 <- best.models[[1]]
      pars.1 <- colnames(red.Cmatrix)[best.eqns[[1]]]
      par.est <- summary(result.1)$coefficients[, 1]
      par.est <- par.est[2:length(par.est)]
      se.est <- summary(result.1)$coefficients[, 2]
      se.est <- se.est[2:length(se.est)]
      max.val <- max(par.est + se.est)
      min.val <- min(par.est - se.est)
      if(max.val < 0){
        max.val <- 0
      }
      if(min.val > 0){
        min.val <- 0
      }
      mp <- barplot(par.est, names.arg=pars.1, 
                    main = "Single Model Means and Cond. SE",
                    ylim = c(minval - .4 * minval, maxval + .4 * maxval), 
                    cex.axis=cex.axis, cex.names=cex.names, cex.main=cex.main)   
      high.se <- par.est + se.est
      low.se <- par.est - se.est
      segments(mp, high.se, mp, low.se, lwd=3)
      segments(mp - 0.1, high.se, mp + 0.1, high.se, lwd=3)
      segments(mp - 0.1, low.se, mp + 0.1, low.se, lwd=3)
      results.1 <- matrix(, 2, length(pars.1))
      results.1[1, ] <- par.est
      results.1[2, ] <- se.est
      colnames(results.1) <- pars.1
      row.names(results.1) <- c("Estimate", "Cond. SE")
    }
  } 
  ## prepare the results to be returned to the user
  final.results <- list()
  mod.names <- list()
  foo2 <- colnames(red.Cmatrix)
  mod.ind <- as.numeric(names(mod.results))
  for(i in 1:length(mod.results)){
    mod.names[[i]] <- paste(foo2[eqns[[mod.ind[i]]]], sep="", collapse=", ")
  }
  names(mod.results) <- mod.names
  if(length(mod.results) > max.models){
    mod.results <- mod.results[order(waic, decreasing = T)[1:max.models]]
    daicc <- daicc[order(daicc, decreasing =F)[1:max.models]]
  }
  final.results[[1]] <- mod.results
  final.results[[2]] <- results
  final.results[[3]] <- daicc
  final.results[[4]] <- new.vars
  names(final.results) <- c("models", "estimates", "daicc", "varimp")
  class(final.results) <- "genarch"
  
  # lets reset peoples graphics paramters so they make sinces again
  par(old.par)
  return(final.results)
}