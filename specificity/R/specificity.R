
dum <- function(Data, Key, Focus, Model, After){
    Data$.allTraits <- rowSums(Data[, Key$start:Key$end][, Key$scale %in% Focus])
    mud <- as.call(Model)
    mud$data <- bquote(Data)
    mud <- eval(mud)
    After(mud)
    }

scramble <- function(Key, Focus, Shuffle) switch(Shuffle, 
    none = Key, 
    inclusive = within(Key, scale <- sample(scale)),
    exclusive = within(Key, scale <- (function(K, F) {
        LK <- length(K)
        NI <- sum(K %in% F)
        K2 <- rep(NA, LK)
        sb <- (1:LK)[! K %in% F] 
        K2[sample(sb, NI)] <- F
        K2
        })(Key$scale, Focus))
        )

specLm <- function(Formula, Data, Key, Focus, Shuffle = "none", R = 1000)
  dum(
    Data, 
    scramble(Key, Focus, Shuffle), 
    Focus,
    Model = list(fun = bquote(lm), formula = Formula),
    After = function(..){
      if(is.factor(..$model[,2])) warning("The independent variable of interest is categorical and therefore the estimates may not be appropriate (plese use specificityEta2() instead)")
      coef(summary(..))[2,]        
    } 
    )

specificityLm <- function(Formula, Data, Key, Shuffle = "exclusive", R = 1000){
  if(is.null(Key$names)) Key$names = c(paste("Trait.no.", (1:length(unique(Key$scale))), sep=""))
  uk <- unique(Key$scale)
  to <- system.time(observed <- lapply(uk, function(..) specLm(Formula, Data, Key, Focus = ..)))
  tr <- system.time(random <- lapply(uk, function(..) replicate(R, specLm(Formula, Data, Key, Focus = .., Shuffle))))
  names(observed) <- names(random) <- Key$names
  out <- list(observed = structure(observed, timing = to), random = structure(random, timing = tr), call=match.call(), key=Key, nsims = R, time = tr)
  class(out) <- "specificity"
  out
}

specGlm <- function(Formula, Data, Key, Focus, Shuffle = "none", Family="binomial", R = 1000)
    dum(
        Data, 
        scramble(Key, Focus, Shuffle), 
        Focus,
        Model = list(fun = bquote(glm), formula = Formula, family=Family),
        After = function(..){
          if(is.factor(..$model[,2])) warning("The independent variable of interest is categorical and therefore the estimates may not be appropriate (plese use specificityEta2() instead)")
          coef(summary(..))[2,]        
        } 
        )

specificityGlm <- function(Formula, Data, Key, Shuffle = "exclusive", Family="binomial", R = 1000){
    if(is.null(Key$names)) Key$names = c(paste("Trait.no.", (1:length(unique(Key$scale))), sep=""))
    uk <- unique(Key$scale)
    to <- system.time(observed <- lapply(uk, function(..) specGlm(Formula, Data, Key, Focus = .., Family=Family)))
    tr <- system.time(random <- lapply(uk, function(..) replicate(R, specGlm(Formula, Data, Key, Focus = .., Shuffle, Family=Family))))
    names(observed) <- names(random) <- Key$names
    out <- list(observed = structure(observed, timing = to), random = structure(random, timing = tr), call=match.call(), key=Key, nsims = R, time = tr)
    class(out) <- "specificity"
    out
    }

eta2 = function(x) {
  a = data.frame(Anova(x, type=3))
  part.eta2 = (a[,1] / (a[,1] + a[nrow(a),1]))[-(nrow(a))]
  result = data.frame(Part.Eta.Sq = round(part.eta2, 3), p.value = round(a[1:(nrow(a)-1),4],3))
  rownames(result) = rownames(a)[-(nrow(a))]  
  return(result)
}

specEta2 <- function(Formula, Data, Key, Focus, Shuffle = "none", R = 1000)
  
  dum(
    Data, 
    scramble(Key, Focus, Shuffle), 
    Focus,
    Model = list(fun = bquote(lm), formula = Formula),
    After = function(..) unlist(eta2(..)[2,]) )

specificityEta2 <- function(Formula, Data, Key, Shuffle = "exclusive", R = 1000){
  if(is.null(Key$names)) Key$names = c(paste("Trait.no.", (1:length(unique(Key$scale))), sep=""))
  uk <- unique(Key$scale)
  to <- system.time(observed <- lapply(uk, function(..) specEta2(Formula, Data, Key, Focus = ..)))
  tr <- system.time(random <- lapply(uk, function(..) replicate(R, specEta2(Formula, Data, Key, Focus = .., Shuffle))))
  names(observed) <- names(random) <- Key$names
  out <- list(observed = structure(observed, timing = to), random = structure(random, timing = tr), call=match.call(), key=Key, nsims = R, time = tr)
  class(out) <- "specificity"
  out
}

summary.specificity <- function(object, ...){
    
    #
    R <- object$nsims
    
    #
    true.results <- round(as.data.frame(t(as.data.frame(object$observed))),3)

    rand.est.list <- lapply(object$random, function(..) as.data.frame(..)[1,])
   
    n.scores <- length(unique(object$key$scale)) 
    rand.results.all.traits <- as.data.frame(matrix(ncol=R, nrow=n.scores))
    for(w in 1:n.scores) rand.results.all.traits[w,] <- rand.est.list[[w]][1,] 
    rownames(rand.results.all.traits) <- object$key$names
  
    # 
    rand.results.all.traits$true.beta <- true.results[,1]
    less      <- function(object) 1 - ( sum(object[length(object)] <= object[-(length(object))] ) / R )
    more      <- function(object) 1 - ( sum(object[length(object)] >= object[-(length(object))] ) / R )
  
    true.results$Spec <- apply( rand.results.all.traits, 1, less ) 
    true.results$Spec[true.results[,1] < 0] <- apply( rand.results.all.traits, 1, more )[true.results[,1] < 0]
  
    rand.results.all.traits$true.beta <- NULL
    true.results$Adj.Est <- round(true.results[,1] - rowMeans(rand.results.all.traits), 3)
    true.results$Adj.Est[ true.results[,4] < 0.05 & true.results$adj.eff*true.results[,1] < 0 ] <- "*"
  
    #
    rand.result.mean <- round(rowMeans(as.data.frame(lapply(rand.est.list, rowMeans))),3) 

    #
    total.result <- list(true.results=true.results, 
                      rand.result.mean=rand.result.mean, 
                      rand.results.all.traits=rand.results.all.traits,
                      number.of.sims = R,
                      call = object$call,
                      rand.analyses.time = object$time
                      )
    class(total.result) <- "summary.specificity"
    total.result
    }

print.summary.specificity <- function(x, ...){
    cat("\nUnivariate associations along with specificity estimates and adjusted effect sizes:\n\n")
    print(x$true.results) 
    cat("\n Mean random association: ")
    cat(x$rand.result.mean, "\n\n") 
    }

## end
