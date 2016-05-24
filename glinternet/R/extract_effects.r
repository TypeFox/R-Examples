extract_effects = function(betahat, activeSet, numLevels){

  make_z = function(levels){
    z1 = matrix(0, prod(levels), levels[1])
    for (i in 1:levels[2]){
      z1[((i-1)*levels[1]+1):(i*levels[1]), ] = diag(levels[1])
    }
    z2 = matrix(0, prod(levels), levels[2])
    for (i in 1:levels[2]){
      z2[((i-1)*levels[1]+1):(i*levels[1]), i] = 1
    }
    list(z1=z1, z2=z2)
  }
  
  extract_interaction = function(beta, levels){
    #case 1: cat-cat
    if (min(levels) > 1){
      z = make_z(levels)
      mainEffects = sapply(z, function(x) lm(beta ~ x - 1)$fitted.values)
      result = beta - apply(mainEffects, 1, sum)
    }
    #case 2: cat-cont
    else if (any(levels > 1)){
      L = max(levels)
      result = beta[(L+1):(2*L)] - mean(beta[(L+1):(2*L)])
    }
    #case 3: cont-cont
    else {
      result = beta[3]
    }
    result
  }

  extract_main_effect = function(beta, levels){
    if (length(levels) > 1){
      #case 1: cat-cat
      if (min(levels) > 1){
        z = make_z(levels)
        result = lapply(z, function(x) lm(beta ~ x - 1)$coefficients)
      }
      #case 2: cat-cont
      else if (any(levels > 1)){
        L = max(levels)
        result = matrix(0, L, 2)
        result[, 1] = beta[1:L]
        result[, 2] = mean(beta[(L+1):(2*L)])
      }
      #case 3: cont-cont
      else {
        result = beta[1:2]
      }
    }
    else {
      result = beta
    }
    result
  }

  group_sizes = function(levels){
    if (length(levels) == 1) return (levels)
    #cat-cat
    if (min(levels) > 1) return (prod(levels))
    #cont-cont
    if(max(levels) == 1) return (3)
    #cat-cont
    return (2*max(levels))
  }

  catRef = which(numLevels > 1)
  contRef = which(numLevels == 1)
  
  #get the main effects and interactions from activeSet
  mainEffects = list(cat=NULL, cont=NULL)
  mainEffectsCoef = list(cat=list(), cont=list())
  interactions = list(catcat=NULL, contcont=NULL, catcont=NULL)
  interactionsCoef = list(catcat=list(), contcont=list(), catcont=list())
  offset = 1
  if (!is.null(activeSet$cat)){#categorical
    mainEffects$cat = c(mainEffects$cat, activeSet$cat)
    sizes = sapply(activeSet$cat, function(x) group_sizes(numLevels[catRef[x]]))
    for (i in 1:length(sizes)){
      mainEffectsCoef$cat = c(mainEffectsCoef$cat, list(betahat[(offset+1):(offset+sizes[i])]))
      offset = offset + sizes[i]
    }
  }
  if (!is.null(activeSet$cont)){#continuous
    mainEffects$cont = c(mainEffects$cont, activeSet$cont)
    sizes = rep(1, length(activeSet$cont))
    for (i in 1:length(sizes)){
      mainEffectsCoef$cont = c(mainEffectsCoef$cont, list(betahat[(offset+1):(offset+sizes[i])]))
      offset = offset + sizes[i]
    }
  }
  if (!is.null(activeSet$catcat)){#categorical-categorical
    interactions$catcat = activeSet$catcat
    for (i in 1:nrow(activeSet$catcat)){
      newMain = activeSet$catcat[i, ]
      levels = numLevels[catRef[newMain]]
      size = group_sizes(levels)
      newCoef = extract_main_effect(betahat[(offset+1):(offset+size)], levels)
      for (j in 1:2){
        if (newMain[j] %in% mainEffects$cat){#main effect already present
          idx = which(mainEffects$cat %in% newMain[j])
          mainEffectsCoef$cat[[idx]] = mainEffectsCoef$cat[[idx]] + newCoef[[j]]
        }
        else {#not present
          mainEffects$cat = c(mainEffects$cat, newMain[j])
          mainEffectsCoef$cat = c(mainEffectsCoef$cat, list(newCoef[[j]]))
        }
      }
      interactionCoef = matrix(extract_interaction(betahat[(offset+1):(offset+size)], levels), nrow=levels[1])
      rownames(interactionCoef) = paste(rep("cat",levels[1]), newMain[1], paste("_",0:(levels[1]-1),sep=""), sep="")
      colnames(interactionCoef) = paste(rep("cat",levels[2]), newMain[2], paste("_",0:(levels[2]-1),sep=""), sep="")
      interactionsCoef$catcat = c(interactionsCoef$catcat, list(interactionCoef))
      offset = offset + size
    }
  }
  if (!is.null(activeSet$contcont)){#continuous-continuous
    interactions$contcont = activeSet$contcont
    for (i in 1:nrow(activeSet$contcont)){
      newMain = activeSet$contcont[i, ]
      size = 3
      newCoef = extract_main_effect(betahat[(offset+1):(offset+size)], c(1,1))
      for (j in 1:2){
        if (newMain[j] %in% mainEffects$cont){
          idx = which(mainEffects$cont %in% newMain[j])
          mainEffectsCoef$cont[[idx]] = mainEffectsCoef$cont[[idx]] + newCoef[j]
        }
        else {
          mainEffects$cont = c(mainEffects$cont, newMain[j])
          mainEffectsCoef$cont = c(mainEffectsCoef$cont, list(newCoef[j]))
        }
      }
      interactionsCoef$contcont = c(interactionsCoef$contcont, list(extract_interaction(betahat[(offset+1):(offset+size)], c(1,1))))
      offset = offset + size
    }
  }
  if (!is.null(activeSet$catcont)){#categorical-continuous
    interactions$catcont = activeSet$catcont
    for (i in 1:nrow(activeSet$catcont)){
      newMain = activeSet$catcont[i, ]
      levels = numLevels[c(catRef[newMain[1]], contRef[newMain[2]])]
      size = group_sizes(levels)
      newCoef = extract_main_effect(betahat[(offset+1):(offset+size)], levels)
      if (newMain[1] %in% mainEffects$cat){
        idx = which(mainEffects$cat %in% newMain[1])
        mainEffectsCoef$cat[[idx]] = mainEffectsCoef$cat[[idx]] + newCoef[, 1]
      }
      else {
        mainEffects$cat = c(mainEffects$cat, newMain[1])
        mainEffectsCoef$cat = c(mainEffectsCoef$cat, list(newCoef[, 1]))
      }
      if (newMain[2] %in% mainEffects$cont){
        idx = which(mainEffects$cont %in% newMain[2])
        mainEffectsCoef$cont[[idx]] = mainEffectsCoef$cont[[idx]] + newCoef[1, 2]
      }
      else {
        mainEffects$cont = c(mainEffects$cont, newMain[2])
        mainEffectsCoef$cont = c(mainEffectsCoef$cont, list(newCoef[1, 2]))
      }
      interactionsCoef$catcont = c(interactionsCoef$catcont, list(extract_interaction(betahat[(offset+1):(offset+size)], levels)))
      offset = offset + size
    }
  }

  list(mainEffects=mainEffects, mainEffectsCoef=mainEffectsCoef, interactions=interactions, interactionsCoef=interactionsCoef)
}
