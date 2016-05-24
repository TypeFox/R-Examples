##' 'glmClust' cluster longitudinal data (trajectories) using
##' the likelihood as a metric of distance, it also deals with
##' multiples covariates with different effects using the generalised
##' linear model 'glm'.
##'
##' 
##' 'glmClust' implements an ECM (esperance classification maximisation)
##' type algorithm which assigns the trajectories to the cluster
##' maximising the likelihood. The procedure is repeated until no change
##' in the partitions or no sufficient increase in the likelihood is possible.\cr
##'
##' 
##' 'glmClust' also deals with multiple covariates with different level
##' effects, different  in each cluster and/or identical for
##' all of them.\cr
##'
##' 
##' The introduction of covariates is possible thanks to 'glm' which
##' fits a generalised linear model and take into account the type of
##' the response (normal, binomial, Poisson ...etc) and the link
##' function. \cr
##'
##' 
##' Several parameters of 'glmClust' are in common with 'glm', like the
##' \code{formula} which requires a particular attention by specifying
##' the covariates with a cluster effect, for e.g. \code{clust(T1+T2+..+Tn)},
##' the covariates with an identical effect in each cluster are specified
##' with the keyword \bold{pop}, for e.g. \code{pop(X1+X2+..+Xn)}, 
##' note that these last covariates are optional. \cr
##' The data are in the long format and no missing values are allowed.\cr
##'
##' 
##' In the parametric case (\code{timeParametric = TRUE})  multiples
##' covariates are allowed, in the non-parametric case only one covariate
##' is allowed.\cr
##'
##' 
##' The algorithm depends greatly on the starting condition,  which
##' is obtained by randomly affecting the trajectories to the
##' clusters unless the user introduce his own partition.
##' To obtain better results it is desirable to run
##' the algorithm several times from different starting points,
##' therefore it is preferable to use \link{kmlCov} which runs the
##' algorithm several times with different number of clusters. 
##'
##' 
##' At the end of the algorithm, an object of class \linkS4class{GlmCluster}
##' is returned and contains information about the affectation of the  
##' trajectories, the proportions, the convergence, ...etc. The main trajectories
##' can be simply visualised by \code{plot(my_GlmCluster_Object)}.
#################################################################'
##' @title Clustering  longitudinal data
##' @aliases glmClust
##' 
##' @usage glmClust(formula, data, ident, timeVar, nClust,  family = 'gaussian',
##' effectVar = '', weights = rep(1,nrow(data)), affUser, timeParametric = TRUE,
##' separateSampling = TRUE, max_itr = 100, verbose = TRUE)
##' @param formula A symbolic description of the model. In the parametric case
##' we write for example 'y ~ clust(time+time2) + pop(sex)', here 'time' and 'time2'
##' will have a different effect according to the cluster, the 'sex' effect is the
##' same for all the clusters. In the non-parametric case only one covariate is allowed.
##' @param data  A [data.frame] in long format (no missing values) which means
##' that each line corresponds to one measure of the observed phenomenon,
##' and one individual may have multiple measures (lines) identified by an
##' identity column. In the non-parametric case the totality of patients must have all
##' the measurements at fixed times.
##' @param nClust The number of clusters, between 2 and 26.
##' @param ident Name of the column identity in the data.
##' @param timeVar Name of the 'time' column  in the data.
##' @param family  A description of the error distribution and link function to be used
##' in the model, by default 'gaussian'. This can be a character string naming a family
##' function, a family function or the result of a call to a family function.  (See
##' \link{family} for more details of family functions).
##' @param effectVar Name of the effect specified or not in the formula is has
##' level cluster effect or not (optional), note that this parameter
##' is useful for the function \link{plot}
##' @param weights Vector of 'prior weights' to be used in the fitting process,
##' by default the weights are equal to one.
##' @param affUser Initial affectation of the individuals in a [data.frame] format,
##' if missing the individuals are randomly assigned to the clusters so it is optional .
##' @param timeParametric By default [TRUE] thus parametric on the time. If [FALSE]
##' then only one covariate is allowed in the formula and the algorithm used is the k-means.
##' @param separateSampling By default [TRUE] it means that the proportions of the clusters
##' are supposed equal in the classification step, the log-likelihood maximised at each step
##' of the algorithm is \eqn{\sum_{k=1}^{K}\sum_{y_i \in P_k} \log(f(y_i, \theta_k))}{},
##' otherwise the proportions of clusters are taken into account and the log-likelihood
##' is \eqn{latex}{\sum_{k=1}^{K}\sum_{y_i \in P_k} \log(\lambda_{k}f(y_i, \theta_k))}.
##' @param max_itr The maximum number of iterations fixed at 100.
##' @param verbose Print the output in the console.
##' @return An object of class \linkS4class{GlmCluster}.
##' @export
##' @seealso \link{kmlCov}
##' @examples
##' data(artifdata)
##' res <- glmClust(formula = Y ~ clust(time + time2 + time3) + pop(treatTime),
##' data = artifdata, ident = 'id', timeVar = 'time', effectVar = 'treatment', nClust = 4)
##' # the trajectories with indices 0 indicate the ones with a normal treatment, 1 indicate a high dose
##' # the color indicates the clusters
##' # the proportions are in the table above the diagram
##' plot(res)
 
glmClust <- function(formula, data,  ident, timeVar, nClust, family = 'gaussian',
                     effectVar = '', weights = rep(1,nrow(data)) , affUser,
                     timeParametric = TRUE, separateSampling = TRUE, max_itr = 100, verbose = TRUE) {
  if (missing(formula) || !is.language(formula) || missing(ident) ||  missing(data) || missing(nClust) || missing(timeVar)) {
    cat(' ### Missing argumetn(s) ! ###\n') #
    return(NULL)
  }
  if (!(ident %in% names(data))) { cat(" Wrong name of column identity !\n"); return(NULL) } # cherche la colonne identite
  if (all.vars(formula)[2] == '.') {
    formula <- rwFormula(formula, names(data), ident)
  }
  ## On verifie que la famille est reconnue, code de la fonction 'glm'
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (timeParametric) {  # on teste si on est dans le cas parametrique
    formulaS <- seperateFormula(formula) # on separe les effets differents par cluster des effets
    formula <- formulaS[[1]]

    if (b_fix_eff <- (length(formulaS)==2)) { # on test si on a un effet commun
      fix_eff <- formulaS[[2]]
    } else {  }
  } else {
    ## Ajouter une instruction pour empecher les nouvelles covariables
    if (length(all.vars(formula))>3) {
      cat('Only one covariate is allowed in the non-paramtric cas.\n')
      return(NULL)
    }
    b_fix_eff <- FALSE
    time_dif <- unique(data[,all.vars(formula)[2]])
    for (i in 1:length(time_dif)) {
      eval(parse(text=paste0( all.vars(formula)[2] , i,'<-', deparse(quote(majIndica(data[,all.vars(formula)[2]], time_dif[i])))))) # time_1 ...time_n
    }
    formula <- paste0(all.vars(formula)[1],'~', paste0(all.vars(formula)[2], 1:length(time_dif), collapse='+'))
    formula <- eval( parse(text = formula))
  }

  nomClust <- 'G'                       # par defaut je le fixe a `G'
  covar <- all.vars(formula)[-1]        # recupere les covariables
  expCov <- addIndic(covar, nClust, timeParametric, nomClust = nomClust)
  exp_fix <- ifelse(b_fix_eff , paste0('+', as.character(fix_eff)[length(as.character(fix_eff))]), '')
  if (b_fix_eff ) { cov <- c(covar, all.vars(fix_eff)) }  else { cov <- covar }

  exp <- paste0(all.vars(formula)[1],'~', expCov, exp_fix, '-1')
  exp_fin <- eval(parse(text=exp))

  nObs <- nrow(data)

  id_data <- unique(data[, ident]) # renvoie les identifiants des individus
  nTraj <- length(id_data)         # nombre d'individus
  nMesInd <- tapply(data[, ident], data[,ident], length) # nb mesure pour chaque individu

  ## la ligne qui suit est facultative dans le cas ou les donnees sont triees par individus par contre
  ## si elles sont triees par numero de mesure alors ce script est INDISPENSABLE
  ## voir meme si elle sont dispersees
  o_ind <- order(order(data[, ident]))

  ## on teste si l'utilisateur a entre une partition de depart, ATTENTION la partition doit etre un data.frame (pour l'instant)
  if (missing(affUser)) {
    aff_init <- aff_cour <- affect_rand(nTraj, nClust)
    aff_init_long <- aff_cour_long <- aff_new_long <- rep(aff_cour, nMesInd)[o_ind]
  } else {
    if (nrow(affUser) == nObs) {
      aff_init_long <- aff_cour_long <- aff_new_long <- affUser[, , drop=T]
      aff_init <- aff_cour <- aff_cour_long[nMesInd]
    } else if ( nrow(affUser)== nTraj) {
      aff_init <- aff_cour <- affUser[, , drop=T]
      aff_init_long <- aff_cour_long <- aff_new_long <- rep(aff_cour, nMesInd)[o_ind]
    }
  }

  for (i in 1:nClust ) {
    eval(parse(text=paste0( nomClust, letters[i],'<-', deparse(quote(majIndica(aff_cour_long, i)))))) # on creee les indicatrices d'appartenance aux clusters
  }
  model.glm <- model.glm.init <- glm(terms(exp_fin, keep.order = TRUE), data = data, family = family, weights = weights) # mod.glm.init n'est pas utilise ) X (, de plus terms() ne sert a rien
  aff_new <- rep(0,nTraj)              # affecte iter courant, new aff
  itrCour <- 1
  bPasChan <- FALSE # indique si il y a eu un changement d'affecttation entre 2 iter
  l.lik <- rep(NA, max_itr+1) # stocke la log vraisemblance aux iterations
  l.lik[1] <- logLik(model.glm)-1 # j'enleve 1 pour passer le test de tolerance
  l.lik[2] <- l.lik[1]+1
  if(!separateSampling) {
    l.lik <- l.lik-(log(nObs)*sum(aff_cour_long)-sum(aff_cour_long*log(aff_cour_long)))
  } else {  }

  tol <- 1e-5 # si la difference de logVrais est <=tol alors on s'arrete, penser a mettre ce parametre en option # A mettre en parametre ds une autre fonction
  pred <- matrix(0, nrow=nObs, ncol=nClust) # Contient les valeurs predites dans chaque colonne qui correspond a un cluster
  convergence <- TRUE

  nom_coef_clust <- list()    # vecteur avec les noms des coefficients
  pred_all <- list()   # contien l'expression pour predire les valeurs
  if (timeParametric) {
    if (b_fix_eff) {
      ## avec effet fixe
      for (itrClust in 1:nClust) {
        nom_coef_clust[[itrClust]] <- eval(parse(text = getNomCoef(covar = covar, cov_fix = all.vars(fix_eff),
                                                   nomClust = nomClust, itrClust = itrClust)))
        pred_all[[itrClust]] <- parse(text = predict_clust(cov = cov, nomCoef = nom_coef_clust[[itrClust]],
                                        deparse(quote(model.glm))))
      }
    } else {
      ## sans effet fixe
      for (itrClust in 1:nClust) {
        nom_coef_clust[[itrClust]] <- eval(parse(text = getNomCoef(covar = covar, nomClust = nomClust,
                                                   itrClust = itrClust)))
        pred_all[[itrClust]] <- parse(text = predict_clust(cov = cov, nomCoef = nom_coef_clust[[itrClust]],
                                        deparse(quote(model.glm))))
      }
    }
  }
  ## ===========================  Debut de ECM  ===========================  #####
  while(!bPasChan && itrCour<max_itr && (l.lik[itrCour+1]-l.lik[itrCour])>tol) {
    itrCour <- itrCour+1
    bPasChan <- TRUE
#######################################
    ## Calcul des valeurs predites
    formula.o <- formula              # on stocke la formule de depart

    if (separateSampling) {
      prop <- rep(1, nClust) # ici on suppose que les clusters ont les mm probas a priori
    } else {
      ## prop <- table(aff_cour_long)  # dans ce cas on prend en compte le nombre de mesures, un cluster avec beaucoup de mesures est plus attractif
      prop <- table(aff_cour) # ici on prend en compte le nombre d'individus par cluster, ainsi un cluster contenant beaucoup d'individus sera plus attractif
    }

    if (timeParametric) {
      ##        yy <- lapply(pred_all, eval)
      for (itrClust in 1:nClust) {
        pred[, itrClust] <- eval(pred_all[[itrClust]])
      }
      ## Calcul de la vraisemblance et affectation des individus

      for(itrTraj in 1:length(id_data)) {
        indTraj <- which(data[,ident] == id_data[itrTraj])
        likelihood <- rep(NA,nClust)   #  en fait la log vraisemblance
        for(itrClust in 1:nClust) {
          ## likelihood[itrClust] <- log_lik(y = data[indTraj, all.vars(formula)[1]], n = rep(1, nMesInd[itrTraj]), mu = pred[indTraj, itrClust], wt = weights[indTraj],  family = family, nparam = mod.glm$rank, disp_mod = summary(mod.glm)[[14]]) # remplacer la formule du dessous par celle du dessus
          likelihood[itrClust] <- log_lik(y = data[indTraj, all.vars(formula)[1]], n = rep(1, length(indTraj)),
                                          mu = pred[indTraj, itrClust], wt = weights[indTraj],  family = family,
                                          nparam = model.glm$rank, disp_mod = summary(model.glm)[[14]])
        }
        aff_new[itrTraj] <- which.max(likelihood*prop)
      }
    } else {
      moy.clust <- t(matrix(coefficients(model.glm), nrow = length(time_dif) )) # les coefficients sont les moy a chaque pas de temps
      for(itrTraj in 1:length(id_data)) {
        indTraj <- which(data[,ident] == id_data[itrTraj])
        dist.moy <- c()
        dist.moy <- apply( moy.clust, 1, function(x) dist(rbind(x, data[indTraj, all.vars(formula)[1] ])))
        aff_new[itrTraj] <- which.min(dist.moy*prop)
      }
    }

    if ( !identical(aff_new,aff_cour) && length(table(aff_new))==length(prop) ) { # 1st stopping condition
      aff_new_long <- rep(aff_new, nMesInd)[o_ind]
      for (i in 1:nClust ) {
        eval(parse(text=paste(nomClust, letters[i],'<-', deparse(quote(majIndica(aff_new_long, i))), sep =''))) # on cree les indicatrices
      }
      model.glm.cour <- glm( terms(model.glm$formula, keep.order = TRUE), data = data, family = family) # MAJ du model

      l.lik[itrCour+1] <- logLik(model.glm.cour)-ifelse(separateSampling, 0, log(nObs)*sum(prop)-sum(prop*log(prop))) # calcul de la vraisemblance

      if ( (l.lik[itrCour+1]-l.lik[itrCour])>tol ) { # 2nd stopping condition
        bPasChan <- FALSE
        aff_cour <- aff_new             # mise a jour des clusters
        aff_cour_long <- aff_new_long
        model.glm <- model.glm.cour
      }
    } else if (length(table(aff_new))<length(prop)) {
      cat('Algorithme did not met convergence\n')
      bPasChan <- TRUE
    }
    if (verbose) { cat('.') } else { }
  }
  if (verbose) { cat(itrCour-1, 'Iterations\n') } else { }
###========================= Fin de ECM
  l.lik <- l.lik[!is.na(l.lik)]
###  return(list(affInit = aff_init, affCour = aff_cour, aff_cour_long=aff_cour_long , nIter = (itrCour-1), ln.likelihood = l.lik[!is.na(l.lik)], model.glm = model.glm ))
  if (effectVar == '') {
    effect <- NULL
  } else {
    effect <- data[, effectVar]
  }
  if (timeParametric) {
    for_ggplot <- data.frame(X = data[, timeVar],Y = model.glm$fitted.values)
  } else {
    for_ggplot <- data.frame(X = unique(data[, timeVar]), Y = coefficients(model.glm))
  }
  
  return(new(Class='GlmCluster', formula=formula.o, nClust=nClust, ident=ident,
             timeVar=timeVar, time = data[, timeVar, drop=T], effectVar=effectVar,
             effect=effect, model.glm=model.glm, timeParametric=timeParametric,
             partition=aff_cour, partition.long=aff_cour_long,
             proportions=round(table(aff_cour)/length(aff_cour), digits = 5),
             criteria=matrix(c(l.lik[length(l.lik)], AIC(model.glm), BIC(model.glm)),
               dimnames=list('',c('log-class-likelihood', 'AIC', 'BIC')), nrow=1),
             converge=new(Class='Converge', nIter=itrCour-1, convergence=TRUE),
             nIter=itrCour-1, for_ggplot = for_ggplot))
}



