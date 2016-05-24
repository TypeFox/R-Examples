
##' Separate the covariates from a 'formula' with a cluster effect from the ones with
##' an identical effect in each cluster if provided.
##'
##' 
##' Given a 'formula' of the form \code{Y ~ clust(T1 + T2 + ...) + pop(X1 + X2 + ...)} or just
##' \code{Y ~ clust(T1 + T2 + ...)}, it returns a list of two or one 'formula' of the form \code{Y ~ T1 + T2 + ...} and \code{ ~ X1 + X2 + ...} if provided.\cr
##' The first element of the list correspond to the covariates with a different effects
##' corresponding to the cluster, the 2nd correspond to covariates having an identical
##' effect in each cluster.\cr
##' In the non-parametric case only \bold{one} covariate is allowed.
##'
##' @title Separate the covariates in a formula
##' @usage seperateFormula(formula)
##' @aliases seperateFormula
##' @param formula  A symbolic description of the model. In the parametric case we write for
##' example \code{ y ~ clust(time+time2) + pop(sex)}, here 'time' and 'time2' will have a different effect according to the cluster, the 'sex' effect is the same for all the clusters. In the non-parametric case only one covariate is allowed.
##' @return A list of 1 or 2 [formula].
##' @include functions4glmClust.R
##' @export
##' @note Meant to be used internally.
seperateFormula <- function(formula) {
  ## test if we have different effects cluster and/or identical in each cluster
  if (any(grepl('clust\\(', deparse(formula) )) && any(grepl('\\+ pop\\(', deparse(formula)))) {
    formula_ch <- unlist(strsplit(deparse(formula), '\\) \\+ p?o?p?u?l?a?t?i?o?n?\\('))
    formula_ch[1] <- sub('clust\\(', '', formula_ch[1])
    formula_ch[2] <- paste('~ ', sub('\\)', '', formula_ch[2]))
    form <- c(formula(formula_ch[1]), formula(formula_ch[2]))
  } else if( any(grepl('clust\\(', deparse(formula)))) {
    formula_ch <- sub('clust\\(', '' ,deparse(formula))
    formula_ch <- sub('\\)', '', formula_ch)
    form <- list(formula(formula_ch))
  } else {
    stop('The formula must be of the form : \'y ~ clust(A+B+..) + pop(C+D+..)\', \'y ~ clust(A+B+..)\' or \'y ~ clust(A)\'\n')
  }
  return(form)
}

##' Rewrite a given formula with all the covariates, so do note have to write them all.
##' @title Rewrite the formula with all the covariates
##' @aliases rwFormula
##' @param formula : An object of type 'formula' of the form \code{y ~ .}
##' @param col.names : Name of the columns in the data.
##' @param ident : Name of the identity column in the data.
##' @return : A 'formula' with all covariates.
##' @export
##' @note Meant to be used internally.
rwFormula <-function(formula, col.names, ident) {
  var <- as.character(formula)[2] # y
  cov <- paste(col.names[-charmatch(c(var, ident), col.names)], sep='', collapse = '+')
  formula <- as.formula(paste(var, cov, sep = '~'))
  return(formula)
}

##'  Write the new formula given the covariates with a level cluster effect, the number of
##' clusters and the type parametric on time or not.
##'
##' Given the covariates and the number of clusters, it returns a character string which will be
##' converted inside \code{\link{glmClust}} to a formula to represent the covariates with different
##' cluster effect.
##' @usage addIndic(covar, nClust, parametric = TRUE, nomClust = 'G')
##' @title Create the new formula with the indicator covariates
##' @aliases addIndic
##' @param covar A 'character' vector of covariates.
##' @param nClust Number of clusters.
##' @param parametric If [TRUE] it means we are parametric on time.
##' @param nomClust The beginning of the name of the indicator covariates 'Ga, Gb, .., etc.
##' @return A character string which will be used as a 'formula'.
##' @export
##' @note Meant to be used internally.
addIndic <- function(covar, nClust, parametric = TRUE, nomClust = 'G') {
  exp <- ''
  for (i in letters[1:nClust]) {
      exp <- paste0(exp,'+', paste0(nomClust, i, ':', covar, collapse = '+' ))
  }
  exp_fin <- sub('\\+', '', exp)
  if (parametric) {
      exp_fin <- paste0(paste0(nomClust, letters[1:nClust], collapse = '+'), '+', exp_fin)
  } else {  }
  return(exp_fin)
}

##' This function creates and return a vector containing the name of the coefficients associated
##' to the current cluster.
##'
##' Given the name of the covariates and the number of the current cluster, it constructs a
##' vector used to retrieve the coefficients from a 'glm' object, these coefficients 
##' are used to calculate the predicted values of the current cluster.
##'
##' @title Get the name of the coefficients in the 'glm' object according to the current cluster
##' @usage getNomCoef(covar, cov_fix, nomClust, itrClust, parametric = TRUE)
##' @aliases getNomCoef
##' @param covar A vector of [character] indicating the covariates with a levec cluster effect.
##' @param cov_fix A vector of [character] indicating the covariates with the same effect in
##' each cluster.
##' @param nomClust The beginning of the name of the undicator covariates, by default 'G[letters]'.
##' @param itrClust The number of the current cluster.
##' @param parametric By default [TRUE] for parametric on time.
##' @return A vector of [character] giving the name of the coefficients associated to awith a
##' given cluster.
##' @export
##' @note Meant to be used internally.
getNomCoef <- function(covar, cov_fix, nomClust, itrClust, parametric = TRUE ) {
  if (parametric) {   covar <- c('', covar)   } else { } # dans ce cas on ajoute le terme correspondant a l'intercept dans chaque cluster
  exp_cov <- paste0("'", nomClust, letters[itrClust] ,':', covar, "'", collapse =',')
  exp_cov <- sub('\\:', '', exp_cov)
  if (!missing(cov_fix)) {
    exp_fix <- paste0("'", cov_fix, "'", collapse =',')
    exp_cov <- paste(exp_cov, exp_fix, sep = ',')
  }
  exp_cov <- paste0('c(', exp_cov, ')')
  return(exp_cov)
}

##' Given the covariates and the name of the coefficients corresponding to a given cluster,
##' the function construct a character string, which will be used to calculate the predicted values.
##'
##' To calculate the predicted values in each cluster, we need the values of the covariates in
##' the data and the right coefficients in the 'glm' object. To do this we construct an expression
##' which will be evaluated inside \code{\link{glmClust}}.
##'
##' @title Creates a character string expression to calculate the predicted values
##' @usage predict_clust(cov, nomCoef, nom_model)
##' @aliases predict_clust
##' @param cov Name of the covariates.
##' @param nomCoef Name of the coefficients.
##' @param nom_model Name of the glm model.
##' @return A character string of the expression of the predicted values of a given cluster.
##' @export
##' @note Meant to be used internally.
predict_clust <- function(cov, nomCoef, nom_model) {
  exp_coef <- paste0('coefficients(', nom_model, ')[', "'", nomCoef[-1],"'", ']')
  exp_data <- paste0('data[,', "'", cov , "'", ']')
  exp <- paste(exp_coef, exp_data, sep = '*', collapse = '+')
  exp <- paste0('coefficients(', nom_model, ')[',"'", nomCoef[1],"'" ,']+', exp)
  return(exp)
}


##' The log-likelihood is calculated with taking into account the type of data ('gaussian',
##' 'binomial', ... etc) and the link function.
##'
##' This function calculates the log-likelihood for the exponential
##' \code{\link{family}}, it uses the 'AIC' function to realise this operatin.
##'
##' @aliases log_lik
##' @title Calculate the log-likelihood
##' @usage log_lik(y, n, mu, wt, family, nparam, disp_mod)
##' @param y Observed values.
##' @param n Vector of '1's and same length  as y.
##' @param mu Predicted values.
##' @param wt Weights.
##' @param family An object of class \code{\link{family}}.
##' @param nparam Number of parameters of the model.
##' @param disp_mod Dispersion of the 'glm' model.
##' @return The log-likelihood of an individual (trajectory).
##' @export
##' @note Meant to be used internally.
log_lik <- function(y, n, mu, wt, family, nparam, disp_mod) { # je peux supprimer
### nparam : devra etre corrige mais ce n'est pas grave puisque l'on est a une constante pres
### log(like) = - AIC/2 + k*n.param # k=2 par defaut
    if( family$family == 'binomial' ) {
        mu <- family$linkinv(mu)
        if(!family$validmu(mu)) {
            cat(' ### Invalid predicted valus (mu) ###\n')
            return(NULL)
        }else {
###            return(-family()$aic(y = y, n = n, mu = family()$linkinv(mu), wt = wt, dev = sum(family()$dev(y = y, mu = mu , wt = wt)))/2+2*(nparam+1)) # ) X (
            return(-family$aic(y = y, n = n, mu = mu, wt = wt, dev = family$dev(y = y, mu = mu , wt = wt))/2+2*(nparam+1)) # )X(
        }
    } else if (family$family == 'gaussian') {
        return(sum(dnorm(x = y, mean = family$linkinv(mu), sd = sqrt(disp_mod), log = TRUE)))
    } else {
        return(-family$aic(y = y, n = n, mu = mu, wt = wt, dev = family$dev(y = y, mu = mu , wt = wt))/2+2*(nparam+1)) # par la suite il faudra remplacer 2 par k
    }
}

##' Affect randomly the individuals to the clusters.
##'
##' Affect randomly the individuals to the clusters providing no empty clusters.
##'
##' @title Affect randomly the individuals to the clusters
##' @usage affect_rand(nObs, nClust)
##' @aliases affect_rand
##' @param nObs Number of observations.
##' @param nClust  Number of clusters.
##' @return  A [vector] of length 'nObs' containing the affectation to the clusters.
##' @export
##' @note Meant to be used internally.
affect_rand <- function(nObs, nClust) {
    if (nClust < 2 || nClust > nObs) {
        cat('Incorrect number of clusters and/or observations. \n ')
        return(NULL)
    } else {
        aff_b <- FALSE
        aff_obs <- rep(0, nObs)
        while (!aff_b) {
            aff_obs <- replicate(nObs, sample(nClust, size = 1))
            if (sum(unique(aff_obs))==sum(1:nClust)) { aff_b = TRUE }
            else { }
        }
    }
    return(aff_obs)
}


##' Calculate and return an indicator vector.
##'
##' @title Calculate an indicator vector
##' @usage majIndica(aff_obs, itrClust)
##' @aliases majIndica
##' @param aff_obs Vector of [integer].
##' @param itrClust Number of current cluster.
##' @return An indicator vector of the belonging to a cluster.
##' @export
##' @note Meant to be used internally.
majIndica <- function(aff_obs, itrClust) {
  return(ifelse(aff_obs==itrClust, 1, 0))
}



