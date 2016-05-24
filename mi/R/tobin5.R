# Part of the mi package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# This is superceded by the tobin5() function below
tobin5 <-
  function(mdf, y, f = NULL) {
    if(!is(mdf, "missing_data.frame")) stop("'mdf' must be a 'missing_data.frame'")
    if(!is.character(y)) stop("'y' must be a character string")
    if(length(y) != 1) stop("'y' must have length one")
    if(!(y %in% colnames(mdf))) stop("'y' must be a variable in 'mdf'")
    
    y <- mdf@variables[[y]]
    NAs <- is.na(y)
    to_drop <- mdf@index[[y@variable_name]]
    X <- mdf@X[,-to_drop]
    
    probit <- bayesglm.fit(X, y = NAs, family = binomial(link = "probit"))
    class(probit) <- c("bayesglm", "glm", "lm")
    gamma <- coef(probit)
    IMR_0 <- dnorm(-fitted(probit)) / pnorm(-fitted(probit))
    IMR_1 <- dnorm( fitted(probit)) / pnorm( fitted(probit))
    
    if(is.null(f)) {
      mark <- colnames(mdf@X)[!grepl("^missing_", colnames(mdf@X))][-1]
      mark <- mark[mark != y@variable_name]
      f <- paste(mark, collapse = " + ")
      f <- paste(y@variable_name, " ~ ", f, " + IMR", sep = "")
      f <- as.formula(f)
    }
    else if(!is(f, "formula")) stop("'f' must be 'NULL' or a formula")
    
    df <- as.data.frame(cbind(mdf@X, IMR = IMR_0))
    if(is(y, "continuous")) {
      model_0 <- bayesglm(f, family = gaussian, data = df, subset = !NAs)
    }
    else stop("only continuous dependent variables are supported at the moment")
    
    df <- as.data.frame(cbind(mdf@X, IMR = IMR_1)) 
    if(is(y, "continuous")) {
      model_1 <- bayesglm(f, family = gaussian, data = df, subset = NAs)
    }
    
    se_0 <- model_0$dispersion
    se_1 <- model_1$dispersion
    
    delta_0 <- IMR_0^2 - fitted(probit) * IMR_0
    delta_1 <- IMR_1^2 + fitted(probit) * IMR_1
    
    betaL_0 <- coef(model_0)
    betaL_0 <- betaL_0[length(betaL_0)]
    betaL_1 <- coef(model_1)
    betaL_1 <- betaL_1[length(betaL_1)]
    
    sigma_0 <- sqrt(se_0^2 + (betaL_0 * delta_0)^2)
    sigma_1 <- sqrt(se_1^2 + (betaL_1 * delta_1)^2)
    
    rho_0 <- -betaL_0 / sigma_0
    rho_1 <-  betaL_1 / sigma_1
    
    ## FIXME: correct vcov(model_0) and vcov(model_1) now
    
    return(list(probit = probit, model_0 = model_0, model_1 = model_1,
                rho_0 = rho_0, rho_1 = rho_1))
  }

tobin5 <-
  function(imputations, y, f = NULL) {
    if(!is(imputations, "mi")) stop("'imputations' must be a 'mi' object")
    if(!is.character(y)) stop("'y' must be a character string")
    if(length(y) != 1) stop("'y' must have length one")
    if(!(y %in% colnames(imputations))) stop("'y' must be a variable in 'imputations'")
    
    dfs <- complete(imputations)
    mdf <- imputations@data[[1]]
    to_drop <- mdf@index[[y@variable_name]]
    cn <- colnames(mdf@X[,-to_drop])[-1]
    f1 <- paste(cn, collapse = " + ")
    NAs <- is.na(mdf@variables[[y]])
    if(paste("missing", y, sep = "_") %in% colnames(mdf@X)) {
      f1 <- paste(paste("missing", y, sep = "_"), "~", f1)
    }
    else for(i in seq_along(dfs)) {
      dfs[[i]] <- cbind(dfs[[i]], NAs)
      colnames(dfs[[i]]) <- c(colnames(dfs[[i]]), paste("missing", y, sep = "_"))
    }
    f1 <- as.formula(f1)
    probit <- pool(f1, data = dfs, family = binomial(link = "probit"))
    
    gamma <- sapply(probit@models, coef)
    Pr <- sapply(probit@models, fitted)
    IMR_0 <- apply(Pr, 2, FUN = function(p) dnorm(-p) / pnorm(-p))
    IMR_1 <- apply(Pr, 2, FUN = function(p) dnorm( p) / pnorm( p))
    
    if(is.null(f)) {
      mark <- colnames(mdf@X)[!grepl("^missing_", colnames(mdf@X))][-1]
      mark <- mark[mark != y]
      f <- paste(mark, collapse = " + ")
      f <- paste(y@variable_name, " ~ ", f, " + IMR", sep = "")
      f <- as.formula(f)
    }
    else if(!is(f, "formula")) stop("'f' must be 'NULL' or a formula")
    
    if(!is(mdf@variables[[y]], "continuous")) {
      stop("only continuous dependent variables are supported at the moment")
    }
    
    for(i in seq_along(dfs)) dfs[[i]]$IMR <- IMR_0[,i]
    model_0 <- pool(f, data = dfs, family = gaussian, subset = !NAs)
    
    for(i in seq_along(dfs)) dfs[[i]]$IMR <- IMR_1[,i]
    model_1 <- pool(f, data = dfs, family = gaussian, subset = NAs)
    
    se_0 <- sapply(model_0@models, FUN = function(m) m$dispersion)
    se_1 <- sapply(model_1@models, FUN = function(m) m$dispersion)
    
    delta_0 <- IMR_0^2 - Pr * IMR_0
    delta_1 <- IMR_1^2 + Pr * IMR_1
    
    betaL_0 <- sapply(model_0@models, coef)
    betaL_0 <- betaL_0[nrow(betaL_0)]
    betaL_1 <- sapply(model_1@models, coef)
    betaL_1 <- betaL_1[nrow(betaL_1)]
    
    sigma_0 <- sqrt(se_0^2 + sweep(delta_0, 2, betaL_0, FUN = "*")^2)
    sigma_1 <- sqrt(se_1^2 + sweep(delta_1, 2, betaL_1, FUN = "*")^2)
    
    rho_0 <- -betaL_0 / sigma_0
    rho_1 <-  betaL_1 / sigma_1
    
    ## FIXME: correct vcov(model_0) and vcov(model_1) now
    
    return(list(probit = probit, model_0 = model_0, model_1 = model_1,
                rho_0 = rho_0, rho_1 = rho_1))
  }
