
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



### !!!!!!!!!!!! ACHTUNG !!!!!!!!!!!! TREND als cov-fct muss
### noch programmiert werden !!!

# RFsimulate:  Not implemented yet: If \code{model} is a formula or of class
#    \command{\dQuote{\link{RFformula}}},
#    the corresponding linear mixed model of the type
 #   \deqn{response = W*b + Z*u + e} is simulated

##   source("~/R/RF/RandomFields/R/MLES.R")

## Printlevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information 

## jetzt nur noch global naturalscaling (ja / nen)
## spaeter eine unktion schreibbar, die den naturscaling umwandelt;
##   im prinzipt CMbuild, aber ruechwaers mit 1/newscale und eingefuegt
##   in eventuell schon vorhandene $ operatoren


#Beim paper lesen im Zug nach Muenchen heute morgen ist mir eine Referenz zu einem R Paket "mlegp: Maximum likelihood estimates of Gaussian processes" aufgefallen. Ist Dir aber sicher schon bekannt! 

#  stop("")
  # problem: natscale; im moment 2x implementiert, 1x mal ueber
  # scale/aniso (user) und einmal gedoppelt -- irgendwas muss raus

## LSQ variogram fuer trend = const.
## kann verbessert werden, insb. fuer fixed effects, aber auch eingeschraenkt
## fuer random effects -> BA/MA


## REML fehlt

## users.guess muss in eine List von meheren Vorschlaegen umgewandelt werden !!! Und dann muss RFfit recursiver call mit allen bisherigen Werden laufen !!


## NAs in data mit mixed model grundsaetzlich nicht kombinierbar !
## NAs in data mit trend (derzeit) nicht kombinierbar


## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


###################################
## to do !!! Mixed Model Equations !!! ##
###################################

######################################################################
## kleine Helfer:
######################################################################
  

model2string <- function(model) {
  ans <- paste(model[[1]], "(", sep="")
  if (length(model) > 1) {
    n <- names(model)
    for (i in 2:length(model)) {
      ans <- paste(ans, if (i>2) ",", n[i], "=",
                   if (is.list(model[[i]])) model2string(model[[i]])
                   else model[[i]])
              }              
  }
  return(paste(ans, ")", sep=""))
}
   

#OneTo <- function(n) return(if (n < 1) NULL else 1:n)
OneTo <- function(n) return(if (length(n) > 1) stop("invalid end of for loop") else if (n < 1) NULL else 1:n)




######################################################################
##  Transform S3 fit result into S4 
######################################################################

list2RMmodelFit <- function(x, isRFsp=FALSE,
                            coords, gridTopology, data.RFparams, T) {
  ## final transformations ....
  
  stopifnot(is.list(x),
            all(c("model", "likelihood", "residuals") %in% names(x)))

  
  if (isRFsp) {
    stopifnot(!missing(coords) &&
              !missing(gridTopology) &&
              !missing(data.RFparams))

    
    ## convert residuals to RFsp class    
    err <-
      try({
        lres <- length(x$residuals)
        if (lres > 0) {
          for (i in 1:lres) {
            gT <- if (length(gridTopology) < i) NULL else gridTopology[[i]]
            co <- if (length(coords)<i) NULL else coords[[i]]            
            if (!is.null(x$residuals[[i]])) {
                  x$residuals[[i]]<-
                  conventional2RFspDataFrame(data=x$residuals[[i]],
                                             coords=co,
                                             gridTopology=gT,
                                             n=data.RFparams[[i]]$n,
                                             vdim=data.RFparams[[i]]$vdim,
                                             T = T,
                                             vdim_close_together=FALSE)
              }
          }
        }
      }
          , silent=TRUE)
    
   if(class(err)=="try-error")
      warning(paste("residuals could not be coerced to class'RFsp';",
                    err))
  } # isRFsp

  if (!is.null(x$trend)) stop("'trend' as list currently does not work anymore")

  return(new(ZF_MODELEXT,
             list2RMmodel(x$model),
             formel = x$formel,
             variab = x$variab,
             param = x$param,
             globalvariance = x$globalvariance,
             covariat = x$covariat,
             likelihood = x$likelihood,
             AIC = x$AIC,
             AICc = x$AICc,
             BIC = x$BIC,
             residuals = x$residuals))
}
  

######################################################################
##  lower, upper, users.model must match 'model'. This is checked here.
######################################################################
GetValuesAtNA <- function(NAmodel, valuemodel, x=NULL, skipchecks, ...) {
  stopifnot(is.list(x), !is.list(x[[1]]))
  aux.reg <- MODEL_AUX
  
  info <- neu <- list()
  models <- list(NAmodel, PrepareModel2(valuemodel, ...))
   
  for (m in 1:2) {    
    info[[m]] <- .Call("SetAndGetModelLikeli", aux.reg,
                       list("Cov",ExpliciteGauss(models[[m]])),
                            x, PACKAGE="RandomFields")
    neu[[m]] <- GetModel(register=aux.reg, modus=GETMODEL_DEL_MLE,
                         spConform=FALSE, do.notreturnparam=TRUE)
  }

  NAs <- sum(info[[1]]$NAs)
  ret <- .Call("Take2ndAtNaOf1st", aux.reg, list("Cov", neu[[1]]),
               list("Cov", neu[[2]]), x$spatialdim, x$Zeit, x$spatialdim, 
               NAs, as.logical(skipchecks), PACKAGE="RandomFields")

  
  
  return(ret)
}




######################################################################
##    'optimal' parscale for 'optim'
######################################################################
ParScale <- function(optim.control, current, lower, upper) {
  if (!is.null(parscale <- optim.control$parscale)) {
    if (!is.numeric(parscale)) {
      stop("non numeric parscale not allowed yet")
    } # else is.numeric: nothing to do
  } else parscale <- rep(NA, length(lower))

  if (!missing(current)) {
    idx <- is.finite(parscale)
    parscale[!idx] <- current[!idx]
    parscale[idx] <- sqrt(parscale[idx] * current[idx]) # geom mittel
  }
  
  if (any(idx <- !is.finite(parscale) | parscale == 0)) {
    parscale[idx] <- pmax(abs(lower[idx]), abs(upper[idx])) / 10 # siehe fit_scale_ratio unten
    idx <- idx & (lower * upper > 0)
    parscale[idx] <- sqrt(lower[idx] * upper[idx])
    stopifnot(all(is.finite(parscale)))
  }
  return(parscale)
}




######################################################################
##   used for detecting splitting directions
######################################################################
vary.variables <- function(variab, lower, upper) {
  n.var <- length(lower)
  w.value <- runif(n.var, 3-0.5, 3+0.5) # weight
  n.modivar <- 5

  idx <- lower <= 0
  modilow <- modiupp <- numeric(n.var)
  modilow[idx] <- (lower[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modilow[!idx] <-(lower[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  modiupp[idx] <- (upper[idx] + w.value[idx] * variab[idx]) / (w.value[idx] + 1)
  modiupp[!idx] <-(upper[!idx]*variab[!idx]^w.value[!idx])^(1/(w.value[!idx]+1))
  
  modivar <- matrix(nrow=n.var, ncol=n.modivar)
  for (i in 1:n.var) {
    modivar[i,] <- seq(modilow[i], modiupp[i], length=n.modivar)
  }
  return(modivar)
}


vary.x <- function(rangex) {
  ## letzter Eintrag ist 0
  npt <- 5
  totdim <- ncol(rangex)
  x <- matrix(0,ncol=npt+1, nrow=totdim)
  for (i in 1:totdim) {
    x[i, 1:npt] <- runif(npt, rangex[1, i], rangex[2, i])
  }
  x
}



ModelSplitXT <- function(splitReg, info.cov, trafo, variab,
                         lower, upper, rangex, modelinfo, model,
                         p.proj=1:length(variab), v.proj=1,
                         report.base) {
  
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  xdimOZ <- modelinfo$xdimOZ
  dist.given <- modelinfo$dist.given

  if (ts.xdim == 1) return(list(p.proj=p.proj,
        x.proj = TRUE,
        v.proj=v.proj))
  lp <- length(variab[p.proj]) # not length(p.proj) sicherheitshalber
  statiso <- info.cov$trans.inv && info.cov$isotropic
   
  abs.tol <- 0

  truely.dep <- dep <- matrix(FALSE, nrow=lp, ncol=ts.xdim)
  varyx <- vary.x(rangex=rangex)
  varyy <- vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)
    
  for (d in 1:ts.xdim) {
    x <- varyx
    x[-d, ] <- 0
    if (dist.given) y <- NULL
    else {
      y <- varyy
      y[-d, ] <- 0
    }

    S <- double(ncol(x) * vdim^2)
    for (i in 1:lp) {

      for (m0 in 1:ncol(modivar)) {
        var <- modivar[, m0]
      
        for (m1 in 1:ncol(modivar)) {
          var[p.proj[i]] <- modivar[i, m1]

          .C("PutValuesAtNA", as.integer(splitReg), as.double(trafo(var)),
             PACKAGE="RandomFields")          
          .Call("CovLoc", splitReg, x, y, xdimOZ, ncol(x), S,
                PACKAGE="RandomFields")
          base::dim(S) <- c(ncol(x), vdim, vdim)
         
          if (any(is.na(S)))
            stop("Model too complex to split it. Please set split=FALSE")

          ## Bei x_j=0, j!=i, gibt es Auswirkungen hinsichtlich
          ## des Parameters auf die Kov-Fkt?
          ## d.h. aendern sich die kov-Werte bei variablem Parameter,
          ## obwohl die x_j = 0 sind?
          if (m1==1) Sstore <- S[, v.proj, v.proj, drop=FALSE]
          else {
            difference <- S[,v.proj, v.proj, drop=FALSE] - Sstore
            if (any(abs(difference) > abs.tol)) {
              dep[i, d] <- TRUE

              ## ist die Kovarianzfunktion der Bauart C_1(x_1) + v C_2(x_2),
              ## und d=1, so bildet v eine Konstante bzg. C_1(x_1)        

              ## i.e.,
              ## note: s_1 C_1(x) + s_2 C_2(t)
              ## then s_2 dep(ends) on x, but not truely
              truely.dep[i, d] <- truely.dep[i, d] || 
                 !(all(apply(difference, 2:3, function(x) all(diff(x) == 0))))
            }
          }
        }
      }
    }

  }

  modelsplit <- NULL
  untreated_x <- rep(TRUE, ts.xdim)
  untrtd_p <- rep(TRUE, lp)
  
  for (d in 1:ts.xdim) {
    if (!untreated_x[d]) next
    rd <- dep[, d]
    if (!any(rd)) next
    same.class <- which(apply(dep == rd, 2, all))
    untreated_x[same.class] <- FALSE

    ## untreated_x wird false falls rd und der Rest keine Abhaengigkeiten
    ## mehr gemein haben
    untreated_x[d] <- !all(!rd | !dep[, -same.class, drop=FALSE])
 
    ## falls es Parameter gibt, die nur zur Koordinate d gehoeren,
    ## werden diese geschaetzt, zusammen mit allen anderen Paremetern,
    ## die zur Koordinate d (und zu anderen Koordinaten) gehoeren 
    trtd <- apply(rd & !dep[, -same.class, drop=FALSE], 1, all)
    if (any(trtd)) {
      untrtd_p <- untrtd_p & !trtd
      report <- "separable"
      modelsplit[[length(modelsplit) + 1]] <-
        list(p.proj = p.proj[dep[, d]], 
             v.proj = v.proj,
             x.proj=as.integer(same.class),
             report = (if (missing(report.base)) report else
                       paste(report.base, report, sep=" : "))
             )
    }
  }

  # !truely.dep, die moegicherweise rausgeflogen waren
  nontruely <- apply(!truely.dep & dep, 1, any)
  if (any(nontruely)) {
    
    if (sum(nontruely) > 1) stop("more than 1 parameter is detected that is not truly dependent -- please use split=FALSE")
      ## Problem ist hier die nicht Eindeutigkeit bei
    ## C = v_1 C_1(x_1) + v_2 C_2(X_2) + v_2 C_3(x_3)
    ## da nun v_2 und v_3 als nontruely auftauchen und somit
    ## dass auf x_1 projezierte Modell nicht mehr identifizierbar ist.
    ## Letztendlich muessen die beiden (oder mehrere) zusammengefasst werden
    ## und in recursive.estimation, use.new.bounds adaequat beruecksichtigt
    ## werden 
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj[nontruely])
    modelsplit[[l + 2]] <-
      list(p.proj = p.proj[nontruely], v.proj = v.proj,
           x.proj=which(apply(nontruely & dep, 2, any))
           ) ## oder vielleicht doch TRUE oder durchzaehlen fuer
    ##                     verschiedene Ebenen
    untrtd_p <- untrtd_p & !nontruely
  }

  ## Parameter die verwoben sind:
  if (any(untrtd_p & !nontruely)) {
    l <- length(modelsplit)
    report <- "simple space-time"
    modelsplit[[l + 1]] <- list(use.new.bounds = which(untrtd_p))
    modelsplit[[l + 2]] <-
      list(p.proj=which(untrtd_p), v.proj=v.proj,
           x.proj=which(apply(untrtd_p & dep, 2, any)),
           report = (if (missing(report.base)) report else
                     paste(report.base, report, sep=" : ")))
  }

  ## komplett, falls notwendig
  if ((any(nontruely) || any(untrtd_p)) && !all(untrtd_p) && !all(nontruely)) {
    l <- length(modelsplit)
    modelsplit[[l + 1]] <- list(use.new.bounds=p.proj)
    modelsplit[[l + 2]] <- list(p.proj=p.proj, x.proj=1:ts.xdim, v.proj=v.proj)
  }

  if (length(modelsplit) == 1) {
    modelsplit <- modelsplit[[1]]
    modelsplit$report <- NULL
  }

  return(modelsplit)
}


ModelSplit <- function(splitReg, info.cov, trafo, variab,
                       lower, upper, rangex, modelinfo, model) {
  vdim <- modelinfo$vdim
  tsdim <- modelinfo$tsdim
  ts.xdim <- modelinfo$ts.xdim
  dist.given <- modelinfo$dist.given
  xdimOZ <- modelinfo$xdimOZ
  refined <- modelinfo$refined

  abs.tol <- 0
  restrictive <- FALSE

  lp <- length(variab)
  statiso <- info.cov$trans.inv && info.cov$isotropic
 # if (xor(is.dist, statiso)) stop("mismatch in ModelSplit -- contact author")

  if (vdim == 1) {
    if (statiso) return(NULL)
    modelsplit <- ModelSplitXT(splitReg=splitReg, info.cov=info.cov,
                               trafo=trafo, variab=variab, lower=lower,
                               upper=upper, rangex=rangex,
                               modelinfo=modelinfo, model=model)

    if (is.list(modelsplit[[1]])) {
      l <- length(modelsplit)
      modelsplit[[l+1]] <- list(use.new.bounds=1:lp)
      modelsplit[[l+2]] <-
        list(p.proj=1:lp,
             x.proj=if (ts.xdim==1 && tsdim > 1) TRUE else 1:ts.xdim,
             v.proj=1:vdim)
      return(modelsplit)
    } else return(NULL)
  }
  
  
  overlap <- dep <- matrix(FALSE, nrow=lp, ncol=vdim)

  x <- vary.x(rangex=rangex)
  S <- double(ncol(x) * vdim^2)
  if (!dist.given) y <- vary.x(rangex=rangex)

  modivar <- vary.variables(variab=variab, lower=lower, upper=upper)

  for (i in 1:lp) {

    for (m0 in 1:ncol(modivar)) {
      var <- modivar[, m0]
      
      for (m1 in 1:ncol(modivar)) {
        var[i] <- modivar[i, m1]
       .C("PutValuesAtNA", splitReg, trafo(var), PACKAGE="RandomFields")
       .Call("CovLoc", splitReg, x, if (dist.given) NULL else y,
              xdimOZ, ncol(x), S, PACKAGE="RandomFields")
        base::dim(S) <- c(ncol(x), vdim, vdim)
        if (any(is.na(S)))
          stop("model too complex to split it. Please set split=FALSE")

         if (m1==1) Sstore <- S
        else {
          for (n in 1:vdim)
            dep[i, n] <-
              dep[i, n] | any(abs(S[,n,n] - Sstore[,n,n]) > abs.tol)
        }
      }
    }
  }

  modelsplit <- list()
  untreated_v <- logical(vdim)
  report <- "indep. multivariate"
  
  for (v in 1:vdim) {
    d1 <- dep[, v]
    if (!any(d1)) next
    overlap[, v] <- apply(dep[, -v, drop=FALSE] & d1, 1, any)
    untreated_v[v] <- if (restrictive) !any(overlap[, v])
             else (sum(overlap[, v]) < sum(d1))
    if (untreated_v[v]) {
      idx <- length(modelsplit) +1
       modelsplit[[idx]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj=which(d1), v.proj = v, report.base=report)
      modelsplit[[idx]]$report <- report
    }
  }
 
  ## bislang oben nur auf univariat heruntergebrochen
  ## man koennte noch auf bivariate runterbrechen
  ## dies aber erst irgendwann;
  ## wuerde ich auch zu weiteren report-Ebenenen 2,3, etc. fuehren
  
  if (any(untreated_v)) {    
    overlapping <- apply(overlap[, untreated_v, drop=FALSE], 1, any)
    depending <- apply(dep[, untreated_v, drop=FALSE], 1, any)

    if (refined && any(depending) && any(overlapping)) {
      idx <- length(modelsplit)
      depend <- which(depending | overlapping)
      notok <- which(!depending & !overlapping)
      modelsplit[[idx + 1]] <- list(use.new.bounds=depend, fix=notok)
      report <- "simple multivariate"
      modelsplit[[idx + 2]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo,
                     variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = depend, v.proj = 1:vdim, report.base=report)
      modelsplit[[idx+2]]$report <- report
    } else notok <- which(!depending | overlapping)
 
    if (length(notok) > 0) {
      idx <- length(modelsplit) + 1
      modelsplit[[idx]] <-
        ModelSplitXT(splitReg=splitReg, info.cov=info.cov, trafo=trafo, variab=variab,
                     lower=lower, upper=upper, rangex=rangex,
                     modelinfo=modelinfo, model=model,
                     p.proj = notok, v.proj = 1:vdim)
    }
  }

  l <- length(modelsplit)
  modelsplit[[l+1]] <- list(use.new.bounds=1:lp)
  modelsplit[[l+2]] <-
    list(p.proj=1:lp,
         x.proj=if (ts.xdim==1 && tsdim > 1) TRUE else 1:ts.xdim,
         v.proj=1:vdim)

  modelsplit[[1]]$all.p <- dep

  return(modelsplit) ## return(NULL) # 11.9.12
}



recurs.estim <- function(split, level, splitReg, Z = Z,
                         lower, upper, users.lower, users.upper, guess,  
                         lsq.methods, mle.methods,
                         optim.control,
                         transform, trafo,
                         spConform,
                         practicalrange,
                         printlevel,
                         fit,
                         minmax,
                         sdvar
                         ) {
  M <- if (length(mle.methods) >= 1) mle.methods[length(mle.methods)]
       else lsq.methods[length(lsq.methods)]

  
  w <- 0.5
   sets <- length(Z$data)

  if (printlevel >= PL_FCTN_DETAILS) Print(split) #
  submodels <- NULL
  submodels_n <- 0
  fixed <- FALSE

  for (s in 1:length(split)) {
    
    sp <- split[[s]]    
    if (!is.null(p <- sp$use.new.bounds)) {
      if (printlevel >= PL_STRUCTURE) {
        cat("    calculating new lower and upper bounds for the parameters ",
            paste(p, collapse=", "), "\n", seq="")
      }
       
      if (!is.null(sp$fix)) {
        fix.zero <- (minmax[sp$fix, MINMAX_PMIN] <= 0) &
          (minmax[sp$fix, MINMAX_PMAX] >=0)
        fix.one <- (minmax[sp$fix, MINMAX_PMIN] <= 1) &
          (minmax[sp$fix, MINMAX_PMAX] >= 1)
        if (!all(fix.zero | fix.one))
          stop("Some parameters could not be fixed. Set 'split_refined = FALSE")
        fix.one <- sp$fix[fix.one & !fix.zero] ### zero has priority
        fix.zero <- sp$fix[fix.zero]
        
        guess[fix.one] <- 1
        guess[fix.zero] <- 0
        
        fixed <- TRUE ## info for the next split[[]] that some variables
        ##               have been fixed. The values must be reported in the
        ##               intermediate results for AIC and logratio test
      }
      next
    }

   
    if (is.list(sp[[1]])) {
     
      res <- recurs.estim(split=sp, level=level+1, splitReg=splitReg,
                          Z = Z,
                          lower=lower,
                          upper=upper,
                          users.lower=users.lower,
                          users.upper=users.upper,
                          guess= guess,  
                          lsq.methods=lsq.methods,
                          mle.methods=mle.methods,
                          optim.control=optim.control,
                          transform=transform,
                          trafo=trafo,
                          spConform = spConform,
                          practicalrange = practicalrange,
                          printlevel=printlevel,
                          fit = fit,
                          minmax = minmax,
                          sdvar = sdvar
                          )
      if (!is.null(sp$report)) {
        submodels_n <- submodels_n + 1
        if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
        sumodels[[submodels_n]] <- res
      }
    } else { # !is.list(sp[[1]])
      if (printlevel>=PL_RECURSIVE) {
        cat("splitting (",
            paste(rep(". ", level), collapse=""), format(s, width=2),
            ") : x-coord=", paste(sp$x, collapse=","),
            "; compon.=", paste(sp$v, collapse=","), sep="")
        cat( "; parameters=", paste(sp$p, collapse=", "), sep="")
        cat(" ")
      }
      guess <- pmax(lower, pmin(upper, guess, na.rm=TRUE), na.rm=TRUE)

      .C("PutValuesAtNAnoInit",
         splitReg, trafo(guess), PACKAGE="RandomFields", NAOK=TRUE) # ok
      old.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE ,
                            do.notreturnparam=TRUE)
      
      guess[sp$p.proj] <- NA
      .C("PutValuesAtNAnoInit",
         splitReg, trafo(guess), PACKAGE="RandomFields", NAOK=TRUE) # ok
      new.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
                            spConform=FALSE, do.notreturnparam=TRUE)

      if (!all(is.na(lower))) {
        .C("PutValuesAtNAnoInit",
           splitReg, trafo(lower),PACKAGE="RandomFields", NAOK=TRUE) # ok
        lower.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
                                spConform=FALSE, do.notreturnparam=TRUE)
      } else lower.model <- NULL
      
     if (!all(is.na(upper))) {
       .C("PutValuesAtNAnoInit",
          splitReg,  trafo(upper), PACKAGE="RandomFields", NAOK=TRUE) # ok
       upper.model <- GetModel(register=splitReg, modus=GETMODEL_DEL_MLE,
                               spConform=FALSE, do.notreturnparam=TRUE)
     } else upper.model <- NULL
      
      if ((new.vdim <- length(sp$v.proj)) < Z$vdim) {
        m <- matrix(0, nrow=new.vdim, ncol=Z$vdim)
        for (j in 1:new.vdim) m[j, sp$v.proj[j]] <- 1
        new.model   <- list("M", M=m, new.model)
        lower.model <- list("M", M=m, lower.model)
        upper.model <- list("M", M=m, upper.model)
        old.model   <- list("M", M=m, old.model)
        new.data <- list()
        vproj <- rep(FALSE, Z$vdim)
        vproj[sp$v.proj] <- TRUE

        for (i in 1:sets) {
          new.data[[i]] <- Z$data[[i]][ , vproj, drop=FALSE]
        }
      } else {
        new.data <- Z$data
        vlen <- length(sp$v.proj)
        for (i in 1:length(new.data)) {
          base::dim(new.data[[i]]) <-
            c(Z$coord[[i]]$restotal, vlen,
              length(new.data[[i]]) / (Z$coord[[i]]$restotal * vlen)) 
        }
      }

      
      x.proj <- sp$x.proj
      ignored.x <- if (is.logical(x.proj)) integer(0) else (1:Z$tsdim)[-x.proj]
      if (Z$Zeit) {
        stopifnot(!is.logical(x.proj))
        ignore.T <- all(x.proj != Z$tsdim)
        x.proj <- x.proj[x.proj != Z$tsdim]
        ignored.x <- ignored.x[ignored.x != Z$tsdim]
        lenT <- sapply(Z$coord, function(x) x$T[3])
      } else {
        lenT <- rep(1, length(new.data))
      }
  
      if (!Z$dist.given && !is.logical(sp$x.proj) &&
          length(sp$x.proj) < Z$tsdim)  {
        if (Z$coord[[1]]$grid) {
          new.x <- Z$coord    
          for (i in 1:length(new.data)) {
            len <- new.x[[i]]$x[3, ]
            d <- base::dim(new.data[[i]])
            new.d <-  c(len, d[-1])
            base::dim(new.data[[i]]) <- new.d
            new.data[[i]] <-
              aperm(new.data[[i]], c(x.proj, length(new.d) + (-1:0), ignored.x))
            repet <- d[3] * prod(len[ignored.x]) ## = repet * ignored
            base::dim(new.data[[i]]) <- c(prod(len[x.proj]), d[2], repet)
            new.x[[i]]$x[, ignored.x] <- c(0, 1, 1)
            new.x[[i]]$restotal <- prod(len)
          }
        } else { ## not grid
          dummy <- new.data
          new.x <- new.data <- list()

          xlen <- sapply(Z$coord, function(x) nrow(x$x))

#          Print(dummy, xlen, ignored.x)
          
          for (i in 1:length(dummy)) {
            d <- base::dim(dummy[[i]])
            base::dim(dummy[[i]]) <-
              c(xlen[i], lenT[i], length(dummy[[i]]) / (xlen[i] * lenT[i])) 
            xyz <- Z$coord[[i]]$x
            ## problem, das hier geloest wird: wenn coordinaten auf 0 gesetzt
            ## werden, muessen die Punkte in diesesn Koordinaten uebereinstimmen.
            ## deshalb werden untermengen gebildet, wo dies der Fall ist
            while(length(dummy[[i]]) > 0) {
              slot <- xyz[1, ignored.x]              
              xyz.ignored <- xyz[-1, ignored.x, drop=FALSE]
              if (length(xyz.ignored) > 0) {
                m <- colMeans(xyz.ignored)
                m[m == 0] <- 1 ## arbitrary value -- components are all 0 anyway

                #Print(xyz.ignored, slot, m, ignored.x, xyz, dummy[[i]])
              
                idx <- colSums(abs(t(xyz.ignored) - slot) / m)
                idx <- c(TRUE, idx < min(idx) * 5) ## take also the first one
              } else idx <- TRUE
              last <- length(new.x) + 1
              new.x[[last]] <- Z$coord[[i]]
              new.x[[last]]$x <- xyz[idx, , drop=FALSE]
              new.x[[last]]$x[, ignored.x] <- 0
              new.x[[last]]$restotal <- new.x[[last]]$l <- nrow(new.x[[last]]$x)
              new.data[[last]] <- dummy[[i]][idx, , ,drop=FALSE]
              dummy[[i]] <- dummy[[i]][-idx, , , drop=FALSE]
              xyz <- xyz[-idx, , drop=FALSE]
            }
          }
        }
        
      } else {
        new.x <- Z$coord
      }

 
      if (!is.null(transform)) {
        isna <- is.na(transform[[2]](guess))
        idx <- transform[[1]][isna]
        stopifnot(max(sp$p.proj) <= length(guess))
                
        f <- function(p) {
          q <- guess
          q[sp$p.proj] <- p
          z <- (transform[[2]](q))[isna]
          z     
        }
      
        
        new.transform <- list(idx, f)
      } else new.transform <- NULL

      general_spConform <- if (s<length(split)) FALSE else spConform

      if (Z$Zeit && ignore.T) {        
        for (i in 1:length(new.x)) new.x[[1]]$T <- double(0)
      }

#     if (length(dim(new.data[[1]])) > 2)
 #       dim(new.data[[1]]) <- c(dim(new.data[[1]])[1], prod(dim(new.data[[1]])[-1]))
#       Print(new.model, new.x, new.data, split, s, sp,
#             length(dim(new.data[[1]]))); 
#      stopifnot(length(dim(new.data[[1]])) == 2)
 #     readline()
 #     stopifnot(s < 1)

      Z1 <- StandardizeData(model=new.model, x=new.x, data=new.data,
                           RFopt=RFoptions())

      res <-
        rffit.gauss(Z=Z1,
                    lower= if (!all(is.na(lower))) lower.model,
                    upper= if (!all(is.na(upper))) upper.model,
                    mle.methods=mle.methods,
                    lsq.methods=lsq.methods,                        
                    users.guess=old.model,  
                    optim.control=optim.control,
                    transform=new.transform,
                    recall = TRUE,
                    fit.split = FALSE,
                    fit.factr = if (s<length(split)) fit$factr_recall
                    else fit$factr,
                    fit.pgtol = if (s<length(split)) fit$pgtol_recall
                    else fit$pgtol,
                    general.practicalrange = practicalrange,
                    general.spConform = general_spConform,
                    sdvar = if (length(sp$v.proj) > 1)
                    sdvar[sp$p.proj, sp$v.proj, drop=FALSE]
                    )

      

      if (!is.null(sp$report)) {
        submodels_n <- submodels_n + 1
        if (general_spConform) {
          res@p.proj <- as.integer(sp$p.proj)
          res@v.proj <- as.integer(sp$v.proj)
          res@x.proj <- sp$x.proj
          res@report <- sp$report
          res@true.tsdim <- as.integer(Z1$tsdim)
          res@true.vdim <- as.integer(Z1$vdim)
          if (fixed) res@fixed <- list(zero=fix.zero, one=fix.one)
        } else {
          res$p.proj <- as.integer(sp$p.proj)
          res$v.proj <- as.integer(sp$v.proj)
          res$x.proj <- sp$x.proj
          res$report <- sp$report
          res$true.tsdim <- as.integer(Z1$tsdim)
          res$true.vdim <- as.integer(Z1$vdim)
          if (fixed) res$fixed <- list(zero=fix.zero, one=fix.one)
        }
        submodels[[submodels_n]] <- res
      }
        
     

      table <- if (general_spConform) res@table else res$table  
      np <- length(sp$p.proj)
      
      lower[sp$p.proj] <- pmin(na.rm=TRUE, lower[sp$p.proj],
                               table[[M]][(np + 1) : (2 * np)])
      upper[sp$p.proj] <- pmax(na.rm=TRUE, upper[sp$p.proj],
                               table[[M]][(2 * np + 1) : (3 * np)])

      
      if (s==length(split)) {
        if (submodels_n >= 1) {
          if (general_spConform) res@submodels <- submodels
          else res$submodels <- submodels
        }
        return(res)
      }
    } # else, (is.list(sp[[1]])) 
 
    guess[sp$p.proj] <- res$table[[M]][1:length(sp$p.proj)]    
    fixed <- FALSE
     
  } # for s in split

  return(res)
} # fit$split




######################################################################
##  main function for fitting Gaussian processes
######################################################################

rffit.gauss <- function(Z, lower=NULL, upper=NULL,
                        mle.methods=MLMETHODS,
                        lsq.methods= LSQMETHODS,
                        ## "internal" : name should not be changed;
                        ## should always be last method!
                        users.guess=NULL,  
                        optim.control=NULL,
                        transform=NULL,
                        recall = FALSE,
                        sdvar = NULL,
                        ...) {



  ## ACHTUNG: durch rffit.gauss werden neu gesetzt:
  ##    practicalrange
  ##

  ## according to MLE.cc, SetAndGetModelInfo


####################################################################
###                      Preludium                               ###
####################################################################
 # internal.examples_reduced = TRUE,
 
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  if (!recall) {
    if (!exists(".Random.seed")) runif(1)
    old.seed <- .Random.seed
    set.seed(0)
    on.exit(set.seed(old.seed), add = TRUE)
  }
  

  if (RFopt$general$modus_operandi == MODENAMES[normal + 1] &&
      RFopt$internal$warn_normal_mode) {
    RFoptions(internal.warn_normal_mode = FALSE)
    message("The modus_operandi='", MODENAMES[normal + 1], "' is save, but slow. If you like the MLE running\nfaster (at the price of being less precise!) choose modus='", MODENAMES[easygoing + 1],
            "' or\neven modus='", MODENAMES[sloppy + 1], "'.")
  }
  
  show.error.message <- TRUE # options
  if (!show.error.message) warning("show.error.message = FALSE")

  
  save.options <- options()
  on.exit(options(save.options), add=TRUE)
  
  general <- RFopt$general
  pch <- general$pch
  printlevel <- general$printlevel
  if (printlevel < PL_IMPORTANT) pch <- ""
   
  fit <- RFopt$fit
  bins <- fit$bins
  nphi <- fit$nphi
  ntheta <- fit$ntheta
  ntime <- fit$ntime
  optimiser <- fit$optimiser

  if (general$practicalrange && fit$use_naturalscaling)
    stop("practicalrange must be FALSE if fit$use_naturalscaling=TRUE")
  if (fit$use_naturalscaling) RFoptions(general.practicalrange = 3)
 
  if (printlevel>=PL_STRUCTURE) cat("\nfunction defintions...\n")
    ##all the following save.* are used for debugging only
  silent <- printlevel <  PL_STRUCTURE   # optimize

  detailpch <- if (pch=="") "" else '+'

### definitions from 
 
  LiliReg <- MODEL_MLE ## for calculating likelihood; either it is
  ##                      not overwritten by 'split' or it is not used
  ##                      anymore later on
  COVreg <- MODEL_LSQ  ## for calculating variogram and getting back models
  split.reg <- MODEL_SPLIT ## for splitting



  
######################################################################
###                function definitions                            ###
######################################################################
    
  # Rcpp need gradient
  nloptr <- NULL ## just to avoid the warning of R CMD check
  if (optimiser != "optim") { # to do : welche optimierer laufen noch;
    ##                          neuen Optimierer einbinden
    #stop("currently unavailable feature")
    optim.call <- optimiser
    if (optim.call == "minqa") optim.call <- "bobyqa"
    else if (optim.call =="pso") optim.call <- "psoptim"
    if (requireNamespace(optimiser, quietly = TRUE)) {
      optim.call <- do.call("::", list(optimiser, optim.call))
    } else {
      stop("to use '", optimiser, "' its package must be installed")
    }
  } else {
    optim.call <- optim
  }
  
  ## optim : standard
  ## optimx: slower, but better 
  ## soma  : extremely slow; sometimes better, sometimes worse
  ## nloptr: viele algorithmen, z.T. gut
  ## GenSA : extrem langsam
  ## minqa : gut, langsamer
  ## pso   : extrem langsam
  ## DEoptim: langsam, aber interessant; leider Dauerausgabe



  nice_modelinfo <- function(minmax) {
    minmax <- minmax[, -MINMAX_NAN, drop=FALSE]
    minmax <- as.data.frame(minmax)
    minmax$type <- TYPEOF_PARAM_NAMES[minmax$type + 1]
    minmax
  }

  INVDIAGHESS <- function(par, fn, control=NULL, ...) {
    if (length(par) == 0) return(list(NULL, NULL))
    if (length(idx <- which("algorithm" == names(control))) > 0)
      control <- control[-idx];
    
    oH <- try(optimHess(par=par, fn=fn, control=control), silent=TRUE)
    
    if (class(oH) != "try-error") {
      zaehler <- 1
      while (zaehler <= 3) {
        e <- eigen(oH)
        values <- Re(e$values) ## some imaginary small values included in eigen
        vectors <- Re(e$vectors)
        if (all(values != 0)) {
          ## solve scheitert am Eigenwert 0
          invH <- vectors %*% (1/values * t(vectors))
          diagH <- -diag(invH)        
          if (any(diagH < 0)) {
            diagH[diagH<0] <- Inf
          }
          return(list(hessian=oH, sd=sqrt(diagH)))
        }
        oH <- oH + rnorm(length(oH), 0, 1e-7)
      }
    }
    if (any(bayes)) return(list(hessian=oH, sd=rep(NA,length(par))))
    else stop("The Hessian matrix is strange and not invertable")
  }

  OPTIM <- function(par, fn, lower, upper, control, optimiser, silent) {
    try(switch(optimiser,
               "optim" = {
               if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
               optim.call(par=par, fn=fn, lower=lower, upper=upper,
                          control=control, method ="L-BFGS-B")
               },
               "optimx" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=fn, lower=lower, upper=upper,
                            control=control, method ="L-BFGS-B")
               },
               "soma" = {
                 optim.call(function(...) - fn(...),
                            bounds=list(min=lower, max=upper),
                            options=list(), strategy="all2one")
               },
               "nloptr" = {
                 if (length(control$xtol_rel) == 0) control$xtol_rel <- 1e-4
                 optim.call(x0=par, eval_f=fn, lb=lower, ub=upper,
                            opts= control[-pmatch(c("parscale", "fnscale"),
                              names(control))] )
               },
               "GenSA" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=function(...) -fn(...),
                            lower=lower, upper=upper,
                            control=control[-pmatch(c("parscale", "fnscale"),
                              names(control))])
               },
               "minqa" =  {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=function(...) -fn(...),
                            lower=lower, upper=upper,
                            control=control[-pmatch(c("parscale", "fnscale"),
                              names(control))])
               },
               "pso" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 optim.call(par=par, fn=fn, lower=lower, upper=upper,
                            control=control[-pmatch(c("fnscale"),
                              names(control))])
               },
               "DEoptim" = {
                 if (length(idx <- which("algorithm" == names(control))) > 0)
                   control <- control[-idx];
                 control <- control[-pmatch(c("parscale", "fnscale"),
                                            names(control))]
                 if (length(control)==0)
                   optim.call(fn=fn, lower=lower,upper=upper)
                 else
                   optim.call(fn=fn, lower=lower,upper=upper, control=control)
               }), silent=silent)
   }

  show <- function(nr, M, OPT, PARAM)
    cat("\n ", M, ", ", switch(nr, "start", "grid ", "re-do"), ": value=",
        format(OPT, dig=6), ", param=", format(PARAM, dig=2), sep="")

  WarningMessage <- function (variab, LB, UB, txt) {
    cat("Note:", txt, ": forbidden values -- if there are too many warnings",
        "try narrower lower and upper bounds for the variables: (",
        paste(variab, collapse=","), ") not in [(",
        paste(LB, collapse=", "),  ") ; (",
        paste(UB, collapse=", "), ")]\n")
  }

  LSQsettings <- function(M) {    
    assign("LSQ.SELF.WEIGHING", M=="self", envir=ENVIR)
    if (!LSQ.SELF.WEIGHING) {
      assign("LSQ.WEIGHTS", weights[[M]], envir=ENVIR)
      if (globalvariance)
        assign("LSQ.BINNEDSQUARE",
               sum(binned.variogram^2 * LSQ.WEIGHTS, na.rm=TRUE),
               envir=ENVIR)
    }
  }

   
  LStarget <- function(variab) {
    variab <- variab + 0## unbedingt einfuegen, da bei R Fehler
    ##                     der Referenzierung !! 16.2.10
    if (printlevel>PL_FCTN_DETAILS) Print(LSMIN, format(variab, dig=20))#
    
  #  if (n.variab==0) return(NA)		##trivial case
    param <- as.double(trafo(variab))

    if (any((variab<LSQLB) | (variab>LSQUB))) {
      ## for safety -- should not happen, older versions of the optimiser
      ## did not stick precisely to the given bounds
      ## 13.12.03 still happens ...
      if (printlevel>=PL_STRUCTURE)
        WarningMessage(variab, LSQLB, LSQUB, "LSQ")      
      assign("BEYOND", BEYOND + 1, envir=ENVIR)
      variab0 <- pmax(LSQLB, pmin(LSQUB, variab))
      penalty <- sum(variab0 - variab)^2
      save <- list(LSMIN, LSPARAM, LSVARIAB)
      assign("LSMIN", +Inf, envir=ENVIR)
      res <- LStarget(variab0)
      res <- res + penalty * (1 + abs(res))
      if (res <= save[[1]]) {              
        assign("LSMIN", save[[1]], envir=ENVIR)
        assign("LSPARAM", save[[2]], envir=ENVIR)
        assign("LSVARIAB", save[[3]], envir=ENVIR)
      } else assign("LSMIN", res, envir=ENVIR)
      return(res)
    }

   .C("PutValuesAtNA", COVreg, param, PACKAGE="RandomFields")

    model.values <- .Call("VariogramIntern", COVreg, PACKAGE="RandomFields")
    
    if (any(!is.finite(model.values))) {
      if (printlevel>=PL_IMPORTANT) {
        message("LSQ missing values!")
        #       Print(info.cov, model.values)
        #        Print(format(variab, dig=20), format(param, dig=20), LSQLB, LSQUB, model.values); print(minmax)
        ## Print(format(variab, dig=20), format(param, dig=20), LSQLB, LSQUB, model.values,  RFgetModelInfo(register=COVreg, level=3, which="internal")); print(minmax)
        #  stop(""); print("del")        
      }
       return(1E300)
    } 

    if (LSQ.SELF.WEIGHING) {
      ## weights = 1/ model.values^2
        gx <- binned.variogram / model.values
        gx <- gx[is.finite(gx)]
       
        if (length(gx) == 0)  {
          res <- Inf
          #print(binned.variogram)
#          print(model.values)
 #         print(param)
  #        print(globalvariance)
   #       print(variab)
          #print(RFgetModelInfo(lev=3, which="internal"))
          stop("The specified model looks odd.")
        }
        if (globalvariance) {
          bgw <- sum(gx^2)
          g2w <- sum(gx)
          res <- bgw - g2w^2 / length(gx)
        } else {
          res <- sum((gx - 1)^2)
        }
    } else {
      if (globalvariance) {
        ## in this case the calculated variogram model is not the one
        ## to which we should compare, but:
        idx <- is.finite(binned.variogram)
        bgw <- sum(binned.variogram * model.values * LSQ.WEIGHTS, na.rm=TRUE)
        g2w <- sum( (model.values^2 * LSQ.WEIGHTS)[idx] )
        res <- LSQ.BINNEDSQUARE - bgw^2/g2w        
      } else {
#        Print("A")
#        print(binned.variogram);
#        print( model.values); ###
        res <- sum((binned.variogram - model.values)^2 * LSQ.WEIGHTS,
                    na.rm=TRUE)
      }
    }

#    Print(res, LSMIN)

    if (res<=LSMIN) {  
      assign("LSMIN", res, envir=ENVIR)
      assign("LSPARAM", param, envir=ENVIR)
      assign("LSVARIAB", variab, envir=ENVIR)
    }

     
    return(res)
  }

  
  MLtarget <- function(variab) {
    ## new version based on C-Code starting wih 3.0.70
   
    if (n.variab > 0) {
      variab <- variab + 0  ## unbedingt einfuegen, da bei R Fehler der Referenzierung !! 16.2.10
      
      if (printlevel>=PL_FCTN_DETAILS ) {
        Print(format(variab, dig=20)) #
        if (printlevel>=PL_FCTN_SUBDETAILS) print(minmax) #
      }
      
      if (any((variab < MLELB) | (variab > MLEUB))) {
        ## for safety -- should not happen, older versions of the optimiser
        ## did not stick precisely to the given bounds
        ## 23.12.03 : still happens
        if (printlevel>=PL_STRUCTURE)
          WarningMessage(variab, MLELB, MLEUB, "MLE")   
        assign("BEYOND", BEYOND + 1, envir=ENVIR)
        penalty <- variab
        variab <- pmax(MLELB, pmin(MLEUB, variab)) 
        penalty <-  - sum(variab - penalty)^2 ## not the best ....
        save <- list(MLEMAX, MLECOVAR, MLEPARAM, MLEVARIAB)
        assign("MLEMAX", -Inf, envir=ENVIR)
        res <- MLtarget(variab)
        res <- res + penalty * (1+ abs(res))       
        if (res > save[[1]]) {
          assign("MLEMAX", save[[1]], envir=ENVIR)
          assign("MLECOVAR", save[[2]], envir=ENVIR)      
          assign("MLEPARAM", save[[3]], envir=ENVIR)
          assign("MLEVARIAB", save[[4]], envir=ENVIR)        
        } else assign("MLEMAX", res, envir=ENVIR)

        return(res )
      }
 
      param <- as.double(trafo(variab))
      .C("PutValuesAtNA", LiliReg, param, PACKAGE="RandomFields")
      options(show.error.messages = show.error.message)
    } else param <- NULL

    if (printlevel > PL_FCTN_SUBDETAILS) {
      cat("\n\nAufruf von MLtarget\n===================\n")
      Print(RFgetModelInfo(register=LiliReg, level=4, which.submodels="call+user"))#
    }


    ##   print(minmax);    Print(param, variab, RFgetModelInfo(register=LiliReg, level=-1, which.submodels="user.but.once+jump"))
 
    ans <- .Call("EvaluateModel", double(0), LiliReg, PACKAGE="RandomFields") 
    res <- ans[1]
    if (is.na(res)) return(1e300)

    ##  Print(ans, MLEMAX, variab, param, MLELB, MLEUB); print(minmax)
  ##  stopifnot(all(is.finite(MLEUB)))
 
   
    if (res >= MLEMAX) {## >= und nicht > da -Inf, -Inf auftreten kann
                        ## z.B. bei autostart
      assign("MLEMAX", res, envir=ENVIR)
      assign("MLECOVAR", ans[-1], envir=ENVIR)      
      assign("MLEPARAM", param, envir=ENVIR)
      assign("MLEVARIAB", variab, envir=ENVIR)

      #Print("mlemax", param, res, variab)
    }

    #    Print(variab, ans)
    if (printlevel>=PL_FCTN_DETAILS ) Print(ans, MLEMAX)

    return(res)
  } # mltarget


  get.residuals <- function(LiliReg) {
    return( .Call("get_logli_residuals", as.integer(LiliReg)))
  }
  
  get.var.covariat <- function(Variab) {
    max <- MLEMAX
    covar <- MLECOVAR
    param <- MLEPARAM
    variab <- MLEVARIAB
    
    assign("MLEMAX", -Inf, envir=ENVIR)
  
    MLtarget(Variab)
    result <- MLECOVAR

 
    assign("MLEMAX", max, envir=ENVIR)
    assign("MLECOVAR", covar, envir=ENVIR)
    assign("MLEPARAM", param, envir=ENVIR)
    assign("MLEVARIAB", variab, envir=ENVIR)
    return(result)
  }

  
  ## to avoid warning on "no visible binding" we define the following
  ## variable that are used in the local functions:
  ENVIR <- environment()
  LSQ.SELF.WEIGHING <- LSQ.WEIGHTS <- LSQ.BINNEDSQUARE <- 
    DO.REML <- DO.REML1 <- RML.A <- RML.data <-
      REML.CORRECTION <- DO.RML1 <- 
          ML.RESIDUALS <- MLEMAX <- MLEINF <- MLECOVAR <- MLEPARAM <- 
            CROSS.DIST <- CROSS.KRIGE <- CROSS.VAR <- CROSSMODEL <-
                LOGDET <- NULL
  BEYOND <- 0

   
######################################################################
###              End of definitions of local functions             ###
######################################################################


######################################################################
###    Initial settings, coord part I (without model info)         ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("\ninitial settings...\n")

  len <- Z$len
  dist.given <- Z$dist.given
  spatialdim <- Z$spatialdim
  time <-  Z$Zeit
  xdimOZ <- Z$xdimOZ
  
  coord <- Z$coord
  C_coord <- trafo.to.C_CheckXT(coord)
  tsdim <- Z$tsdim
  mindistances <- Z$mindist
  maxdistances <-
    sqrt(sum((if (dist.given) Z$rangex[2, ] else Z$rangex[2, ]-Z$rangex[1,])^2))

  #print(summary(Z$coord[[1]]$x)); Print("dsfas", Z$rangex, dist.given)

################    analyses of orginal model        ###############
##### variables needed for analysis of trend, upper and lower input --
##### user cannot know what the internal represenatation is



  if (printlevel>=PL_STRUCTURE) cat("\nfirst analysis of model  ...\n")
  
  info.cov <- .Call("SetAndGetModelLikeli", LiliReg,
                    list("RFloglikelihood", data = Z$data, Z$model),
                    C_coord, PACKAGE="RandomFields")
 
  ## hier zum ersten mal model verwendet wichtig,
  ## da interne Darstellung abweichen kann. Z.B. dass ein optionaler Parameter
  ## auf einen Standardwert gesetzt wird
  ## -- wichtig fuer u.a. GetValuesAtNA

  ### ACHTUNG: Z hat nicht VAR immer noch auf NA (falls globalvariance!!) ?????

 
  modelinfo <- RFgetModelInfo(register=LiliReg, level=2, spConform=FALSE)
  vdim <- modelinfo$vdim

  #print(info.cov); Print(modelinfo); kkk


  NAs <-  info.cov$NAs ## NAs per model
  trans.inv <- info.cov$trans.inv ## note: only with respect to the
  ##              coordinates, mixed effect koennen andere wirkung haben
  isotropic <- info.cov$isotropic
  
  if (dist.given && (!trans.inv || !isotropic)) 
    stop("only domain and isotropic models go along with distances")
  
  minmax <- info.cov$minmax # 4 Spalten: 1:min, 2:max, 3:type, 4:is.nan, not na
  bayes <- as.logical(minmax[, MINMAX_BAYES])
  
  if (printlevel >= PL_SUBIMPORTANT + recall) print(minmax) #

# Print(info.cov); print(minmax); xxxx
  
  stopifnot(sum(NAs) == nrow(minmax))
  n.param <- nrow(minmax)
  ptype <- minmax[, MINMAX_TYPE]
  diag.idx <- which(ptype == DIAGPARAM)
  if (length(diag.idx)>0) minmax[diag.idx[1], MINMAX_PMIN] <- fit$min_diag
  effect <- info.cov$effect

  if (!any(effect == RemainingError))
    stop("there must be an error component in the model")
  if (Z$matrix.indep.of.x.assumed && !info.cov$matrix.indep.of.x)
    stop("x-coordinates are neither given by 'x' nor by 'distances' nor by 'data',\n  but the model seem to require them")
   
  eff <- rep(FALSE, nrow(minmax))  ## lokale Hilfsvariable
  csNAs <- cumsum(c(0, NAs))
  
  if (any(effect >= RandomEffect & effect <= SpVarEffect))
    stop("mixed effects currently not programmed")
   for (k in 1:length(effect)) {
    if (effect[k] >= RandomEffect && effect[k] <= SpVarEffect && NAs[k] > 0)
      eff[(csNAs[k]+1) : csNAs[k]] <- TRUE
  } 
  minmax[ptype == MIXEDVAR & eff, MINMAX_PMIN] <- fit$minmixedvar
  minmax[ptype == MIXEDVAR & eff, MINMAX_PMAX] <- fit$maxmixedvar
  rm("eff")

 
  ## anyFixedEffect <- any(effect == FixedEffect | effect == FixedTrendEffect)
  anyFixedEffect <- any(effect == FixedTrendEffect)
 
  ts.xdim <- xdimOZ + time
  sets <- length(Z$data)
  repet <- as.integer(dummy <- sapply(Z$data, function(x) ncol(x) / vdim[1]))

#  Print(repet, Z, vdim)
  
  stopifnot(all(repet == dummy))
   
  N <- S <- Sq <- 0 
  for (i in 1:sets) {
    N <- N + rowSums(matrix(colSums(!is.na(Z$data[[i]])), ncol=repet[i]))
    S  <- S  + rowSums(matrix(colSums(Z$data[[i]], na.rm=TRUE),
                              ncol=repet[i]), na.rm=TRUE)
    Sq <- Sq + rowSums(matrix(colSums(Z$data[[i]]^2, na.rm=TRUE),
                              ncol=repet[i]), na.rm=TRUE)
  }
  if (vdim == 1) {
    N <- sum(N)
    S <- sum(S)
    Sq <- sum(Sq)
  }
  mean.data <- S / N
  var.data <- Sq / N - mean.data^2
 


######################################################################
## model specific settings, upper & lower
######################################################################
  if (printlevel>=PL_STRUCTURE) cat("\nupper & lower...\n")
  
  trafoidx <- trafo <- NULL
  if (!is.null(transform)) {
    if (!recall) message("Note: if 'transform' is used, 'anisoT' within 'RMS' should be used instead of 'Aniso', if off-diagonal elements of an anisotropy matrix are estimated")
    if (trafo <- is.list(transform) && length(transform)>0) {
      trafoidx <- transform[[1]]
      if (length(transform)==2 && nrow(minmax) == length(trafoidx) &&
          is.logical(trafoidx)) {        
        trafo <- transform[[2]]
      }
    }
  }

  if (!is.null(trafo) && !is.function(trafo)) {
    #if (length(trafoidx) != nrow(minmax) || !is.numeric(trafoidx))   
    .C("PutValuesAtNA", LiliReg, as.double((1:nrow(minmax)) * 1.11),
       PACKAGE="RandomFields", NAOK=TRUE) # ok
    model_with_NAs_replaced_by_ranks <-
      GetModel(register=LiliReg, modus=GETMODEL_DEL_MLE,
               which.submodels="user.but.once+jump")
    cat("transform must be a list of two elements:\n* A vector V of logicals whose length equals the number of identified NAs in the model (see below).\n* A function taking a vector whose length equals sum(V). The output\n  is the original model where the NAs are replaced by real numbers.\n\nThe currently identified NAs are:\n")
    print(minmax[, c(MINMAX_PMIN, MINMAX_PMAX)]) #
    cat("\nwithin the following model (NA positions are given by n.nn) :\n")
    str(model_with_NAs_replaced_by_ranks, vec.len=20) #
    cat("\nHowever, 'RFfit' was called by\ntransform =")
    str(transform) #
    cat("\n")
    if (length(trafoidx) > nrow(minmax))
      cat("Note that the parameters of the trend are optimised analytically, hence they may not be considered, here.\n")
    return(NULL)
  }


#######################   upper,lower,user     ########################
  if (printlevel>=PL_STRUCTURE) cat("\nlower and upper ...\n")
  
  users.lower <- users.upper <- NULL
  if (!is.null(lower)) {
    if (is.numeric(lower)) {
      users.lower <- lower
      if (length(users.lower) != n.param)
        stop("number of parameters of 'lower' does not match the model. Better use the model definition also for 'lower'.")
    } else {
      users.lower <- try(GetValuesAtNA(NAmodel=Z$model, valuemodel=lower,
                                       x = C_coord[[1]], skipchecks=!is.null(trafo),... ))     
      if (!is.numeric(users.lower)) {
        if (printlevel>=PL_IMPORTANT) Print(Z$model, lower, users.lower)#
        stop("'lower' does not match 'model'")
      }
    }
  }
  if (!is.null(upper)) {
    if (is.numeric(upper)) {
      users.upper <- upper
      if (length(users.upper) != n.param)
        stop("number of parameters of 'upper' does not match the model. Better use the model definition also for 'upper'.")
   } else {      
      users.upper <- try(GetValuesAtNA(NAmodel=Z$model, valuemodel=upper,
                                       x = C_coord[[1]], skipchecks=!is.null(trafo), ...))
      if (!is.numeric(users.upper)) {
        if (printlevel>=PL_IMPORTANT) Print(Z$model, upper, users.upper)#
        stop("'upper' does not match 'model'")
      }#
    }
   }
 

  if(!is.null(users.guess)) {
    if (is.numeric(users.guess)) {
      Users.guess <- as.vector(users.guess)
      if (length(Users.guess) != n.param)
        stop("number of parameters of 'users.guess' does not match the model. Better use the model definition also for 'users.guess'.")
    } else {
      Users.guess <-
        (GetValuesAtNA(NAmodel=Z$model, valuemodel=users.guess,
                          x=C_coord[[1]], skipchecks=!is.null(trafo), ...))   
      if (!is.numeric(Users.guess)) {
        if (printlevel>=PL_IMPORTANT)
          Print(Z$model, PrepareModel2(users.guess, ...), Users.guess) #
        stop("'users.guess' does not match 'model'")
      }
    }

  
    if (any(is.finite(Users.guess))) users.guess <- Users.guess
    else users.guess <- NULL


  }


###########################      transform     #######################
## either given bu users.transform + users.min, users.max
## DIESER TEIL MUSS IMMER HINTER GetValuesAtNA STEHEN
  
  if (printlevel>=PL_STRUCTURE) cat("\ntransform ...\n")
  
  lower <- minmax[, MINMAX_PMIN] 
  upper <- minmax[, MINMAX_PMAX]
  delete.idx <- rep(FALSE, length(lower))
  if (is.null(trafo)) {
    if (any(minmax[, MINMAX_NAN]==1)) ## nan, not na
      stop("NaN only allowed if transform is given.")
    trafo <- function(x) x;
  } else { ## is function
    origNAs <- nrow(minmax)
    delete.idx <- !trafoidx
    if (any(minmax[, MINMAX_NAN]==1)) {
      if (any(minmax[delete.idx, MINMAX_NAN] != 1) ||
          any(minmax[delete.idx, MINMAX_NAN] == 1)) 
        stop("NaNs do not match logical vector of transform")
    }

    try(z <- trafo(lower[!delete.idx]))    
    if (!is.numeric(z) || !all(is.finite(z)) || !is.vector(z))
      stop("The transformation does not return a vector of finite numerical values.")
    if (length(z) != length(lower))
      stop("\n'transform' returns a vector of length ",
           length(z), ", but one of length ", origNAs, " is expected.",
           if (length(z) > origNAs) " Note that the parameters of the trend are optimised analytically, hence they may not be considered, here.",
           " Call 'RFfit' with `transform=list()' to get more information on the parameters of the model.\n\n")
  }
  

  ## Achtung which, da upper,lower etc um box-cox-Variable verlaengert
  ## werden koennten !
  SDVAR.IDX <- ptype == SDPARAM | ptype == VARPARAM | ptype == NUGGETVAR
  SIGN.VAR.IDX <- ptype == SIGNEDVARPARAM
  SIGN.SD.IDX <- ptype == SIGNEDSDPARAM
  ALL.SDVAR <- SDVAR.IDX | SIGN.VAR.IDX | SIGN.SD.IDX

  
  MINMAX_COL <- 9
  MINMAX_ROW <-  10
  if (is.null(sdvar)) {
    sdvar <- matrix(FALSE, nrow=n.param, ncol=vdim)
    for (i in 1:vdim) {
      sdvar[ , i] <- (minmax[, MINMAX_COLS] == i |
                      minmax[, MINMAX_ROWS] == i) & SDVAR.IDX
    }    
  } else stopifnot(all(rowSums(sdvar[SDVAR.IDX, ]) >= 1))
  SCALE.IDX <- ptype == SCALEPARAM  ## large capitals 
  var.idx <- which(ptype == VARPARAM)
  sd.idx <- which(ptype == SDPARAM)
  nugget.idx <- which(ptype == NUGGETVAR)
  MIXED.IDX <- which(ptype == MIXEDVAR)
  mixed.idx <-  which(ptype == MIXEDVAR)

  if (vdim ==1) {
    varmin <- varmax <- rep(var.data, n.param)
  } else {
    ## sdvar : matrix of indices
  
    varmax <- apply(sdvar, 1, function(x) if (any(x)) max(var.data[x]) else NA)
    varmin <- apply(sdvar, 1, function(x) if (any(x)) min(var.data[x]) else NA)
    if (printlevel >= PL_IMPORTANT && !recall) {
      if (mean(abs(mean.data)) != 0 &&
          any(log(sd(mean.data) / mean(abs(mean.data))) > 1.5))
        message("Are the average values of the components rather different? If so, it might be\n worth thinking of standardising the values before calling RFfit.\n")
      else if (any(abs(log(var.data / var.data[1])) > 2.0))
        message("The standard deviations of the components are rather different. It might be\n better to standardise the components of the data before calling RFfit.\n")
    }
  } # vdim > 1
 

#################################################################
##############     prepare constants in S, X,etc      ###########
#################################################################
  if (printlevel>=PL_STRUCTURE) cat("\ndistances and data...")

##############         distances              #################
## to do: distances auf C berechnen falls vorteilhaft!


 
##############         Coordinates & data    #################
    ## note: the direct C call needs matrix where points are given column-wise
    ##       whereas the R function CovarianceFct need them row-wise,
    ##                   except for fctcall==CovarianceMatrix
  

  
  ## to do: missing values -- should be distinguished between
  ## lots of missing and not that many?
  ## balanced <- rep(TRUE, sets)
  

  if (vdim>1 && printlevel>=PL_IMPORTANT && !recall)
    message("Due to the covariance model a ", vdim,
        "-variate random field is expected. Therefore, \nthe data matrix",
        " is assumed to consist of ", repet,
        " independent measurements for\neach point.",
        " Each realisation is given as the entries of ", vdim,
        " consecutive \ncolumns.")

 
##############      find upper and lower bounds      #################
  if (printlevel>=PL_STRUCTURE) cat("\nbounds...")
 
  txt <- "lower and upper are both lists or vectors of the same length or NULL"
  lengthmismatch <- "lengths of bound vectors do not match model"
  structuremismatch <- "structures of bounds do not match the model"

  varnames <- minmax.names <- attr(minmax, "dimnames")[[1]]
   
  ## autostart will give the starting values for LSQ
    ## appears if trafo is given. Then better do not search for
    ## automatic bounds

  autostart <- numeric(length(lower))
  neg <- lower <= 0
  autostart[neg] <-  0.5 * (lower[neg] + upper[neg])
  autostart[!neg] <- sqrt(lower[!neg]*upper[!neg])

  idx <- which(SDVAR.IDX)
 
  if (length(idx) > 0) {
    ## lower bound of first model is treated differently!
    ## so the "main" model should be given first!             !!!!!
     
    
    ## lower[idx] <- 0
    ## first.idx <- nugget.idx
    ## if (is.null(first.idx)) first.idx <- var.idx
    ## if (is.null(first.idx)) first.idx <- sd.idx
     lower[idx] <- varmin[idx] / fit$lowerbound_var_factor / length(idx)
   #  lower[idx] <- 0 ## ??
   
    if (fit$lowerbound_var_factor == Inf && length(idx)>1) {
      idx2 <- idx[if (is.null(users.guess)) length(idx) else
                 which(users.guess == max(users.guess, na.rm=TRUE))[1] ]
      lower[idx2] <- varmin[idx] / 1e8
    }
    
    upper[idx] <- varmax[idx] * fit$upperbound_var_factor
    autostart[idx] <- sqrt(varmin[idx] * varmax[idx]) /length(idx)
    
    if (length(sd.idx) > 0) {
      lower[sd.idx] <- sqrt(lower[sd.idx])
      upper[sd.idx] <- sqrt(upper[sd.idx])
      autostart[sd.idx] <- sqrt(autostart[sd.idx])
    }       
  }

  if (any(SIGN.VAR.IDX)) {
    lower[SIGN.VAR.IDX] <-
      -(upper[SIGN.VAR.IDX]<- varmax[SIGN.VAR.IDX] * fit$upperbound_var_factor);
    autostart[SIGN.VAR.IDX] <- 0 
  }
  
  if (any(SIGN.SD.IDX)) {
    lower[SIGN.SD.IDX] <-
      -(upper[SIGN.SD.IDX]<-sqrt(varmax[SIGN.SD.IDX]*fit$upperbound_var_factor))
    autostart[SIGN.SD.IDX] <- 0 
  }
 
#

 
  lb.s.ls.f <- fit$lowerbound_scale_ls_factor
  up.s.f <- fit$upperbound_scale_factor

  #Print(lb.s.ls.f, up.s.f, mindistances, maxdistances)
  
  if (ZF_EARTHCOORD_NAMES[1] %in% coord[[1]]$coordunits) {
    if ("km" %in% coord[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 40 # 40 = approx 7000 / 180
      up.s.f <- up.s.f * 40
    } else if ("miles" %in% coord[[1]]$new_coordunits) {
      lb.s.ls.f <- lb.s.ls.f / 20 # 20 = approx 7000 / 180
      up.s.f <- up.s.f * 20      
    } else stop("unknown coordinate transformation")
  }
  if (any(idx <- ptype == DIAGPARAM)) {
    lower[idx] <- 1 / (up.s.f * maxdistances)
    upper[idx] <- lb.s.ls.f / mindistances
    autostart[idx] <- 8 / (maxdistances + 7 * mindistances)
  }

#  Print("A", lower, upper, any(idx <- ptype == ANISOPARAM),
#        any(idx <- ptype == DIAGPARAM), any(idx <- ptype == ANISOPARAM),
#        lb.s.ls.f, up.s.f)
#  print(minmax)
 
  if (any(idx <- ptype == ANISOPARAM)) {
    if (is.null(trafo))
      warning("The algorithms RandomFields transpose the matrix Aniso to aniso -- this may cause problems when applying transform to the anisotropy parameters. To be safe, use only the parameter anisoT in RMfit.")
    lower[idx] <- -lb.s.ls.f / mindistances
    autostart[idx] <- 0
  }

  if (any(SCALE.IDX)) {
    idx <- which(SCALE.IDX)
    lower[idx] <- mindistances / lb.s.ls.f
    upper[idx] <- maxdistances * up.s.f
    autostart[idx] <- (maxdistances + 7 * mindistances) / 8      
   }

  #Print("Z", lower, upper)

  
###########################        split       #######################
  if (printlevel>=PL_STRUCTURE) cat("\nsplit...")
  if (fit$split > 0 && length(autostart)>=fit$split) {
    stopifnot(fit$split > 1)
    
    new.param <- if (is.null(users.guess)) autostart else users.guess


### Achtung!! delete.idx darf davor nur fuer trafo gesetzt werden!!

    stopifnot(spatialdim == tsdim - time)
    
    if (length(C_coord[[1]]$y) > 0) stop("x/y mismatch -- pls contact author")
    splitxy <- matrix(as.double(1:(spatialdim^2)), ncol=spatialdim)
    split_l <- spatialdim
    if (dist.given) {
      if (xdimOZ == 1) {
        splitxy <- t(as.vector(dist(splitxy)))
        split_l = length(splitxy)
      } else stop("vector-valued distances currently not allowed")
    }

  #  Print(info.cov, modelinfo, Z$dist.given, Z$xdimOZ, C_coord[[1]]);dfdsfs

       
    splitcoord <- list(x=splitxy,                 #0
                       y=if (!dist.given) splitxy else double(0),#1
                       T=C_coord[[1]]$T,          #2
                       grid = FALSE,              #3
                       spatialdim = spatialdim,   #4
                       Zeit = C_coord[[1]]$Zeit,  #5
                       dist.given = dist.given,   #6
                       restotal = spatialdim,     #7
                       l = spatialdim,            #8
                       coordunits = C_coord[[1]]$coordunits,
                       new_coordunits = C_coord[[1]]$new_coordunits)

    ##    Print(C_coord, splitcoord); ooo

#    Print("BBB")
 #   Print(Z$model)
#    Print(splitcoord); 

    .Call("SetAndGetModelLikeli", split.reg, list("Cov", Z$model),
          splitcoord,
          PACKAGE="RandomFields")
    rm("splitcoord")
    
    stopifnot(ncol(Z$rangex) == ts.xdim)
    keep <- !delete.idx
    split <- try(ModelSplit(splitReg=split.reg, info.cov=info.cov, trafo=trafo,
                            variab=new.param[keep],
                            lower=lower[keep], upper=upper[keep],
                            rangex = Z$rangex,
                            ## ts.xdim != tsdim falls distances (x Zeit)
                            modelinfo=list(ts.xdim=ts.xdim, tsdim=tsdim,
                                xdimOZ = xdimOZ, vdim=vdim,
                                dist.given=dist.given,
                                refined = fit$split_refined),
                            model=Z$model))
    #Print("BA", split)

    if (is(split, "try-error")) {
      message("Splitting failed (", split[[1]], "). \nSo, standard optimization is tried")
    } else {
      if (printlevel>=PL_STRUCTURE) cat("\nsplitted...")  

      if (length(split) == 1)
        stop("split does not work in this case. Use split=FALSE")
      if (length(split) > 1) {
        if (printlevel>=PL_STRUCTURE) {
          cat("\n")
        }
        
        if (printlevel >= PL_RECURSIVE) Print(split) #
        
        
        return(recurs.estim(split=split, level=0,  splitReg=split.reg,
                            Z = Z,
                            lower= rep(NA, sum(keep)),
                            upper= rep(NA, sum(keep)),
                            users.lower = {
                              if (is.null(users.lower)) rep(-Inf,sum(keep))
                              else users.lower[keep]
                            },
                            users.upper = {
                              if (is.null(users.upper)) rep(Inf, sum(keep))
                              else users.upper[keep]
                            },
                            guess=new.param[keep], # setzt default werte
                            lsq.methods = LSQMETHODS,
                            mle.methods=mle.methods,
                            optim.control=optim.control,
                            transform=transform,
                            trafo=trafo,
                            spConform = general$spConform,
                            practicalrange = general$practicalrange,
                            printlevel=printlevel,
                            minmax=minmax,
                           fit = RFopt$fit,
                            sdvar = apply(split[[1]]$all.p, 2,
                                function(x) {x[!ALL.SDVAR] <- FALSE; x})
                            ## sdvar in spaltenrichtung vdim, in zeilenrichtung
                            ## die parameter
                            ))
      }
    }
  }

  delete.idx <- which(delete.idx) ## !!
  
 
######################################################################
###                                                                ###
###   check which parameters are NA -- only old style is allowed   ###
###                                                                ###
###     certain combinations of NA allow for faster algorithms     ###
###                                                                ###
###     !is.na(sill) needs special treatment, hence must be        ###
###     identified --- currently deleted; see before 3.0.70        ###
###                                                                ###
###                                                                ###
###     scaling method must be identified                          ###
###                                                                ###
###     some autostart values are calculated                       ###
###                                                                ###
###                                                                ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("\nauto...")  
  
  ssm <-  nonugget <- novariance <- FALSE
  var.global <- var.idx
  
  if (is.null(users.lower)) {
    users.lower <- rep(-Inf, length(lower))
  } else {
    idx <- !is.na(users.lower)
    lower[idx] <- users.lower[idx]
    users.lower[!idx] <- -Inf
  }
  if (is.null(users.upper)) {
    users.upper <- rep(Inf, length(upper))
  } else {
    idx <- !is.na(users.upper)
    upper[idx] <- users.upper[idx]
    users.upper[!idx] <- Inf
  }

  bounds <- minmax[, c(MINMAX_MIN, MINMAX_MAX), drop=FALSE]
  if (any(bayes)) {
    lower[bayes] <- pmax(lower[bayes], minmax[bayes, MINMAX_PMIN])
    upper[bayes] <- pmin(upper[bayes], minmax[bayes, MINMAX_PMAX])    
    bounds[bayes, 1] <- pmax(bounds[bayes, 1], minmax[bayes, MINMAX_PMIN])
    bounds[bayes, 2] <- pmin(bounds[bayes, 2], minmax[bayes, MINMAX_PMAX])
    users.lower[bayes] <- pmax(users.lower[bayes], minmax[bayes, MINMAX_PMIN])
    users.upper[bayes] <- pmin(users.upper[bayes], minmax[bayes, MINMAX_PMAX])
  }

  bounds <- apply(abs(bounds), 1, max)
  paramnames <- varnames
  param.lower <- lower
  param.upper <- upper
  if (length(delete.idx)>0) {
    upper <- upper[-delete.idx]
    lower <- lower[-delete.idx]
    users.lower <- users.lower[-delete.idx]
    users.upper <- users.upper[-delete.idx]
    autostart <-autostart[-delete.idx]
    varnames <- varnames[-delete.idx]
    SCALE.IDX <- SCALE.IDX[-delete.idx]
    SDVAR.IDX <- SDVAR.IDX[-delete.idx]
    ptype <- ptype[-delete.idx]
    bounds <- bounds[-delete.idx]    
  }


  if (any(autostart<lower) || any(autostart>upper)) {
    if (printlevel >= PL_ERRORS)
      Print(cbind(lower, autostart, upper)) #, orig.lower,orig.upper)#
    autostart <- pmin(upper, pmax(lower, autostart))
  }

  
  if (is.null(likeli.info <- .Call("get_likeliinfo", LiliReg)))
    stop("bug in likelihood. Please inform author.")
  n.covariat <- likeli.info$betas
  betanames <- likeli.info$betanames
  globalvariance <- likeli.info$estimate_variance
  sum.not.isna.data <- likeli.info$sum_not_isna_data

  n.variab <- length(lower)

  if (any(idx <- lower >= upper)) {
    lu <- cbind(lower=lower, upper=upper, idx)  
    stop(paste("Some lower bounds for the parameters are greater than ",
               "or equal to the upper bound\n",
               paste(collapse="\n ", dimnames(lu)[[1]], ":",
                     apply(lu, 1, function(x)
                           paste("\tlower=", signif(x[1]),
                                 ",\tupper=", signif(x[2]),
                                 if (x[3]) "  \tnot ok!", sep=""))
                     )
               ))
  }


  fill.in <-  trafo(autostart)
  if (length(MIXED.IDX) > 0) {  ## Aufruf .C(Cov*Loc wird vor optim benoetigt
    ## somit muessen die NA's (jedoch nur) in mixed.var gesetzt sein
    autostart[MIXED.IDX] <- 1.0
  }
#  .C("PutValuesAtNA", R e g, trafo(autostart), PACKAGE="RandomFields")   

  if (n.variab == 0) {
    if (globalvariance) {
      if (RFopt$internal$warn_onlyvar) {
        message("Only the variance has to be estimated (except for some parameter in the linear model part). Note that 'RFlikelihood' does already this job and is much simpler.")
        RFoptions(warn_onlyvar = FALSE)
      }
    } else {
      #Print(n.covariat, options()$warn)
       curwarn <- options()$warn
      options(warn = 1)
      if (n.covariat) warning("No genuine variable has to be estimated (except some parameter in the linear model part). Note that 'RFlikelihood' does already the job.")
      else warning("No genuine variable has to be estimated. You should use 'RFlikelihood' instead.")
      options(warn = curwarn)
    }
  }

  
######################################################################
######################################################################
###                     Estimation part itself                     ###
######################################################################
######################################################################

  ## check optim.control 
  ## parscale will give the magnitude of the parameters to be eliminated
  ##     passed to optim/optimise so that the optimiser eliminates
  ##     values around 1 ##to do: irgendwo in paper die wichtigkeit beschreiben?

  parscale <- ParScale(optim.control, lower=lower, upper=upper) 
  fit.fnscale <- optim.control$fnscale

  
  if (length(optim.control)>0) {
    opt.control <- optim.control
    stopifnot(is.list(opt.control))
    forbidden.param <- c("parscale", "fnscale", "algorithm")
    ## fnscale=-1 turns the problem into a maximisation problem, see below
    forbidden <- which(!is.na(pmatch(names(opt.control), forbidden.param)))
    forbidden.opt <- opt.control[forbidden]  
    if (length(forbidden) > 0)  opt.control <- opt.control[-forbidden]
   } else {
    opt.control <- list()
  }

  if (length(fit$algorithm) > 0 && fit$algorithm != "")
      opt.control$algorithm <- fit$algorithm
  if (length(optim.control)>0) {
    if (length(forbidden.opt$algorithm) > 0)
      opt.control$algorithm <- forbidden.opt$algorithm
  }

  if (fit$optimiser=="optim") {
    if (length(opt.control$pgtol)==0) #not given by user
      opt.control$pgtol <- fit$pgtol
   if (length(opt.control$factr)==0) #not given by user
      opt.control$factr <- fit$factr
  }


###################  preparation  ################
  if (printlevel>=PL_STRUCTURE) cat("\npreparing fitting...")
  ## methods
  formals <- formals()
  allprimmeth <- c("autostart", "users.guess")
  nlsqinternal <- 3 ## cross checked after definition of weights below
  alllsqmeth <- c(LSQMETHODS[-length(LSQMETHODS)],
                  paste("internal", 1:nlsqinternal, sep=""))
              
  allmlemeth <- eval(formals$mle.methods)
  if (length(allmlemeth) != 1 || allmlemeth != "ml")
    stop("reml currently not programmed")

  
  allcrossmeth <- NULL
  allmethods <- c(allprimmeth, alllsqmeth, allmlemeth, allcrossmeth)

  ## how preceding methods have been considered ?
  ## note cm is used again at the very end when error checking
  cm <- cumsum(c(0, length(allprimmeth), length(alllsqmeth),
                     length(allmlemeth), length(allcrossmeth)))
  cm <- cbind(cm[-length(cm)] + 1, cm[-1])
  cm <- apply(cm, 1, function(x) x[1] : x[2])
  names(cm) <- c("prim", "lsq", "mle", "cross")


   methodprevto <-
    if (fit$only_users) {
      list(lsq="users.guess",mle="users.guess",cross="users.guess")
    } else list(lsq=c(cm$prim),
              mle=c(cm$prim, cm$lsq),
              cross=c(cm$prim, cm$lsq, cm$cross)
              )

 
  ## index (start, end) to the various categories of
  ## information to be stored
  IDX <- function(name) { idx <- tblidx[[name]]; idx[1]:idx[2]}
  tblidx <- cumsum(c(0,
                     n.variab, # variables used in algorithm
                     n.variab, # their lower bounds
                     n.variab, # ... and upper bounds
                     n.variab,  # sd of variabs
                     n.param,# param values to be eliminated
                     rep(1, length(allmethods) - length(allprimmeth)),#method
                     ##                                                 score
                     1, 1, 1, ## AIC, AICc, BIC,
                     as.integer(globalvariance),
                     n.covariat, #
                     # whether there has been a global variance to estimated
                     # coeff to eliminated for covariates, i.e.
                     ##           trend parameters
                     n.param ## sd of params
                    ))

 
  tblidx <- rbind(tblidx[-length(tblidx)] + 1, tblidx[-1])
  idx <- tblidx[1, ] > tblidx[2, ]
  tblidx[, idx] <- 0 
 
   
  dimnames(tblidx) <- list(c("start", "end"),                           
                           c("variab", "lower", "upper", "sdvariab",
                             "param", 
                             allmethods[-1:-length(allprimmeth)],
                             "AIC", "AICc", "BIC",
                             "glbl.var",
                             "covariat",
                             "sdparam"
                             ##,  "lowbeta", "upbeta", only used for
                             ## cross-validation
                             ))
  if (tblidx[2, "covariat"] != 0) tblidx[2, "glbl.var"] <- tblidx[2, "covariat"]
  if (tblidx[1, "glbl.var"] == 0) tblidx[1, "glbl.var"] <- tblidx[1, "covariat"]

  maxtblidx <- max(tblidx)
  tblidx <- data.frame(tblidx)

   ## table of all information; col:various method; row:information to method
  tablenames <-
     c(if (n.variab > 0) paste("v", varnames, sep=":"),        
       if (n.variab > 0) paste("lb", varnames, sep=":"),
       if (n.variab > 0) paste("ub", varnames, sep=":"),
       if (n.variab > 0) paste("sdv", varnames, sep=":"),
      minmax.names,
##       if (nrow(minmax) > 0) paste("p", 1:nrow(minmax), sep=":")
  ##                       else minmax.names,
      allmethods[-1:-length(allprimmeth)],
       "AIC", "AICc", "BIC",
      if (globalvariance) "glbl.var",  
      betanames,
      ## do not try to join the next two lines, since both
      ## varnames and betanames may contain nonsense if
      ## n.variab==0 and n.covariat==0, respectively
      if (n.variab > 0)  paste("sd", minmax.names, sep=":")
       )

#  Print(tablenames, allmethods, nrow=maxtblidx, ncol=length(allmethods), tblidx)
  
  param.table <- data.frame(matrix(NA, nrow=maxtblidx, ncol=length(allmethods),
                                   dimnames=list(tablenames, allmethods)))


  fit.fnscale <- if (is.null(fit.fnscale)) rep(NA, length(allmethods)) else
                  -abs(fit.fnscale)


#############################################################
## end preparation; remaining part is elimination  ###########
#############################################################

##################################################
###############    PRIMITIVE METHODS   ###########
##################################################
 
  MLELB <- LSQLB <- lower
  MLEUB <- LSQUB <- upper



#  Print(LSQLB,LSQUB); ooo


 
  ##****************    autostart    *****************
  if (printlevel>=PL_STRUCTURE) cat("\nautostart...")
  M <- "autostart"
  primMethods <- M
  default.param <- param.table[[M]][IDX("variab")] <- autostart


#  Print(autostart, trafo(autostart), param.table[[M]][IDX("param")]); 
  param.table[[M]][IDX("param")] <- trafo(autostart) 
  MLEVARIAB <- autostart
  param.table[[M]][IDX("glbl.var")] <- get.var.covariat(autostart)



  ## ****************    user's guess    *****************
  if (!is.null(users.guess)) {
    M <- "users.guess"
    primMethods <- c(primMethods, M)
    if (length(delete.idx) > 0) users.guess <- users.guess[-delete.idx]
    idx <- users.guess < lower | users.guess > upper
    if (any(idx)) {
      if (recall) users.guess <- NULL
      else if (general$modus_operandi ==  MODENAMES[careless + 1]) {
        lower <- pmin(lower, users.guess)
        upper <- pmax(upper, users.guess)
        MLELB <- LSQLB <- lower
        MLEUB <- LSQUB <- upper
     } else {
        m <- cbind(lower, users.guess, upper, idx)
        dimnames(m) <- list(rep("", length(lower)),
                            c("lower", "user", "upper", "outside"))
        cat("\n")
        print(m) ## nicht loeschen!
        
        stop("not all users.guesses within bounds\n change values of `lower' and `upper' or \nthose of the `lowerbound*.factor's and `upperbound*.factor's")
      }
    }

    if (!is.null(users.guess)) {
      param.table[[M]][IDX("variab")] <- users.guess
      param.table[[M]][IDX("param")] <- trafo(users.guess)
      MLEVARIAB <- users.guess
      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(users.guess)
    }
  }


##################################################
################### Empirical Variogram   ########################

 
##################################################
################### Empirical Variogram   ########################

  ## see above for the trafo definitions
  ##
  ## zuerst regression fit fuer variogram,
  ## dann schaetzung der Parameter, dann berechnung der covariates
  ## T.o.D.o.: kann verbessert werden durch einschluss der Trendschaetzung
  ## Xges auf RFsimulate basieren

 
 lsqMethods <- NULL
  ev <- list()
  if (length(lsq.methods) == 0) {  # not trans.inv
    if (!is.null(lsq.methods)) warning("submethods are not allowed")
  } else { # trans.inv
    if (printlevel>=PL_STRUCTURE) cat("\nempirical variogram ...\n")
    sd <- vector("list", vdim)
    index.bv <- NULL

    
   # str(coord); kkk

    
    if (vdim == 1 && !dist.given) {
       residuals <- .Call("simple_residuals", LiliReg) 
    } else {
      residuals <- list()
      for (i in 1:sets) {
        residuals[[i]] <- Z$data[[i]]

      #  Print(Z, residuals[[i]], Z$coord[[i]]$restotal, vdim, Z$repetitions[i])
        
        base::dim(residuals[[i]]) <-
          c(Z$coord[[i]]$restotal, vdim, Z$repetitions[i])
      }
    }

    bin <-
        if (length(bins)>1) bins
        else c(-1, seq(0, fit$bin_dist_factor * maxdistances, len=bins+1))
    
    binned.variogram <- NULL

    #Print(vdim, dist.given)
    
    for (j in 1:vdim) {
      if (!dist.given) {        
        ev <-
          RFempiricalvariogram(coord,
                               data= if (vdim == 1) residuals else
                               lapply(residuals, function(x) x[, j, ]),
                               bin=bin,
                               phi=if ((spatialdim>=2) && !isotropic) nphi,
                               theta=if ((spatialdim>=3) && !isotropic) ntheta,
                               deltaT=if (Z$Zeit) ntime,
                               spConform=FALSE, boxcox=c(Inf, Inf))
      } else {
        if (Z$xdimOZ != 1) stop("Distance vectors are not allowed.")
        n.bin <- vario2 <- vario <- rep(0, length(bin))
        for (i in 1:sets) {
          W <- residuals[[i]][, j, ]
          lc <- Z$len[i]
          rep <- 2 * Z$repetitions[i]
          k <- 1
          for (g in 1:(lc-1)) {
            for (h in (g+1):lc) {
              idx <- sum(coord[[i]]$x[k] > bin) ## distances
              dW <- W[g, ] - W[h, ]
              n.bin[idx] <- n.bin[idx] + sum(!is.na(dW))
              vario[idx] <- vario[idx] + sum(dW^2, na.rm=TRUE)
              vario2[idx] <- vario2[idx] + sum(dW^4, na.rm=TRUE)
              k <- k + 1
            }
          }
        }
        n.bin <- n.bin * 2
        vario <- vario / n.bin ## Achtung! n.bin ist bereits gedoppelt        
        sdvario <-
          sqrt(pmax(0, vario2 / n.bin / 2.0 - vario^2)) ## numerische Fehler
        sdvario[1] <- vario[1] <- 0
        n.bin[1] <- sum(Z$len)
        centers <- 0.5 * (bin[-1] + bin[-length(bin)])
        centers[1] <- 0
        ev <- list(centers=centers,
                   emp.vario=vario[-length(n.bin)],
                   sd=sdvario[-length(n.bin)],
                   n.bin=n.bin[-length(n.bin)])
      }
      n.bin <- ev$n.bin
      sd[[j]] <- ev$sd # j:vdim; ev contains the sets
      if (is.null(binned.variogram))
        binned.variogram <- array(dim=c(length(ev$emp.vario), vdim, vdim))
      binned.variogram[, j, j] <- ev$emp.vario
    } # j in 1:vdim

#Print("OK", binned.variogram, bin)
    
    if (!(sum(binned.variogram, na.rm=TRUE) > 0) )
      stop("not more than 1 value in empirical variogram that is not NA; check values of bins and bin_dist_factor")

#Print("OK 2")
    

    bin.centers <- as.matrix(ev$centers)

    if (!is.null(ev$phi)) {
      if (spatialdim<2) stop("x dimension is less than two, but phi is given") 
      bin.centers <- cbind(as.vector(outer(bin.centers, cos(ev$phi))),
                           as.vector(outer(bin.centers, sin(ev$phi))))
    }
    if (!is.null(ev$theta)) {
      if (spatialdim<3)
        stop("x dimension is less than three, but theta is given") 
      if (ncol(bin.centers)==1) bin.centers <- cbind(bin.centers, 0)
      bin.centers <- cbind(as.vector(outer(bin.centers[, 1],
                                           cos(ev$theta))),
                           as.vector(outer(bin.centers[, 2],
                                           cos(ev$theta))),
                           rep(sin(ev$theta), each=nrow(bin.centers)))
    } else {
      
      ##  warning("must be ncol()")      
      if (ncol(bin.centers) < spatialdim) { # dimension of bincenter vector
        ##                       smaller than dimension of location space      
        bin.centers <- 
          cbind(bin.centers, matrix(0, nrow=nrow(bin.centers),
                                    ncol=spatialdim - ncol(bin.centers)
                                    ))
      }
    }
#Print("OK 3")
    


    ## es muessen beim direkten C-aufruf die componenten der Punkte
    ## hintereinander kommen (siehe auch variable coord, Xdistance). Deshalb t()
    ## aus den vdim sd's muss ein Gewicht gemacht werden
    
    if (length(sd) > 0)  { 
      evsd <- sapply(sd, function(x) x)^2    
      evsd <-
        if (is.matrix(evsd)) rowSums(evsd/rowSums(is.finite(evsd)), na.rm=TRUE)
        else evsd[!is.finite(evsd)] <- 0    
      evsd[evsd==0] <- 10 * sum(evsd, na.rm=TRUE) ## == "infinity"
      evsd <- as.double(evsd)
    } else evsd <- 1

    bins             <- length(n.bin)
    binned.n         <- as.integer(n.bin)

    weights <- cbind(NA,                      # self
                     rep(1, bins),            # plain 
                     sqrt(binned.n),          # sqrt(#)
                     1 / evsd,                # sd^-1
                     sqrt(bins:1 * as.double(binned.n)), # internal1 # kann sonst
                     ##                   fehler verursachen, da integer overflow
                     bins:1,                  # internal2
                     sqrt(bins:1)             # internal3
                     )
    stopifnot(ncol(weights)==length(alllsqmeth))
    dimnames(weights) <- list(NULL, alllsqmeth)
    weights <- data.frame(weights)

    ##################################################
    #######################  LSQ  ####################
    ##***********   elimination part itself   **********     
    ## find a good initial value for MLE using weighted least squares
    ## and binned variogram
    ##
    ## background: if the number of observations (and the observation
    ## field) tends to infinity then any least square algorithm should
    ## yield the same result as MLE
    ## so the hope is that for a finite number of points the least squares
    ## find an acceptable initial values

    ## advantage of the following way is that the for-loop is run through
    ## in an ordered sense -- this might be useful in case partial results
    ## are reused
    if (printlevel>=PL_STRUCTURE) cat("\nelimination part...")

    # Print("CCC")
    .Call("SetAndGetModelLikeli", COVreg, list("Cov", Z$model),
           C_CheckXT(x = bin.centers, T = ev$T),        
           PACKAGE="RandomFields")
      
    
    LSMIN <- Inf
    lsqMethods <- LSQMETHODS[pmatch(lsq.methods, LSQMETHODS)]
    if (!is.null(lsqMethods) &&
        any(is.na(lsqMethods))) stop("not all lsq.methods could be matched")
    if ("internal" %in% lsqMethods)
      lsqMethods <- c(lsqMethods, paste("internal", 1:nlsqinternal, sep=""))
    
    firstoptim <- TRUE
    for (M in c(lsqMethods[1], alllsqmeth)) {
     
      if (!(M %in% lsqMethods)) next;
      if (printlevel>=PL_STRUCTURE) cat("\n", M) else cat(pch)

      param.table[[M]][IDX("variab")] <- default.param
      
      LSQsettings(M)
     
      LSMIN <- Inf ## must be before next "if (n.variab==0)"
      LSPARAM <- LSVARIAB <- NA 
      ##      if (n.variab == 0) {
 #       warning("trivial case may cause problems")
 #     } else {
      param.table[[M]][IDX("lower")] <- LSQLB
      param.table[[M]][IDX("upper")] <- LSQUB
      options(show.error.messages = show.error.message) ##

      if (n.variab == 0) {
        LStarget(param.table[IDX("variab"), methodprevto$lsq[1]])
      } else {
        min <- Inf
        min.variab <- NULL
        for (i in methodprevto$lsq) { ## ! -- the parts that change if
          ##                               this part is copied for other method
          if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {
            value <- LStarget(variab) ##

             if (is.finite(value)) {
              param.table[tblidx[[M]][1], i] <- value
              if (value < min) {
                min.variab <- variab
                min <- value
              } else {
                param.table[tblidx[[M]][1], i] <- NaN
                next
              }
            }
          }
        }

        stopifnot(min==LSMIN) ## check
        if (min.variab != LSVARIAB && printlevel > PL_SUBIMPORTANT) {
          Print("optimal variables, directly returned and by optim, differ",#
                min.variab, LSVARIAB)
        }
          
        
        fnscale <- if (is.null(fit.fnscale) || is.na(fit.fnscale[M]))
          min else fit.fnscale[M]

        lsq.optim.control <-
          c(opt.control, list(parscale=parscale, fnscale=fnscale))

        #Print(LSVARIAB, lower, upper, autostart, info.cov, parscale); 
        #options(warn=2); Print(LSQLB, LSQUB); print(minmax);
        ## stop("")

        OPTIM(LSVARIAB, LStarget, lower = LSQLB, upper = LSQUB,
              control=lsq.optim.control,optimiser=optimiser,silent=silent)

      } # n.variab > 0
      options(show.error.messages = show.error.message)  
      ## side effect: minimum so far is in LSMIN and LSPARAM
      ## even if the algorithm finally fails


      if (is.finite(LSMIN)) {
        param.table[[M]][tblidx[[M]][1]] <- LSMIN
        param.table[[M]][IDX("variab")] <- LSVARIAB
        param.table[[M]][IDX("param")] <- LSPARAM

        ps <- abs(LSVARIAB)
        zero <- ps == 0
        parscale[!zero] <- ps[!zero]
        
      } else {
        param.table[[M]] <- if (n.variab==0) NA else NaN
      }

      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(LSVARIAB)

      firstoptim <- FALSE
    } # for M   
  } # trans.inv


  if (length(alllsqmeth) > 0) {
    ps <- matrix(NA, nrow=n.variab, ncol=length(alllsqmeth))
    for (iM in 1:length(alllsqmeth)) {
      M <- alllsqmeth[iM]

      if (!is.na(param.table[[M]][tblidx[[M]][1]])) {
        ps[ , iM] <- abs(param.table[[M]][IDX("variab")])
      }
    }
    ps <- apply(ps, 1, median, na.rm=TRUE)
    parscale <- ParScale(optim.control, ps, lower, upper)
  }


##################################################
### optional parameter grid for MLE and CROSS  ###
  if (printlevel>=PL_STRUCTURE) cat("\nmle param...")
  
  
  idx <- IDX("variab")
  glblvar.idx <- IDX("glbl.var")

  
  gridmax <- as.matrix(param.table[idx, cm$lsq])
  if (!any(is.finite(gridmax))) gridmax <- param.table[idx, , drop=FALSE]
   
  gridmin <- apply(gridmax, 1, min, na.rm=TRUE)
  gridmax <- apply(gridmax, 1, max, na.rm=TRUE)

  gridbound <- lower
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  abase <- 0.25
  a <- is.na(gridmin[idx]) * (1-abase) + abase
  ## maybe there have not been any lsq elimination; then a=1
  gridmin[idx] <- (1-a) * gridmin[idx] + a * gridbound[idx]
  stopifnot(all(gridmin >= lower))
  gridbound <- upper
  gridbound[!is.finite(gridbound)] <- NA
  idx <- !is.na(gridbound)
  a <- is.na(gridmax[idx]) * (1-abase) + abase
  gridmax[idx] <- (1-a) * gridmax[idx] + a * gridbound[idx]
  stopifnot(all(gridmax <= upper))

 

##################################################
###################   MLE    #####################


  if (printlevel>=PL_STRUCTURE) cat("\nMLE XXX...")
  mleMethods <- (if (is.null(mle.methods)) NULL else
                 allmlemeth[pmatch(mle.methods, allmlemeth)])

  
  if ("reml" %in% mleMethods && n.covariat == 0)
    mleMethods <- c("ml")# to do, "reml", "rml")

  ## lowerbound_scale_ls_factor <  lowerbound_scale_factor, usually
  ## LS optimisation should not run to a boundary (what often happens
  ## for the scale) since a boundary value is usually a bad initial
  ## value for MLE (heuristic statement). Therefore a small
  ## lowerbound_scale_ls_factor is used for LS optimisation.
  ## For MLE elimination we should include the true value of the scale;
  ## so the bounds must be larger. Here lower[SCALE] is corrected
  ## to be suitable for MLE elimination
  
  if (any(MLELB > MLEUB))
    stop("the users lower and upper bounds are too restricitve")

  ## fnscale <- -1 : maximisation
  for (M in c(allmlemeth)) {
 
    if (!(M %in% mleMethods)) next;
    if (printlevel>=PL_STRUCTURE ) cat("\n", M) else cat(pch)
    param.table[[M]][IDX("variab")] <- default.param    
    
    if (M!="ml" && !anyFixedEffect) { ## also reml nicht nochmal rechnen...
      param.table[[M]] <- param.table[["ml"]]
      ## ML-value is now REML value:
      param.table[[M]][tblidx[[M]][1]] <- param.table[[M]][tblidx[["ml"]][1]]
      next
    }
     
    MLEMAX <- -Inf ## must be before next "if (nMLEINDEX==0)"
    MLEVARIAB <- Inf ## nachfolgende for-schleife setzt MLEVARIAB
    MLEPARAM <- NA
    onborderline <- FALSE
    if (length(MLELB) == 0) { ## n.variab == 0
      MLtarget(NULL)
    } else {
      param.table[[M]][IDX("lower")] <- MLELB
      param.table[[M]][IDX("upper")] <- MLEUB
      options(show.error.messages = show.error.message) ##
      max <- -Inf
      
      for (i in methodprevto$mle) { ## ! -- the parts that change if
        ##                             this part is copied for other methods
        ## should mle be included when M=reml?
        ## same for lsq methods as well: should previous result be included?

        if (!any(is.na(variab <- param.table[IDX("variab"), i]))) {

          value <- MLtarget(variab) ## !
          if (is.finite(value)) {
            param.table[tblidx[[M]][1], i] <- value
            if (value > max) {
              max.variab <- variab
              max <- value
            }
          } else {
            param.table[tblidx[[M]][1], i] <- NaN
            next
          }
        }
      }

   
      fnscale <-
        if (is.null(fit.fnscale) || is.na(fit.fnscale[M]))
          -max(abs(max), 0.1) else fit.fnscale[M]

      mle.optim.control <-
        c(opt.control, list(parscale=parscale, fnscale=fnscale))


      
      stopifnot(length(parscale)==0 || length(parscale) == length(MLEVARIAB))
      
      MLEINF <- FALSE

      if (fit$critical < 2) {
        OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper=MLEUB,
              control=mle.optim.control, optimiser=optimiser, silent=silent)
        
        if (MLEINF) {
          if (printlevel>=PL_STRUCTURE )
            Print("MLEINF", MLEVARIAB, MLEMAX) else cat("#") #
          OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper=MLEUB,
                control=mle.optim.control,optimiser=optimiser,silent=silent)
          if (printlevel>=PL_STRUCTURE )
            Print("MLEINF new", MLEVARIAB, MLEMAX) #
        }
              
        options(show.error.messages = TRUE) ##
        mindistance <-
          pmax(fit$minbounddistance, fit$minboundreldist * abs(MLEVARIAB))
        
        onborderline <- 
          (abs(MLEVARIAB - MLELB) <
           pmax(mindistance,  ## absolute difference
                fit$minboundreldist * abs(MLELB) ## relative difference
                )) |
        (abs(MLEVARIAB - MLEUB) <
         pmax(mindistance, fit$minboundreldist * abs(MLEUB)))
      }
    } # length(MLELB) > 0
  
    if (printlevel>=PL_STRUCTURE )
      Print("mle first round", MLEVARIAB, MLEPARAM, MLEMAX) #
    
    if (!is.finite(MLEMAX)) {
      if (printlevel>=PL_IMPORTANT ) message(M, ": MLtarget I failed.")
      param.table[[M]] <- MLEPARAM <- NaN
      variab <- MLELB ## to call for onborderline
      ml.residuals <- NA
    } else {
      param.table[[M]][tblidx[[M]][1]] <- MLEMAX
      param.table[[M]][IDX("variab")] <- MLEVARIAB

      stopifnot(all(MLEVARIAB >= MLELB & MLEVARIAB <= MLEUB))
      
      param.table[[M]][IDX("param")] <- MLEPARAM
      param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)     
      ml.residuals <- ML.RESIDUALS

      
 

      if (FALSE) Print(recall, fit$critical>=2 || any(onborderline),#
             !fit$only_users , fit$critical >= 0)
      
      if ((fit$critical>=2  || any(onborderline))
          && !fit$only_users && fit$critical >= 0) {
        ## if the MLE result is close to the border, it usually means that
        ## the algorithm has failed, especially because of a bad starting
        ## value (least squares do not always give a good starting point,helas)
        ## so the brutal method:
        ## calculate the MLE values on a grid and start the optimization with
        ## the best grid point. Again, there is the believe that the
        ## least square give at least a hint what a good grid is
     
        if (fit$critical == 0) {
          MLEgridmin <- gridmin
          MLEgridmax <- gridmax
          
          if (any(is.na(MLEgridmin)) || any(is.na(MLEgridmax))) {
            if (printlevel >= PL_SUBIMPORTANT ) {
              Print(cbind(MLELB, variab, MLEUB, onborderline), #
                    MLEgridmin, MLEgridmax)
            }
            warning(paste(M, "converged to a boundary value -- ",
                          "better performance might be obtained",
                          "when allowing for more lsq.methods"))
          } else {
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(1, M, MLEMAX, MLEVARIAB) else cat(detailpch)
            MLEgridlength <-
              max(3, round(fit$approximate_functioncalls^(1/n.variab)))
            ## grid is given by the extremes of the LS results
            ## so, therefore we should examine above at least 4 different sets
            ## of weights wichtig: gridmin/max basiert auf den reduzierten
            ## Variablen 
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-2) # grid starts
                                        # bit outside
            MLEgridmin <- pmax(MLEgridmin - runif(length(MLEgridmin)) * step/2,
                               MLELB)   # the extremes of LS
            MLEgridmax <- pmin(MLEgridmax + runif(length(MLEgridmax)) * step/2,
                               MLEUB)
            step <- (MLEgridmax - MLEgridmin) / (MLEgridlength-1)
            
            startingvalues <- vector("list", length(step))
            for (i in 1:length(step)) {
              startingvalues[[i]] <-
                MLEgridmin[i] + step[i] * 0:(MLEgridlength-1)
            }
            
            startingvalues <- do.call(base::expand.grid, startingvalues)
            
            limit <- 10 * fit$approximate_functioncalls
            if ((rn <- nrow(startingvalues)) > limit) {              
              if (printlevel>=PL_STRUCTURE)
                cat("using only a random subset of the", rn, "grid points")
              rand <- runif(rn)
              startingvalues <-
                startingvalues[rand < quantile(rand, limit / rn), , drop=FALSE]
              gc()
            }
            
            MLEMAX <- -Inf
            
            apply(startingvalues, 1,
                  function(x) {
                  try(MLtarget(x), silent=silent)
                  ## try-result:
                  ##     if (!is.numeric(z) && abs(z)<1e10) cat(paste(x, sep=","), "\n") else cat(".\n")
                  ##  if (MLEINF) stop("stop")
                })
            
            if (printlevel>=PL_FCTN_SUBDETAILS)
              Print("mle grid search", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            ## side effect:Maximum is in MLEMAX!
            ##                             and optimal parameter is in MLEVARIAB
            if (printlevel>=PL_FCTN_SUBDETAILS)
              show(2, M, MLEMAX, MLEVARIAB)
            
            cat(detailpch)
            options(show.error.messages = show.error.message) ##

            OPTIM(MLEVARIAB, MLtarget, lower = MLELB, upper = MLEUB,
                  control=mle.optim.control, optimiser=optimiser,
                  silent=silent)
            options(show.error.messages = TRUE) ##
            if (!is.finite(MLEMAX) &&(printlevel>=PL_IMPORTANT))
              message("MLtarget II failed.\n")
            ## do not check anymore whether there had been convergence or not.
            ## just take the best of the two strategies (initial value given by
            ## LS, initial value given by a grid), and be happy.
            if (printlevel>=PL_FCTN_SUBDETAILS )
              show(3, M, MLEMAX, MLEVARIAB) else 
            if (printlevel>=PL_STRUCTURE )
              Print("mle second round", MLEVARIAB, MLEPARAM, MLEMAX) #
            
            if (is.finite(MLEMAX) && MLEMAX >=param.table[[M]][tblidx[[M]][1]]){
              param.table[[M]][tblidx[[M]][1]] <- MLEMAX
              param.table[[M]][IDX("variab")] <- MLEVARIAB
              stopifnot(all(MLEVARIAB >= MLELB & MLEVARIAB <= MLEUB))
              param.table[[M]][IDX("param")] <- MLEPARAM
              param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)
              ml.residuals <- ML.RESIDUALS
            }
          } # (is.na(MLEgridmin[1]))
        } else { # fit$critical > 0  

          critical <- ptype == CRITICALPARAM;
          if (fit$critical>=3 || !any(critical)) critical <- rep(TRUE, n.variab)
          ncrit <- as.integer(fit$n_crit^(1/sum(critical)))
          if (!is.finite(ncrit)) ncrit <- 2
          if (ncrit^sum(critical) > 100 && printlevel>=PL_IMPORTANT)
            message("The optimisation may last a pretty long time!")
          
          if (ncrit > 1 || fit$critical>=2) {
            if (!is.null(transform)) {
              stop("if 'transform' is given, 'critical' must be '0'")
            }
            w <- 0.5
            lowlist <- upplist <- list()
            for (i in 1:n.variab) {
              if (critical[i]) {               
                newparam <- seq(if (SDVAR.IDX[i]) 0 else MLELB[i], MLEUB[i],
                                len=ncrit+1)
                lowlist[[i]] <- newparam[-length(newparam)]
                upplist[[i]] <- newparam[-1]
              } else {
                lowlist[[i]] <- MLELB[i]
                upplist[[i]] <- MLEUB[i]
              }
            }
            lowlist <- as.matrix(do.call(base::expand.grid, lowlist))
            upplist <- as.matrix(do.call(base::expand.grid, upplist))
            
            orig.MLEVARIAB <- MLEVARIAB
            b.idx <- is.finite(bounds)

            stopifnot(length(lowlist) > 1)
          
            for (i in 1:nrow(lowlist)) {
              cat(detailpch)

              new.lower.vector <- trafo(lowlist[i, ])
              new.upper.vector <- trafo(upplist[i, ])
              new.bounds <- TRUE
                        
              new.parscale <- guess <- fnscale <- NULL   
              first.passage <- TRUE
              while (TRUE) {
                if (new.bounds) {
                  new.bounds <- FALSE
                  .C("PutValuesAtNA", LiliReg, as.double(new.lower.vector),
                     PACKAGE="RandomFields")
                  new.lower <- GetModel(register=LiliReg,modus=GETMODEL_DEL_MLE,
                                        which.submodels = "user.but.once+jump",
                                        spConform=FALSE, do.notreturnparam=TRUE)

                  .C("PutValuesAtNA", LiliReg, as.double(new.upper.vector),
                     PACKAGE="RandomFields")
                  new.upper <- GetModel(register=LiliReg,modus=GETMODEL_DEL_MLE,
                                        which.submodels="user.but.once+jump",
                                        spConform=FALSE, do.notreturnparam=TRUE)
                }
               
                if (printlevel>=PL_STRUCTURE) cat("\nrecall rffit...\n")

                res <-
                  rffit.gauss(Z= Z,
                              lower=new.lower,
                              upper=new.upper,
                              mle.methods="ml",
                              lsq.methods=lsq.methods,
                              users.guess = guess,
                              optim.control=c(opt.control,
                                fnscale=list(fnscale),
                                parscale=list(new.parscale)),
                              transform=transform,
                              recall = TRUE,
                              general.pch = if (pch == "") "" else ":",
                              general.practicalrange = general$practicalrange,
                              general.spConform = FALSE,
                              fit.critical = -1,
                              fit.split =FALSE)
                
 
                guess <- res$table[[M]][IDX("param")]
                ## scale_ratio <- abs(log(abs(guess[-delete.idx]/new.parscale)))

                stopifnot(length(users.lower) + length(delete.idx)
                          == length(new.lower.vector))
                
                if (!all(guess >= new.lower.vector & guess <= new.upper.vector)
                    || !all(new.lower.vector[-delete.idx] > users.lower &
                            new.upper.vector[-delete.idx] < users.upper)) {
                  new.lower.vector <- pmin(new.lower.vector, guess)
                  new.upper.vector <- pmax(new.upper.vector, guess)
                  new.lower.vector[-delete.idx] <-
                    pmax(new.lower.vector[-delete.idx], users.lower)
                  new.upper.vector[-delete.idx] <-
                    pmin(new.upper.vector[-delete.idx], users.upper)
                  guess <- pmax(new.lower.vector, pmin(new.upper.vector, guess))
                  new.bounds <- TRUE
                }

                
                likelihood <- res$table[[M]][tblidx[[M]][1]]

                if (is.finite(likelihood)) {
                  if (likelihood > param.table[[M]][tblidx[[M]][1]]) {
                    if (printlevel > PL_RECURSIVE && !first.passage)
                      cat("parscale: table improved by ",
                          likelihood -  param.table[[M]][tblidx[[M]][1]], "\n") 
                    param.table[[M]][tblidx[[M]][1]] <- MLEMAX <- likelihood
                    for (j in c("variab", "param", "glbl.var"))
                      param.table[[M]][IDX(j)] <- res$table[[M]][IDX(j)]
                    ml.residuals <- res$ml$residuals
                  }
                }
                
                if (first.passage) {
                  old.likelihood <- likelihood
                } else {
                  if (printlevel > PL_RECURSIVE) {
                    if (likelihood > old.likelihood) {                      
                      cat("parscale: mle improved by", likelihood-old.likelihood,
                          "\n")
                    }
                  }
                  break;
                }
                            
                ## urspruengliches Abbruchkriterium, das nicht gut fkt:
                ##value.ratio <- abs(log (abs(likelihood) / abs(MLEMAX)))
                ##outside <- (scale_ratio > fit$scale_ratio)  & MLEVARIAB != 0
                ##outside <- outside & (!b.idx | new.parscale >
                ##                      exp(fit$scale_ratio) * 
                ##                      (w * abs(guess) + (1-w) * bounds))
                ## if (!any(outside) && value.ratio <= fit$scale_ratio) {
               ##   break
               ## }

                
                zero <- new.parscale == 0
                new.parscale[zero] <-
                  pmax(abs(lower[zero]), abs(upper[zero])) / 10                
                new.parscale[b.idx] <-
                  w * new.parscale[b.idx] + (1-w) * bounds[b.idx]
                #new.parscale <- abs(guess[-delete.idx])

                restable <- as.matrix(res$table)
                used <- which(!apply(is.na(restable), 2, all)[-1:-2]) # ohne auto,user
                names.tbl <- names(tblidx)
                start <- which(names.tbl %in% alllsqmeth)              
                start <- if (length(start) == 0) 0 else min(start) - 1
                if (start>0) names.tbl <- names.tbl[-1:-start]
                 
                fnscale <- numeric(length(used))                
                for (j in 1:length(used)) {
                  fnscale[j] <- abs(restable[tblidx[[start + used[j]]][1],
                                              2 + used[j]])
                  if (names.tbl[used[j]] == "ml")
                    fnscale[j] <- -max(0.1, fnscale[j])
                }

                .C("PutValuesAtNA", LiliReg, guess,PACKAGE="RandomFields",
                   NAOK=TRUE) # ok
                guess <- GetModel(register=LiliReg, modus=GETMODEL_DEL_MLE,
                                  which.submodels = "user.but.once+jump",
                                  spConform=FALSE, do.notreturnparam=TRUE)
                first.passage <- FALSE
              } # while true
            } # for 1:nrow(lowlist)            
          } # ncrit  > 1          
        } # fit$critical > 0
      }
   
       if (!fit$only_users && (fit$reoptimise || any(SDVAR.IDX)) &&
          fit$critical >= 0 && fit$critical <= 1){
        if (pch!="") cat("$")
        new.parscale <- abs(revariab <- param.table[[M]][IDX("variab")])
        new.parscale[new.parscale == 0] <- 1e-3
        fnscale <- -max(0.1, abs(MLEMAX <- param.table[[M]][tblidx[[M]][1]]))
      
        old.control<- if (exists("mle.optim.control")) mle.optim.control else NA
        mle.optim.control <- c(opt.control,
                               list(parscale=new.parscale, fnscale=fnscale))
        low <- MLELB
        
        if (any(SDVAR.IDX)) low[SDVAR.IDX] <- pmax(0, users.lower[SDVAR.IDX])


        OPTIM(revariab, MLtarget, lower = low, upper=MLEUB,
              control=mle.optim.control, optimiser=optimiser, silent=silent)

        if (MLEMAX >= param.table[[M]][tblidx[[M]][1]]) {   
          if (printlevel > PL_SUBIMPORTANT)
            Print(old.control, fnscale, revariab, #
                  param.table[[M]][IDX("param")],
                  MLEMAX, MLEVARIAB,  MLEPARAM, MLELB, MLEUB)
          param.table[[M]][tblidx[[M]][1]] <- MLEMAX
          param.table[[M]][IDX("variab")] <- MLEVARIAB         
          param.table[[M]][IDX("param")] <- MLEPARAM
          param.table[[M]][IDX("glbl.var")] <- get.var.covariat(MLEVARIAB)
        }
      }
    } # is.finite(MLEMAX)     
  } ## 



######################################################################
###                     error checks                               ###
######################################################################
  ## if rather close to nugget and nugget==fixed, do exception handling.
  ## This part is active only if
  ## scale_max_relative.factor < lowerbound_scale_ls_factor
  ## By default it is not, i.e., the following if-condition
  ## will "always" be FALSE.
  
  if (globalvariance && length(nugget.idx) == 0 && any(SCALE.IDX)) {
    idx <- IDX("variab")
    alllsqscales <- param.table[idx, cm$lsq][SCALE.IDX, ]
   
    if (any(alllsqscales < mindistances/fit$scale_max_relative_factor,
            na.rm=TRUE))
      warning(paste(sep="",
                    "Chosen model seems to be inappropriate!\n Probably a (larger) nugget effect should be considered")
              )
  }


######################################################################
###                   prepare lower and upper                      ###
######################################################################

  ## if the covariance functions use natural scaling, just
  ## correct the final output by GNS$natscale
  ## (GNS$natscale==1 if no rescaling was used)
  ##

  lower <- trafo(lower)
  upper <- trafo(upper)
  idx <- lower == upper
  lower[idx] = upper[idx] = NA
  .C("PutValuesAtNA", LiliReg, as.double(lower), PACKAGE="RandomFields",
     NAOK=TRUE) # ok
  lower <- GetModel(register=LiliReg, modus=GETMODEL_SOLVE_MLE,
                    which.submodels = "user.but.once+jump")
  .C("PutValuesAtNA", LiliReg, as.double(upper), PACKAGE="RandomFields",
     NAOK=TRUE) # ok
  upper <- GetModel(register=LiliReg, modus=GETMODEL_SOLVE_MLE,
                    which.submodels = "user.but.once+jump")
 

######################################################################
###                   prepare models for returning                 ###
######################################################################

  if (printlevel>=PL_STRUCTURE) cat("preparing for returning ...\n");

  idxCovar <- IDX("covariat")
  idxpar <- IDX("param") # henceforth
  if (globalvariance) idxvar <- IDX("glbl.var")
  res <- values.res <- list()
  nparam <- as.integer(n.variab + n.covariat) ## not n.param
  AICconst <-
    2 * nparam + 2 * nparam * (nparam + 1) / (sum.not.isna.data - nparam - 1) 

  allMethods <- c(primMethods, lsqMethods, mleMethods)
  for (i in OneTo(length(allmethods))) {
    if (!(allmethods[i] %in% allMethods)) next
    
    if (is.na(param.table[1, i]) ## && !is.nan(param.table[1,i]) ?? unklar
        && (length(idxpar)!=1 || idxpar != 0)) next ## result with NAs
    
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##     calculate all target values for all optimal parameters     +++
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    if (printlevel>= PL_STRUCTURE) cat("calculating ...\n")
    p <- param.table[IDX("variab"), i]      
    for (M in alllsqmeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (M %in% lsqMethods) {          
        LSQsettings(M)
        param.table[tblidx[[M]][1], i] <- LStarget(p)
      }
    } 
    
    for (M in allmlemeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (is.na(cur[1]) && !is.nan(cur[1]) && M %in% mleMethods) {
        param.table[tblidx[[M]][1], i] <- MLtarget(p)
      }

      if (allmethods[i] == M) {
        param.table[[M]][IDX("sdvariab")] <-
          INVDIAGHESS(param.table[[M]][IDX("variab")], MLtarget,
                      control=mle.optim.control)$sd
      }
    }


    for (M in allcrossmeth) {
      cur <- param.table[tblidx[[M]][1], i]
      if (is.na(cur) && !is.nan(cur) && M %in% NULL) { # crossMethods) {
        stop("not programmed ")
          ##  crosssettings(M) ## uncomment
        variab <- p
        if (n.covariat > 0) {
          variab <- c(variab, param.table[idxCovar, i])
        }
        param.table[tblidx[[M]][1], i] <- stop("") # crosstarget(variab)
      }
    }
  
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## die nachfolgenden Zeilen aendern die Tabellenwerte -- somit
    ## muessen diese Zeilen unmittelbar vor dem return stehen !!!
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (fit$use_naturalscaling && any(SCALE.IDX)) {
      param <- as.double(param.table[idxpar, i]  + 0.0) ## + 0.0 MUSS STEHEN
      ## register muss mit cov initialisiert sein, sonst
      ## wird laueft er in cov->key rein!
      .C("PutValuesAtNA", LiliReg, param, PACKAGE="RandomFields")
      .C("expliciteDollarMLE", LiliReg, param, PACKAGE="RandomFields")
      param.table[idxpar, i]  <- param
    }

    #Print(fit$use_naturalscaling, any(SCALE.IDX), globalvariance,
     #     as.double(param.table[idxpar, i] )); 


    

    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##  models with NAs filled
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (printlevel>= PL_STRUCTURE)cat("with nas filled ...\n")

    .C("PutValuesAtNA", LiliReg, as.double(param.table[idxpar, i] ),
       PACKAGE="RandomFields")
    if (globalvariance) .C("PutGlblVar", as.integer(LiliReg),
                           as.double(param.table[glblvar.idx, i]))
    modelres <- GetModel(register = LiliReg, modus = GETMODEL_SOLVE_MLE,
                        spConform = if (is(Z$orig.model, "formula")) 2
                         else general$spConform,
                         solve_random = TRUE,
                         which.submodels = "user.but.once+jump"
                         )
    if (is(Z$orig.model, "formula")) {
      VarNames <- extractVarNames(Z$orig.model)      
      VarNames <-
        if (length(VarNames) == 0) ""
        else paste("c(", paste(VarNames, collapse=","), ")")
      txt <- paste(VarNames, "~", model2string(modelres))
      formel <- as.formula(txt)
    } else formel <- NULL
    
    #    Print(allmethods[i], modelres, i, allmethods); 
    
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ##   AIC etc for MLE
    ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (printlevel>= PL_STRUCTURE) cat("AIC...\n")
    M <- allmethods[i]
 #   if (M == "ml")  {
      ## jetzt muss RFlikelihood-initialisierung sein!

    lilihood <- .Call("EvaluateModel", double(0), LiliReg,
                      PACKAGE="RandomFields")
    residu <- get.residuals(LiliReg)
    likelihood <- param.table[tblidx[["ml"]][1], i]
    
    if (abs(lilihood[1] - likelihood) > 1e-7 * (abs(likelihood) + 1)) {
      #stop("The two ways of calculating the likelihood do not match: ",
     #      lilihood[1], " != ", likelihood)
    }
    AIC <- 2 * nparam  - 2 * likelihood
    AICc <- AICconst - 2 * likelihood
    BIC <- log(sum.not.isna.data) * nparam - 2 * likelihood      
    param.table[tblidx[["AIC"]][1], i] <- AIC
    param.table[tblidx[["AICc"]][1], i] <- AICc
    param.table[tblidx[["BIC"]][1], i] <- BIC
  #  }

    
    if (n.covariat > 0) {
      c.table <- rbind(param.table[IDX("covariat"), i], NA)
      dimnames(c.table) <- list(NULL, betanames) 
    } else c.table <- NULL
    v.table <- rbind(param.table[IDX("variab"), i],
                     param.table[IDX("sdvariab"), i],
                     param.table[IDX("lower"), i],
                     param.table[IDX("upper"), i]
                     )

     
    dimnames(v.table) <- list(c("value", "sd", "lower", "upper"), varnames)
    p.table <- if (n.variab > 0) rbind(param.table[IDX("param"), i], NA)
               else matrix(ncol=0, nrow=2)
    dimnames(p.table) <- list(c("value", "sd"), paramnames)
    
    res[[M]] <-
      list(model=modelres,
           formel = formel,
           variab = v.table,
           param = p.table,
           covariat = c.table,
           globalvariance=if (globalvariance) param.table[IDX("glbl.var"),i][1],
           hessian = NULL,
           likelihood = if (TRUE || M == "ml") likelihood else as.integer(NA),
           AIC = if (TRUE || M == "ml")  AIC else as.integer(NA),
           AICc = if (TRUE || M == "ml")  AICc else as.integer(NA),
           BIC = if (TRUE || M == "ml")  BIC else as.integer(NA),
           residuals = if (TRUE || M == "ml") residu else as.integer(NA)
           )
    
    class(res[[M]]) <- "RM_modelFit"
  }
  
  if (pch!="" && !recall) cat("\n")
                     
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ## Hess-Matrix for the parameters (for the variables see above)
  ## and further error checking
  ##
  ## next lines must be the very last call as standards are overwritten
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   if (printlevel>= PL_STRUCTURE) cat("hessian ...\n")

  ## prepare MLtarget for getting 'par3am' directly
  trafo <- function(x) x
  MLELB <- param.lower
  MLEUB <- param.upper

  for (i in OneTo(length(mleMethods))) {
    p <- param.table[IDX("param"), i] ## variab ist oben!
    if (any(p < MLELB | p > MLEUB)) 
      stop(paste("estimated parameters outside bounds --",
                 if (is.null(transform)) "please inform author"
                 else "'transform' is likely not correct"))
    
    if (!is.na(p[1])) {        
      M <- mleMethods[i]     
      cur <- param.table[tblidx[["ml"]][1], i]
      if (!is.na(cur)) { ## kann bei autostart + Bayes auftreten
        if (abs(cur - MLtarget(p)) > 1e-10 && printlevel >= PL_IMPORTANT)
          message("likelihoods are calculated in two different ways leading",
                  " to the values ", cur, " and ",  MLtarget(p),
                  ", which are not the same. Please contact author.\n")
        H <- INVDIAGHESS(p, MLtarget)
        res[[M]]$param[2, 1:n.param] <- inv <- H$sd
        res[[M]]$hessian <- H$hessian
        param.table[IDX("sdparam"), i] <- inv
      }
    }
  }

  
  if (printlevel >= PL_IMPORTANT && BEYOND > 100)
    cat("\nNote: There are ",
        if (BEYOND > 1000)"very strong" else
        if (BEYOND > 250) "strong" else "some",
        " indications that the model might be overparametrised\nor that the bounds for the variables are too wide. Try narrower lower and\nupper bounds for the variables in the latter case. One of the critical\nparameters is 'lowerbound_var_factor' whose value (", fit$lowerbound_var_factor, ") might be reduced.\n", sep="")
  
  if (printlevel>= PL_STRUCTURE) cat("S3 / S4 ...\n")


#  print(param.table)
  
  L <- list(ev = ev,
            table=param.table,
            n.variab = as.integer(n.variab),
            n.param = as.integer(n.param),
            n.covariates = as.integer(n.covariat),
            deleted = as.integer(delete.idx),
            lowerbounds=lower,
            upperbounds=upper,
            transform = transform,
            
            number.of.data= as.integer(sum.not.isna.data),
            number.of.parameters = nparam,
            modelinfo = nice_modelinfo(minmax),               
            p.proj = integer(0),
            v.proj = integer(0),
            x.proj = TRUE,
            fixed = NULL,
            true.tsdim = as.integer(tsdim),
            true.vdim = as.integer(vdim),
            report = "",
            submodels = NULL)


#  Print(res)
  
  if (!general$spConform) {
    Res <- c(L,
             list(vario = "'$vario' is defunctioned. Use '$ml' instead of '$value$ml'!" ## must be the very last of the parameter list!
                  ),
             res)
    class(Res) <-  "RF_fit"
    return(Res)
  } else {

    #str(res)
    
     res2 <- lapply(res, FUN=list2RMmodelFit, isRFsp=Z$isRFsp,
                   coord=Z$RFsp.coord, gridTopology=Z$gridTopology,
                   data.RFparams=Z$data.RFparams, coord[[1]]$T)
     return(rffit.gauss2sp(res2, L, Z))
  }
}

rffit.gauss2sp <- function(res2, L, Z) {
  #Print(res2, L, Z, L$ev, Z$varunits, res2) ; 
  if (length(L$ev) > 0) L$ev <- list2RFempVariog(L$ev)
  L$lowerbounds <- list2RMmodel(L$lower)
  L$upperbounds <- list2RMmodel(L$upper)
  if (is.null(L$transform)) L$transform <- list()
  do.call.par <- c(list(Class = "RFfit",
                        Z=Z,
                        coordunits=Z$coordunits,
                        varunits= if (length(L$ev) > 0) L$ev@varunits
                        else ""
                       ),
                   L,
                   res2)

#  str(do.call.par, max.level=2)
  
  z <- do.call(methods::new, args=do.call.par)
  z
}

  
