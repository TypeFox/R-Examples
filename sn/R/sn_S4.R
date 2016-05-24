#  file sn/R/sn_S4.R (S4 methods and classes)
#  This file is a component of the package 'sn' for R
#  copyright (C) 1997-2014 Adelchi Azzalini
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#---------
setClass("SECdistrUv",
   representation(family="character", dp="numeric", name="character"),
   validity=function(object){
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     np <- 3 + as.numeric(object@family %in% c("ST","ESN"))
     if(length(object@dp) != np) return(FALSE) 
     if(object@dp[2] <= 0) return(FALSE)
     TRUE
   }
)

setClass("summary.SECdistrUv",
   representation(family="character", dp="numeric", name="character",
     cp="numeric", cp.type="character", aux="list"),
   validity=function(object){
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     np <- 3 + as.numeric(object@family %in% c("ST","ESN"))
     if(length(object@dp) != np) return(FALSE) 
     if(object@dp[2] <= 0) return(FALSE)
     # if(length(object@op) != length(object@dp)) return(FALSE)
     if(length(object@cp) != length(object@dp)) return(FALSE)
     TRUE
   }
)

setClass("SECdistrMv",
   representation(family="character", dp="list", name="character", 
                  compNames="character"),
   validity=function(object){
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     np <- 3 + as.numeric(object@family %in% c("ST","ESN"))
     dp <- object@dp
     if(mode(unlist(dp)) != "numeric") return(FALSE)
     if(length(dp) != np) return(FALSE) 
     d <- length(dp[[3]])
     Omega <- dp[[2]]
     if(length(dp[[1]]) != d | any(dim(Omega) != c(d,d))) return(FALSE)
     if(any(Omega != t(Omega))) {message("non-symmetric Omega"); return(FALSE)}
     if(any(eigen(Omega, symmetric=TRUE, only.values=TRUE)$values <= 0)) { 
        message("Omega not positive-definite")
        return(FALSE)}
     if(object@family == "ST") { if(dp[[4]] <= 0) return(FALSE) }
     if(length(object@compNames) != d) return(FALSE)
     if(length(object@name) != 1) return(FALSE)
     TRUE
   }
)

setClass("summary.SECdistrMv",
   representation(family="character", dp="list", name="character", 
     compNames="character",     # op="list", 
     cp="list", cp.type="character", aux="list"),
   validity=function(object){
     family <- object@family
     if(!(family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     np <- 3 + as.numeric(family %in% c("ST","ESN"))
     dp <- object@dp
     if(mode(unlist(dp)) != "numeric") return(FALSE)
     if(length(dp) != np) return(FALSE) 
     d <- length(dp[[3]])
     if(length(dp[[1]]) != d | any(dim(dp[[2]]) != c(d,d))) return(FALSE)
     if(family == "ST") { if(dp[[4]] <= 0) return(FALSE) }
     if(length(object@compNames) != d) return(FALSE)
     if(length(object@name) != 1) return(FALSE)
     if(length(object@cp) != length(object@dp)) return(FALSE)
     # if(length(object@op) != length(object@dp)) return(FALSE)
     TRUE
   }
)


setMethod("show", "SECdistrUv",
  function(object){
    if(object@name != "") 
      cat("Probability distribution of variable '", object@name, "'\n", sep="")
    cat("Skew-elliptically contoured distribution of univariate family", 
      object@family,"\nDirect parameters:\n")
    print(object@dp)
  }
)

setMethod("show","SECdistrMv",
  function(object){
    if(object@name != "")  
      cat("Probability distribution of variable '", object@name, "'\n", sep="")
    dp <- object@dp
    attr(dp[[2]],"dimnames") <- 
         list(paste("Omega[", object@compNames, ",]", sep=""), NULL)
    cat("Skew-elliptically contoured distribution of ", length(dp[[3]]),
        "-dimensional family ", object@family,"\nDirect parameters:\n", sep="")
    out <- rbind(xi=dp[[1]], Omega=dp[[2]], alpha=dp[[3]])
    colnames(out) <- object@compNames
    print(out)
    if(object@family=="ST") cat("nu", "=", dp[[4]], "\n")
    if(object@family=="ESN") cat("tau", "=", dp[[4]], "\n")
  }
)
#

#--------------------

setMethod("show", "summary.SECdistrUv",
  function(object){
    obj <- object
    if(obj@name != "")  
      cat("Probability distribution of variable '", obj@name, "'\n", sep="")
    cat("\nSkew-elliptical distribution of univariate family", obj@family,"\n")
    cat("\nDirect parameters (DP):\n")
    print(c("", format(obj@dp)), quote=FALSE)
    # cat("\nOriginal parameters (OP):\n")
    # print(c("", format(obj@op)), quote=FALSE)
    cp <- obj@cp
    note <- if(obj@cp.type == "proper") NULL else ", type=pseudo-CP" 
    cat(paste("\nCentred parameters (CP)", note, ":\n", sep=""))
    print(c("", format(cp)), quote=FALSE)
    cat("\nAuxiliary quantities:\n")
    print(c("", format(c(delta=obj@aux$delta, mode=obj@aux$mode))), quote=FALSE)
    cat("\nQuantiles:\n")
    q <- obj@aux$quantiles
    q0 <- c("q", format(q))
    names(q0) <- c("p", names(q))
    print(q0,  quote=FALSE)
    measures <- rbind(obj@aux$std.cum, obj@aux$q.measures)
    cat("\nMeasures of skewness and kurtosis:\n ")
    attr(measures, "dimnames") <- list(
       c("  std cumulants", "  quantile-based"), c("skewness", "kurtosis"))
    print(measures)
    }
)

setMethod("show","summary.SECdistrMv",
  function(object){
    obj <- object
    #------ DP
    dp <- obj@dp
    if(obj@name != "") cat("Probability distribution of",obj@name,"\n")
    cat("Skew-elliptically contoured distribution of ", length(dp[[3]]),
        "-dimensional family ", obj@family,"\n", sep="")
    cat("\nDirect parameters (DP):\n")
    attr(dp[[2]], "dimnames") <- 
         list(paste("  Omega[", obj@compNames, ",]", sep=""),NULL)
    out.dp <- rbind("  xi"=dp[[1]], omega=dp[[2]],"  alpha"=dp[[3]])
    colnames(out.dp) <- obj@compNames
    print(out.dp)
    if(length(dp) > 3){
      extra <- unlist(dp[-(1:3)])
      names(extra) <- paste("  ",names(dp[-(1:3)]), sep="")
      # print(extra)
      for(j in 1:length(extra)) cat(names(extra)[j], "=", extra[j], "\n")
      }
    #------ OP
    if(FALSE) {
    op <- obj@op  
    cat("\nOriginal parameters (OP):\n")
    attr(op[[2]], "dimnames") <- 
         list(paste("  Psi[", obj@compNames, ",]", sep=""),NULL)
    out.op <- rbind("  xi"=op[[1]], "  psi"=op[[2]],"  lambda"=op[[3]])
    colnames(out.op) <- obj@compNames
    print(out.op)
    if(length(op) > 3){
      extra <- unlist(op[-(1:3)])
      names(extra) <- paste("  ",names(op[-(1:3)]), sep="")
      # print(extra)
      for(j in 1:length(extra)) cat(names(extra)[j], "=", extra[j], "\n")
      }  
    }
    #------ CP  
    cp <- obj@cp
    note <- if(obj@cp.type == "proper") NULL else ", type=pseudo-CP" 
    cat("\nCentred parameters (CP)", note, ":\n", sep="")
    attr(cp[[2]], "dimnames") <- 
       list(paste("  var.cov[", obj@compNames, ",]", sep=""),NULL)
    out.cp <- rbind("  mean"=cp[[1]], cp[[2]], "  gamma1"=cp[[3]])
    colnames(out.cp) <- obj@compNames
    print(out.cp)
    if(length(cp) > 3) {
      extra <- unlist(cp[-(1:3)])
      names(extra) <- paste("  ", names(cp[-(1:3)]), sep="")
      for(j in 1:length(extra)) cat(names(extra)[j], "=", extra[j], "\n")
      }
    aux <- obj@aux
    out.aux <- rbind("  delta" = aux$delta, "  mode" = aux$mode) 
        #"  lambda"=aux$lambda, 
    colnames(out.aux) <- obj@compNames
    cat("\nAuxiliary quantities:\n")
    print(out.aux)
    cat("\nGlobal quantities:\n")
    cat("  alpha* =", format(aux$alpha.star), 
        ", delta* =", format(aux$delta.star), "\n")
    mardia <- obj@aux$mardia
    cat("  Mardia's measures: gamma1M = ", format(mardia[1]),
        ", gamma2M = ", format(mardia[2]),"\n", sep="")
    invisible()
    }
)
 
setClass("selm",
   representation(call="call", family="character",  logL="numeric", 
     method="character",
     param="list", param.var="list",  size="vector", fixed.param="vector",
     residuals.dp="numeric", fitted.values.dp="numeric", control="list", 
     input="list", opt.method="list"),
   validity=function(object){
     if(class(object) != "selm") return(FALSE)
     if(!is.numeric(object@logL)) return(FALSE)
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     if(!is.vector(object@param$dp)) return(FALSE)
     TRUE
   }
)

setMethod("coef", "selm", coef.selm)

setMethod("logLik", "selm", 
   function(object){
     logL <- slot(object,"logL")
     attr(logL, "df") <- as.numeric(slot(object, "size")["n.param"])    
     class(logL) <- "logLik"
     return(logL)
     }
   )
   
   
setMethod("vcov", "selm", function(object, param.type="CP") {
    vcov <- slot(object, "param.var")[[tolower(param.type)]]
    if(is.null(vcov) & tolower(param.type) == "cp") {
        message("CP not defined, consider param.type='DP' or 'pseudo-CP'")
        return(NULL)}
    vcov}
   )
  
setMethod("show", "selm",
  function(object){
    # cat("Object: ", deparse(substitute(obj)),"\n")
    cat("Object class:", class(object), "\n") 
    cat("Call: ")
    print(object@call)
    cat("Number of observations:", object@size["n.obs"], "\n")
    if(!is.null(slot(object,"input")$weights))
      cat("Weighted number of observations:", object@size["nw.obs"], "\n")
    cat("Number of covariates:", object@size["p"], "(includes constant term)\n")
    cat("Number of parameters:", object@size["n.param"], "\n")
    cat("Family:", slot(object,"family"),"\n")
    fixed <- slot(object, "param")$fixed
    if(length(fixed) > 0) { 
      fixed.char <- paste(names(fixed), format(fixed), sep=" = ", collapse=", ")
      cat("Fixed parameters:", fixed.char, "\n") }  
    method <- slot(object, "method") 
    u <- if(length(method)==1) NULL else 
         paste(", penalty function:", method[2])
    cat("Estimation method: ", method[1], u, "\n", sep="")
    logL.name <- paste(if(method[1]=="MLE") "Log" else "Penalized log",
        "likelihood:", sep="-")
    cat(logL.name, format(object@logL, nsmall=2),"\n")
    if(object@param$boundary) 
      cat("Estimates on/near the boundary of the parameter space\n")
    invisible(object)  
  }
)


#----------------------------------------------------------

setClass("summary.selm",
   representation(call="call", family="character", logL="numeric",
     method="character",
     param.type="character",   param.table="matrix", param.fixed="list",
     resid="numeric",  control="list", aux="list",
     size="vector", boundary="logical", note="character"),
   validity=function(object){
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     TRUE
   }
)

#----------------------------------------------------------

setClass("mselm",
   representation(call="call", family="character", logL="numeric", 
     method="character", param="list", param.var="list", size="vector",
     residuals.dp="matrix", fitted.values.dp="matrix", control="list",
     input="list", opt.method="list"),
   validity=function(object){
     if(class(object) != "mselm") return(FALSE)
     if(!is.numeric(object@logL)) return(FALSE)
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     if(!is.list(object@param$dp)) return(FALSE)
     TRUE
   }
)

setMethod("coef", "mselm", coef.mselm)  

setMethod("logLik", "mselm",
   function(object){
     logL <- slot(object,"logL")
     attr(logL, "df") <- as.numeric(slot(object, "size")["n.param"])
     class(logL) <- "logLik"
     return(logL)
     }
   )   

setMethod("vcov", "mselm", function(object, param.type="CP") {
    vcov <- slot(object,"param.var")[[tolower(param.type)]]
    if(is.null(vcov) & tolower(param.type) == "cp") {
        message("CP not defined, consider param.type='DP' or 'pseudo-CP'")
        return(NULL)}
    vcov}
   )

setMethod("show", "mselm",
  function(object){
   cat("Object class:", class(object), "\n") 
    cat("Call: ")
    print(object@call)
    cat("Number of observations:", object@size["n.obs"], "\n")
    if(!is.null(slot(object,"input")$weights))
      cat("Weighted number of observations:", object@size["nw.obs"], "\n")
    cat("Dimension of the response:", object@size["d"], "\n")
    cat("Number of covariates:", object@size["p"], "(includes constant term)\n")
    cat("Number of parameters:", object@size["n.param"], "\n")
    cat("Family:", slot(object, "family"),"\n")
    fixed <- slot(object,"param")$fixed
    if(length(fixed) > 0) { 
      fixed.char <- paste(names(fixed), format(fixed), sep=" = ", collapse=", ")
      cat("Fixed parameters:", fixed.char, "\n") }
    method <- slot(object, "method")
    u <- if(length(method) == 1) NULL else 
            paste(", penalty function:", method[2])   
    cat("Estimation method: ", method[1], u, "\n", sep="")
    logL.name <- paste(if(method[1]=="MLE") "Log" else "Penalized log",
        "likelihood:", sep="-")
    cat(logL.name, format(object@logL, nsmall=2),"\n")
    if(object@param$boundary) 
      cat("Estimates on/near the boundary of the parameter space\n")
    invisible(object)
  }
)

#----------------------------------
setClass("summary.mselm",
   representation(call="call", family="character", logL="numeric",
     method="character",
     param.type="character",  param.fixed="list",  resid="matrix",
     coef.tables="list", scatter="list", slant="list", tail="list",
     control="list", aux="list", size="vector", boundary="logical"),
   validity=function(object) {
     if(!(object@family %in% c("SN","ST","SC","ESN"))) return(FALSE)
     TRUE
   }
)

setMethod("mean", signature(x="SECdistrUv"), mean.SECdistrUv)
setMethod("mean", signature(x="SECdistrMv"), mean.SECdistrMv)
setMethod("sd", signature(x="SECdistrUv"), sd.SECdistrUv)
setMethod("vcov", signature(object="SECdistrMv"), vcov.SECdistrMv)

setMethod("plot", signature(x="SECdistrUv", y="missing"), plot.SECdistrUv)
setMethod("plot", signature(x="SECdistrMv", y="missing"), plot.SECdistrMv)
setMethod("plot", signature(x="selm"), plot.selm) # y="missing" not required?
setMethod("plot", signature(x="mselm"), plot.mselm)

setMethod("show", signature(object="summary.selm"),  print.summary.selm)
setMethod("show", signature(object="summary.mselm"),  print.summary.mselm)

setMethod("summary", signature(object="SECdistrUv"), summary.SECdistrUv)
setMethod("summary", signature(object="SECdistrMv"), summary.SECdistrMv)
setMethod("summary", signature(object="selm"), summary.selm)
setMethod("summary", signature(object="mselm"), summary.mselm)

setMethod("fitted", signature(object="selm"), fitted.selm)
setMethod("fitted", signature(object="mselm"), fitted.mselm)

setMethod("residuals", signature(object="selm"), residuals.selm)
setMethod("residuals", signature(object="mselm"), residuals.mselm)
 
