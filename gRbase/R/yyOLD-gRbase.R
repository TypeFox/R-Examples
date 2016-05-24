
## ## #################################################################
## ##
## ## Add and drop edges
## ##
## ## #################################################################

#dropEdge <- function(object, name.1, name.2) UseMethod("dropEdge")
dropEdge.gModel <-
          function(object,name.1,name.2) {

            ## cat("Drop:",name.1,name.2,"\n",sep=" ")
            ## edit hllm formula
            form <- formula(object)
            listform <- readf(form[2])
            new.form <- .delete.edge(listform,c(name.1,name.2))

            form <- paste("~",showf(new.form))
            formula(object) <- as.formula(form)

            if (inherits(object,"gRfit"))
              object <- fit(object)

            return(object)
          }


#addEdge <- function(object, name.1, name.2) UseMethod("addEdge")
addEdge.gModel <-
          function(object,name.1,name.2) {

            new.object <- object
            ## edit hllm formula
            form <- formula(object)
            listform <- readf(form[2])
            new.form <- .add.edge(listform,c(name.1,name.2))
            form <- paste("~",showf(new.form))
            formula(new.object) <- as.formula(form)

            if (inherits(new.object,"gRfit"))
              new.object <- fit(new.object)

            return(new.object)
          }
















## dropVertex <- function(object, name) UseMethod("dropVertex")
## dropVertex.gModel <-
##   function(object,name) {
##             ## edit hllm formula
##             form <- formula(object)
##             listform <- readf(form[2])

##             ## delete 'name' from generators
##             new.form <- lapply(listform,setdiff,name)
##             form <- paste("~",showf(new.form))
##             formula(object) <- as.formula(form)

##             if (inherits(object,"gRfit"))
##               object <- fit(object)

##             return(object)
##           }


## addVertex <- function(object, name) UseMethod("addVertex")
## addVertex.gModel <-
##           function(object,name) {
##             ## edit formula
##             form <- formula(object)
##             listform <- readf(form[2])
##             listform[[length(listform)+1]] <- name
##             form <- paste("~",showf(listform))
##             formula(object) <- as.formula(form)
##             if (inherits(object,"gRfit"))
##               object <- fit(object)

##             return(object)
##             }



ggm <- function(formula=~.^1, gmData, marginal){
  value <- processFormula(formula,gmData, marginal,"Continuous")
  value$gmData <- gmData
  class(value) <- c("ggm","gModel")
  return(value)
}

fit.ggm <- function(object, ...){
  Ydf  <- observations(object$gmData)
  nobs <- nrow(Ydf)
  gc <- object$numformula
  Ymat <- as.matrix(Ydf)
  Smat   <- cov(Ymat)*(nobs-1)/nobs
  ipsFit <- ips(gc,Smat)
  fit      <- outfun( ipsFit$MLE, Smat,nrow(Ydf))
  fit$n    <- nobs
  fit$mean <- apply(Ymat,2,mean)
  fit$df   <- length(which(fit$part==0))/2
  fit$iterations <- ipsFit$iterations
  value<-object
  value$fit <- fit
  class(value) <- c("gRfit", "ggm",class(object))
  return(value)
}

## Partial correlation matrix ##
## computes partial correlation matrix for covariance matrix ##
partial.corr.matrix <- function(S){
  A <- solve(S)
  temp <- diag(1/sqrt(diag(A)))
  temp <- zapsmall(-temp%*%A%*%temp)
  diag(temp) <- 1
  return(temp)
}

## Output function ##
outfun <- function(Sigma, S, n){
  return(list(Sigma=round(Sigma,3),
              eigenvalues=eigen(Sigma)[[1]],
              correlation=cov2cor(Sigma),###corr.matrix(Sigma),
              partial.correlations=partial.corr.matrix(Sigma),
              loglik=ell(Sigma,S,n)))
}




## UNDIRECTED GRAPHS ###

## cliques must be a list with all cliques ##
## the components of this list must be ##
## vectors enumerating the vertices in a clique ##
ips <- function(cliques, S){
  if(!is.matrix(S)){
    return("Second argument is not a matrix!")
  }
  if(dim(S)[1]!=dim(S)[2]){
    return("Second argument is not a square matrix!")
  }
  if(min(eigen(S)[[1]])<=0){
    return("Second argument is not a positive definite matrix!")
  }
  start <- diag(diag(S)) # starting value
  p <- dim(S)[1] # dimensionality
  K <- solve(start)
  i <- 0
  if(length(cliques)==1){
    return(list(MLE=S, iterations=1))
  }
  my.complement <- function(C) return(setdiff(1:p,C))
  cliques.complements <- lapply(cliques, my.complement)
  repeat{
    K.old <- K
    i <- i+1
    for(j in 1:length(cliques)){
      C <- cliques[[j]]
      notC <- cliques.complements[[j]]
      K[C,C] <- solve( S[C,C] ) +
        K[C,notC]%*%solve(K[notC,notC])%*%K[notC,C]
    }
    if(sum(abs(K.old-K)) < 1e-10) break
  }
  return(list(MLE=solve(K), iterations=i))
}

globalVariables(c("rawdata", "loglm.formula"))

##gmData.R ---
##Author          : Claus Dethlefsen
##Created On      : Mon May 02 09:34:40 2005
##Last Modified By:
##Last Modified On:
##Update Count    : 0
##Status          : Unknown, Use with caution!
##

### Some generic functions

"latent.gmData" <- function(x){attr(x,"latent")}
"latent" <- function(x) UseMethod("latent")

"latent<-.gmData" <- function(tmp,value){attr(tmp,"latent")<-value; return(tmp)}
"latent<-" <- function(tmp,value) UseMethod("latent<-")

valueLabels.gmData<- function(x) attr(x,"valueLabels")
valueLabels       <- function(x) UseMethod("valueLabels")

"valueLabels<-.gmData"<- function(tmp,value){attr(tmp,"valueLabels")<-value; return(tmp)}
"valueLabels<-"       <- function(tmp,value) UseMethod("valueLabels<-")

observations.gmData <- function(x) attr(x,"observations")
observations    <- function(x) UseMethod("observations")
obs             <- function(x) UseMethod("observations")

"observations<-.gmData"<- function(tmp,value){attr(tmp,"observations")<-value; return(tmp)}
"observations<-"       <- function(tmp,value)UseMethod("observations<-")

## "description.gmData" <- function(x){attr(x,"description")}
## "description" <- function(x) UseMethod("description")

"description<-.gmData" <- function(tmp,value){attr(tmp,"description")<-value; return(tmp)}
"description<-" <- function(tmp,value) UseMethod("description<-")

"varTypes.gmData" <- function(x){structure(x$varTypes, .Names=varNames(x))}
"varTypes" <- function(x) UseMethod("varTypes")

"varTypes<-.gmData" <- function(tmp,value){ tmp$varTypes <-value; return(tmp)}
"varTypes<-" <- function(tmp,value) UseMethod("varTypes<-")

varNames.gmData <- function(x)as.vector(x$varNames)
varNames <- function(x)UseMethod("varNames")

"varNames<-.gmData" <- function(tmp,value){ tmp$varNames <-value; return(tmp)}
"varNames<-" <- function(tmp,value) UseMethod("varNames<-")

nLevels.gmData <- function(x)structure(as.vector(x$nLevels), .Names=varNames(x))
nLevels <- function(x)UseMethod("nLevels")

"nLevels<-.gmData" <- function(tmp,value){ tmp$nLevels <-value; return(tmp)}
"nLevels<-" <- function(tmp,value) UseMethod("nLevels<-")


shortNames.gmData <- function(x)structure(as.vector(x$shortNames), .Names=varNames(x))
shortNames <- function(x)UseMethod("shortNames")

"shortNames<-.gmData" <- function(tmp,value){ tmp$shortNames <-value; return(tmp)}
"shortNames<-" <- function(tmp,value) UseMethod("shortNames<-")

dataOrigin.gmData   <- function(x) attr(x,"dataOrigin")[1]
dataOrigin   <- function(x)UseMethod("dataOrigin")

"ordinal"           <- function(tmp) UseMethod("ordinal")
"ordinal<-"         <- function(tmp,value) UseMethod("ordinal<-")

"ordinal.gmData" <- function(tmp)attr(tmp,"ordinal")

"ordinal<-.gmData" <- function(tmp,value){
  varTypes(tmp)[match(value, varNames(tmp))]<-"Ordinal"
  return(tmp)}

"nominal"           <- function(tmp) UseMethod("nominal")
"nominal<-"         <- function(tmp,value) UseMethod("nominal<-")


"nominal.gmData" <- function(tmp){
  varNames(tmp)["Discrete"==varTypes(tmp)]
}

"nominal<-.gmData" <- function(tmp,value){
  varTypes(tmp)[match(value, varNames(tmp))]<-"Discrete"
  return(tmp)}






##################################################################################
as.gmData       <- function(from) UseMethod("as.gmData")
##################################################################################

print.gmData  <- function(x, ...){
  xx<-attr(x,"description")
  if (!is.null(xx))
    cat("Description:", xx, "\n")
  print.data.frame(x);
  ##cat("Data origin:     ", .dataOrigin(x),"\n")
  if (!is.null(latent(x)))
    cat ("Latent variables:", paste(latent(x),collapse=' '), "\n")
  if (!is.null(valueLabels(x)))
  cat("To see the values of the factors use the 'valueLabels' function\n")
  if (!is.null(observations(x)))
  cat("To see the data use the 'observations' function\n")
  return(invisible(x))
}

# summary.gmData  <- function(object, ...){
#   print(table(object$varTypes))
#   if (!is.null(observations(object))) {
#     cat("\nObservation summary:\n")
#     print(summary(obs(object)))
#   }
#   invisible(object)
# }



summary.gmData <- function(object, ...){
  print(object)
  mapply(function(xx,ll){
    cat("Factor:", ll, "\n Levels:", paste(xx,sep=' '),"\n")
  }, valueLabels(object),names(valueLabels(object)))
  return(invisible(object))

}






#### ##############################################################

# newgmData <- function(varNames,
#                    varTypes=rep(validVarTypes()[1],length(varNames)),
#                    nLevels=NA,
#                    latent=NA,
#                    valueLabels=NULL,
#                    observations=NULL,
#                    description=NA,
#                    shortNames=c(letters,LETTERS)
#                    ){
#   value <- data.frame(varNames, abbreviate(varNames,1),row.names=NULL)

#   names(value) <- c("varNames","shortNames")
#   value$varTypes <- factor(varTypes,levels=validVarTypes())
#   value$nLevels  <- nLevels

#   obsclass <- class(observations)
#   class(value) <- c("gmData","data.frame")


#   attr(value,"valueLabels")    <- valueLabels
#   attr(value,"latent")         <- latent
#   attr(value,"description")    <- description
#   attr(value,"observations")   <- observations
#   ##switch(class(data),
#   switch(obsclass,
#          "table"=     { attr(value,"dataOrigin")     <- "table"      },
#          "data.frame"={ attr(value,"dataOrigin")     <- "data.frame" },
#          NULL=        { attr(value,"dataOrigin")     <- "table"      })
#   return(value)
# }



newgmData <-
function (varNames,
          varTypes = rep(validVarTypes()[1], length(varNames)),
          nLevels  = NULL,
          latent   = NULL,
          valueLabels  = NULL,
          observations = NULL,
          description  = NULL,
          shortNames   = NULL)
{

  cl <- match.call()

  .is.subset <- function(x,y){
    setequal(intersect(x,y),x)
  }

  .simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
  }

  ## Find good short names...
  ##
  if (is.null(shortNames)){
    nam <- varNames
    nama  <- abbreviate(nam,1)
    nc    <- nchar(nama)
      rest  <- setdiff(c(letters,LETTERS), nama[nc==1])
    if (length(which(nc>1)) <= length(rest))
      nama[nc>1]<- rest[1:length(which(nc>1))]
  } else {
    nama <- shortNames
  }

  value <- data.frame(varNames, nama, row.names = NULL)
  names(value) <- c("varNames", "shortNames")

  ## Deal with abbreviated varTypes
  ##
  varTypes <- sapply(varTypes, .simpleCap)
  varTypes <- sapply(varTypes, match.arg, choices=validVarTypes(), several.ok = FALSE)

  value$varTypes <- factor(varTypes, levels = validVarTypes())
  discidx        <- which("Discrete"==varTypes | "Ordinal"==varTypes)
  aa             <- rep(NA, length(varNames))

  ## If valueLabels=c(1,2,3) turn into list(c(1,2,3))
  ##
  if (!is.null(valueLabels) & !is.list(valueLabels))
    valueLabels <- list(valueLabels)

  if (is.null(nLevels) & is.null(valueLabels)){
    ## If neither nLevels or valueLabels are given; make all
    ## categorical variables binary
    ##
    aa[discidx]   <-  2
    nLevels <- aa
  }

  if (!is.null(valueLabels)){
    ## If valueLabels are given, infer nLevels from these; recycle if necessary...
    ##
    if (!.is.subset(varNames[discidx], names(valueLabels))){
      vl <- rep(valueLabels, length(discidx))
      valueLabels <- vl[1:length(discidx)]
      names(valueLabels) <- varNames[discidx]
    }
    uu            <- valueLabels[varNames[discidx]]
    uu            <- sapply(uu,length)
    aa[discidx]   <- uu
    value$nLevels <- unlist(aa)
  } else {
    ## Use nLevels as given; recycle if necessary
    ## Infer valueLabels from these
    v <- nLevels[discidx]
    v <- v[!is.na(v)]
    if (length(v)==0)
      v <- 2
    v <- rep(v, length(discidx))
    v <- v[discidx]
    aa[discidx]   <-  v
    value$nLevels <- unlist(aa)
    uu <- varNames[discidx]
    ##v  <<- v

    valueLabels <- mapply(function(nn,vv){paste(nn,1:vv,sep='')},uu,v, SIMPLIFY=FALSE)

  }


  class(value) <- c("gmData", "data.frame")
  attr(value, "valueLabels")  <- valueLabels
  attr(value, "latent")       <- latent
  attr(value, "description")  <- description
  attr(value, "observations") <- observations
  attr(value, "dataOrigin")   <- class(observations)

  obsclass <- class(observations)

  if (is.null(obsclass)){
    attr(value, "dataOrigin") <- NULL
  } else {
    if(is.element("table", obsclass))
      attr(value, "dataOrigin") <- c("table",setdiff(obsclass, "table"))
    else{
      if(is.element("data.frame", obsclass))
        attr(value, "dataOrigin") <- c("data.frame", setdiff(obsclass, "data.frame"))
      else
        attr(value, "dataOrigin") <- "other"
    }
  }



  return(value)
}







#     switch(class(observations),
#     table = {
#         attr(value, "dataOrigin") <- "table"
#     }, data.frame = {
#         attr(value, "dataOrigin") <- "data.frame"
#     }, "NULL" = {
#         attr(value, "dataOrigin") <- "table"
#     })

validVarTypes <- function() {c("Discrete","Ordinal","Continuous")}


## ####################################################################
## Convert data.frame into gmData

as.gmData.data.frame <- function(from){
  fact   <- unlist(lapply(1:ncol(from), function(j)
                          is.factor(from[,j])))
  Types <- rep(validVarTypes()[3],length(fact))
  Types[fact] <- validVarTypes()[1]

  levels <- unlist(lapply(1:ncol(from),
                          function(j)
                          {
                            if(is.factor(from[,j]))
                              length(levels(from[,j]))
                            else NA}
                          )
                   )

  if (length(which(fact))>0){
    vallabels <- list()
    for (j in which(fact)){
      vallabels <- c(vallabels, list(levels(from[,j])))
    }
    names(vallabels) <- names(from[which(fact)])
  } else {
    vallabels <- list()
  }

  newgmData(
      varNames=names(from),
      varTypes=Types,
      nLevels=levels,
      valueLabels=vallabels,
      observations=from
 )
}



## ####################################################################
## Convert table into gmData

as.gmData.table <- function(from){
  counts <- as.vector(from)
  dn     <- dimnames(from)
  name   <- names(lapply(dn,function(x)names(x)))
  dim    <- unlist(lapply(dn,length))
  newgmData(
         varNames=name,
         varTypes=rep("Discrete",length(name)),
         nLevels=dim,
         valueLabels=dn,
         observations=from
         )
}


## ####################################################################
## Convert array into gmData

as.gmData.array <- function(from){
  res <- as.gmData(as.table(from))
  observations(res) <- from
  res
}

##                              -*- Mode: Ess -*-
##gModel.R ---
##Author          : Claus Dethlefsen
##Created On      : Mon May 02 09:35:24 2005
##Last Modified By:
##Last Modified On:
##Update Count    : 0
##Status          : Unknown, Use with caution!
##


gModel <- function(formula, gmData){
  value <- list(formula=formula, gmData=gmData)
  class(value) <- "gModel"
  return(value)
}

"formula<-.gModel" <- function(tmp,value){tmp$formula<-value; return(tmp)}
"formula<-" <- function(tmp,value) UseMethod("formula<-")

"gmData.gModel" <- function(x){x$gmData}
"gmData" <- function(x) UseMethod("gmData")

"gmData<-.gModel" <- function(tmp,value){tmp$gmData<-value; return(tmp)}
"gmData<-" <- function(tmp,value) UseMethod("gmData<-")

print.gModel <- function(x, ...){
  cat("Model information (gRbase)\n")
  cat(" Class:   ", paste(class(x),collapse=' <- '),"\n")
  cat(" Formula: ", paste(paste(x$formula),collapse=''), "\n")
}

## necessary for use with dynamicGraph
## setOldClass("gModel")

##gRfit.R ---
##Author          : Claus Dethlefsen
##Created On      : Mon May 02 09:35:51 2005
##Last Modified By:
##Last Modified On:
##Update Count    : 0
##Status          : Unknown, Use with caution!
##



"getFit.gRfit" <- function(x){x$fit}
"getFit" <- function(x) UseMethod("getFit")

"getFit<-.gRfit" <- function(tmp,value){ tmp$fit <-value; return(tmp)}
"getFit<-" <- function(tmp,value) UseMethod("getFit<-")

print.gRfit <- function(x, ...){
  print.gModel(x)
  cat("Fit information (gRbase)\n")
  cat("   logL", deviance(getFit(x)), "df", x$fit$df,"\n")
}

summary.gRfit <- function(object,...)
  summary(getFit(object))



hllm <- function(formula = ~.^1,  gmData, marginal){
  stop("function 'hllm' from gRbase is defunct. Please use the gRim package for hierarchical log-linear models.")

##   value <- processFormula(formula, gmData, marginal,"Discrete")
##   value$gmData <- gmData
##   class(value) <- c("hllm","gModel")
##   return(value)
}

fit.hllm <- function(object,engine="loglm",...){
  stop("function 'fit.hllm' from gRbase is defunct. Please use the gRim package for hierarchical log-linear models.")
##   rawdata <- observations(object$gmData)
##   if (is.data.frame(rawdata)){
##     rawdata <- xtabs(~., rawdata)
##   }
##   value <- object
##   mimform <- processFormula(formula(object),gmData(object),type="Discrete")
##                                         #cat("processFormula done...\n")
##   switch(engine,
##          "loglm"={
##            loglm.formula <- mimform$formula
##            val <- loglm(loglm.formula, data=rawdata)
##            val$call$data <- rawdata
##            class(value) <- c("gRfit","loglm", class(object))
##          },
##          {stop("Only engine 'loglm' currently implemented...")
##         })
##   value$fit <- val
##   return(value)
}

stepwise.hllm <-    function (object, ...)  {
  stop("function 'stepwise.hllm' from gRbase is defunct. Please use the gRim package for hierarchical log-linear models.")
##   if (!exists("rawdata",envir=.GlobalEnv)&
##       !exists("loglm.formula",envir=.GlobalEnv)) {
##     assign("rawdata",observations(gmData(object)),.GlobalEnv)
##     assign("loglm.formula",formula(object),.GlobalEnv)

##     if (!inherits(object,"gRfit"))
##       object <- fit(object)

##     res <- step(getFit(object),...)
##     gRobj <- object
##     formula(gRobj) <- formula(res$formula)
##     getFit(gRobj) <- res
##     rm(rawdata,loglm.formula,envir=.GlobalEnv)
##     return(gRobj)
##   }
##   else
##     stop("You will have to move/rename rawdata and loglm.formula from .GlobalEnv\n")
}




# loglmSHD <- function (formula, data, subset, na.action, ...)
# {
#     .call <- match.call()
#     if (missing(data) || inherits(data, "data.frame")) {
#         m <- match.call(expand.dots = FALSE)
#         m$... <- NULL
#         m[[1]] <- as.name("model.frame")
#         data <- eval.parent(m)
#         .formula <- as.formula(attr(data, "terms"))
#     }
#     else {
#       trms <- attr(data, "terms") <- terms(formula <- denumerate(formula))
#       .formula <- renumerate(as.formula(trms))
#     }
#     loglm1(formula, data, ..., .call = .call, .formula = .formula)
# }



## .fit.hllm <- function(m,engine="loglm",...){

##   rawdata <- observations(m$gmData)
##   if (is.data.frame(rawdata)){
##     rawdata <- xtabs(~., rawdata)
##   }
##   value <- m
##   mimform <- processFormula(formula(m),gmData(m),type="Discrete")
##   switch(engine,
##          "loglm"={
##            mimformula <- mimform$mimformula
##            loglm.formula <- formula(paste("~",mimformula))
##            ##val <- loglm(loglm.formula, rawdata)
##            val <- loglm(loglm.formula, data=rawdata)
##            val$call$data <- rawdata
##            class(value) <- c("gRfit","loglm", class(m))
##          },
##          {stop("Only engine 'loglm' currently implemented...")
##         })
##   value$fit <- val
##   return(value)
## }










# ## example of how to extende the user menu in dynamicGraph
# UserMenus <-
#   list(
#        MainUser =
#        list(label = "Stepwise",
#               command = function(object, ...)
#             stepwise(object,...)
#               )
#           )


# for use with dynamicGraph
##setOldClass("hllm")
##setIs("hllm","gModel")

update.gModel <- function(object, addedge=NULL, dropedge=NULL, ...){

  if (!is.null(addedge))
    object <- addEdge.gModel(object, addedge[1], addedge[2])

  if (!is.null(dropedge))
    object <- dropEdge.gModel(object, dropedge[1], dropedge[2])

  return(object)


}

