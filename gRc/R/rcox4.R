
## RCOX models; internal representation (slot intRep)
## eccI[[i]]: Colour class as px2 matrix (lowest index always in 1. column)
## vccU[[i]]: Colour class as px2 matrix
## eccV[[i]]: Colour class as vector; each element is a pair coded as a 6-digit numer
## vccV[[i]]: Colour class as vector; 

rcox <- function(gm=NULL, vcc=NULL, ecc=NULL, 
                 type   = c('rcon','rcor'),
                 method = "ipm",
                 fit=TRUE, data=NULL, S=NULL, n=NULL, Kstart=NULL,
                 control=list(),
                 details=1,
                 trace=0){

  type           <- match.arg(type)
  method         <- match.arg(method,c("ipm", "scoring", "matching", "ipms", "hybrid1"))
  modelSpec      <- list(gm=gm, vcc=vcc, ecc=ecc)

  ####
  #### Control list
  ####
  con <- list(
              vccfit   = TRUE,
              eccfit   = TRUE,
              vcov     = "inf", #{if(method=="scoring" || method=="ipm") "inf" else NULL},
              nboot    = 100, 
              maxouter = 500,
              maxinner = 10,
              logL     = FALSE,
              logLeps  = 1e-6,
              deltaeps = 1e-3,
              hybrid1switch = 10,
              short    = FALSE
              )
  con[(namc <- names(control))] <- control

  ## List representation (only temporary)
  gmN   <- formula2names(gm)
  vccN  <- formula2names(vcc)
  eccN  <- formula2names(ecc)
  usedVars <- unique(unlist(c(gmN, vccN, eccN)))

  ####
  #### Data representation
  ####
  if (trace>=2)
    cat("..Building data representation\n")

  dataRep  <- .buildDataRepresentation (data, S, n, usedVars, type, trace)
  nodes    <- dataRep$nodes

    if (trace>=2)
    cat("..Building data representation - done\n")

  ####
  #### Standard representation (remove gm part)
  ####
  if (trace>=2)
    cat("..Building standard representation\n")

  stdRep   <- .buildStandardRepresentation (gmN, vccN, eccN,
                                            dataNames=dataRep$dataNames, trace)


  ## List representation of model (with names of variables)
  ##
  vccN     <- .addccnames(stdRep$vccN, type="vcc")

  eccN     <- .addccnames(stdRep$eccN, type="ecc")

  if (trace>=2)
    cat("..Building internal representation\n")

  ####
  #### Internal representation (matrices and vectors)
  ####
  intRep <- .buildInternalRepresentation(vccN=vccN, eccN=eccN,
                                         dataNames=dataRep$dataNames, trace)

  ####
  #### Create model object
  ####
  ans <- structure(list(vcc     = vccN,                        
                        ecc     = eccN,
                        ###dim     = length(vccN)+length(eccN),
                        nodes   = nodes,
                        intRep  = intRep,
                        dataRep = dataRep,
                        Kstart  = Kstart,
                        type    = type,
                        method  = method,
                        trace   = trace,
                        control = con
                        ),
                   class=c(type, "rcox"))

  ####
  #### Fit if required
  ####
  if (fit){
    ans$fitInfo <- fit(ans, method=method, trace=trace, returnModel=FALSE)
  }  
  return(ans)
}

print.rcox <- function(x, ...){

  cat(toupper(getSlot(x,"type")), "model: ")

  if (!is.null(x$fitInfo)){
    cat("logL=", x$fitInfo$logL, " dimension=", dimension(x),
        " method=", x$method, " time=", fitInfo(x,"time"),sep='')
  }
  cat("\n")

  if (!(x$control$short)){
    xcc <- getvcc(x)
    xcc <- getecc(x)
    
    xf  <- names2formula(xcc)
    xs  <- formula2string(xf)
    cat("vcc: ",paste(unlist(xs),collapse=", "),"\n")
    
    if (length(xcc)){
      xf  <- names2formula(xcc)
      xs  <- formula2string(xf)
      cat("ecc: ",paste(unlist(xs),collapse=", "),"\n")
    }
  }

}


.buildDataRepresentation <- function(data=NULL, S=NULL, n=NULL, nodes, type="rcon",
                                     trace=2){
  if (is.null(data) & is.null(S)){
    stop("No data given...\n")
  }

  if (!is.null(data)){    
    dataNames <- names(data)  
    S         <- cov(data)
    n         <- nrow(data)##+1
  } else {
    dataNames <- colnames(S)
  }

  nodes     <- dataNames[sort(match(nodes, dataNames))]

  S   <-  S[nodes, nodes]  
  ans <- list(S=S, n=n, dataNames=rownames(S), nodes=nodes)

  return(ans)
}

.buildStandardRepresentation <- function(gmN, vccN, eccN, dataNames, trace=2){

  ## Get from formulas to lists
  ##
  t0 <- proc.time()
  
  ## Get vertices/edges from gm-spec
  ##
  idx         <- unlistPrim(lapply(gmN, length))
  gmNvertices <- lapply(uniquePrim(unlistPrim(gmN)),list) ## lapply(gmN[idx==1], list)
  x           <- unlist(lapply(gmN[idx>1], names2pairs),recursive=FALSE)
  gmNedges    <- lapply(x, list)

  ## Make standard representation
  ##
  eccN <- c(eccN, gmNedges)
  uuu  <- unlistPrim(eccN)  
  uuu  <- lapply(uuu, as.list)
  vccN <- uniquePrim(c(uuu, vccN, gmNvertices))

  vccI <- names2indices(vccN, dataNames, matrix=FALSE)
  eccI <- names2indices(eccN, dataNames, matrix=FALSE)

  #vccI <<- vccI
  #ri   <- .redundant.index(vccI)

  xxxx2 <- lapply(vccI, function(x3)
                   {z<-do.call("rbind",x3); z[,1]<-z[,1]*10000; rowSumsPrim(z)})
  ri <- which(removeRedundant(xxxx2, index=TRUE)>0)

  vccN <- vccN[ri]

  #ri   <- .redundant.index(eccI)
  if (length(eccI)){
    xxxx2 <- lapply(eccI, function(x3)
                    {z<-do.call("rbind",x3); z[,1]<-z[,1]*10000; rowSumsPrim(z)})
    ri <- which(removeRedundant(xxxx2, index=TRUE)>0)
    eccN <- eccN[ri]
  }
  
  
  varNames <- uniquePrim(unlistPrim(c(vccN,eccN)))

  ans <- list(vccN=vccN, eccN=eccN, varNames=varNames)
  return(ans)
}


.redundant.index <- function(xxxx){
  
  if (length(xxxx)==0)
    return(FALSE) ## Not really intuitive, but it works...
  else
    if (length(xxxx)==1)
      return(TRUE)
  
  xxxx2 <- lapply(xxxx, function(x3)
                  {z<-do.call("rbind",x3); z[,1]<-z[,1]*10000; rowSums(z)})
  ind <- rep(TRUE, length(xxxx2))
  for (i in 1:length(xxxx2)){
    xi <- xxxx2[[i]]
    for (j in 1:length(xxxx2)){
      if (i != j){
        xj <- xxxx2[[j]]      
        if (subsetof(xj,xi) && ind[i])
          ind[j] <- FALSE      
      }
    }
  }
  ind
}


.buildInternalRepresentation <- function(vccN, eccN, dataNames, trace=2){

  vccI <- names2indices(vccN, dataNames, matrix=TRUE)
  eccI <- names2indices(eccN, dataNames, matrix=TRUE)

  eccV <- indices2vectors(eccI)
  vccV <- indices2vectors(vccI)
  
  ans <- list(vccI=vccI, eccI=eccI, vccV=vccV, eccV=eccV)
  return(invisible(ans))
}


