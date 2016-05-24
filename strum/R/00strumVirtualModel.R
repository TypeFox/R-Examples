#==============================================================================
# File: strumVirtualModel.R
#
# Author: Nathan Morris
#
# Notes: strumVirtualModel class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumVirtualModel class
#------------------------------------------------------------------------------

setClass("strumVirtualModel",
         representation(
           varList          = "data.frame",
           formulas         = "character",
           allRandomEffects = "character",
           paramNames       = "character",
           ascertainment    = "ANY",
           E                = "list",
           Z                = "list",
           L                = "function",
           B                = "function",
           Gs               = "function",
           Gm               = "function",
           thToThB          = "numeric"),
         representation("VIRTUAL"),
         prototype = list(
           varList          = data.frame(), #names(varList)=c("name","type","obs","covariate","exogen")
           formulas         = character(),
           allRandomEffects = character(),
           paramNames       = character(),
           ascertainment    = NULL,
           E                = list(),
           Z                = list(),
           L                = function(onefamilydf) return(TRUE),
           B                = function(onefamilydf) return(TRUE),
           Gs               = function(onefamilydf) return(TRUE),
           Gm               = function(onefamilydf) return(TRUE),
           thToThB          = numeric())
) 

#------------------------------------------------------------------------------
# 'covariatePhenotypes' accessor functions:
#------------------------------------------------------------------------------
setGeneric('covariatePhenotypes', function(object) standardGeneric('covariatePhenotypes'))
setMethod('covariatePhenotypes', 'strumVirtualModel',
          function(object)
          {
            return(as.character(object@varList$name[object@varList$covariate & object@varList$distribution!="mendelian"]))
          }
)

#------------------------------------------------------------------------------
# 'varList' accessor functions:
#------------------------------------------------------------------------------
setGeneric('varList', function(object) standardGeneric('varList'))
setMethod('varList', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@varList)
          }
)

setGeneric('varList<-', function(object,value) standardGeneric('varList<-'))
setMethod('varList<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@varList <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'formulas' accessor functions:
#------------------------------------------------------------------------------
setGeneric('formulas', function(object) standardGeneric('formulas'))
setMethod('formulas', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@formulas)
          }
)

setGeneric('formulas<-', function(object,value) standardGeneric('formulas<-'))
setMethod('formulas<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@formulas <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'allRandomEffects' accessor functions:
#------------------------------------------------------------------------------
setGeneric('allRandomEffects', function(object) standardGeneric('allRandomEffects'))
setMethod('allRandomEffects', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@allRandomEffects)
          }
)

setGeneric('allRandomEffects<-', function(object,value) standardGeneric('allRandomEffects<-'))
setMethod('allRandomEffects<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@allRandomEffects <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'paramNames' accessor functions:
#------------------------------------------------------------------------------
setGeneric('paramNames', function(object) standardGeneric('paramNames'))
setMethod('paramNames', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@paramNames)
          }
)

setGeneric('paramNames<-', function(object,value) standardGeneric('paramNames<-'))
setMethod('paramNames<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@paramNames <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'ascertainment' accessor functions:
#------------------------------------------------------------------------------
setGeneric('ascertainment', function(object) standardGeneric('ascertainment'))
setMethod('ascertainment', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@ascertainment)
          }
)

setGeneric('ascertainment<-', function(object,value) standardGeneric('ascertainment<-'))
setMethod('ascertainment<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@ascertainment <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumVirtualModel"),
          function(object) 
          {
            .showModel(object, "strumVirtualModel")
          }
)

#------------------------------------------------------------------------------
# plot generic functions
#------------------------------------------------------------------------------
setMethod("plot", "strumVirtualModel",
          function(x, y, name="strumVirtualModel", toFile=TRUE, fileType="dot", ...) 
          {
            if( missing(y) )
              y = "dot"

            .plotModel(x, layoutType=y, name=name, toFile=toFile, fileType=fileType, ...)
          }
)

#------------------------------------------------------------------------------
# Common function to show Model class 
#------------------------------------------------------------------------------
.showModel = function(object, name) 
{
  cat("\nBasic properties of the model:\n")
    .printInfoLine("Model Class", name, 40)

  if( is.null(object@ascertainment) )
    .printInfoLine("Ascertainment", "FALSE", 40)
  else
    .printInfoLine("Ascertainment", "TRUE", 40)

  cat("\nList of all variables:\n")
  myvarList = object@varList
  myrowNames = myvarList$name
  myvarList = myvarList[,names(myvarList)!="endogenous" & names(myvarList)!="name" ]
  mynames = names(myvarList)
  mynames[mynames=="variableLevels"] = "Level"
  mynames[mynames=="exogenous"]      = "Exogen."
  mynames[mynames=="distribution"]   = "Distr."
  mynames[mynames=="parameter"]      = "Param."
  mynames = paste(toupper(substring(mynames, 1,1)), substring(mynames, 2), sep="")
  names(myvarList) = mynames
  print(myvarList, row.names = myrowNames)

  cat("\nModel formulas:\n")
  cat(object@formulas)
}

#------------------------------------------------------------------------------
# Common function to plot Model class 
#------------------------------------------------------------------------------
.plotModel = function(model, layoutType, name, toFile, fileType,
                      rankDir, nodeColor, fontColor, nodeFont, edgeColor, size, ...)
{
  B  = model@B(NA)
  L  = model@L(NA)
  Gm = model@Gm(NA)
  Gs = model@Gs(NA)
  Z  = lapply(model@Z, function(v) v(NA))
  E  = lapply(model@E, function(v) v(NA))

  varList = model@varList
  eitaNames = varList$name[varList$inEita]
  yNames    = varList$name[varList$inY]
  covNames  = c("[intercept]",varList$name[varList$covariate])
  reNames   = model@allRandomEffects

  Bij = .getParamIndex(B)
  Lij = .getParamIndex(L)
  
  Gsij = .getParamIndex(matrix(Gs[,-1], nrow=nrow(Gs)))
  Gmij = .getParamIndex(matrix(Gm[,-1], nrow=nrow(Gm)))
  
  frB = .getNodeNames(Bij, eitaNames, 2, 0)
  toB = .getNodeNames(Bij, eitaNames, 1, 0)
  
  frL = .getNodeNames(Lij, eitaNames, 2, 0)
  toL = .getNodeNames(Lij, yNames,    1, 0)
  
  frGs = .getNodeNames(Gsij, covNames, 2, 1)
  toGs = .getNodeNames(Gsij, eitaNames, 1, 0)
  
  frGm = .getNodeNames(Gmij, covNames, 2, 1)
  toGm = .getNodeNames(Gmij, yNames, 1, 0)

  ftZ = .getVCNodeEdge(Z, eitaNames, "Z", reNames)
  ftE = .getVCNodeEdge(E, yNames, "E", reNames)

  frZ = ftZ$fr; frE = ftE$fr; frZCov = ftZ$frCov; frECov = ftE$frCov;
  toZ = ftZ$to; toE = ftE$to; toZCov = ftZ$toCov; toECov = ftE$toCov;
  
  frL = frL[!(frL %in% yNames)]
  toL = toL[!(toL %in% eitaNames)]
  ey = eitaNames %in% yNames
  ye = yNames %in% eitaNames
  yNames = yNames[!ye]
  #eitaNames = eitaNames[!ey]

  fr = c(frB, frL, frGs, frGm, frZ, frE, frZCov, frECov)
  to = c(toB, toL, toGs, toGm, toZ, toE, toZCov, toECov)

  g = ftM2graphNEL(cbind(fr, to))

  # global attributes
  #
  if( missing(rankDir) )
    rankDir = "TB"

  if( missing(nodeColor) )
    nodeColor = "cornflowerblue"

  if( missing(fontColor) )
    fontColor = "white"

  if( missing(nodeFont) )
    nodeFont = "helvetica"

  if( missing(edgeColor) )
    edgeColor = "gray30"

  if( missing(size) )
    size = "5,5"
  
  aN = list(color=nodeColor,
            fillcolor=nodeColor,
            style="filled",
            fontcolor=fontColor,
            fixedsize="false")
  aE = list(color=edgeColor, arrowhead="normal")
  aG = list(rankdir=rankDir, fontname=nodeFont, size=size)

  # node attributes
  #
  nAt = list()

  zNodeLabel = ftZ$vcNodeL
  eNodeLabel = ftE$vcNodeL
  names(zNodeLabel) = frZ
  names(eNodeLabel) = frE

  nAt$label = c(zNodeLabel, eNodeLabel)

  eiShape = rep("ellipse", length(eitaNames))
  ycShape = rep("box",     length(yNames)+length(covNames[-1]))
  vcShape = rep("circle",  length(frZ)+length(frE))
  eiShape[which(ey)] = "box"

  vcStyle   = rep("solid",  length(frZ)+length(frE))
  vcFicolor = rep("transparent",  length(frZ)+length(frE))
  vcFocolor = rep("black",  length(frZ)+length(frE))
  vcWidth   = rep("0.3", length(frZ)+length(frE))
  vcHeight  = rep("0.3", length(frZ)+length(frE))
  vcFixed   = rep("true", length(frZ)+length(frE))

  names(eiShape) = eitaNames
  names(ycShape) = c(yNames, covNames[-1])
  names(vcShape) = names(vcStyle) = names(vcFicolor) = names(vcFocolor) = names(vcWidth) = names(vcHeight) = names(vcFixed) = c(frZ, frE)

  nAt$shape     = c(eiShape, ycShape, vcShape)
  nAt$style     = c(vcStyle)
  nAt$fillcolor = c(vcFicolor)
  nAt$fontcolor = c(vcFocolor)
  nAt$width     = vcWidth
  nAt$height    = vcHeight
  nAt$fixedsize = vcFixed

  # edge attributes
  #
  eAt = list()

  bNames  = if( length(frB) > 0 ) paste(frB, "~", toB, sep="") else c()
  lNames  = if( length(frL) > 0 ) paste(frL, "~", toL, sep="") else c()
  gsNames = if( length(frGs) > 0 ) paste(frGs, "~", toGs, sep="") else c()
  gmNames = if( length(frGm) > 0 ) paste(frGm, "~", toGm, sep="") else c()
  zVarNames = if( length(frZ) > 0 ) paste(frZ, "~", toZ, sep="") else c()
  eVarNames = if( length(frE) > 0 ) paste(frE, "~", toE, sep="") else c()
  zCovNames = if( length(frZCov) > 0 ) paste(frZCov, "~", toZCov, sep="") else c()
  eCovNames = if( length(frECov) > 0 ) paste(frECov, "~", toECov, sep="") else c()

  allW  = rep(1000, length(fr))
  allAh = rep("normal", length(fr))
  
  bColor  = rep("maroon", length(bNames))
  lColor  = rep("black", length(lNames))
  gsColor = rep("navy", length(gsNames))
  gmColor = rep("navy", length(gmNames))
  zVarColor = rep(nodeColor, length(zVarNames))
  eVarColor = rep(nodeColor, length(eVarNames))
  zCovColor = rep(nodeColor, length(zCovNames))
  eCovColor = rep(nodeColor, length(eCovNames))
  
  bLabel  = B[Bij]; bLabel[is.na(bLabel)] = "&beta;"
  lLabel  = L[Lij]; lLabel[is.na(lLabel)] = "&lambda;"
  gsLabel = matrix(Gs[,-1], nrow=nrow(Gs))[Gsij]; gsLabel[is.na(gsLabel)] = "&gamma;"
  gmLabel = matrix(Gm[,-1], nrow=nrow(Gm))[Gmij]; gmLabel[is.na(gmLabel)] = "&gamma;"
  zVarLabel = ftZ$vEdgeL
  zCovLabel = ftZ$cEdgeL
  eVarLabel = ftE$vEdgeL
  eCovLabel = ftE$cEdgeL

  names(allW)    = c(bNames,lNames,gsNames,gmNames,zVarNames,eVarNames,zCovNames,eCovNames)
  names(allAh)   = c(bNames,lNames,gsNames,gmNames,zVarNames,eVarNames,zCovNames,eCovNames)
  names(bColor)  = names(bLabel)  = bNames
  names(lColor)  = names(lLabel)  = lNames
  names(gsColor) = names(gsLabel) = gsNames
  names(gmColor) = names(gmLabel) = gmNames
  names(zVarColor) = names(zVarLabel) = zVarNames
  names(eVarColor) = names(eVarLabel) = eVarNames
  names(zCovColor) = names(zCovLabel) = zCovNames
  names(eCovColor) = names(eCovLabel) = eCovNames
  
  eAt$weight = allW
  eAt$arrowhead = allAh
  eAt$color = c(zVarColor, eVarColor, zCovColor, eCovColor)
  
  # subgraph
  #
  sg1 <- subGraph(unique(frE), g)
  sg2 <- subGraph(unique(toE), g)
  sg3 <- subGraph(eitaNames, g)
  sg4 <- subGraph(unique(frZCov), g)
  
  subGList <- vector(mode="list", length=2)
  subGList[[1]] <- list(graph=sg1, cluster=FALSE, attrs=c(rank="max"))
  subGList[[2]] <- list(graph=sg2, cluster=FALSE, attrs=c(rank="same"))
  subGList[[3]] <- list(graph=sg4, cluster=FALSE, attrs=c(rank="same"))
  #subGList[[4]] <- list(graph=sg4, cluster=FALSE, attrs=c(rank="same"))
  
  invisible(plot(g, y=layoutType, ..., name=name, main=name,
                 attrs = list(node = aN, edge = aE, graph = aG),
                 nodeAttrs = nAt,
                 edgeAttrs = eAt,
                 subGList=subGList))
    
  if( toFile == TRUE )
  {
    eL = c(bLabel, lLabel, gsLabel, zVarLabel, eVarLabel, zCovLabel, eCovLabel)
    if( length(which(names(eL)=="")) )
      eL = eL[-which(names(eL)=="")]
    if( length(which(is.na(names(eL)))) )
      eL = eL[-which(is.na(names(eL)))]
    
    eAt$label = eL
    
    rg = agopen(g, name=name,
                attrs=list(node = aN, edge = aE, graph = aG, cluster=list()),
                nodeAttrs=nAt, edgeAttrs=eAt, subGList=subGList)
    
    invisible(toFile(rg, filename=paste(name, "_", layoutType, ".", fileType, sep=""),
                     layoutType=layoutType, fileType=fileType))
  }
}

.getParamIndex = function(M)
{
  Mna = which(is.na(M), arr.ind=T)
  M1  = which(M!=0, arr.ind=T)
  Mij = rbind(Mna, M1)

  return(Mij)
}

.getNodeNames = function(Mij, nNames, col, offset)
{
  ftN = vector(mode = "character", length = 0)

  if( length(Mij) > 0 )
    ftN = nNames[Mij[,col] + offset]

  return(ftN)
}

.addVC = function(Mv, nNames, vc, reName, ftVC)
{
  Mij = .getParamIndex(Mv)
  
  if( length(Mij) > 0 )
  {
    vcPre   = paste(vc, reName, sep="")
    vcNames = paste(vcPre, sort(unique(c(Mij))), sep="")

    frVar = paste(vcPre, Mij[,2], sep="")
    toVar = .getNodeNames(Mij, nNames, 1, 0)

    vELabel = Mv[Mij]
    
    cov = which(apply(Mij, 1, function(m) {return(m[1] != m[2])}))
    
    if( length(cov) > 0 )
    {
      frVar = frVar[-cov]
      toVar = toVar[-cov]

      frCov = vcNames[Mij[cov,][,2]]
      toCov = vcNames[Mij[cov,][,1]]

      vELabel = Mv[Mij[-cov,]]
      cELabel = Mv[Mij[cov,]]

      vELabel[is.na(vELabel)] = ""
      cELabel[is.na(cELabel)] = ""

      ftVC$frCov = c(ftVC$frCov, frCov)
      ftVC$toCov = c(ftVC$toCov, toCov)
      ftVC$cEdgeL = c(ftVC$cEdgeL, cELabel)
    }

    vcNLabel = rep(reName, length(frVar))

    ftVC$fr = c(ftVC$fr, frVar)
    ftVC$to = c(ftVC$to, toVar)

    ftVC$vcNodeL = c(ftVC$vcNodeL, vcNLabel)
    ftVC$vEdgeL  = c(ftVC$vEdgeL, vELabel)
  }
  
  return(ftVC)
}

.getVCNodeEdge = function(M, nNames, vc, reNames)
{
  fr = vector(mode = "character", length = 0)
  to = vector(mode = "character", length = 0)
  
  frCov = vector(mode = "character", length = 0)
  toCov = vector(mode = "character", length = 0)
  
  vcNodeL = vector(mode = "character", length = 0)
  vEdgeL  = vector(mode = "character", length = 0)
  cEdgeL  = vector(mode = "character", length = 0)

  ftVC = list(fr=fr, to=to, frCov=frCov, toCov=toCov,
              vcNodeL=vcNodeL, vEdgeL=vEdgeL, cEdgeL=cEdgeL)

  for( v in seq_along(reNames) )
  {
    ftVC = .addVC(M[[v]], nNames, vc, reNames[v], ftVC)
  }
  
  return(ftVC)
}
