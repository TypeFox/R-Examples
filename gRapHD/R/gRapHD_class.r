################################################################################
#   This file is part of gRapHD R package.
#
#   gRapHD R package
#   Copyright (C) 2009 Gabriel Coelho Goncalves de Abreu, Rodrigo Labouriau,
#   and David Edwards
#
#   gRapHD R package is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the Free
#   Software Foundation, either version 3 of the License, or any later version.
#
#   gRapHD R package program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

setClass("gRapHD",
         representation(edges="matrix",
                        homog="logical",
                        minForest="integer",
                        numCat="integer",
                        numP="integer",
                        p="integer",
                        stat.minForest="character",
                        stat.stepw="character",
                        stat.user="character",
                        statSeq="numeric",
                        stepw="integer",
                        userDef="integer",
                        vertNames="character"),
         prototype(edges=matrix(integer(0),ncol=2),
                   homog=TRUE,
                   minForest=integer(2),
                   numCat=integer(0),
                   numP=integer(0),
                   p=integer(1),
                   stat.minForest=character(0),
                   stat.stepw=character(0),
                   stat.user=character(0),
                   statSeq=numeric(0),
                   stepw=integer(2),
                   userDef=integer(2),
                   vertNames=character(0)),package="gRapHD")

################################################################################
# Define the function used for initializing a new gRapHD2 object
initialize.gRapHD <- function(.Object,edges=matrix(integer(0),ncol=2),
                                      homog=TRUE,
                                      minForest=integer(2),
                                      numCat=integer(0),
                                      numP=integer(0),
                                      p=integer(1),
                                      stat.minForest=character(0),
                                      stat.stepw=character(0),
                                      stat.user=character(0),
                                      statSeq=numeric(0),
                                      stepw=integer(2),
                                      userDef=integer(2),
                                      vertNames=character(0))
{
  zeros <- function(x,ch,n)
  {
    if (as.integer(log10(x)) >= n)
      stop("x is longer than n.")
    result <- paste(paste(rep(ch,n-(as.integer(log10(x))+1)),sep="",collapse=""),x,sep="")
    return(result)
  }

  .Object@edges <- edges
  .Object@homog <- homog
  .Object@minForest <- as.integer(minForest)
  .Object@numCat <- as.integer(numCat)
  .Object@numP <- as.integer(numP)
  .Object@stat.minForest <- stat.minForest
  .Object@stat.stepw <- stat.stepw
  .Object@stat.user <- stat.user
  .Object@statSeq <- statSeq
  .Object@stepw <- as.integer(stepw)
  .Object@userDef <- as.integer(userDef)
  .Object@vertNames <- vertNames
  .Object@p <- as.integer(max(c(length(unique(as.vector(.Object@edges))),p)))

  if (length(.Object@edges)==0)
  {
    .Object@numP <- integer(0)
    .Object@statSeq <- numeric(0)
  }
  else
    if (NCOL(.Object@edges)==2) # matrix of edges
    {
      n <- NROW(.Object@edges)
      .Object@numP <- rep(.Object@numP,length.out=n)
      .Object@statSeq <- rep(.Object@statSeq,length.out=n)
      a <- .Object@edges[,1]-.Object@edges[,2]
      .Object@edges[a>0,] <- .Object@edges[a>0,c(2,1)] # sort rows
      if (0%in%a)
      {
        warning("Self-loops removed.")
        .Object@edges <- .Object@edges[-which(a==0),]
        .Object@numP <- .Object@numP[-which(a==0)]
        .Object@statSeq <- .Object@statSeq[-which(a==0)]
      }
      z0 <- (.Object@edges[,1]-1)*.Object@p-(.Object@edges[,1]-1)*.Object@edges[,1]/2+.Object@edges[,2]-.Object@edges[,1]
      z1 <- unique(z0)
      if (length(z1)<n) # there are multiple edges
      {
        warning("Multiple edges removed.")
        .Object@edges <- .Object@edges[match(z1,z0),]
        .Object@numP <- .Object@numP[match(z1,z0)]
        .Object@statSeq <- .Object@statSeq[match(z1,z0)]
      }
      .Object@userDef[1] <- as.integer(1)
      .Object@userDef[2] <- as.integer(nrow(.Object@edges))
    }
    else
      stop("Object must have dimension (k,2).")

  if (length(.Object@vertNames)<.Object@p)
    .Object@vertNames <- c(.Object@vertNames,
                           paste("X",apply(matrix((length(.Object@vertNames)+(length(.Object@vertNames)==0)):.Object@p,1),
                                           2,zeros,ch="0",n=as.integer(log10(.Object@p)+1)),sep=""))
  else
    if (length(.Object@vertNames)>.Object@p)
      .Object@vertNames <- .Object@vertNames[1:.Object@p]

  if (length(.Object@numCat)<.Object@p)
    .Object@numCat <- as.integer(c(.Object@numCat,rep(0,.Object@p-length(.Object@numCat))))
  else
    if (length(.Object@numCat)>.Object@p)
      .Object@numCat <- .Object@numCat[1:.Object@p]

  if (!all(.Object@numCat>=0))
    stop("numCat must be greater or equal zero.")

  return(.Object)
}
################################################################################

setAs(from="matrix",to="gRapHD",def=matrix.gRapHD)
setAs(from="gRapHD",to="graphNEL",def=gRapHD.graphNEL)
setAs(from="graphNEL",to="gRapHD",def=graphNEL.gRapHD)
setMethod("initialize",signature(.Object="gRapHD"),initialize.gRapHD)
setMethod("show",signature(object="gRapHD"),show.gRapHD)
setGeneric("plot",function(x,...)standardGeneric("plot"))
setMethod("plot",signature(x="gRapHD"),plot.gRapHD)
setMethod("print",signature(x="gRapHD"),print.gRapHD)
setMethod("summary",signature(object="gRapHD"),summary.gRapHD)
