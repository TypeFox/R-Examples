# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

Data <- function(object,
                 inclZeroWRES=FALSE,
                 onlyfirst=FALSE,
                 subset=NULL) {

  data <- object@Data
  if(!inclZeroWRES) {
    data <- data[data[,xvardef("wres",object)]!=0,]
  }

  if(onlyfirst) {
    data <- data[!duplicated(data[,xvardef("id",object)]),]
  }
  
  if(!is.null(subset)) {
    #attach(data)
    #on.exit(detach(data))
    
                                        # data <- data[eval(parse(text=subset)),]
    #browser()
    #data <- data[eval(parse(text=paste("data$", subset))),] # fix subsets 22/3/06
    data<-with(data,data[eval(parse(text=subset)),])
    

    if(dim(data)[1]==0) return(NULL)
  }

  return(data)

}

"Data<-" <- function(object,quiet=TRUE,keep.structure=F,value) {

  if(is.null(value)) {
    return(NULL)
  }

  object@Data     <- value
  totab.nam  <- names(value)
                                        # cat(totab.nam, "\n")

  ## define the important variables in the dataset
  singdef <-
    c("id","idlab","idv","dv","pred","ipred","iwres","wres","res","cwres")

  ## Check if all variable definitions exist in the data
  for(v in singdef) {
    if(is.null(xvardef(v,object))){
      next
    }
    if(is.na(pmatch(xvardef(v,object), totab.nam))) {
      object@Prefs@Xvardef[[v]] <- NULL

      ## if DV undefined check what 4th from last in list files says
      ## must make sure no append is not called
      if(v=="dv"){
        tab.names.list <-
          c(".sdtab.names.tmp",".patab.names.tmp",".catab.names.tmp",".cotab.names.tmp")
        append.list <- c("PRED","RES","WRES")
        dv.found <- FALSE
        for(vv in tab.names.list){
          if(!dv.found){
            if (file.exists(vv)){
              tmp.list <- scan(vv, what="list", quiet=T)
              tmp.length <- length(tmp.list)
              if (all(tmp.list[(tmp.length-2):tmp.length]==append.list)){
                object@Prefs@Xvardef[[v]] <- tmp.list[[tmp.length-3]]
                dv.found <- TRUE
              }
            }
          }
        }
      }

      ## IF TIME isn't in tables is TAD?
      if(v=="idv"){
        if(!is.na(pmatch("TAD",totab.nam))){
          object@Prefs@Xvardef[[v]] <- "TAD"
        } else {
          if(!is.na(pmatch("CP",totab.nam))){
            object@Prefs@Xvardef[[v]] <- "CP"
          }
        }
      }
    }
  }

  
  ## sort out covariates and parameters
  ## look at the difference between the names in singdef
  ## and the names in *.tmp files
  if (file.exists(".patab.names.tmp")) {
    xparms <- setdiff(scan(".patab.names.tmp", what="list", quiet=T), sapply(singdef, xvardef, object))
    ##unlink("patab.names.txt")
  } else {
    xparms <- NULL
  }
  
  if (file.exists(".catab.names.tmp")) {
    xcatcovs <- setdiff(scan(".catab.names.tmp", what="list", quiet=T), sapply(singdef, xvardef, object))
    ##unlink("catab.names.txt")
  } else {
    xcatcovs <- NULL
  }
  
  if (file.exists(".cotab.names.tmp")) {
    xconcovs <- setdiff(scan(".cotab.names.tmp", what="list", quiet=T), sapply(singdef, xvardef, object))
    ##unlink("cotab.names.txt")
  } else {
    xconcovs <- NULL
  }
  
  if (!is.null(xcatcovs)) {
    if (!is.null(xconcovs)) {
      xcovs <- c(xcatcovs, xconcovs)
    } else {
      xcovs <- xcatcovs
    }
  } else {
    if (!is.null(xconcovs)) {
      xcovs <- xconcovs
    } else {
      xcovs <- NULL
    }    
  }

  ## check for categorical variables
  if(!keep.structure){
    for(v in totab.nam) {
      if (!is.na(match("TIME",v))) next
      if (!is.na(match("TAD",v))) next
      if (!is.na(match("CP",v))) next
      if (!is.na(match("PRED",v))) next
      if (!is.na(match("IPRE",v))) next
      if (!is.na(match("IPRED",v))) next
      num <- length(unique(object@Data[,v]))
      if(is.null(xvardef("dv",object))){
        if(num <= object@Prefs@Cat.levels) {
          if(!is.factor(object@Data[,v])){
            if(!quiet){
              cat("\n  Inferring that ",v," is categorical (",num," levels).\n",sep="")
              cat("  Transforming",v,"from continuous to categorical\n",sep=" ")
            }
            object@Data[,v] <- as.factor(object@Data[,v])
          }
        } else {
          if(is.factor(object@Data[,v])){
            if(!quiet){
              cat("\n  Inferring that",v,"is continuous.\n",sep=" ")
              cat("  Transforming",v,"from categorical to continuous\n",sep=" ")
            }
            object@Data[,v] <- as.numeric(levels(object@Data[,v]))[object@Data[,v]]
          }
        }
      } else {
        if(is.na(match(xvardef("dv",object), v))) {
          if(num <= object@Prefs@Cat.levels) {
            if(!is.factor(object@Data[,v])){
              if(!quiet){
                cat("\n  Inferring that ",v," is categorical (",num," levels).\n",sep="")
                cat("  Transforming",v,"from continuous to categorical\n",sep=" ")
              }
              object@Data[,v] <- as.factor(object@Data[,v])
            }
          } else {
            if(is.factor(object@Data[,v])){
              if(!quiet){
                cat("\n  Inferring that",v,"is continuous.\n",sep=" ")
                cat("  Transforming",v,"from categorical to continuous\n",sep=" ")
              }
              object@Data[,v] <- as.numeric(levels(object@Data[,v]))[object@Data[,v]]
            }
          }
        } else {
          if(num <= object@Prefs@DV.Cat.levels) {
            if(!is.factor(object@Data[,v])){
              if(!quiet){
                cat("\n  Inferring that DV is categorical (",num," levels).\n",sep="")
                cat("  Transforming DV from continuous to categorical\n")
              }
              object@Data[,v] <- as.factor(object@Data[,v])
            }
          } else {
            if(is.factor(object@Data[,v])){
              if(!quiet){
                cat("\n  Inferring that DV is continuous.\n",sep="")
                cat("  Transforming DV from categorical to continuous\n")
              }
              object@Data[,v] <- as.numeric(levels(object@Data[,v]))[object@Data[,v]]
            }
          }
        }
      }
    }
    if(is.factor(object@Data[,xvardef("dv",object)])){
      for (v in totab.nam){
        num <- length(unique(object@Data[,v]))
        if (!is.na(match("PRED",v)) | !is.na(match("IPRE",v)) | !is.na(match("IPRED",v))) {
          if(!quiet){
            cat("\n  Inferring that ",v," is categorical (",num," levels).\n",sep="")
            cat("  Transforming",v,"from continuous to categorical\n",sep=" ")
          }
          object@Data[,v] <- as.factor(object@Data[,v])
        }
      }
    } else {
      for (v in totab.nam){
        if (!is.na(match("PRED",v)) | !is.na(match("IPRE",v)) | !is.na(match("IPRED",v))) {
          if(is.factor(object@Data[,v])){
            if(!quiet){
              cat("\n  Inferring that",v,"is continuous.\n",sep=" ")
              cat("  Transforming",v,"from categorical to continuous\n",sep=" ")
            }
            object@Data[,v] <- as.numeric(levels(object@Data[,v]))[object@Data[,v]]
          }
        }
      }
    }
  }
  
  for(v in xcatcovs) {
    object@Data[,v] <- as.factor(object@Data[,v])
  }

  
  if (is.null(xparms)) {
    ep <- c()
    for(v in xvardef("parms",object)) {
      if(!is.na(pmatch(v, totab.nam))) {
        ep <- c(v,ep)
      }
    }
    object@Prefs@Xvardef[["parms"]] <- ep
  } else {
    object@Prefs@Xvardef[["parms"]] <- xparms
  }

  ep <- c()
  for(v in xvardef("tvparms",object)) {
    if(!is.na(pmatch(v, totab.nam))) {
      ep <- c(v,ep)
    }
  }
  object@Prefs@Xvardef[["tvparms"]] <- ep

  ep <- c()
  ranpar.loc <- grep("^ETA",totab.nam)
  if (length(ranpar.loc)!=0){
    ep <- totab.nam[ranpar.loc]
  }
  ##   for(v in xvardef("ranpar",object)) {
  ##     if(!is.na(pmatch(v, totab.nam))) {
  ##       ep <- c(v,ep)
  ##     }
  ##   }
  if ((is.vector(ep)) || (is.vector(ep))) {
    object@Prefs@Xvardef[["ranpar"]] <- sort(ep)
  } else {
    object@Prefs@Xvardef[["ranpar"]] <- ep
  }

  
  ##   ep <- c()
  ##   for(v in xvardef("cat.cov",object)) {
  ##     if(!is.na(pmatch(v, totab.nam))) {
  ##       ep <- c(v,ep)
  ##       object@Data[,v] <- as.factor(object@Data[,v])
  ##     }
  ##   }
  ##   object@Prefs@Xvardef[["cat.cov"]] <- ep
  
  if (is.null(xcovs)) {
    ep <- c()
    for(v in xvardef("covariates",object)) {
      if(!is.na(pmatch(v, totab.nam))) {
        ep <- c(v,ep)
      }
    }
    object@Prefs@Xvardef[["covariates"]] <- ep
  } else {
    object@Prefs@Xvardef[["covariates"]] <- xcovs
  }

  ## Fix the labels
  for(v in totab.nam) {
    if(is.null(xlabel(v,object))) {
      object@Prefs@Labels[[v]] <- v
    }
  }

  return(object)
}
