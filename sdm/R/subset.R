# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2015
# Version 1.0
# Licence GPL v3

.deleteRecords <- function(x,id) {
  species <- .getSpeciesRecords(x,id)
  species <- colnames(species)[2:ncol(species)]
  for (sp in species) {
    if (x@species[[sp]]@type == 'Presence-Absence') {
      w <- which(x@species[[sp]]@presence %in% id)
      if (length(w) > 0) x@species[[sp]]@presence <- x@species[[sp]]@presence[-w]
      if (length(x@species[[sp]]@presence) == 0) {
        x@species[[sp]]@type <- 'Absence-Only!'
        x@species[[sp]]@presence <- NULL
      }
      #-
      w <- which(x@species[[sp]]@absence %in% id)
      if (length(w) > 0) x@species[[sp]]@absence <- x@species[[sp]]@absence[-w]
      if (length(x@species[[sp]]@absence) == 0) {
        x@species[[sp]]@type <- 'Presence-Only'
        x@species[[sp]]@absence <- NULL
      } 
    } else if (x@species[[sp]]@type == 'Presence-Only') {
      w <- which(x@species[[sp]]@presence %in% id)
      if (length(w) > 0) x@species[[sp]]@presence <- x@species[[sp]]@presence[-w]
      if (length(x@species[[sp]]@presence) == 0) {
        x@species[[sp]]@type <- ''
        x@species[[sp]]@presence <- NULL
        if (is.null(x@species[[sp]]@background)) warning('no records remained!')
      }
    } else if (x@species[[sp]]@type %in% c('Abundance','Abundance-constant!')) {
      w <- which(x@species[[sp]]@abundance$rID %in% id)
      if (length(w) > 0) x@species[[sp]]@abundance <- x@species[[sp]]@abundance[-w,]
      if (nrow(x@species[[sp]]@abundance) == 0) {
        x@species[[sp]]@type <- ''
        x@species[[sp]]@abundance <- NULL
        warning('no records remained!')
      }
    } else if (x@species[[sp]]@type == 'Absence-Only!') {
      w <- which(x@species[[sp]]@absence %in% id)
      if (length(w) > 0) x@species[[sp]]@absence <- x@species[[sp]]@absence[-w]
      if (length(x@species[[sp]]@absence) == 0) {
        x@species[[sp]]@type <- ''
        x@species[[sp]]@absence <- NULL
        if (is.null(x@species[[sp]]@background)) warning('no records remained!')
      }
    } else if (x@species[[sp]]@type == 'Multinomial') {
      w <- which(x@species[[sp]]@Multinomial$rID %in% id)
      if (length(w) > 0) x@species[[sp]]@Multinomial <- x@species[[sp]]@Multinomial[-w,]
      if (nrow(x@species[[sp]]@Multinomial) == 0) {
        x@species[[sp]]@type <- ''
        x@species[[sp]]@Multinomial <- NULL
        warning('no records remained!')
      }
    }
    
    if (!is.null(x@species[[sp]]@background)) {
      w <- which(x@species[[sp]]@background %in% id)
      if (length(w) > 0) x@species[[sp]]@background <- x@species[[sp]]@background[-w]
      if (length(x@species[[sp]]@background) == 0) {
        x@species[[sp]]@background <- NULL
      }
    }
  }
  #----
  w <- which(x@features$rID %in% id)
  if (length(w) > 0) x@features <- x@features[-w,]
  if (nrow(x@features) == 0) {
    x@features <- NULL
    x@features.name <- NULL
    x@factors <- NULL
  }
  #-----
  if (!is.null(x@info)) {
    if (!is.null(x@info@info)) {
      w <- which(x@info@info$rID %in% id)
      if (length(w) > 0) x@info@info <- x@info@info[-w,]
      if (nrow(x@info@info) == 0) x@info@info <- NULL
    }
    if (!is.null(x@info@time)) {
      w <- which(x@info@time$rID %in% id)
      if (length(w) > 0) x@info@time <- x@info@time[-w,]
      if (nrow(x@info@time) == 0) x@info@time <- NULL
    }
    if (!is.null(x@info@coords)) {
      w <- which(x@info@coords[,1] %in% id)
      if (length(w) > 0) x@info@coords <- x@info@coords[-w,]
      if (nrow(x@info@coords) == 0) x@info@coords <- NULL
    }
  }
  #--------
  grp <- .getGroupNames(x)
  for (g in grp) {
    for (l in as.character(x@groups[[g]]@values[,2])) {
      w <- which(x@groups[[g]]@indices[[l]] %in% id)
      if (length(w) > 0) x@groups[[g]]@indices[[l]] <- x@groups[[g]]@indices[[l]][-w]
      if (length(x@groups[[g]]@indices[[l]]) == 0) {
        x@groups[[g]]@indices <- x@groups[[g]]@indices[-which(names(x@groups[[g]]@indices) == l)]
        x@groups[[g]]@values <- x@groups[[g]]@values[x@groups[[g]]@values$values != l,]
        if (nrow(x@groups[[g]]@values) == 0) x@groups <- x@groups[-which(names(x@groups) == g)]
      }
    }
  }
  x
}
#------

.subsetRecords <- function(x,id) {
  # x is sdmdata
  species <- .getSpeciesRecords(x,id)
  species <- colnames(species)[2:ncol(species)]
  for (sp in species) {
    if (x@species[[sp]]@type == 'Presence-Absence') {
      w <- which(x@species[[sp]]@presence %in% id)
      if (length(w) > 0) x@species[[sp]]@presence <- x@species[[sp]]@presence[w]
      else x@species[[sp]]@presence <- NULL
      
      w <- which(x@species[[sp]]@absence %in% id)
      if (length(w) > 0) x@species[[sp]]@absence <- x@species[[sp]]@absence[w]
      else x@species[[sp]]@absence <- NULL
      
    } else if (x@species[[sp]]@type == 'Presence-Only') {
      w <- which(x@species[[sp]]@presence %in% id)
      if (length(w) > 0) x@species[[sp]]@presence <- x@species[[sp]]@presence[w]
      else x@species[[sp]]@presence <- NULL
      
    } else if (x@species[[sp]]@type %in% c('Abundance','Abundance-constant!')) {
      w <- which(x@species[[sp]]@abundance$rID %in% id)
      if (length(w) > 0) x@species[[sp]]@abundance <- x@species[[sp]]@abundance[w,]
    } else if (x@species[[sp]]@type == 'Absence-Only!') {
      w <- which(x@species[[sp]]@absence %in% id)
      if (length(w) > 0) x@species[[sp]]@absence <- x@species[[sp]]@presence[w]
      else x@species[[sp]]@absence <- NULL
    } else if (x@species[[sp]]@type == 'Multinomial') {
      w <- which(x@species[[sp]]@Multinomial$rID %in% id)
      if (length(w) > 0) x@species[[sp]]@Multinomial <- x@species[[sp]]@Multinomial[w,]
      else x@species[[sp]]@Multinomial <- NULL
    }
    
    if (!is.null(x@species[[sp]]@background)) {
      w <- which(x@species[[sp]]@background %in% id)
      if (length(w) > 0) x@species[[sp]]@background <- x@species[[sp]]@background[w]
      else x@species[[sp]]@background <- NULL
    }
  }
  #----
  w <- which(x@features$rID %in% id)
  if (length(w) > 0) x@features <- x@features[w,]
  else x@features <- NULL
  #-----
  if (!is.null(x@info)) {
    if (!is.null(x@info@info)) {
      w <- which(x@info@info$rID %in% id)
      if (length(w) > 0) x@info@info <- x@info@info[w,]
      else x@info@info <- NULL
    }
    if (!is.null(x@info@time)) {
      w <- which(x@info@time$rID %in% id)
      if (length(w) > 0) x@info@time <- x@info@time[w,]
      else x@info@time <- NULL
    }
    if (!is.null(x@info@coords)) {
      w <- which(x@info@coords[,1] %in% id)
      if (length(w) > 0) x@info@coords <- x@info@coords[w,]
      else x@info@coords <- NULL
    }
  }
  #--------
  grp <- .getGroupNames(x)
  for (g in grp) {
    for (l in as.character(x@groups[[g]]@values[,2])) {
      w <- which(x@groups[[g]]@indices[[l]] %in% id)
      if (length(w) > 0) x@groups[[g]]@indices[[l]] <- x@groups[[g]]@indices[[l]][w]
      else x@groups[[g]]@indices[[l]] <- c()
      if (length(x@groups[[g]]@indices[[l]]) == 0) {
        x@groups[[g]]@values <- x@groups[[g]]@values[x@groups[[g]]@values$values != l,]
      }
    }
    if (nrow(x@groups[[g]]@values) == 0) x@groups <- x@groups[-which(names(x@groups) == g)]
  }
  x
}
#---------
.getDataFame <- function(x,id=NULL,sp=NULL,grp=NULL,time=NULL) {
  if (is.null(id)) id <- .getIndex(x,sp=sp,groups=grp,time=time)
  f <- NULL
  if (length(id) > 0) {
    f <- .getSdmDataFrame(x,ind = id)
    if (!is.null(x@info@coords)) f <- cbind(f,data.frame(x@info@coords[x@info@coords[,1] %in% id,2:3]))
    if (!is.null(x@info@time)) f <- cbind(f,x@info@time[x@info@time[,1] %in% id,2:ncol(x@info@time)])
    if (!is.null(x@info@info)) f <- cbind(f,x@info@info[x@info@info[,1] %in% id,2:ncol(x@info@info)])
  }
  f
}


#-------
setMethod("[", c("sdmdata"),
          function(x, i, ...,drop=FALSE) {
            if (missing(i)) stop('i is missing!')
            if(drop) .getDataFame(x,i)
            else .subsetRecords(x,i)
          }
)
#--------

setMethod("[[", c("sdmModels","ANY","ANY"),
          function(x,i,drop=TRUE,...) {
            if ( missing(i)) { stop('you must provide an index') }
            mi <- .getModel.info(x,w=i)
            if (nrow(mi) == 0) stop('the specified index in i are not within the range of model IDs')
            
            mi$newRID <- mi$replicationID
            if (drop) mi$newID <- 1:nrow(mi)
            
            
            sp <- as.character(unique(mi[,2]))
            smo <- new('sdmModels',setting=x@setting,data=x@data)
            
            
            
            for (s in sp) {
              mj <- mi[which(mi[,2] == s),]
              m <- as.character(unique(mj[,3]))
              r <- unique(mj[,5])
              for (ri in seq_along(r)) mi[mi[,2] == s & mi[,5] == r[ri],10] <- ri
              smo@replicates[[s]] <- x@replicates[[s]][r]
              smo@models[[s]] <- list()
              
              for (mo in m) {
                smo@models[[s]][[mo]] <- list()
                mk <- mj[which(mj[,3]  == mo),]
                if (drop) {
                  for (id in mk[,11]) {
                    smo@models[[s]][[mo]][[as.character(id)]] <- x@models[[s]][[mo]][[as.character(id)]]
                  }
                } else {
                  for (id in mk[,1]) {
                    smo@models[[s]][[mo]][[as.character(id)]] <- x@models[[s]][[mo]][[as.character(id)]]
                  }
                }
              }
            }
            
            mi[,5] <- mi[,10]
            if (drop) {
              for (j in 1:nrow(mi)) {
                x@models[[mi[j,2]]][[mi[j,3]]][[as.character(mi[j,1])]]@mID <- mi[j,11]
              }
              mi[,1] <- mi[,11]
            }
            smo@run.info <- mi[,1:9]
            smo
          }
)



if (!isGeneric("subset")) {
  setGeneric("subset", function(x, ...)
    standardGeneric("subset"))
}

setMethod("subset","sdmModels",
          function(x, subset, drop=TRUE, ...) {
            if (missing(subset)) { stop('you must provide an index') }
            x[[subset,drop=drop,...]]
          })

