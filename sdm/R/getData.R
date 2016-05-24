
.getFeatureRecords <- function(d,ind=NULL,n=NULL) {
  if (is.null(ind)) ind <- .getIndex(d)
  if (is.null(n)) n <- d@features.name
  else {
    if (!all(n %in% d@features.name)) stop('some feature names do not exist!')
  }
  n <- c('rID',n)
  #------
  w <- match(ind,d@features$rID)
  d@features[w,n]
}



.getSpeciesRecords <- function(d,ind=NULL,sp=NULL,grp=NULL,time=NULL) {
  if (is.null(ind)) ind <- .getIndex(d,sp=sp,groups = grp,time = time)
  else {
    ind <- ind[ind %in% .getIndex(d,sp=sp,groups = grp,time = time)]
  }
  if (is.null(sp)) sp <- d@species.names
  else {
    if (!all(sp %in% d@species.names)) stop('the species name(s) does not exist')
  }
  o <- data.frame(rID=ind)
  for (s in sp) {
    if (!is.null(d@species[[s]]@abundance)) {
      w <- match(ind,d@species[[s]]@abundance$rID)
      o[,s] <- d@species[[s]]@abundance$abundance[w]
    } else if (!is.null(d@species[[s]]@Multinomial)) {
      w <- match(ind,d@species[[s]]@Multinomial$rID)
      o[,s] <- d@species[[s]]@Multinomial$name[w]
    } else {
      if (!is.null(d@species[[s]]@presence)) {
        w <- which(ind %in% d@species[[s]]@presence)
        if (length(w) > 0) o[w,s] <- 1
      }
      
      if (!is.null(d@species[[s]]@absence)) {
        w <- which(ind %in% d@species[[s]]@absence)
        if (length(w) > 0) o[w,s] <- 0
      }
      
      if (!is.null(d@species[[s]]@background)) {
        w <- which(ind %in% d@species[[s]]@background)
        if (length(w) > 0) o[w,s] <- 0
      }
    }
  }
  o
}
#------
.getSdmDataFrame <- function(d,ind=NULL,sp=NULL,grp=NULL,time=NULL,preds=NULL) {
  o1 <- .getSpeciesRecords(d,ind,sp,grp,time)
  o2 <- .getFeatureRecords(d,o1$rID,preds)[,-1,drop=FALSE]
  cbind(o1,o2)
}
#--------
.getData.sdmMatrix <- function(formula,data,normalize=FALSE) {
  if (normalize) {
    sp <- as.character(formula[[2]])
    if (sp %in% colnames(data)) data <- .normalize(data,except=sp)
    else data <- .normalize(data)
  }
  if (length(formula) == 3) formula <- as.formula(paste('~',deparse(formula[[3]])))
  model.matrix(formula,data)[,-1,drop=FALSE]
}
.getData.sdmY <- function(formula,data) {
  data[,as.character(formula[[2]])]
}
#-------
#, row.names = NULL, optional = FALSE
setMethod('as.data.frame', signature(x='sdmdata'), 
          function(x, ...) {
            .getSdmDataFrame(x,...)
          }
)

