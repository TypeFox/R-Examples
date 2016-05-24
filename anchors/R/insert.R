#######################################################################
##
## Function: anchors()
## Author  : Olivia Lau and Jonathan Wand <wand(at)stanford.edu>
## Created : 2006-12-01 (OL)
##  
## DESCRIPTION: Function for inserting rank values in anchors object back to a dataset
##
## OUTPUT:   data.frame
##
## INPUT:
##   data  : dataframe
##   object: anchors.rank object
##   vnames: (option) column names for B and C in new dataset
##           format: list(B = c(Bs,Be), C=c(Cs,Ce))
##   overwrite : replace existing columns with same names in vnames?
##
## MODIFIED:
##    2007-01-19 : JW
##    - uses structure of 'anchors' output
##    - corrects error checking on column names
##    - ensures that data is a dataframe
##    - does matching both ways on rownames C->data, data->C
##
##    2007-09-01 : JW
##    - added usage for B
## 
##    2008-04-20 : JW
##    - allow single call for both B and C inserts
##
#######################################################################
insert <- function(data,
                   obj, 
                   vnames, overwrite = FALSE,
                   debug=0) {

  ## only anchors objs make sense..
  if (!(class(obj) %in% c("anchors.rank","anchors.rank.type")))
    stop("Second argument, obj, must be of class 'anchors.rank' or 'anchors.rank.type' \n")
  
  if (class(obj) == "anchors.rank" && class(obj$rank) == "anchors.rank.type")
    obj <- obj$rank
  
  ## add or modify format of names
  if (missing(vnames)) 
    vnames <- colnames(obj$span)

  data <- insert.anchors.rank.type( data, obj$span, overwrite, vnames , debug=debug)

  return(data)
}

insert.anchors.rank.type <- function(data,obj,overwrite,vnames,debug=0) {
  ## Error checks
  if (is.null(vnames) || length(vnames) != ncol(obj))
    stop(paste("\n",  "Error in 'vnames'.  Cannot use:",  paste( vnames, collapse = ","),  "\n"))

  idx <- vnames %in% names(data)
  if (any(idx) & !overwrite)
    stop(paste("\n",paste(vnames[idx], collapse = " and "), "matches column name(s) in data.  Specify different 'vnames' or choose 'overwrite = TRUE'"))

  ## data must be a list or dataframe
  if ( !(typeof(data) == c("list") ) )
    data <- as.data.frame(data)

  ## create NA default
  data[[vnames[1]]] <- data[[vnames[2]]] <- rep(NA,NROW(data))

  ## does matching both ways on rownames C->data, data->C
  didx <- rownames(data) %in% rownames(obj) 
  cidx <- rownames(obj) %in% rownames(data)

  if (debug > 0) {
    print(rownames(data))
    print(rownames(obj))
#    print(didx)
#    print(cidx)
  }
  
  if (sum(didx) != sum(cidx)) stop("insert: mismatch in rowname matches")
  
  data[[vnames[1]]][didx] <- obj[cidx,1]
  data[[vnames[2]]][didx] <- obj[cidx,2]
  
  return(data)
}
