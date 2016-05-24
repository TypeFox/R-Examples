#'
#' Extracting important values for ecological Inference
#'
#' @description
#' Extracting important values for calculation of the ecological Inference with the
#'       \code{\link[eiwild]{runMBayes}}-function
#' 
#' @param form \code{formula} in this format 
#'          \code{cbind(column_1,column_2, ...,column_c)~cbind(row_1,row_2,...,row_r))}
#' @param aggr \code{data.frame} with aggregate data.
#'          One district per line and one column giving one ID per district. (see Details)
#' @param indi \code{data.frame} with individual data. 
#'          One district per line and one column giving one ID per district. (see Details)
#'          If no individual data are present it defaults to \code{NULL}
#' @param IDCols vector of length 2 giving the column-names or numbers of ID column
#' 
#' @details
#' \code{indi} is a \var{districts x [(r*c)+1]} \code{data.frame} containing one district per line. 
#' One column gives the ID of the districts which will be connection to the ID column in the \code{aggr}-data.frame.
#'
#' For example a 2x3 ecological Inference problem with \code{formula}
#' \code{cbind(col1,col2,col3) ~ cbind(row1,row2)}
#' will have the  row format :
#' \code{[ID, row1.col1, row1.col2, row1.col3, row2.col1, row2.col2, row2.col3]}
#' 
#' It is important that the \code{formula} names correspond to the exact column number in the \var{indi}-data.frame.
#' 
#' @return \code{list} with components needed for the Metropolis algorithm in \code{\link[eiwild]{runMBayes}}
#' 
#' @seealso
#' \code{\link[eiwild]{runMBayes}}, \code{\link[coda]{mcmc}}
#' \code{\link[eiwild]{tuneVars}}
#' 
#' @examples
#' \dontrun{
#' # loading some fake election data
#' data(topleveldat)
#' form <- cbind(CSU_2, SPD_2, LINK_2, GRUN_2) ~ cbind(CSU_1, SPD_1, Link_1)
#' 
#' conv <- convertEiData(form=form, aggr=aggr, indi=indi, IDCols=c("ID", "ID"))
#' 
#' str(conv)
#' }
#' 
#' @export


convertEiData <- function(form, aggr, indi=NULL, IDCols=c("ID")){
  
  aggr <- as.data.frame(aggr)
  
  if(IDCols[1] %in% colnames(aggr)){
    aggrIDs <- as.character(aggr[,IDCols[1]])
  } else {
    stop("name given in \"IDCols[1]\" doesn't exist in aggr-data.frame!", call.=FALSE)
  }
  
  modfr <- model.frame(form, aggr)
  rowdf <- modfr[[2]]
  coldf <- modfr[[1]]
  # alles geht über rownames
  rownames(coldf) <- rownames(rowdf) <- aggrIDs
  
  if(!all(rowSums(rowdf)==rowSums(coldf)))
    stop(paste0("aggr: row and column margins of districts ",
                paste(which(rowSums(rowdf)!=rowSums(coldf)),collapse=", "),
                " are not equal!"), call.=FALSE)
  
  c <- ncol(coldf)#
  r <- ncol(rowdf)#
  
  indiNames <- sprintf("%s.%s",
                       rep(colnames(rowdf),each=length(colnames(coldf))), #rep rows each c-times
                       rep(colnames(coldf),times=length(colnames(rowdf))))#duplicate cols r-times
  
  if(is.null(indi)){ #create fake individual data for the algorithm to work
    indi <- as.data.frame(matrix(0,1,1+r*c))
    indi[1,1] <- aggrIDs[1]
    names(indi) <- c("ID",indiNames)
    IDCols <- c(IDCols[1], "ID")
  } else{
    indi <- as.data.frame(indi)
  }
  
  if( IDCols[2] %in% colnames(indi) ){
    indiIDs <- as.character(indi[,IDCols[2]])
  } else{
    stop("name given in \"IDCols[2]\" doesn't exist in indi-data.frame!", call.=FALSE)
  }
  
  if(!( all(indiIDs %in% aggrIDs) ))
    stop("Missing individual data District IDs in aggregate data IDs!", call.=FALSE)
  
  indPrec <- nrow(indi)#
  aggPrec <- nrow(aggr)#
  N <- rowSums(coldf)
  names(N) <- aggrIDs
  
  if(!all(indiNames %in% colnames(indi)))
    stop("Colnames of indi-data have to be of type row1.col1, row1.col2,... !", .call=FALSE)
  inditmp <- indi[,indiNames]
  
  if(!( ncol(inditmp)==r*c ))
    stop("indi: individual data doesn't have \"r x c\" (=", r*c, ") data columns!!", call.=FALSE)
  
  ###### create M_r and Z_c
  Mdf <- matrix(NA,indPrec,r)
  rownames(Mdf) <- indiIDs
  colnames(Mdf) <- paste0("M_",colnames(rowdf))
  rowBegSeq <- seq(1,r*c, by=c) #sequence of row beginnings
  rowEndSeq <- rowBegSeq+(c-1) #sequence of row endings
  #    inditmp[,rowBegSeq]
  #    inditmp[,rowEndSeq]
  
  
  for(i in 1:r)
    Mdf[,i]<- rowSums(inditmp[, rowBegSeq[i]:rowEndSeq[i] ])
  if(nrow(Mdf)==1){
    M<-sum(Mdf[,1:r])
  } else{
    M <- rowSums(Mdf)
  }
  if(any(N[indiIDs]<M))
    stop("precinct with ID ",
         paste(indiIDs[N[indiIDs]<M],collapse=" "),
         " bigger individual than aggregate data!", call.=FALSE)
  if(any(rowdf[indiIDs,]<Mdf)){
    stop("precinct with ID ",
         paste(indiIDs[ceiling(which(t(rowdf[indiIDs,]<Mdf))/r)], collapse=" & "),
         " bigger individual than aggregate data in rows!", call.=FALSE)
  }
  #   Mdf
  #   M
  
  Zdf <- matrix(NA, indPrec, c)
  rownames(Zdf) <- indiIDs
  colnames(Zdf) <- paste0("Z_", colnames(coldf))
  for(i in 1:c)
    Zdf[,i] <- rowSums(inditmp[, rowBegSeq+ (i-1)]) # - 1 damit es bei null anfängt
  if(any(coldf[indiIDs,]<Zdf)){
    stop("precinct with ID ",
         paste(indiIDs[ceiling(which(t(coldf[indiIDs,]<Zdf))/c)], collapse=" & "),
         " bigger individual than aggregate data in columns!", call.=FALSE)
  }
  
  #coldf[indiIDs,][10,]<Zdf[10,]
  
  
  #   Zdf
  if(!( all(rowSums(Zdf)==M) ))
    stop("indi: row and colum margins are not equal", call.=FALSE)
  
  ###### create Zrc für agg
  Zrc <- array(0,c(r,c,aggPrec))
  dimnames(Zrc) <- list(colnames(rowdf), colnames(coldf), aggrIDs)
  for(i in 1:r) #jede reihe der aggbezirke wird mit den jeweiligen Individualdaten gefüllt. Rest bleibt 0
    Zrc[i,,indiIDs] <- t(inditmp[, rowBegSeq[i]:rowEndSeq[i] ])
  
  ###### create X_r
  Xdf <- matrix(NA, aggPrec,r)
  rownames(Xdf) <- aggrIDs
  colnames(Xdf) <- paste0("X_",colnames(rowdf))
  Xdf[,] <- rowdf
  Xdf[rownames(Mdf),] <- Xdf[rownames(Mdf),]-Mdf #N_r-M_r
  Nmod <- N
  Nmod[names(M)] <- Nmod[names(M)]-M # N - M
  Xdf <- Xdf/Nmod # (N_r-M_r)/(N-M)
  Xdf[which(is.nan(Xdf))]<-0
  #   Xdf
  
  ###### create T_c - Z_c
  Tdf <- matrix(NA, aggPrec,c)
  rownames(Tdf) <- aggrIDs
  colnames(Tdf) <- paste0("T_",colnames(coldf))
  Tdf[,] <- coldf
  Tdf[rownames(Zdf),] <- Tdf[rownames(Zdf),]-Zdf
  # Tdf
  
  ret <- list(formula=form,
              rowdf=rowdf, coldf=coldf,
              M=M, Mdf=Mdf, Zdf=Zdf, Zrc=Zrc, # indi components
              N=N, Xdf=Xdf, Tdf=Tdf) # aggr components
  
  return(ret) 
}

