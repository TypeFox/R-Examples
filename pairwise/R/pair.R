#' @title Rasch Item Parameter (Main Function) 
#' @export pair
#' @description This is the (new) main function for calculation of the item parameter for the dichotomous Rasch Model (Rasch, 1960) and its extension for polytomous items (thurstonian thresholds) according to the Partial Credit Model (Masters, 1982), using a generalization of the pairwise comparison algorithm (Choppin, 1968, 1985; Wright & Masters, 1982). The number of (response) categories may vary accross items. 
#' Missing values up to an high amount in data are allowed, as long as items are proper linked together.
#' 
#'@details Parameter calculation is based on the construction of a paired comparison matrix M\emph{nicjc} with entries f\emph{icjc} representing the number of respondents who answered to item \emph{i} in category \emph{c} and to item \emph{j} in category \emph{c-1} widening Choppin's (1968, 1985) conditional pairwise algorithm to polytomous item response formats. This algorithm is simply realized by matrix multiplication.
#'
#' To avoid numerical problems with off diagonal zero's when constructing the pairwise comparison matrix M\emph{nij}, powers of the M\emph{nicjc} matrix, can be used (Choppin, 1968, 1985). Using powers \emph{k} of M\emph{nicjc} - argument \code{pot=TRUE} (default), replaces the results of the direct comparisons between \emph{i} and \emph{j} with the sum of the indirect comparisons of \emph{i} and \emph{j} through an intermediate \emph{k}.
#' 
#'In general, it is recommended to use the argument with default value \code{pot=TRUE}.
#'
#'For a graphic representation of the item 'estimates' the plotting S3 method \code{\link{plot.pair}} is available. For plotting the item category probabilities the function \code{\link{catprob}} can be used.  
#'
#'    
#' @param daten a \code{data.frame} or \code{matrix} with optionaly named colums (names of items), potentially with missing values, comprising polytomous or dichotomous (or mixted category numbers) responses of \code{n} respondents (rows) on \code{k} items (colums) coded starting with 0 for lowest category to \emph{m}-1 for highest category, with \emph{m} beeing a vector (with length k) with the number of categories for the respective item.
#' @param m an integer (will be recycled to a vector of length k) or a vector giving the number of response categories for all items - by default \code{(m = NULL)}, \code{m} is calculated from data, assuming that every response category is at least once present in data. For \emph{'sparse' data} it is \emph{strongly recomended} to explicitly \emph{define the number of categories} by defining this argument.
#' @param pot either a logical or an integer  >= 2 defining the power to compute of the pairwise comparison matrix. If TRUE (default) a power of three of the pairwise comparison matrix is used for further calculations. If FALSE no powers are computed.
#' @param zerocor logical, if TRUE (default) unobserved combinations (1-0, 0-1) in data for each pair of items are given a frequency of one conf. proposal by Alexandrowicz (2011, p.373).
#' @param ccf logical with default \code{ccf=FALSE} to perform normal item parameter calculation, if set to \code{ccf=TRUE} just the conditional item (category) frequencies are returned.
#' @param ... additional parameters passed through.
#' @return A (list) object of class \code{"pair"} containing the item category thresholds and difficulties sigma, also called item location.
#' @exportClass pair
#' @references Alexandrowicz, R. W. (2011). 'GANZ RASCH': A Free Software for Categorical Data Analysis. \emph{Social Science Computer Review, 30}(3), 369-379.
#' @references Choppin, B. (1968). Item Bank using Samplefree Calibration. \emph{Nature, 219}(5156), 870-872.
#' @references Choppin, B. (1985). A fully conditional estimation procedure for Rasch model parameters. \emph{Evaluation in Education, 9}(1), 29-42.
#' @references Masters, G. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika, 47}(2), 149–174.
#' @references Rasch, G. (1960). \emph{Probabilistic models for some intelligence and attainment tests.} Copenhagen: Danmarks pædagogiske Institut.
#' @references Wright, B. D., & Masters, G. N. (1982). \emph{Rating Scale Analysis.} Chicago: MESA Press.
#' 
#' @examples data(bfiN) # loading example data set
#' # calculating itemparameters for 5 neuroticism items with 6 answer categories (0-5).
#' neuro_itempar<-pair(daten = bfiN, m = 6) 
#' summary(neuro_itempar) 
#' summary(neuro_itempar, sortdif=TRUE) # ordered by difficulty 
#' # plotting threshold profiles for 5 neuroticism items.
#' plot(neuro_itempar) 
#' plot(neuro_itempar, sortdif=TRUE) # plotting ordered by difficulty 
#' ################ with unequal number of categories 
#' data(sim200x3)
#' res<-pair(sim200x3)
#' summary(res)
#' plot(res)

pair <- function(daten, m=NULL, pot=TRUE, zerocor=TRUE, ccf=FALSE, ...){
  
  fuargs <- list(daten=daten, m=m, pot=pot, zerocor=zerocor, ccf=ccf, ...=...)
  #### some internal functions ------
  dataprep1<-function(X){
  X<-as.matrix(X)
  if(length(colnames(X))==0){
    Iname<-nchar(paste(dim(X)[2]))
    colnames(X)<-paste("I",formatC(1:dim(X)[2], width = Iname, format = "d", flag = "0"),sep="") 
    cat("no item names found in data" ,"\n", "items are named", colnames(X)[1], "(first item) to",  colnames(X)[dim(X)[2]],"(last item)","\n")
  }
  if(length(rownames(X))==0){
    Pname<-nchar(paste(dim(X)[1]))
    rownames(X)<-paste("P",formatC(1:dim(X)[1], width = Pname, format = "d", flag = "0"),sep="")  
    cat("no person names (IDs) found in data" ,"\n", "persons are named", rownames(X)[1], "(first row) to",  rownames(X)[dim(X)[1]],"(last row)","\n")
  }
  return(X)
  }
  
  mat.mult<-function(A,B,na.rm=TRUE){
  # new version based on remarks by A. Robitzsch
  A<-as.matrix(A);B<-as.matrix(B)
  if(dim(A)[2]!=dim(B)[1]){stop("not conformable Matrices A B")}
  A[ is.na(A) ] <- 0
  B[ is.na(B) ] <- 0
  erg <- A %*% B
  return(erg)
  }
  
  katdat<-function(X,mVector, INFO=FALSE){
  #if(length(mVector)==0){mVector<-mV(X)}
  nitem<-dim(X)[2]
  XList<-as.list(data.frame((X)))
  foo<-function(q,w){ lapply(q,function(x){ (x==w)*1 })    }
  katdatList<-mapply(foo, (lapply(as.list(mVector-1),function(x){c(0:x[1])})), XList, SIMPLIFY = F )
  katdat<-do.call(cbind, lapply(katdatList,function(x){(do.call(cbind,x))}))
  colnames(katdat)<-unlist(mapply(paste, mapply(rep,colnames(X),mVector,SIMPLIFY = F), ( lapply(as.list(mVector-1),function(x){c(0:x[1])})  ) , MoreArgs=list(sep = "."), SIMPLIFY = F))
  rownames(katdat)<-rownames(X)
  if(INFO==FALSE){return(katdat)}
  if(INFO==TRUE){
    erg<-list(katdat=katdat, mVector=mVector)
    class(erg)<-"katdat"
    return(erg)
  }
  }
    
  Fmat<-function(X,INFO=FALSE,ipsative=TRUE){
    # der normale Fall X ist ein objekt der Klasse "katdat"
  stopifnot(class(X)=="katdat")
  mVector<-X$mVector; X<-X$katdat
  condf <- mat.mult(t(X),(X))
  spaltenraus <- cumsum(mVector)
  zeilenraus <- cumsum(mVector)-((mVector)-1)
  fmat <- condf[-zeilenraus,-spaltenraus]
  ### setzten der ipsativen item vergleiche auf 0
  if(ipsative==FALSE){
    bisnull <- cumsum((mVector)-1)
    vonnull <- (cumsum((mVector)-1) ) - ((mVector)-2)
    y<-do.call(c,mapply(function(x,y){rep(x,each=y) },mapply(seq,vonnull,bisnull,SIMPLIFY = FALSE),mVector-1,SIMPLIFY = FALSE))
    x<-do.call(c,mapply(function(x,y){rep(x,times=y) },mapply(seq,vonnull,bisnull,SIMPLIFY = FALSE),mVector-1,SIMPLIFY = FALSE))
    for (i in 1:length(y)){fmat[(y[i]),(x[i])]<-0}
  }
  ####
  if(INFO==FALSE){return(fmat)}
  if(INFO==TRUE){
    erg<-list(fmat=fmat, mVector=mVector)
    class(erg)<-"Fmat"
  return(erg)
  }
  }
  
#   cube<-function(M){
#   Mp<-M%*%M%*%M
#   return(Mp)
#   }
  
  cube<-function(M, POTT=3){
    Mp <- M
    for (i in 1:(POTT-1)){
      Mp<-Mp%*%M
    }
    return(Mp)
  }
    
  Dmat<-function(X){
  dmat <- t(X) / X
  return (dmat) 
  }
#### compute arbitrary powers of matrix 
pow <- function (x, k) 
{  n <- nrow(x)
   # return(pow(solve(x), -k))
   x.k <- x
   i <- 2
   while (i <= k) {
     x.k <- x.k %*% x
     i <- i + 1
   }
   return(x.k)
}
  #### END internal functions ------
  
  #### start data preparation ------
  d<-dataprep1(daten)
  
  ### some category checks for m ----
  if(length(m)==0){if (all(apply(d,2,function(x){min(x,na.rm=TRUE)})== 0) == FALSE){stop("item categories must start with 0") }}
  if (length(m)==0){m<-apply(d,2,function(x){max(x,na.rm=TRUE)+1}); comment(m) <- "estimated from data"}
  if (length(m)==1){m<-rep(m,dim(d)[2]); comment(m) <- "user defined (recycled)"}else{comment(m) <- "user defined"} 
  if (any (m < apply(d,2,function(x){max(x,na.rm=TRUE)+1}))){stop("some items in data have more categories than defined in m","\n","max item categories in data are: ",paste(apply(d,2,function(x){max(x,na.rm=TRUE)+1}),collapse=" , "),"\n", "but m was defined: ",paste(m,collapse=" , ")  )}

  #### start itemparameter calculation ------
  if (ccf==FALSE){
    
  f <- Fmat(katdat(d, mVector=m, INFO=T),INFO=T,ipsative=FALSE)
  mVector<-f$mVector
  # powers of nij ?
  if(class(pot)=="logical"){
  if(pot==FALSE){nij <- f$fmat}
  if(pot==TRUE){nij <- cube(f$fmat) }
  }
  if(class(pot)=="numeric"){
    if(pot<2){stop("pot must be >= than 2")}
    else{
      nij <- pow(f$fmat,pot)
    }
  }
   
  # new zero correction as an option --> conf. R. Alexandrowicz(2011) p.373 
  if(zerocor==TRUE){nij[nij==0]<-1}
  
  # calculation thresholds / dificulties
  tau<-rowMeans(log(Dmat(nij)), na.rm = TRUE)
  # formatieren der ausgabe
  vo <- (cumsum((mVector)-1) ) - ((mVector)-2)
  bi <- cumsum((mVector)-1)
  er1<-mapply(function(v,b){ tau[v:b]},vo,bi,SIMPLIFY = FALSE)
  names(er1) <- colnames(d)

### ungleiche oder gleiche kategoriezahl -- sb als liste, thresholds ggf mit NA als matrix
  sigma <- sapply(er1, mean)
  threshold <- er1
  ### berechnung der sb; sb als liste!!!!
  sb<-(lapply(threshold,function(x){cumsum(x)})) # richtig: umrechnung threshold (tau) in sb !!!!!!!!!
  if(length(unique(mVector))!=1){
    maxLen <- max(sapply(threshold, length))
    # create a new list with elements padded out with NAs
    newthreshold <- lapply(threshold, function(.ele){c(.ele, rep(NA, maxLen))[1:maxLen]})
    threshold <- do.call(rbind, newthreshold)
    colnames(threshold) <- c(1:dim(threshold)[2])
  }
  if(length(unique(mVector))==1){
    threshold<-do.call(rbind, threshold)
    colnames(threshold) <- c(1:(unique(mVector)-1))
    # rownames(threshold) <- colnames(d)
  }

  erg<-list(threshold=threshold,sigma=sigma,sb=sb,resp=d, fuargs=fuargs, m=m)
  class(erg) <- c("pair","list")
  
  #invisible(erg)
   return(erg)

  }

  if (ccf==TRUE){
  
  f <- Fmat(katdat(d, mVector=m, INFO=T),INFO=T,ipsative=FALSE)
  mVector<-f$mVector
  # powers of nij ?
  if(class(pot)=="logical"){
    if(pot==FALSE){nij <- f$fmat}
    if(pot==TRUE){nij <- cube(f$fmat) }
  }
  if(class(pot)=="numeric"){
    if(pot<2){stop("pot must be >= than 2")}
    else{
      nij <- cube(f$fmat, POTT=pot)
    }
  }
  
  # new zero correction as an option --> conf. R. Alexandrowicz(2011) p.373 
  if(zerocor==TRUE){nij[nij==0]<-1}
  
  return(nij)
  }

}
