#' @title Graphical Model Check
#' @export grm
#' @description This function makes the basic calculations for the graphical model check for dicho- or polytomous item response formats. It is more or less a wraper function, internaly calling the function \code{\link{pairSE}}. Several splitting options are available (see arguments).    
#' 
#'@details The data is splitted in two or more subsamples and then item thresholds, the parameter (Sigma) and their standard errors (SE) for the items according the PCM  are calculated for each subsample. Additional arguments (see description of function \code{\link{pairSE}}) for parameter calculation are passed through. 
#'
#'WARNING: When using data based on booklet designs with systematically missing values (by design) you have to ensure that in each of the booklet the maximum raw value to reach is equal while using the raw value as splitting criterion.
#' 
#' @section A note on standard errors: Estimation of standard errors is done by repeated calculation of item parameters for subsamples of the given data. This procedure is mainly controlled by the arguments \code{nsample} and \code{size} (see arguments). With regard to calculation time, the argument \code{nsample} is the 'time killer'. On the other hand, things (estimation of standard errors) will not necessarily get better when choosing large values for \code{nsample}. For example choosing \code{nsample=400} will only result in minimal change for standard error estimation in comparison to (\code{nsample=30}) which is the default setting (see examples). 
#'
#' @param daten a data.frame or matrix with optionaly named colums (names of items), potentially with missing values, comprising polytomous or dichotomous (or mixted category numbers) responses of \code{n} respondents (rows) on \code{k} items (colums) coded starting with 0 for lowest category to \emph{m}-1 for highest category, with \emph{m} beeing a vector (with length k) with the number of categories for the respective item.
#' @param m an integer (will be recycled to a vector of length k) or a vector giving the number of response categories for all items - by default \code{m = NULL}, \code{m} is calculated from data, assuming that every response category is at least once present in data. For sparse data it is strongly recomended to explicitly define the number of categories by defining this argument.
#' 
#' @param split Specifies the splitting criterion. Basically there are three different options available - each with several modes - which are controlled by passing the corresponding character expression to the argument. 
#' 
#' 1) Using the rawscore for splitting into subsamples with the following modes: \code{split = "median"} median raw score split - high score group and low score group; \code{split = "mean"} mean raw score split - high score group and low score group.
#' 
#' 2) Dividing the persons in \code{daten} into subsamples with equal size by random allocation with the following modes: \code{split = "random"} (which is equivalent to \code{split = "random.2"}) divides persons into two subsamples with equal size. In general the number of desired subsamples must be expressed after the dot in the character expression - e.g. \code{split = "random.6"} divides persons into 6 subsamples (with equal size) by random allocation etc. 
#' 
#' 3) The third option is using a manifest variable as a splitting criterion. In this case a vector with the same length as number of cases in \code{daten} must be passed to the argument grouping the data into subsamples. This vector should be coded as \code{"factor"} or a \code{"numeric"} integer vector with min = 1.
#' 
#' @param splitseed numeric, used for \code{set.seed(splitseed)} for random splitting - see argument \code{split}.
#' 
#' @param verbose logical, if \code{verbose = TRUE} (default) a message about subsampling is sent to console when calculating standard errors.
#' @param ... additional arguments \code{nsample}, \code{size}, \code{seed}, \code{pot} for caling \code{\link{pairSE}} are passed through - see description for \code{\link{pairSE}}.
#' 
#' @return A (list) object of class \code{"grm"} containing the item difficulty parameter sigma and their standard errors for two or more subsamples.
#' @exportClass grm
#' @references description of function \code{\link{pairSE}}\code{{pairwise}}.
#' 
#' @examples data(bfiN) # loading example data set
#' 
#' data(bfi_cov) # loading covariates to bfiN data set
#' 
#' # calculating itemparameters and SE for two random allocated subsamples
#' grm_gen <- grm(daten=bfiN, split = bfi_cov$gender)
#' summary(grm_gen)
#' plot(grm_gen)
#' 
#' grm_med <- grm(daten=bfiN, split = "median")
#' summary(grm_med)
#' plot(grm_med)
#' 
#' grm_ran<-grm(daten=bfiN, split = "random") 
#' 
#' summary(grm_ran)
#' 
#' # some examples for plotting options
#' # plotting item difficulties for two subsamples against each other 
#' # with elipses for a CI = 95% .
#' plot(grm_ran) 
#' 
#' # using triangles as plotting pattern
#' plot(grm_ran,pch=2) 
#' 
#' #plotting without CI ellipses
#' plot(grm_ran,ci=0,pch=2) 
#' 
#' # plotting with item names
#' plot(grm_ran,itemNames=TRUE) 
#' 
#' # Changing the size of the item names
#' plot(grm_ran,itemNames=TRUE, cex.names = 1.3)
#' 
#' # Changing the color of the CI ellipses
#' plot(grm_ran,itemNames=TRUE, cex.names = .8, col.error="green")
#' 
#' ###### example from details section 'Some Notes on Standard Errors' ########
#' grm_def<-grm(daten=bfiN, split = "random",splitseed=13)
#' plot(grm_def)
#' ######
#' grm_400<-grm(daten=bfiN, split = "random", splitseed=13 ,nsample=400)
#' plot(grm_400) 
#' 
#' 
############## funktions beginn ########################################################

grm<-function(daten, m=NULL, split="random", splitseed="no", verbose=FALSE, ...){ 
#### abfragen der teilungskriterien und teiler vorbereiten
  teil <- split  # übergabe an internes argument
if(!(length(teil) > 1)) {  
  if(teil=="no"){
    teiler<-rep(1,dim(daten)[1])
    #OK
  }
  if(teil=="random"){
    if (class(splitseed)=="numeric"){set.seed(splitseed)}
    teiler<-as.numeric(cut(sample(1:(dim(daten)[1])),2))
    #OK
    }
  if(nchar(teil)>6){
    nteil<-as.numeric(unlist(strsplit(teil,".",TRUE))[2]) 
    if (class(splitseed)=="numeric"){set.seed(splitseed)}
    teiler<-as.numeric(cut(sample(1:(dim(daten)[1])),nteil))
    #OK
  }     
  if(teil=="mean"){
    daten<-as.matrix(daten)
    rscore<-rowSums(daten,na.rm = TRUE)
    teiler<-factor(ifelse(rscore > round(mean(rscore)),"above mean" ,"mean and below" ))
    #OK
  }
  if(teil=="median"){
    daten<-as.matrix(daten)
    rscore<-rowSums(daten,na.rm = TRUE)
    teiler<-factor(ifelse(rscore > median(rscore),"above median" ,"median and below" ))
    #OK
  }
  
}
  
  if((class(teil)=="integer") | (class(teil)=="numeric") | (class(teil)=="factor")){
    #teiler<-daten[,teil]
    if( (dim(daten)[1])!=length(teil) ){stop("length of argument 'split' dose not match with 'data'")}
    teiler<-teil
    #if (class(teiler)=="factor"){teiler<-(as.numeric(teiler))}
    #if (min(teiler!=1)){stop("argument teil is not valid specified")}
    # daten<-daten[,-teil]
    #OK
  }
#### ENDE abfragen der teilungskriterien und teiler vorbereiten  
# vorbereiten des objektes datalist anhand des vectors teiler
subsamp <- names(table(teiler))

datalist<-vector("list",length=length(subsamp)) #vorber. leere datalist   
 for (i in 1:length(datalist)){
   datalist[[i]]<-daten[which(teiler==subsamp[i]),]  #hier die zuordnung der subsamples aus daten
 }
names(datalist) <- paste(subsamp,"sample")
     
     erg <- lapply(datalist, pairSE, m=m, verbose=verbose, ...) 
     class(erg) <- c("grm","list")

  return(erg)
# summary.grm(grm)
}
