#'  Pearson product-moment correlations
#'
#'  R implementation of the SPSS \code{CORRELATIONS} function.
#'
#' @usage xpssCorrelations(x, 
#'                variables = NULL, 
#'                miss = list(alternative = "pairwise", 
#'                          missings = "exclude"), 
#'                print = list(test = "twotail", 
#'                        level = "sig"), 
#'                matrix = NULL, 
#'                statistics = NULL)
#' @param x a (non-empty) data.frame or input data of class "xpss-Frame".
#' @param variables atomic character or character vektor with the name of the variables.   
#' @param miss method which indicates what should happen when the data contain NAs. Default for alternative is 'pairwise', optionally listwise can be used as treatment for missings. The visualisation of the NAs can be specied via the argument missings. Default is 'exclude', optionally 'include' can be chosen to add missings in the statistics.
#' @param print method which indicates what significnace level shall be used. Default significance test is 'twotail', optionally 'onetail' can be chosen. Default significance level is 'sig' to add significance asterisks, optionally 'nosig'.
#' @param matrix exports the correlation matrix with observations, stddevs, means and variable names. Default is NULL.
#' @param statistics method which enumerate the deskriptive statistics. Default is NULL. Optionally 'descriptives', 'xprod' or 'all' can be chosen.
#'
#' @return Returns a matrix of Pearson's r correlation. 
#' 
#' @details \code{xpssCorrelations} produces Pearson product-moment correlations with significance levels and, optionally, univariate statistics, covariances, and cross-product deviations.
#'  
#'
#' @author Benjamin Piest
#' @seealso \code{\link{cor}} \code{\link{cor.test}} \code{\link{rcorr}}
#' @importFrom Hmisc rcorr
#' @examples 
#'
#'data(fromXPSS)
#'
#'
#'
#'xpssCorrelations (fromXPSS, 
#'                  variables =c("V5","V6","V7_2"))
#'
#'
#'
#'xpssCorrelations (fromXPSS, 
#'                  variables =c("V5","V6","V7_2") ,
#'                  miss = list(alternative = "pairwise",
#'                                 missings = "exclude"),
#'                  print = list(test = "onetail",
#'                               level = "sig"),
#'                  statistics="all")
#'
#'
#'
#'xpssCorrelations (fromXPSS, 
#'                  variables =c("V5","V6","V7_2") ,
#'                  miss = list(alternative = "pairwise",
#'                                 missings = "include"),
#'                  print = list(test = "onetail",
#'                               level = "sig"),
#'                  statistics="all")
#' 
#' 
#'
#'xpssCorrelations (fromXPSS, 
#'                  variables =c("V5","V6","V7_2") ,
#'                  miss = list(alternative = "listwise",
#'                                 missings = "exclude"),
#'                  print = list(test = "twotail",
#'                               level = "sig"),
#'                  statistics="all")
#' 
#' 
#'
#'xpssCorrelations (fromXPSS,
#'                  variables =c("V5","V6","V7_2"),
#'                  statistics = "all",
#'                  matrix = paste0(getwd(),"/correlations.txt"))
#' 
#' @export 

xpssCorrelations <- function(x, 
                            variables = NULL,
                            miss = list(alternative = "pairwise",
                                           missings = "exclude"),
                            print = list(test = "twotail",
                                         level = "sig"),
                            matrix = NULL,
                            statistics=NULL){

#--------------------------------------------------------#
#---------------------  Exceptions  ----------------------#
#--------------------------------------------------------#
  
  if(is.null(variables)) stop("Variables are missing")
  if(length(variables)<2) stop("Variables needs two or more objects")
  not <- variables[variables %in% names(x)=="FALSE"]
  if(length(not)==1) stop(paste("Variable",not ,"is not element of x"))
  if(length(not)>1) stop(paste("Variables",paste(not[1:length(not)-1],collapse=", "),"and",paste(not[length(not)]),"are not element of x"))
  
#--------------------------------------------------------#
#------------------ print test --------------------------#
#--------------------------------------------------------#

if(print$test=="twotail"){
  if(miss$missings=="include"){
    x <- xpssMissingValues(x, variables = variables)
  }
  x1 <- x[variables] 
  var.label <- 1:length(x1)
  for(i in 1:length(x1)){
    var.label[i] <- lapply(x1, attributes)[[i]][["variable.label"]]
  }
  names(x1) <- var.label
  var <- as.matrix(x1)
  if(miss$alternative=="listwise"){
    var <- na.omit(var)
  }
  cor <- rcorr(var)
  for(i in 1:ncol(var)){
    cor$n[names(var[i,])[i],names(var[i,])[i]]=NROW(na.omit(var[,i]))
  }
  attributes(cor) <- NULL
  note <- "(2-tailed)"
  sig.name <- paste0("Sig_",print$test) 
}

if(print$test=="onetail"){
  if(miss$missings=="include"){
    x <- xpssMissingValues(x, variables = variables)
  }
  x1 <- x[variables]
  var.label <- 1:length(x1)
  for(i in 1:length(x1)){
    var.label[i] <- lapply(x1, attributes)[[i]][["variable.label"]]
  }
  names(x1) <- var.label
  var <- as.matrix(x1)
  if(miss$alternative=="listwise"){
    var <- na.omit(var)
  }
  cor <- rcorr(var)
  for(i in 1:ncol(var)){
    cor$n[names(var[i,])[i],names(var[i,])[i]]=NROW(na.omit(var[,i]))
  }
  attributes(cor) <- NULL
  cor[[3]] <- cor[[3]]/2
  note <- "(1-tailed)"
  sig.name <- paste0("Sig_",print$test) 
}
  
#--------------------------------------------------------#
#------------------ statistics --------------------------#
#--------------------------------------------------------#

#---- descriptives ----#

if(is.null(statistics)==FALSE){
  if(statistics=="descriptives"){
    mean <- round(apply(var,2,mean,na.rm=T),2)
    sd <- round(apply(var,2,sd,na.rm=T),2)
    N <- round(apply(var, 2, function(x) length(which(!is.na(x)))),2)
    all <- rbind(mean,sd,N)
    all <- t(all)
    ifelse(miss$alternative=="listwise",cor <- list("descriptive_statistics"=all,"Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3)), cor <- list("descriptive_statistics"=all,"Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3),"N"=round(cor[[2]],2)))
    names(cor)[which(names(cor) %in% "sig.name")] <- sig.name
    }
  

  if(statistics=="xprod"){
    cov <- round(cov(var,use="pairwise.complete.obs"),2)
    df <- cor[[2]]-1
    crossproduct <- round(cov(var,use="pairwise.complete.obs")*df,2)
    ifelse(miss$alternative=="listwise",cor <- list("Pearson_Correlation"=round(cor[[1]],2),"sig.name"=round(cor[[3]],2),"Crossproduct"=crossproduct,"Covariance"=cov),cor <- list("Pearson_Correlation"=round(cor[[1]],2),"sig.name"=round(cor[[3]],2),"Crossproduct"=crossproduct,"Covariance"=cov,"N"=round(cor[[2]],2)))
    names(cor)[which(names(cor) %in% "sig.name")] <- sig.name
  }

#---- all ----#

if(statistics=="all"){
  mean <- round(apply(var,2,mean,na.rm=T),2)
  sd <- round(apply(var,2,sd,na.rm=T),2)
  N <- round(apply(var, 2, function(x) length(which(!is.na(x)))),2)
  all <- rbind(mean,sd,N)
  all <- t(all)
  cov <- round(cov(var,use="pairwise.complete.obs"),2)
  df <- cor[[2]]-1
  crossproduct <- round(cov(var,use="pairwise.complete.obs")*df,2)
  ifelse(miss$alternative=="listwise",cor <- list("descriptive_statistics"=all,"Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3),"Crossproduct"=crossproduct,"Covariance"=cov),cor <- list("descriptive_statistics"=all,"Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3),"Crossproduct"=crossproduct,"Covariance"=cov,"N"=round(cor[[2]],2)))
  names(cor)[which(names(cor) %in% "sig.name")] <- sig.name
}}

#-----------------------------------------------------#
#----------------- Correlations only -----------------#
#-----------------------------------------------------#



if(length(cor)==3){
  ifelse(miss$alternative=="listwise",cor <- list("Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3)),cor <- list("Pearson_Correlation"=round(cor[[1]],2),sig.name=round(cor[[3]],3),"N"=round(cor[[2]],2)))
  names(cor)[which(names(cor) %in% "sig.name")] <- sig.name
}

#--------------------------------------------------------#
#-----------------  Print level  ------------------------#
#--------------------------------------------------------#

  for(i in 1:length(cor[[sig.name]])){
    if(is.na(cor[[sig.name]][[i]])==TRUE){
      cor[[sig.name]][[i]]=999
    }}

if(print$level=="sig"){
  sig99 <- FALSE
  sig95 <- FALSE
  for(i in 1:length(cor[[sig.name]])){
    if(cor[["Pearson_Correlation"]][[i]]>=0){ 
      cor[["Pearson_Correlation"]][[i]] <- paste0(" ",cor[["Pearson_Correlation"]][[i]])
    }
    if(cor[[sig.name]][[i]]<=0.01){ 
      cor[["Pearson_Correlation"]][[i]] <- paste0(cor[["Pearson_Correlation"]][[i]],"**")
      sig99 <- TRUE
    
    }
    if(cor[[sig.name]][[i]]>0.01 & cor[[sig.name]][[i]]<=0.05){ 
      cor[["Pearson_Correlation"]][[i]] <- paste0(cor[["Pearson_Correlation"]][[i]],"*")
      sig95 <- TRUE
    }
  }
  
  if(sig99=="TRUE" & sig95=="TRUE"){
    note <- matrix(c(paste(".Correlation is significant at the 0.01 level",note),paste(".Correlation is significant at the 0.05 level",note)),2,1,byrow = TRUE,dimnames=list(c("**","*"),c("")))
    cor <- c(cor,list("Note"=note))
       }
  
    if(sig99=="TRUE" & is.null(cor$Note)){
      note <- matrix(c(paste(".Correlation is significant at the 0.01 level",note)),1,1,byrow = TRUE,dimnames=list(c("**"),c("")))
      cor <- c(cor,list("Note"=note))
    }
  
    if(sig95=="TRUE" & is.null(cor$Note)){
      note <- matrix(c(paste(".Correlation is significant at the 0.05 level",note)),1,1,byrow = TRUE,dimnames=list(c("*"),c("")))
      cor <- c(cor,list("Note"=note))
    }
}

for(i in 1:length(cor[[sig.name]])){
if(cor[[sig.name]][[i]]==999){
  cor[[sig.name]][[i]]=""
}}

#-----------------------------------------------------#
#-----------------  Matrix argument  -----------------#
#-----------------------------------------------------#

if(!is.null(matrix)){
  
  sink(matrix, append=TRUE)
  cat(" \n")
  cat("===================================\n")
  cat("Pearson product-moment correlations\n")
  cat("===================================\n")
  cat(" \n")
  print(noquote(cor))
  sink()
}
cor <- noquote(cor)
return(cor)
}


