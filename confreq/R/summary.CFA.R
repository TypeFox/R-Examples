#' @export summary.CFA
#' @S3method summary CFA
#' @title S3 Summary for CFA
#' @description S3 summary method for object of class\code{"CFA"}
#' @param object object of class\code{"CFA"}
#' @param digits integer rounds the values to the specified number of decimal places, default is \code{digits=3}.
#' @param type character with default \code{type="ex.bin.test"}, to return wether the observed pattern are 'Types', 'Antitypes' or not significant at all. Possible options for \code{type} are \code{"pChi"}, \code{"z.pChi"}, \code{"z.pChi"}, \code{"z.pBin"} and \code{"p.stir"}. 
#' @param sorton sort results of local test by any column. by default the output is not sorted. Other options may be \code{"pat."}, \code{"obs."}, \code{"exp."}, \code{"Type"}, \code{"Chi"}, etc. ...
#' @param decreasing logical. Should the sort be increasing or decreasing? see \code{\link{order}}
#' @param showall logical with default \code{showall = FALSE} to return only significant pattern (types / antitypes).   
#' @param ... other parameters passed trough

#sorted by the p-value selected in argument \code{type}
########################### hier die summary method #class CFA #######################
summary.CFA<-function(object, digits=3, type="z.pChi",sorton=NULL, decreasing=FALSE, showall=FALSE, ...){
  local.test <- object$local.test
  global.test <- object$global.test
    
  #object$bonferroni.alpha
  temp1 <- local.test[,which(names(local.test)==type)] < object$bonferroni.alpha
  temp2 <- ifelse(test=local.test$obs. > local.test$exp., yes="+",no="-") 
  
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="." )   },x=temp1,y=temp2 )
  
  temp3 <- round(as.double(local.test$exp.),10)!=round(as.double(local.test$obs.),10)#(28-04-2015)
  Type <- mapply(FUN=function(x,y){ifelse(test=(x==TRUE),yes=y, no="b" )   },x=temp3,y=Type )#(28-04-2015)
  
  templocal <- data.frame(pat.=local.test[,1],round(local.test[,2:3],digits=digits),Type, round(local.test[,4:11],digits=digits)    ,cor.=local.test[,12],round(local.test[,13:14],digits=digits) )
  #print(templocal)
  # sort output by sorton ...
  if (length(sorton)!=0){
  sorter <- order(templocal[,which(names(templocal)==sorton)],decreasing=decreasing)
  erg <- templocal[sorter,]
  }
  if (length(sorton)==0){
    erg <- templocal
  }
  if (showall==FALSE){
    erg <- erg[erg$Type!=".",]
  }
  
  
  cat("function Call:","\n","-------------","\n","Formula:",paste(object$used.formula,collapse = " "),"\n","Variables:", names(object$variables),"\n","Categories:", object$variables,"\n")
  
  cat("\n","results of global tests:","\n", "-----------------------")
  cat("\n","pearson Chi-square test:","\n")
  print(data.frame(global.test$pearson))
  
  cat("\n","likelihood ratio test:","\n")
  print(data.frame(global.test$likelihood.ratio))
  
  cat("\n","Information Criteria:","\n")
  print(data.frame(global.test$infocrit))
  
  cat("\n","results of local tests:","\n", "-----------------------","\n","Type (+) / Antitype (-) based on:", type, "; Bonferoni adj. alpha:", object$bonferroni.alpha,"\n")
  
  if(any(temp3==FALSE)){ cat("\n","Type (b): blanked out (functional CFA)","\n")}#(28-04-2015)}
  
  return(erg)
}