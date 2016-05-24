#' Paste with parentheses
#' 
#' Paste with parentheses
#' 
#' @param x a vector
#'@keywords helper
#'@export
#'@examples
#'pstprn(c(1,2,3,4,5))
#'pstprn(c("Hello","Hi",2))
pstprn<-function(x){paste(x[1]," (",paste(x[-1],collapse=","),")",sep="")}

#' Round and paste with parentheses
#' 
#' Round and paste with parentheses
#' 
#' @param x a numeic vector
#' @param y integer corresponding to the number of digits to round by
#'@keywords helper
#'@export
#'@examples
#'psthr(c(1.111111,2.2222222,3.333333))
psthr<-function(x,y=2){
  x<-round(x,y)
  pstprn(x)
}
covnm<-function(betanames,call){
  sapply(betanames,function(betaname){
    
    indx<-which(sapply(call,function(cov){charmatch(cov,betaname)})==1)
    if(length(indx)==1) return(call[indx])
    #If one  facorname is a subset of another
    indx2<-which.max(sapply(call[indx],nchar))
    if(length(indx2)==1) return(call[indx[indx2]])
    indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
    if(length(indx3)==1)  return(call[indx[indx2[indx3]]])                      
  })  
}

alleql<-function(x,y){
  !any((x==y)==F)
}


betaindx<-function(x){
  i=1
  out<-1
  result<-NULL
  while(TRUE){
    if(i+1>length(x)){
      result<-c(result,list(out))
      return(result)
    }
    else if(alleql(x[[i+1]],x[[i]])){
      out<-c(out,i+1)
    }
    else{
      result<-c(result,list(out))
      out<-i+1
    }
    i=i+1
  }
}

#' Capitalize a string
#' 
#' Calitalize a string
#' 
#' @param x string
#' @keywords helper
#' @export
cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#'Clean strings for printing
#'
#' Returns strings with . and _ replaced by a space. This is nice when printing column names of your dataframe in a report
#' @param strings vector of strings to give a nice name
#' @keywords helper
#' @export
nicename<-function(strings){
  out<-sapply(strings,function(x){
    x<-chartr(".", " ",x)
    x<-chartr("_", " ",x)
    return(x)})
  return(out)
}

#' Formats p-values
#' 
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits
#' 
#' @param x an integer
#' @export
pvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("<0.001")
  else return(signif(x,2))
}

sanitize <- function(str) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", 
                 result, fixed = TRUE)
  return(result)
}

#' Sanitizes strings to not break LaTeX
#' 
#' Strings with special charaters will break LaTeX if returned 'asis' by knitr. This happens every time we use one of the main reportRx functions. We first sanitize our strings with this function to stop LaTeX from breaking.
#'
#'@param str a vector of strings to sanitize
#'@export
sanitizestr<-function(str){
  as.vector(sapply(str,function(char){sanitize(char)}))
}

#'Bold strings in LaTeX
#'
#'Bold strings in LaTeX.
#'
#'@param strings A vector of strings to bold.
#'@export
lbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste("\\textbf{",x,"}",sep=""))})}

#'Add spaces to strings in LaTeX
#'
#'Add spaces to strings in LaTeX. Returns appends ~~~ before the string
#'
#'@param x string
#'@export
addspace<-function(x){
  paste("~~~",x,sep="")
}
#' Formats p-values for LaTeX
#' 
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to 2 significant digits. Will bold the p-value if it is <= 0.05
#' @param x an integer
#' @export
lpvalue<-function(x){
  if(is.na(x)|class(x)=="character") return(x)
  else if (x<=0.001) return("\\textbf{$<$0.001}")
  else x=signif(x,2)
  if(x<=0.05) return(paste("\\textbf{",x,"}",sep=""))
  else return(x)
}



removedollar<-function(x){
  colnms<-strsplit(x,":")
  indx<-unlist(lapply(colnms,function(colnm) sapply(colnm, function(coln) regexpr("$",coln,fixed=T)[1]+1)))
  if(length(unique(indx))==1){
    if(unique(indx)!=0) x<-unlist(lapply(colnms,function(colnm) paste(substring(colnm,indx[1]),collapse=":")))
  }
  return(x)  
}

modelmatrix<-function(f,data=NULL){
  k<-as.character(f)
  y<-NULL
  if(!length(k)%in%c(2,3)) stop("formula not properly formed")
  if(length(k)==3) {
    f<-as.formula(paste("~",k[2],"+",k[3],sep=""))
    y<-model.matrix(as.formula(paste("~",k[2],sep="")),data)[,-1,drop=F]}
  x<-model.matrix(f,data)[,-1,drop=F]
  colnames(x)<-removedollar(colnames(x))
  if(!is.null(y)){
    return(list(x[,1:ncol(y),drop=F],x[,(ncol(y)+1):ncol(x),drop=F]))
  }else{
    return(x)
  }}

matchcovariate=function(betanames,ucall){
  out=as.vector(sapply(betanames,function(betaname){
    splitbetaname=unlist(strsplit(betaname,":",fixed=T))
    out=sapply(splitbetaname,function(bname){
      indx=which(sapply(ucall,function(cov)charmatch(cov,bname))==1)
      if(length(indx)==1)return(indx)
      #If one  facorname is a subset of another
      indx2<-which.max(sapply(ucall[indx],nchar))
      if(length(indx2)==1) return(indx[indx2])
      indx3<-which(sapply(ucall[indx2],function(c){substr(betaname,1,nchar(c))==c}))
      if(length(indx3)==1)  return(ucall[indx[indx2[indx3]]])  
      return(-1)      
    })
    if(-1 %in% out) return(-1)
    result=0
    n=length(out)
    for(i in 1:length(out)){
      result=result+out[i]*100^(n-1)
      n=n-1
    }
    return(result)}))
  if(-1 %in% out) return(-1)
  return (out)
}
