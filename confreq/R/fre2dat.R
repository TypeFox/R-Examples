#' @title pattern frequency to dataset conversion 
#' @export fre2dat
#' @description Given a (response) pattern frequencies table this function returns a dataset representation of it.
#' @details No details
#' 
#' @param x an object of class "matrix" which is a (response) pattern frequencies table. It is assumed, that the last column of the object \code{x} represents the frequencies of the (respose) patern represented by the other columns in \code{x}.
#' @param fact logical, default is \code{(fact=FALSE)}. If this argument is set to \code{(fact=TRUE)} the result is coerced to a data.frame with factor variables.
#' @param ... additional parameters passed trough. This is an option to assign factor labels to the resulting \code{data.frame} (when setting argument \code{fact=TRUE}) --> see \code{factor} in the \code{base} package and examples.
#' 
#' WARNING using this option will only work correct when all 'pattern' columns (variables) in the frequencies table share the same number of categories    

#' @return An object of class "matrix" or "data.frame" (depending on the argument \code{fact}) containing the dataset representation of the (response) pattern frequencies table give in the argument \code{x}. 
#' @references No references in the moment 
#' @examples #######################################
#' data(LienertLSD)# loading example pattern frequencies table
#' fre2dat(LienertLSD)# coverting it into a (data) matrix
#' # for a matrix without colnames
#' colnames(LienertLSD)<-NULL # first removing the colnames
#' fre2dat(LienertLSD) # conversion with automatic new colnames
#' # requesting a data.frame using factor levels 
#' fre2dat(LienertLSD,fact=TRUE,labels=c("yes","no"))
############### start of function definition ##################
########### pattern frequency to data conversion ##############
################ jhheine at googlemail.com ####################
fre2dat <- function(x,fact=FALSE,...) {
	x <- data.matrix(x)
	lp <- dim(x)[2]
	pp <- dim(x)[1] # könnte zu kurz sein! je nach x
  kat <- apply(x[,1:(lp-1)],2,max) # Anzahl kategorien je variable
	npos<-prod(kat) #anzahl möglicher pattern
  ######## sicherheitshalber richtiges (zeilenweises) sortieren von x
	vars<-1:(lp-1)
	orderlist <- list()
	for(i in 1:length(vars)){
	  orderlist[[i]]<-rank(x[, vars[i]])
	}
	x<-(x[do.call("order", orderlist), ])
	######## ENDE sicherheitshalber richtiges sortieren von x
  
	######## wenn x unvollständig ist
	if (npos!=pp){ 
	  cat("pattern frequencies representation covers not all possible paterns!","\n")
	  ### vervollständigen (sortiert ist ja schon):
	  pos_p<-apply(pos_cfg_cfa(kat, fact = FALSE),1,function(x){paste(x,collapse="")})
	  emp_p<-apply(x[ ,1:(lp-1)],1,function(x){paste(x,collapse="")})
	  fehlend_p<-pos_p[!pos_p%in%emp_p]
	  fehlende<-cbind(matrix(as.numeric(unlist(strsplit(fehlend_p, split=""))),ncol=(lp-1),byrow = TRUE),rep(0,length(fehlend_p)) )
	  x<-rbind(x,fehlende) #fehlende pattern an x anhängen
	  # nun x nochmals sortieren
	  vars<-1:(lp-1)
	  orderlist <- list()
	  for(i in 1:length(vars)){
	    orderlist[[i]]<-rank(x[, vars[i]])
	  }
	  x<-(x[do.call("order", orderlist), ])
	  pp <- dim(x)[1] # neuen wert für pp festsetzten nach vervollständigung
	  cat("missing pattern were added and were given frequencies of zero","\n")
	}
  
  ######## wenn x vollständig ist
  if (npos==pp){    
    fun <- function(x){matrix(rep(c(x[1:(lp-1)]), each=x[lp]), ncol=(lp-1))}
    dat<-do.call(rbind,apply(x,1,fun))
    if(length(colnames(x))!=0){colnames(dat)<-colnames(x)[1:(lp-1)]} # variablennmen aus x uebertragen
    if(length(colnames(x))==0){colnames(dat)<-LETTERS[1:(lp-1)]}
    if(fact==T){dat<-as.data.frame(apply(dat,2,factor,...))
    } 
  }

  return(dat)
}