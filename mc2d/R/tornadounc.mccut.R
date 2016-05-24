#<<BEGIN>>
tornadounc.mccut <- function(mc,output = length(mc), quant=c(0.5,0.75,0.975),use = "all.obs",	method=c("spearman","kendall","pearson"),...)
#ISALIAS tornadounc.mc
#--------------------------------------------
{
	method <- match.arg(method)
 	na.method <- pmatch(use, lesmet <- c("all.obs", "complete.obs", "pairwise.complete.obs"))

  mc <- mc[[which(sapply(mc,inherits,"summary.mccut")==TRUE)]]

  if(length(output) > 1) stop("Only one output permitted")
  if(is.numeric(output)) output <- names(mc)[output]
	if(!(output %in% names(mc)))	stop("Output is not a valid value")

  typen <- sapply(mc,attr,which="type")
  typeout <- typen[output]
	if(typeout!="U" && typeout!="VU")	stop("Output should be a 'U' or a 'VU' node")

  #Select nodes according to the type of the output
  quel <- typen == "U" | (typeout=="VU" & typen=="VU")
  if(sum(quel) < 2) stop("No valid pairs of mcnodes")

  # Select the nodes
  mc <- mc[quel]
  typen <- typen[quel]
  nomn <- names(rapply(mc,function(x)1))


  # Select the column
  yaprob <- length(quant) > 0
  if(yaprob) {nom1 <- paste(quant*100,"%",sep="")
  	nom1[quant==0] <- "Min"
  	nom1[quant==1] <- "Max"} else nom1 <- NULL
  nom1 <- c("mean","sd",nom1)
  mc <- rapply(mc,drop,how="replace")
  mc <- rapply(mc,function(x) x[,nom1],classes="matrix",how="replace")
  mc <- rapply(mc,as.matrix,classes="numeric",how="replace")
  
  ninc <- rapply(mc,function(x) dim(x)[1])[1]

  lesnom <- function(x,nom){
    if(is.list(x)){
      name2 <- names(x)
      if(is.null(name2)) name2 <- paste(".",1:length(x),sep="")
      else name2 <- paste(": ",name2," of variates",sep="")
      nom <- paste(nom,name2,sep="")
      lecol <- colnames(x[[1]])
      return(lapply(1:length(nom),function(x) paste(lecol,nom[x])))
      }
    return(paste(colnames(x),nom))
    }

  noml <- mapply(lesnom,mc,names(mc),SIMPLIFY=FALSE)

  out <- mc[[output]]
  nomout <- noml[[output]]
  inp <- matrix(unlist(mc[names(mc)!= output]),nrow=ninc)
  nomin <- unlist(noml[names(noml)!= output])

   # Which columns are complete
  quelk <-  unique(rapply(mc,function(x) which(apply(x,1,is.na))))

  if(length(quelk) != 0){
    if(na.method==1) stop("NA values. Change 'use' option")
    if(na.method==2){
      out <- lapply(out,"[",,-quelk,,drop=FALSE)    # ! perte structure
      inp <- inp[-quel,]    # ! perte structure
      use <-  "all.obs"}
  }

  if(!is.list(out)) out <- list(out)
  if(!is.list(nomout)) nomout <- list(nomout)
  lescorr <- mapply(function(x) as.matrix(cor(x,inp,method=method,use=use)),out,SIMPLIFY=FALSE)
  lescorr <- lapply(lescorr,"colnames<-",value=nomin)
  lescorr <- mapply("rownames<-",lescorr,value=nomout,SIMPLIFY=FALSE,USE.NAMES=TRUE)

  tc <- list(value = lescorr, output = output, method = method, use = use)
	class(tc) <- "tornadounc"
  return(tc)
	}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
