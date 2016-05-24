.ls_multipar <- list(slice_var="k",
                     hemisphere="both",
                     norm_mu=FALSE,
                     norm_sigma=FALSE,
                     as.logical=FALSE,
                     breaks=50,
                     type.breaks="range",
                     palette="terrain.colors",
                     cex=1,
                     col.NA="lightyellow",
                     pch.NA=8,
                     col.midplane="red",
                     axes=TRUE,
                     window=FALSE,
                     legend=TRUE,
                     mfrow=NULL,
                     mar=rep(1.5,4),
                     mgp=c(2,0.5,0),
                     pty=NULL,
                     asp=1,
                     bg="lightblue",
                     xlab="",
                     ylab="",
                     main=NULL,
                     num.main=TRUE,
                     cex.main=1.5,
                     quantiles.legend=TRUE,
                     digit.legend=3,
                     cex.legend=1.5,
                     mar.legend=c(2,7,2,2),
                     main.legend=NULL,
                     outline.index=FALSE,
                     cex.index=c(1,1,1),
                     pch.index=20:22,
                     col.index=c("red","purple","green"),
                     filter.index="2D_N4",
                     width=1000,
                     height=500,
                     path=NULL,
                     unit="px",
                     res=NA                     
                     )

multipar <- function(...,reinit.multipar=FALSE){
  
  if(reinit.multipar==TRUE){
    
#     unlock <- "unlockBinding"
#     eval(parse(text=paste(
#       unlock,"(\".ls_multipar\", env=environment(multipar))",
#       sep="")))  
    unlockBinding(".ls_multipar", env=environment(multipar)) # \Achanger
    assign(".ls_multipar",
           list(slice_var="k",hemisphere="both",norm_mu=FALSE,norm_sigma=FALSE,as.logical=FALSE,
                breaks=50,type.breaks="range",palette="terrain.colors",cex=1,
                col.NA="lightyellow",pch.NA=8,col.midplane="red",axes=TRUE,
                window=FALSE,legend=TRUE,mfrow=NULL,mar=rep(1.5,4),mgp=c(2,0.5,0),pty=NULL,asp=1,bg="lightblue",
                xlab="",ylab="",main=NULL,num.main=TRUE,cex.main=1.5,
                quantiles.legend=TRUE,digit.legend=3,cex.legend=1.5,mar.legend=c(2,7,2,2),main.legend=NULL,
                outline.index=FALSE,cex.index=c(1,1,1),pch.index=20:22,col.index=c("red","purple","green"),filter.index="2D_N4",
                width=1000,height=500,path=NULL,unit="px",res=NA),               
                envir=environment(multipar)     
    )
    lockBinding(".ls_multipar", env=environment(multipar))
    
    return(invisible(get(".ls_multipar",envir=parent.env(environment()))))
  }
  
  args <- list(...)
  n.args <- length(args)
  
  #### si lecture 
  if(n.args==0){ # retourne tout
    
    value <- selectMultipar()
    single <- FALSE
    
  }else if(all(unlist(lapply(args, is.character)))){ # retourne uniquement les arguments demandes
    
    args <- unlist(args)
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- selectMultipar(args)
    
  } else { #### si ecriture 
    if (length(args) == 1) {single <- TRUE} else {single <- FALSE}
    value <- allocMultipar(args,n.args)
  }
  
  #### export
  if (single){ 
    value <- value[[1L]]
  }
  
  if(!is.null(names(args))){
    invisible(value)
  }else{
    value
  } 
}

selectMultipar <- function(field=NULL){
  
  .ls_multipar <- get(".ls_multipar",envir=environment(multipar))
   
  if(is.null(field)){
    
    return(.ls_multipar)
    
  }else{
    
    if(any(field %in% names(.ls_multipar) == FALSE)){
      stop("selectMultipar : wrong arguments \n",
           "unknown argument : ",paste(field[field %in% names(.ls_multipar) == FALSE],collapse=" ")," \n",
           "available arguments : ",paste(names(.ls_multipar)[names(.ls_multipar) %in% field == FALSE],collapse=" "),"\n")
    }
    
    return(.ls_multipar[field])
  }
  
}

allocMultipar <- function(field,n.args){
  
  names.field <- names(field)
  .ls_multipar <- get(".ls_multipar",envir=parent.env(environment()))
#   .ls_multipar <-  <- get(".ls_multipar",envir=.GlobalEnv)
  
  if(any(names.field %in% names(.ls_multipar) == FALSE)){
    stop("selectMultipar : wrong arguments \n",
         "unknown argument : ",paste(names.field[names.field %in% names(.ls_multipar) == FALSE],collapse=" ")," \n",
         "available arguments : ",paste(names(.ls_multipar)[names(.ls_multipar) %in% field == FALSE],collapse=" "),"\n")
  }
  
  for(iter_field in 1:n.args){
  .ls_multipar[[names.field[iter_field]]] <- field[[iter_field]]
  }
  
  validCharacter(.ls_multipar$slice_var,name="slice_var",type="slice_var",method="allocMultipar")
  validCharacter(.ls_multipar$hemisphere,name="hemisphere",type="hemisphere",method="allocMultipar")
  validCharacter(.ls_multipar$norm_mu,name="norm_mu",type="norm_mu",method="allocMultipar")
  validCharacter(.ls_multipar$norm_sigma,name="norm_sigma",type="norm_sigma",method="allocMultipar")
  validLogical(.ls_multipar$as.logical,name="as.logical",method="allocMultipar")
  validNumeric(.ls_multipar$breaks,name="breaks",integer=TRUE,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$type.breaks,name="type.breaks",type="type.breaks",method="allocMultipar")
  validCharacter(.ls_multipar$palette,name="palette",type="palette",method="allocMultipar")
  validNumeric(.ls_multipar$cex,name="cex",integer=TRUE,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$col.NA,name="col.NA",method="allocMultipar")
  validNumeric(.ls_multipar$pch.NA,name="pch.NA",integer=TRUE,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$col.midplane,name="col.midplane",method="allocMultipar")
  validLogical(.ls_multipar$axes,name="axes",method="allocMultipar")
  validCharacter(.ls_multipar$legend,name="legend",type="legend",method="allocMultipar")
  validNumeric(.ls_multipar$mfrow,name="mfrow",integer=TRUE,length=2,min=0,refuse.NULL=FALSE,method="allocMultipar")
  validNumeric(.ls_multipar$mar,name="mar",length=4,min=0,refuse.NULL=TRUE,method="allocMultipar")
  validNumeric(.ls_multipar$mgp,name="mgp",length=3,min=0,refuse.NULL=TRUE,method="allocMultipar")
  validCharacter(.ls_multipar$pty,name="pty",validValues=c("m","s"),refuse.NULL=FALSE,method="allocMultipar")
  validNumeric(.ls_multipar$asp,name="asp",integer=FALSE,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$bg,name="bg",method="allocMultipar")
  validCharacter(.ls_multipar$xlab,name="xlab",method="allocMultipar")
  validCharacter(.ls_multipar$ylab,name="ylab",method="allocMultipar")
  validCharacter(.ls_multipar$main,name="main",refuse.NULL=FALSE,method="allocMultipar")
  validLogical(.ls_multipar$quantiles.legend,name="quantiles.legend",method="allocMultipar")
  validNumeric(.ls_multipar$digit.legend,name="digit.legend",integer=FALSE,min=0,method="allocMultipar")
  validNumeric(.ls_multipar$cex.legend,name="cex.legend",min=0,method="allocMultipar")
  validNumeric(.ls_multipar$mar.legend,name="mar.legend",integer=TRUE,length=4,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$main.legend,name="main.legend",refuse.NULL=FALSE,method="allocMultipar")
  validNumeric(.ls_multipar$width,name="width",min=0,max=NULL,refuse.NA=TRUE,method="allocMultipar")
  validNumeric(.ls_multipar$height,name="height",min=0,max=NULL,refuse.NA=TRUE,method="allocMultipar")
  validCharacter(.ls_multipar$path,name="path",refuse.NULL=FALSE,method="allocMultipar")
  validPath(.ls_multipar$path,name="path",method="allocMultipar")
  validCharacter(.ls_multipar$unit,name="unit",type="unit",method="allocMultipar")  
  validNumeric(.ls_multipar$res,name="res",refuse.NA=FALSE,method="allocMultipar")
  validLogical(.ls_multipar$outline.index,name="outline.index",method="allocMultipar")
  validNumeric(.ls_multipar$cex.index,name="cex.index",length=3,min=0,method="allocMultipar")
  validNumeric(.ls_multipar$pch.index,name="pch.index",length=3,min=0,method="allocMultipar")
  validCharacter(.ls_multipar$col.index,name="col.index",length=3,method="allocMultipar")
  validCharacter(.ls_multipar$filter.index,name="filter.index",type="Neighborhood",method="allocMultipar")

#   unlock <- "unlockBinding"
#   eval(parse(text=paste(
#     unlock,"(\".ls_multipar\", env=environment(multipar))",
#   sep="")))  
  unlockBinding(".ls_multipar", env=environment(multipar)) # A changer
  assign(".ls_multipar",.ls_multipar,envir=environment(multipar))
  lockBinding(".ls_multipar", env=environment(multipar))
  
  return(field)
}

#### validation functions 
validCharacter <- function(value,name,length=1,validValues=NULL,type=NULL,refuse.NULL=TRUE,method){
    
  ## type
  if(!is.null(type) && type=="hemisphere" && is.null(validValues)){
    validValues <- c("both","left","right","lesion","contralateral")
  }
  
  if(!is.null(type) && type=="legend" && is.null(validValues)){
    validValues <- c(TRUE,FALSE,"only")
    refuse.NULL <- FALSE
  }
  
  if(!is.null(type) && type=="norm_mu" && is.null(validValues)){
    validValues <- c(FALSE,
                     "global","global_1slice","global_3slices",
                     "contralateral","contralateral_1slice","contralateral_3slices",
                     "default_value")
  }
  
  if(!is.null(type) && type=="norm_sigma" && is.null(validValues)){
    validValues <- c(FALSE,
                     "global","global_1slice","global_3slices",
                     "contralateral","contralateral_1slice","contralateral_3slices",
                     "default_value")
  }
  
  if(!is.null(type) && type=="palette" && is.null(validValues)){
    if(length(value)==1){
    validValues <- c("rgb","hsv",
                     "rainbow","grey.colors","heat.colors","terrain.colors","topo.colors","cm.colors")
    }else{
      length  <- NULL
    }
    
  }
  
  if(!is.null(type) && type=="type.breaks" && is.null(validValues)){
    validValues <- c("range","range_center","quantile")
  }
  
  if(!is.null(type) && type=="unit" && is.null(validValues)){
    validValues <- c("px","in","cm","mm")
  }
  
  if(!is.null(type) && type=="window" && is.null(validValues)){
    validValues <- c(TRUE,FALSE,"png","eps","svg","pdf")
    refuse.NULL <- FALSE
  }
  
  if(!is.null(type) && type=="Neighborhood" && is.null(validValues)){
    validValues <- c("2D_N4","2D_N8","3D_N4","3D_N6","3D_N8","3D_N10","3D_N18","3D_N26")
  }
  
  ## tests
  if(!is.null(value) || refuse.NULL==TRUE){
    
    if(!is.null(length) && length(value)!=length){
      stop(method," : wrong specification of \'",name,"\' \n",
           "\'",name,"\' must have length ",length,"  \n",
           "length(",name,") : ",length(value),"\n")
    }
    
    if( any( (is.character(value)==FALSE)*(is.logical(value)==FALSE)>0 ) ){
      stop(method," : wrong specification of \'",name,"\' \n",
           "\'",name,"\' must be a character or logical \n",        
           "is(",name,") : ",paste(is(value),collapse=" "),"\n")
    }
    
    if(!is.null(validValues) && any(value %in% validValues==FALSE)){
      stop(method," : wrong specification of \'",name,"\' \n",
           "\'",name,"\' must not be one of ",if(refuse.NULL==FALSE){"NULL"}," \"",paste(validValues,collapse="\" \""),"\" \n",
           "refused values of \'",name,"\' : \"",paste(value[value %in% validValues==FALSE],collapse=" "),"\"\n")
    }
  }

}

validLogical <- function(value,name,length=1,refuse.NA=TRUE,method){
  
  if(length(value)!=length){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must have length ",length,"  \n",
         "length(",name,") : ",length(value),"\n")
  }
  
  if(is.logical(value)==FALSE){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must be TRUE or FALSE \n",        
         "is(",name,") : ",paste(is(value),collapse=" "),"\n")
  }
  
  if(is.na(value) && refuse.NA==TRUE){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must not be an NA \n")
  }
  
}
  
validNumeric <- function(value,name,integer=FALSE,length=1,min=NULL,max=NULL,refuse.NA=TRUE,refuse.NULL=TRUE,method){
  
  if(refuse.NULL==TRUE || is.null(value)==FALSE){
  if(length(value)!=length){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must have length ",length,"  \n",
         "length(",name,") : ",length(value),"\n")
  }
  
  if(any(is.na(value)) && refuse.NA==TRUE){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must not contain NA \n",
         "index of NA values : ",which(paste(is.na(value),collapse=" ")),"\n")
  }
  if(any( (is.numeric(value)==FALSE)*(is.na(value)==FALSE) )){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must be a numeric \n",        
         "is(",name,") : ",paste(is(value),collapse=" "),"\n")
  }
  
  if(!is.null(min) && any(stats::na.omit(value)<min)){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must be bigger than ",min," \n",        
         "invalid value(s) in ",name," : ",paste(value[stats::na.omit(value)<min],collapse=" "),"\n")
  }
  
  if(!is.null(max) && any(stats::na.omit(value)>max)){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must be smaller than ",max," \n",        
         "invalid value(s) in ",name," : ",paste(value[stats::na.omit(value)>max],collapse=" "),"\n")
  }
  
  if(integer==TRUE && any(value %% 1 >0)){
    stop(method," : wrong specification of \'",name,"\' \n",
         "\'",name,"\' must contain integers not doubles \n",        
         "invalid value(s) in ",name," : ",paste(value[value %% 1 >0],collapse=" "),"\n")
  }
  }
}

validPath <- function(value,name,method){
  
  if(!is.null(value)){
    try_path <- try(setwd(value),silent=TRUE)
    
    if(is(try_path)=="try-error"){
      stop(method," : wrong specification of \'",name,"\' \n",
           "proposed ",name," : \"",value,"\" \n",
           "error : ",try_path[1],
           "current ",name," : ",getwd(),"\n")
    }
    
    
    if(substr(value,start=nchar(value),stop=nchar(value))!="/"){
      warning(method," : possible bad specification of \'",name,"\' \n",
              "\'",name,"\' should end with a fsep (e.g. \"/\") \n",
              "proposed ",name," : ",value,"\n")
    }
  }
}