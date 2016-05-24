#' Explore and visualize models
#' 
#' Explore and visualize \samp{Mlxtran} and \samp{pharmML} models 
#' with the \samp{Mlxplore} software.
#' 
#' See http://simulx.webpopix.org/mlxr/mlxplore/ for more details.
#' @param model a \code{Mlxtran} or \code{PharmML} model
#' @param output a list with fields: 
#' \itemize{
#'   \item \code{name}: a vector of output names
#'   \item \code{time}: a vector of times 
#' }
#' @param parameter a vector of parameters with their names and values
#' @param treatment a list with fields
#' \itemize{
#'   \item \code{time} : a vector of input times,
#'   \item \code{amount} : a scalar or a vector of amounts,
#'   \item \code{rate} : a scalar or a vector of infusion rates (default=\code{Inf}),
#'   \item \code{type} : the type of input (default=1),
#'   \item \code{target} : the target compartment (default=NULL). 
#' }
#' @param group a list with unique field: 
#' \itemize{
#'   \item \code{treatment} : a list,
#' }
#' @importFrom tools file_ext
#' @export
mlxplore <- function(model,parameter=NULL,output=NULL,group=NULL,treatment=NULL)
{
  
  # ########################################################################################  
  #  mlxplore.R is governed by the CeCILL-B license. 
  #  You can  use, modify and/ or redistribute the software under the terms of 
  #  the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
  #  http://www.cecill.info/index.en.html
  #
  #  mlxplore.R was developed by Marc Lavielle and the Inria Popix team for the DDMoRe project. 
  # ########################################################################################  
  initMlxLibrary()
  session <- Sys.getenv("session.mlxplore")
  if (is.na(file.info(session)$isdir))
    stop("You need to provide the path of Mlxplore in the file \"initMlR.R\"")
  
  
  Sys.setenv(LIXOFT_HOME=session)
  model <- tools::file_path_as_absolute(model)
  tmproject <- paste0(dirname(model),"/temp_mlxplore.txt")  
  tmpmodel  <- "temp_model.txt"  
  model_ext <- file_ext(model)
  if(model_ext=="xml"){
    model = pharmml2mlxtran(model)
  }  
  
  str <- "<MODEL>"
  
  if(model==tmpmodel){
    str <-c(str,readLines(tmpmodel))
  }else{
    str <- c(str,paste0("file = '",model,"'"))
  }
  
  if(!is.null(treatment)){
    str <- c(str,"","<DESIGN>","[ADMINISTRATION]")
    if (!is.list(treatment[[1]])){
      nadmin <- 1
      treatment[[1]] = treatment
    } else{
      nadmin <-length(treatment)
    }
    for(k in seq(1,nadmin)){
      admk=treatment[[k]]
      strk2=exstr1(admk)
      strk=paste(paste0('adm',k),'=',strk2)
      str <- c(str, strk)
    }  
    
    if(nadmin>1){
      str <- c(str,"[TREATMENT]")
      strt="trt={"
      for(k in seq(1,nadmin)){
        strt=paste(strt,paste0('adm',k))
        if(k<nadmin)
          strt=paste0(strt,",")
      }
      strt <- paste(strt,"}")
      str <- c(str,strt)
    }
  }
  
  if(!is.null(group)){
    if(!is.null(group[[1]]$treatment)){
      str <- c(str,"","<DESIGN>","[ADMINISTRATION]")
      ng=length(group);
      ik=0;
      for(i in seq(1,ng)){
        treati <- group[[i]]$treatment
        if (!is.null(treati$time)){
          treati <- list(treati)
          group[[i]]$treatment <- treati
        }
        nadmin <-length(treati)
        for(k in seq(1,nadmin)){
          ik <- ik+1
          admk=treati[[k]]
          strk2=exstr1(admk)
          strk=paste(paste0('adm',ik),'=',strk2)
          str <- c(str, strk)
        } 
      }
      str <- c(str,"[TREATMENT]")
      ik=0
      for(i in seq(1,ng)){
        strt <- paste0("trt",i," = { ")
        treati <- group[[i]]$treatment
        nadmin <-length(treati)
        for(k in seq(1,nadmin)){
          ik <- ik+1
          strt=paste(strt,paste0('adm',ik))
          if(k<nadmin)
            strt=paste0(strt,",")
        }
        strt <- paste(strt,"}")
        str <- c(str,strt)
      } 
    }}
  
  str <- c(str,"","<PARAMETER>")
  if (isfield(parameter,"name")){
    p.name=parameter$name
    p.value=parameter$value
  }else{
    p.name <- names(parameter)
    p.value <- parameter
  }
  lpn <- length(p.name)
  for(k in seq(1,lpn)){
    strk <- paste(p.name[k],"=",p.value[k])
    str <- c(str,strk)
  }
  
  str <- c(str,"","<OUTPUT>")
  oname=paste(output$name,collapse=", ")
  stro=paste0("list = {",oname,"}")
  str <- c(str,stro)
  
  grid=output$time
  strg=paste0("grid = ",grid[1],":",grid[2]-grid[1],":",grid[length(grid)])
  str <- c(str,strg)
  
  cat(str,file=tmproject,sep="\n")
  
  #--------------------------------------------------------
  
  str=paste0('"',session,'/lib/mlxPlore" --multiple-windows=true --project=',tmproject)  
   system(str, wait=F, invisible=F)
#   system(str, wait=F)
  
}

exstr1 <- function(s){
  ns <- names(s)
  lns <- length(ns)
  u="{"
  for(k in seq(1,lns)){
    uk=paste(ns[k],"=",s[ns[k]])
    uk=gsub("c(","{",uk,fixed=TRUE)
    uk=gsub(")","}",uk,fixed=TRUE)
    u=paste(u,uk)
    if(k<lns){
      u=paste0(u,",")
    }
  }
  u=paste(u,"}")
  return(u)
}
