#' Convert a Monolix Project  into an executable for the simulator  Simulx 
#' @param project : the name of a Monolix project 
#' @param parameter : string $(NameOfTypeOfParameter), the type of specific parameters to use 
#'                   example: "mode", "mean"...
#' @param group : a list with the number of subjects 
#' @param open : load the R script created if \code{open=TRUE}
#' @param r.data : read the data if \code{r.data=TRUE}
#' @param fim : Fisher information matrix
#' @return  creates a folder projectNameR  containing files : 
#' \itemize{
#'   \item \code{projectName.R} :  executable R code for the simulator,
#'   \item \code{treatment.txt} :  contains the treatment informations,
#'   \item \code{populationParameter.txt} : contains the  population parameters estimated from Monolix,
#'   \item \code{individualParameter.txt} : contains the  individual parameters (mode/mean) estimated from Monolix (if used for the simulation),
#'   \item \code{individualCovariate.txt} : contains the individual covariates,
#'   \item \code{originalId.txt} : contains the original id's when group is used with a different size than the original one,
#'   \item \code{outputi.txt} : contains the output number i informations (time, id),
#'   \item \code{$(NameOfTypeOfParameter)s.txt} : contains the specific parameter used.
#' }       
#' 
#' @examples
#' \dontrun{
#' project.file <- 'monolixRuns/theophylline1_project.mlxtran'  #relative path
#' monolix2simulx(project=project.file,open=TRUE)
#' monolix2simulx(project=project.file,parameter=list("mean",c(a=0, b=0)),open=TRUE)
#' }
#' @export

monolix2simulx <-function(project,parameter=NULL,group=NULL,open=FALSE,r.data=TRUE,fim=NULL)
{ 
  #------- project to be converted into Simulx project
  myOldENVPATH = Sys.getenv('PATH');
  initMlxLibrary()
  session=Sys.getenv("session.simulx")
  Sys.setenv(LIXOFT_HOME=session)
  if  (!is.null(names(group)))
    group <- list(group)
  ans <- processing_monolix(project=project,
                            model=NULL,
                            treatment=NULL,
                            parameter=parameter,
                            output=NULL,
                            group=group,
                            r.data=r.data,
                            fim=fim)
  model         <- ans$model
  treatment     <- ans$treatment
  parameter     <- ans$param
  output        <- ans$output
  group         <- ans$group
  regressor     <- ans$regressor
  occasion      <- ans$occ
  fim           <- ans$fim
  mlxtranpath <- dirname(project)
  mlxtranfile = file_path_sans_ext(basename(project))
  mypath <- getwd()
  Rproject <- file.path(mypath,paste0(mlxtranfile,"_simulx"))
  if(file.exists(Rproject) )
    unlink(Rproject, recursive = TRUE, force = TRUE)
  modelname = basename(model)
  Sys.sleep(0.2)
  dir.create(Rproject, showWarnings = FALSE, recursive = FALSE, mode = "0777")
  file.copy(model, Rproject, overwrite = FALSE)
  file.remove(model)
  model<-file.path(Rproject,modelname)
  
  #configure and write output 
  RprojectPath <- dirname(model)
  mlxtranfile = file_path_sans_ext(basename(project))
  projectExe <- file.path(RprojectPath,paste0(mlxtranfile,".R"))
  cat(paste0("# File generated automatically on ", Sys.time(),"\n \n"), file =projectExe, fill = FALSE, labels = NULL,append = TRUE)
  cat("library(mlxR)  \n \nsetwd(dirname(parent.frame(2)$ofile)) \n\n# model \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  cat(paste0("model<-\"",modelname,"\"\n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  # write  treatment 
  if(!(is.null(treatment)) && length(treatment)>0){ 
    if (!is.null(treatment$value)){
      treat2<-matrix(treatment$value,nrow=nrow(treatment$value),ncol=ncol(treatment$value))
      colnames(treat2)<-treatment$colNames
      treatment <- treat2
    }
    write.table(treatment,file=file.path(Rproject,"/treatment.txt"),row.names=FALSE,quote=FALSE)
    cat("\n# treatment\n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    cat("trt <- read.table(\"treatment.txt\", header = TRUE) \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  }
  
  # write  parameters   
  if(!(is.null(parameter))){  
    param.list <- NULL
    cat("\n# parameters \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    
    if (!is.null(ans$id)){
      outfile = file.path(Rproject,paste0("/originalId.txt"))      
      write.table(ans$id,file=outfile,row.names=FALSE,quote=FALSE)
      cat(paste0("originalId<- read.table('originalId.txt', header=TRUE) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
    }
    
    populationParameter <- parameter[[1]]
    if (!is.null(populationParameter)){
      outfile = file.path(Rproject,paste0("/populationParameter.txt"))      
      cat(paste0("populationParameter<- read.vector('populationParameter.txt') \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
      write.table(populationParameter,file=outfile,col.names=FALSE,quote=FALSE)
      
      if (!is.null(param.list))
        param.list <- paste(param.list,"populationParameter",sep=",")  
      else
        param.list <- "populationParameter"
    } 
    
    individualCovariate <- parameter[[2]]
    if (!is.null(individualCovariate)){
      outfile = file.path(Rproject,paste0("/individualCovariate.txt"))      
      cat(paste0("individualCovariate<- read.table('individualCovariate.txt', header = TRUE) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
      write.table(individualCovariate,file=outfile,row.names=FALSE,quote=FALSE)
      if (!is.null(param.list))
        param.list <- paste(param.list,"individualCovariate",sep=",")  
      else
        param.list <- "individualCovariate"
      i.factor <- which(sapply(individualCovariate[-1], is.factor))
      if (length(i.factor)>0){
        cat(paste0("individualCovariate[,",i.factor+1,"]<- as.factor(individualCovariate[,",i.factor+1,"]) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
      }
    } 
    individualParameter <- parameter[[3]]
    if (!is.null(individualParameter)){
      outfile = file.path(Rproject,paste0("/individualParameter.txt"))      
      cat(paste0("individualParameter<- read.table('individualParameter.txt', header = TRUE) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
      write.table(individualParameter,file=outfile,row.names=FALSE,quote=FALSE)
      if (!is.null(param.list))
        param.list <- paste(param.list,"individualParameter",sep=",")  
      else
        param.list <- "individualParameter"
    } 
    
    param.list <- paste(param.list,sep=",")
    param.str <- paste0("list.param <- list(",param.list,")")
    cat(param.str, file =projectExe, fill = FALSE, labels = NULL, append = TRUE)   
  }
  
  # write f.i.m
  if(!(is.null(fim))) 
    write.table(fim,file=file.path(Rproject,"/fim.txt"),row.names=FALSE,quote=FALSE) 
  
  # write  requested output 
  if(!(is.null(output)))
  {  
    cat("\n# output \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    
    if(length(output)==1)
    {
      out1 <- output[[1]]
      out1.name <- out1$name 
      cat(paste0("name<-\"",out1.name,"\"\n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
      cat(paste0("time<-read.table(\"output.txt\",header=TRUE)\n"),file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
      cat(paste0("out<-list(name=name,time=time) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
      outfile = file.path(Rproject,"/output.txt")
      write.table(out1$time,file=outfile,row.names=FALSE,quote=FALSE) 
    } else {    # many types of output could exist
      for(i in seq(1:length(output))) {
        outi <- output[[i]]
        outi.name <- outi$name
        cat(paste0("name<-\"",outi.name,"\"\n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
        if(is.data.frame(outi$time)) {
          cat(paste0("time<-read.table(\"output",i,".txt\",header=TRUE)\n"),file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
          cat(paste0("out",i,"<-list(name=name,time=time) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
          outfile = paste0(file.path(Rproject,paste0("/output",i)),".txt")
          write.table(outi$time,file=outfile,row.names=FALSE,quote=FALSE) 
        } else {
          cat(paste0("out",i,"<-list(name=name) \n"), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
        }
      }
      
      cat("out<-list(out1", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
      for(i in seq(2,length(output))) {
        cat(paste0(",out",i), file =projectExe, fill = FALSE, labels = NULL, append = TRUE)   
      }
      cat(")\n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    }
  }
  
  # regressor    
  if(!(is.null(regressor))) {  
    cat("\n# regressor \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    outfile = file.path(Rproject,paste0("/regressor.txt"))    
    write.table(regressor,file=outfile,row.names=FALSE,quote=FALSE)
    cat(paste0("regressor <-read.table(\"regressor.txt\", header = TRUE)\n"),file =projectExe, fill = FALSE, labels = NULL, append = TRUE)             
  }
  
  # occasion    
  if(!(is.null(occasion))) {  
    cat("\n# occasion \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
    outfile = file.path(Rproject,paste0("/occasion.txt"))      
    write.table(occasion,file=outfile,row.names=FALSE,quote=FALSE)
    cat(paste0("occasion <-read.table(\"occasion.txt\", header = TRUE)\n"),file =projectExe, fill = FALSE, labels = NULL, append = TRUE)             
  }
  
  # call the simulator
  cat("\n# call the simulator \n", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  cat("res <- simulx(model=model", file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  if(!(is.null(treatment))&& length(treatment)>0)  
    cat(",treatment=trt",file =projectExe, fill = FALSE, labels = NULL, append = TRUE) 
  
  if(!(is.null(parameter)))
    cat(",parameter=list.param",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  if(!(is.null(group))) 
    cat(",group=grp",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  if(!(is.null(output)))
    cat(",output=out",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  if(!(is.null(regressor)))
    cat(",regressor=regressor",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  if(!(is.null(occasion)))
    cat(",varlevel=occasion",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  cat(")\n",file =projectExe, fill = FALSE, labels = NULL, append = TRUE)
  
  Sys.setenv('PATH'=myOldENVPATH);
  if( (Sys.getenv("RSTUDIO")=="1") & (open==TRUE) ) {
    eval(parse(text='file.edit(projectExe)'))
    # file.edit(projectExe) 
    setwd(mypath)
  }
  return(projectExe)
}

