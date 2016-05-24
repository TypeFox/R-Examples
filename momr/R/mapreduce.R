
#-------------------------------------------------------------------
# New parallel computing procedures. To be re-engeniered with an indexing
# approach and architecture of the files.
#-------------------------------------------------------------------

#' \code{splitData} 
#' @title splitData
#' @description This function splits a dataframe object on a given number of equally sized shares
#' and saves them in the disk as RData objects with an incremental name.
#' @author Edi Prifti
#' @param data : the dataframe to be split
#' @param shares : the number of shares, default=30
#' @param name : the name of the split files in the disk preceding the incremental number, default="data_part"
#' @param rows : logical parameter, if default rows=="TRUE" the rows will be split
#' @return does not return anything
splitData <- function(data, shares=30, name="data_part", rows=TRUE){
  if(rows){
    genes_per_share <- ceiling(nrow(data)/shares)
    coords <- as.data.frame(matrix(NA, ncol=2, nrow=length(seq(1,nrow(data),genes_per_share))))
    coords[,1] <- seq(1,nrow(data),genes_per_share)
    coords[,2] <- c(seq(genes_per_share, nrow(data), genes_per_share), nrow(data))
    colnames(coords)<- c("start","stop")
    for(i in 1:nrow(coords)){
      print(paste("Cutting share",i))
      data_part <- data[coords$start[i]:coords$stop[i],]
      save(data_part,file=paste(name,"_",i,".RData",sep=""))
    }
  }else{
    samples_per_share <- ceiling(ncol(data)/shares)
    coords <- as.data.frame(matrix(NA,ncol=2,nrow=length(seq(1,ncol(data),samples_per_share))))
    coords[,1] <- seq(1,ncol(data),samples_per_share)
    coords[,2] <- c(seq(samples_per_share,ncol(data),samples_per_share),ncol(data))
    coords <- coords[-nrow(coords),]
    coords[nrow(coords),2] <- ncol(data)
    colnames(coords)<- c("start","stop")
    for(i in 1:nrow(coords)){
      print(paste("Cutting share",i))
      data_part <- data[,coords$start[i]:coords$stop[i]]
      save(data_part,file=paste(name,"_",i,".RData",sep=""))
    }
  }
}

# modified splitData (momr) April 8th 2015 to obtain "shares"  files and not "shares-1"   elechat
# splitData <- function(data, shares=30,name="data_part",rows=TRUE){
#   if(rows){
#     genes_per_share <- ceiling(nrow(data)/shares)
#     coords <- as.data.frame(matrix(NA,ncol=2,nrow=shares))
#     coords[,1] <- seq(1,nrow(data),genes_per_share)
#     coords[,2] <- c(seq(genes_per_share,nrow(data),genes_per_share),nrow(data))
#     colnames(coords)<- c("start","stop")
#     for(i in 1:nrow(coords)){
#       print(paste("Cutting share",i))
#       data_part <- data[coords$start[i]:coords$stop[i],]
#       save(data_part,file=paste(name,"_",i,".RData",sep=""))
#     }
#   }else{
#     samples_per_share <- ceiling(ncol(data)/shares)
#     coords <- as.data.frame(matrix(NA,ncol=2,nrow=length(seq(1,ncol(data),samples_per_share))))
#     coords[,1] <- seq(1,ncol(data),samples_per_share)
#     coords[,2] <- c(seq(samples_per_share,ncol(data),samples_per_share),ncol(data))
#     colnames(coords)<- c("start","stop")
#     for(i in 1:nrow(coords)){
#       print(paste("Cutting share",i))
#       data_part <- data[,coords$start[i]:coords$stop[i]]
#       save(data_part,file=paste(name,"_",i,".RData",sep=""))
#     }
#   }
# }

#' \code{launchTask} 
#' @title launchTask
#' @description This function distributes the calculations as separate processes in a multi-thread server.
#' @author Edi Prifti
#' @param input : a folder containing elements that are to be processed
#' @param output : a folder where the processed results will be writen
#' @param script : the R script that is to be executed in parallel
#' @return nothing
#' @note the number of processors may be specified. TODO: add the argument
launchTask <- function(input, output, script){
  command <- paste("nohup ls -1 ", input,"/* | xargs -l -i --max-procs=`grep -c '^processor' /proc/cpuinfo` Rscript ",script," {} ",output,"/ &>/dev/null &",sep="")
  system(command)
}


#' \code{watchProgress} 
#' @title watchProgress
#' @description This function prints the progress % of the overall threads that are treated. It is based
#' on the input and output folders by comparing the number of output objects.
#' @author Edi Prifti
#' @param input : a folder containing elements that are to be processed
#' @param output : a folder where the processed results will be writen
#' @return nothing
watchProgress <- function(input, output){
  input.files <- dir(input)
  output.files <- dir(output)
  return(length(output.files)/length(input.files)*100)
}


#' \code{mapResults} 
#' @title mapResults
#' @description This function loads the results once finished and merges them together in one dataframe
#' @author Edi Prifti
#' @param folder : the folder where the results are found
#' @param pattern : the pattern of the split files in the disk preceding the incremental number, default="data_part"
#' @param type : a string c(rows, cols, list) indicating how to merge the data. This should be the same as the one 
#' used in the splitData procedure.
#' @return a merged dataframe : Caution ! The merged results may be shuffeled. It is up to the user to reorder
#' the data accordingly.
mapResults <- function (folder=".", pattern = "_result.rda", type = "rows"){
  if (type %in% c("rows","cols","list")) stop("Please provide a valid mapping type rows,cols,list")
  
  files <- dir(folder, pattern = paste(pattern, sep = "."))
  files.share <- as.numeric(gsub("_result.rda","",files))
  files <- files[order(files.share)] 
  if (type %in% c("rows","cols")){
    res <- c()
    for (i in 1:length(files)) {
      if(i %% 10 == 0) print(paste("i =",i))
      load(paste(folder,files[i],sep="/"))
      data_part_treated <- get(paste("treated"))
      if (type=="rows") {
        res <- rbind(res, data_part_treated)
      } else {
        res <- cbind(res, data_part_treated)
      }
    }  
  }else { # if list
    res <- list()
    for (i in 1:length(files)) {
      if(i %% 10 == 0) print(paste("i =",i))
      load(paste(folder,files[i],sep="/"))
      data_part_treated <- get(paste("treated"))
      res <- c(res,data_part_treated)
    }
  }
  warning("Caution! The merged results may be shuffeled. It is up to the user to reorder the data accordingly")
  return(res)
}


#' \code{deleteData} 
#' @title deleteData
#' @description This function will delete the original and treated split temporary files to celan the workspace.
#' @author Edi Prifti
#' @param name : the name of the split files in the disk preceding the incremental number, default="data_part"
#' @return nothing to be returned
deleteData <- function(name="data_part"){
  files <- dir(pattern=paste(name,sep="_"))
  for(i in 1:length(files)){
    system(paste("rm -f",files[i]))
  }
}

