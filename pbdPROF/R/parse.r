#' @importFrom utils read.table
#' @importFrom stats embed



### Internal parsing methods
parse.prof <- function(x, ...) UseMethod("parse.prof")



### default method (Don't delete/change this.)
parse.prof.default <- function(x, ...){
  stop("Profiler type is not found.")
}



parse.prof.fpmpi <- function(x, ...){
  # For return
  ret <- list()
  
  # Get informative region
  id.start <- grep("^Details for each MPI", x)
  if(length(id.start) == 0){
    stop("Input file format is invalid.")
  } else{
    id.start <- id.start + 5
  }

  id.end <- grep("Summary of target processes for point-to-point", x)
  if(length(id.end) == 0){
    id.end <- max(which(x != ""))
  } else{
    id.end <- id.end - 2
  }

  x.sub <- x[id.start:id.end]

  # Check tab or space is being used.
  id.func <- grep("^\\t", x.sub)
  if(length(id.func) > 0){
    head.string <- "^\\t"
  } else{
    head.string <- " .*"
  }
  
  # Get functions
  id.func <- grep(head.string, x.sub, invert = TRUE)
  
  j <- 0
  for(i in id.func){
    j <- j + 1
    ret[[j]] <- list(Routine = NULL,
                     Calls = NULL,
                     Time = NULL,
                     Data.Sent = NULL,
                     SyncTime = NULL)
    
    # Get MPI routine name
    k <- i
    ret[[j]]$Routine <- gsub("^(.*):", "\\1", x.sub[k])
    
    # Get Calls if any
    k <- k + 1
    if(length(grep(paste(head.string, "Calls", sep = ""), x.sub[k])) > 0){
      tmp <- strsplit(x.sub[k], ":") 
      tmp <- strsplit(tmp[[1]][2], " ") 
      id.tmp <- which(tmp[[1]] != "")
      ret[[j]]$Calls <- as.integer(tmp[[1]][id.tmp[1]])
    }
    
    # Get Time if any
    k <- k + 1
    if(length(grep(paste(head.string, "Time", sep = ""), x.sub[k])) > 0){
      tmp <- strsplit(x.sub[k], ":") 
      tmp <- strsplit(tmp[[1]][2], " ") 
      id.tmp <- which(tmp[[1]] != "")
      ret[[j]]$Time <- as.double(tmp[[1]][id.tmp[1]])
    }
    
    # Get Data Sent if any
    k <- k + 1
    if(length(grep(paste(head.string, "Data Sent", sep = ""), x.sub[k])) > 0){
      tmp <- strsplit(x.sub[k], ":") 
      tmp <- strsplit(tmp[[1]][2], " ") 
      id.tmp <- which(tmp[[1]] != "")
      ret[[j]]$Data.Sent <- as.integer(tmp[[1]][id.tmp[1]])
    }
    
    # Get SyncTime if any
    k <- k + 1
    if(length(grep(paste(head.string, "SyncTime", sep = ""), x.sub[k])) > 0){
      tmp <- strsplit(x.sub[k], ":") 
      tmp <- strsplit(tmp[[1]][2], " ") 
      id.tmp <- which(tmp[[1]] != "")
      ret[[j]]$SyncTime <- as.double(tmp[[1]][id.tmp[1]])
    }
  }

  # Drop empty finction calls.
  new.ret <- ret
  for(j in 1:length(ret)){
    if(ret[[j]]$Routine == ""){
      new.ret[[j]] <- NULL
    }
  }
  
  # Cast return as dataframe
  ret <- parsed_fpmpi_2_df(new.ret)
  
  return( ret )
} # End of parse.prof.fpmpi().



parse.prof.mpip <- function(x, ...){
  #Rscript for profiling mpiP
  
  ret_mpip <- list()

  #lines <- readLines("12-null-null.exe.1.24393.1.mpiP")
  lines <- x 
  #selection the region between ---- and --- putting it in time series space using embed
  regions <- t(t(embed(grep("@---", lines), 2)) + c(-2, 2)) 
  #mapply on set of regions
  ret_mpip <- mapply(function(start,stop) {
    #converting to character without having to worry about spaces and empty lines
    chunk <- paste(lines[start:stop], collapse = "\n")
    #resusing the chunk
    chunk <- gsub("Line Parent_Funct", "Line_Parent_Funct", chunk)
    #making a text file connection to read it as table since it follows pattern
    tc <- textConnection(chunk)
    #reading as table
    df <- read.table(tc, header = T)
    #closing connection
    close(tc)
    #returning the result
    ret_mpip <- df
    return(ret_mpip)
    #arguments of mapply
  }, regions[,2], regions[,1])

  return(ret_mpip)
} # End of parse.prof.mpip().



parse.prof.tau <- function(x, ...){
  stop("parse.prof for TAU is not implemented yet.")
} # End of parse.prof.tau().

