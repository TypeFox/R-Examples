ListDir <- function(dir.name="", dir.path="", print.curl=FALSE, 
                    shared.username=NULL, suppress.Warnings=FALSE,
                    show.hidden=FALSE) {
  # This function lists all files in the 'dir.name' contained in 'dir.path'
  #   A user can also list files shared with them from the shared user.
  #
  # Args:
  #   dir.name: Name of directory
  #   dir.path: Path to current directory
  #   shared.username: String, valid iPlant username with whom the object
  #     is being shared.
  #   print.curl: Prints the associated curl statement
  #   suppress.Warnings: Don't do any error checking (faster)
  #
  # Returns:
  #   Returns nothing unless an error if 'dir.name' does not exist
  if (is.null(shared.username)){
    web <- paste(rplant.env$weblist, rplant.env$user, sep="/")
  } else {
    web <- paste(rplant.env$weblist, shared.username, sep="/")
  }

  if (dir.path == ""){
    web <- paste(web, dir.name, sep="/")
  } else {
    web <- paste(web, dir.path, dir.name, sep="/")
  }
  # Check 'dir.name'
  Check(dir.name, dir.path, suppress.Warnings, shared.username) 

  if (print.curl){
    curl.string <- paste(rplant.env$first, " ", web, sep="")
    print(curl.string)
  }
  Renew()
  tmp <- tryCatch(expr  = fromJSON(getURL(web, curl=rplant.env$curl.call)), 
                  error = function(err) {return(paste(err))})
  if (!suppress.Warnings){Error(tmp)}
  nms <- NULL  #create names vector to parse hiddens
  type <- NULL
  for(i in sequence(length(tmp$result))){
    nms <- c(nms, tmp$result[[i]]$name)
    type <- c(type, tmp$result[[i]]$type)
  }
    
  toIgnore <- grep("^\\.+$", nms) # always ignore home directory
  whichHiddens <- grep("^[.]\\D+", nms)
  if(show.hidden)
    toIgnore <- union(toIgnore, whichHiddens)

  nms <- nms[-toIgnore]
  type <- type[-toIgnore]

  # This portion is necessary to weed out the artifact folders.  It probably
  #   won't be implemented because it is slow.
  newnms <- NULL
  newtype <- NULL
  for (i in 1:length(nms)) {
    path <- paste(dir.path, dir.name, nms[i], sep="/")
    dir.exist <- fromJSON(getURL(url  = paste(rplant.env$webcheck, path, sep=""), curl = rplant.env$curl.call))
    if (dir.exist$status != 'error') {
      newnms <- append(newnms, nms[i])
      newtype <- append(newtype, type[i])
    }
  }

  res <- matrix(nrow=length(new), ncol=2)
  colnames(res) <- c("name", "type")

  for (i in sequence(dim(res)[1])) {
    res[i, 1] <- newnms[i]
    res[i, 2] <- newtype[i]
  }
  return(res)
}
