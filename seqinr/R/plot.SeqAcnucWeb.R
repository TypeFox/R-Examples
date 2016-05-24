################################################################################
#
#                              plot.SeqAcnucWeb
#
################################################################################

plot.SeqAcnucWeb <- function(x, types = getType()$sname, socket = autosocket(), ...){
  verbose <- FALSE # a passer en argument si besoin est
  #
  # Check arguments:
  #
  if(!inherits(x, "SeqAcnucWeb")) stop("Sequence of class SeqAcnucWeb is needed")

  #
  # Save graphical parameters:
  #
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  if(verbose) cat(paste("types:", types, sep = "\n"))
  
  #
  # Get the parent sequence:
  #
  
  GiveMeTheMotherSequence <- paste("me n=", x, sep = "") 
  query(listname = "me", query = GiveMeTheMotherSequence, socket = socket)
  MotherLength <- as.numeric(getLength(get("me", .seqinrEnv)$req[[1]]))
  MotherName <- get("me", .seqinrEnv)$req[[1]]
  if(verbose) cat("\nMotherLength = ", MotherLength)

  #
  # Plot organization:
  #
  par(mar = c(2.1, 0.1, 4.1, 0.1), lend = "square", ljoin = "mitre")
  cx <- c(0, MotherLength)
  cy <- c(0, 1)
  plot(cx, cy, ann = FALSE, type = "n", axes = FALSE)
  axis(1)
  title(main = paste("Physical position of subsequences on the parent sequence",
          MotherName, "(", MotherLength, "bp )", sep=" "))

  #
  # Look for subsequences:
  #
 
  GiveMeAllSubsequences <- paste("fi n=", MotherName, sep = "")
  query(listname = "filles", query = GiveMeAllSubsequences, socket = socket)

  n <- length(types) # number of potential subsequences types
  ispresent <- rep(FALSE, n) # will be TRUE if one or more subsequence of this type is found
  nb <- numeric(n) # count of subsequences for available types
  posi <- vector(mode = "list", length = n) # position of subsequences
    
  for(i in seq_len(n)){
    q <- paste("filles et t=", types[i], sep = "")
    if(verbose) cat("query = ", q, "\n")

    result <- try(query(socket = socket, listname = "tmp", query = q))
    if( inherits(result, "try-error")) next
    if(get("tmp", .seqinrEnv)$nelem == 0) next
    if(is.na(get("tmp", .seqinrEnv)$req[[1]])) next
    if(get("tmp", .seqinrEnv)$req[[1]] == x ) next

    ispresent[i] <- TRUE
      
    u <- lapply(get("tmp", .seqinrEnv)$req, getLocation)
    names(u) <- get("tmp", .seqinrEnv)$req
    nb[i] <- length(u)
    posi[[i]] <- u
  }

  #
  # Draw subsequences:
  #
  posi <- posi[ispresent]
  nb <- nb[ispresent]
  types <- types[ispresent]
  n <- length(types)
  for(i in seq_len(n)){
    for(j in seq_len(length(posi[[i]]))){
      xleft <- posi[[i]][[j]][1]
      ybottom <- (i - 1)/(n + 1)
      xright <- posi[[i]][[j]][2]
      ytop <- i/(n + 1)
      rect(xleft, ybottom, xright, ytop, col = i, border = "black", lend = "square", ljoin = "mitre" )
    }
  }
    
  #
  # Draw legend:
  #
  legend("top", legend = paste(types, "(", nb, ")", sep = ""), fill = seq_len(n), horiz = TRUE, bty = "n")

  resu <- lapply(posi,function(x){lapply(x,unlist)})
  names(resu) <- types

  #
  #  workspace cleanup
  #
  
  rm("me", pos = .seqinrEnv)
  rm("filles", pos = .seqinrEnv)
  rm("tmp", pos = .seqinrEnv)
  
  #
  # Return invisibly the result:
  #
  invisible(resu)
}
