
read.tree2<-function (file = "", text = NULL, tree.names = NULL, skip = 0, 
                      comment.char = "#", keep.multi = FALSE, ...) 
{
  read.tree <-function (file = "", text = NULL, tree.names = NULL, skip = 0, 
                        comment.char = "#", keep.multi = FALSE, ...) 
  {
    unname <- function(treetext) {
      nc <- nchar(treetext)
      tstart <- 1
      while (substr(treetext, tstart, tstart) != "(" && tstart <= 
               nc) tstart <- tstart + 1
      if (tstart > 1) 
        return(c(substr(treetext, 1, tstart - 1), substr(treetext, tstart, nc)))
      return(c("", treetext))
    }
    if (!is.null(text)) {
      if (!is.character(text)) 
        stop("argument `text' must be of mode character")
      tree <- text
    }
    else {
      tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE, 
                   skip = skip, comment.char = comment.char, ...)
    }
    if (identical(tree, character(0))) {
      warning("empty character string.")
      return(NULL)
    }
    tree <- gsub("[ \t]", "", tree)
    tree <- unlist(strsplit(tree, NULL))
    y <- which(tree == ";")
    Ntree <- length(y)
    x <- c(1, y[-Ntree] + 1)
    if (is.na(y[1])) 
      return(NULL)
    STRING <- character(Ntree)
    for (i in 1:Ntree) STRING[i] <- paste(tree[x[i]:y[i]], sep = "", 
                                          collapse = "")
    tmp <- unlist(lapply(STRING, unname))
    tmpnames <- tmp[c(TRUE, FALSE)]
    STRING <- tmp[c(FALSE, TRUE)]
    if (is.null(tree.names) && any(nzchar(tmpnames))) 
      tree.names <- tmpnames
    colon <- grep(":", STRING)
    if (!length(colon)) {
      obj <- lapply(STRING, clado.build)
    }
    else if (length(colon) == Ntree) {
      obj <- lapply(STRING, tree.build)
    }
    else {
      obj <- vector("list", Ntree)
      obj[colon] <- lapply(STRING[colon], tree.build)
      nocolon <- (1:Ntree)[!1:Ntree %in% colon]
      obj[nocolon] <- lapply(STRING[nocolon], clado.build)
    }
    for (i in 1:Ntree) {
      ROOT <- length(obj[[i]]$tip.label) + 1
      if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 
            1) 
        stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading Newick file aborted at tree no.", 
                   i))
    }
    if (Ntree == 1 && !keep.multi) 
      obj <- obj[[1]]
    else {
      if (!is.null(tree.names)) 
        names(obj) <- tree.names
      class(obj) <- "multiPhylo"
    }
    obj
  }
  unname <- function(treetext) {
    nc <- nchar(treetext)
    tstart <- 1
    while (substr(treetext, tstart, tstart) != "(" && tstart <= 
             nc) tstart <- tstart + 1
    if (tstart > 1) 
      return(c(substr(treetext, 1, tstart - 1), substr(treetext, tstart, nc)))
    return(c("", treetext))
  }
  if (!is.null(text)) {
    if (!is.character(text)) 
      stop("argument `text' must be of mode character")
    tree <- text
  }
  else {
    tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE, 
                 skip = skip, comment.char = comment.char, ...)
  }
  if (identical(tree, character(0))) {
    warning("empty character string.")
    return(NULL)
  }
  dollarSignReader <- function(text)
  {
    text <- gsub("[ \t]", "", text)
    text <- paste(text, sep="", collapse ="")
    y<- text
    y<- gsub("\n","", y)
    y<- gsub(" " , "" ,y)
    y<- unlist(strsplit(y, ";"))
    ny <-length(y)
    for(a in 1:ny)
    {
      x<- y[a]
      nx<-nchar(x)
      x<-strsplit(x, "")
      x<- unlist(x)
      x[nx +1] <- ","
      for(i in 1:nx)
      {
        isThereDollar <- 0
        if(x[i] == "(" || x[i] == "," || x[i] == ")")
        {
          for(j in (i+1):(nx)) 
          {
            if(x[j] == ":")
            {
              for(k in i:j)
              {
                if(x[k] == "$")
                {
                  isThereDollar <-1
                }
              }  
              if(isThereDollar==0)
                x[j] = "$mu0:"
              break
            }
          }
          
        }
      }
      x<- paste(x, sep="", collapse ="")
      nx<-nchar(x)
      x<-strsplit(x, "")
      x<- unlist(x)
      for(i in 1:nx)
      {
        if(x[i] == ":")
        {
          for(j in (i+1):(nx+1)) 
          {
            if(x[j] == "," || x[j] == ")"  )
            {
              for(k in i:(j-1))
              {
                x[k] = ""
              }
              break
            } 
          }
          
        }
      }
      x<- x[1:(nx-1)]
      x[length(x)+1] = ";" 
      x<-paste(x , collapse = "", sep = "")
      x<-gsub("\\$" , ":", x)
      y[a] <-x
    }
    y<- paste(y , collapse = "")
    y<- unlist(strsplit(y, ";"))
    ny <-length(y)
    
    tempchar2 <- vector()
    tempvec2 <- vector()
    for(a in 1:ny)
    {
      x<- y[a]
      nx<-nchar(x)
      x<-strsplit(x, "")
      x<- unlist(x)
      x[nx +1] <- ","
      x[nx +1] <- ","
      count = 1
      for(i in 1:nx)
      {
        if(x[i] == ":")
          for(j in (i+1):nx)
          {
            if(x[j] == "," || x[j] ==")" )
            {
              tempchar <- vector()
              tempvec <- vector()
              kk <- j
              for(k in (i+1):(j-1))
              {
          
                if(x[k] == "=")
                {
                  #kk <- k
                }
                if( k <kk)
                {
                tempchar[k] <- x[k]
                }
                if( k >kk)
                tempvec[k] <-x[k]
                x[k] <- ""
              }
              match<- paste( tempchar[!is.na(tempchar)], collapse ="", sep = "")
              matchvec <- paste( tempvec[!is.na(tempvec)], collapse ="", sep = "")

              if(match %in% tempchar2)
              {
                index <- which(tempchar2 ==  match )
                if(!is.na(tempvec2[index]) && !identical(tempvec2[index], matchvec))
                  stop("The values assigned to ", tempchar2[index], "'s do not match.")
                x[i] <- paste(":",index, sep="", collapse ="")  	
              }
              
              if(!match %in% tempchar2)
              {
                tempchar2[count]<- paste( tempchar[!is.na(tempchar)], collapse ="", sep = "")
                tempvec2[count] <- paste( tempvec[!is.na(tempvec)], collapse ="", sep = "")
                x[i] <- paste(":", count, sep="", collapse = "")
                count <- count+1
              }
              break
            }
          }
      }
      x[length(x)] <- ";"
      x<-paste(x , collapse = "", sep = "")
      y[a] <-x
    }
    y<- paste(y , collapse = "")
    k<- list()


    k[[1]] <- y
    k[[2]] <- tempchar2
    k[[3]] <- as.numeric(tempvec2)
    return(k)
  }
  dollarSignRemover <- function(text)
  {
    text <- gsub("[ \t]", "", text)
    text <- paste(text, sep="", collapse ="")
    y<- text
    y<- gsub("\n","", y)
    y<- gsub(" ","", y)
    y<- unlist(strsplit(y, ";"))
    ny <-length(y)
    
    for(a in 1:ny)
    {
      x<- y[a]
      nx<-nchar(x)
      x<-strsplit(x, "")
      x<- unlist(x)
      x[nx +1] <- ":"
      
      for(i in 1:nx)
      {
        if(x[i] == "$")
        {
          for(j in (i+1):(nx+1)) 
          {
            if(x[j] == ":" )
            {
              for(k in i:(j-1))
              {
                x[k] = ""
              }
              break
            } 
          }
        }
      }
      x<- x[1:nx]
      x[length(x)+1] = ";" 
      x<-paste(x , collapse = "", sep = "")
      x<-gsub("\\$" , ":", x)
      y[a] <-x
    }
    y<- paste(y , collapse = "")
    return(y)
  }  
  obj2 <- read.tree(text = dollarSignReader((tree) )[[1]] )  
  tempchars <- dollarSignReader(tree)[[2]]
  tempvalues <- dollarSignReader(tree)[[3]]
  tree<- dollarSignRemover(tree)
  tree <- gsub("[ \t]", "", tree)
  tree <- unlist(strsplit(tree, NULL))
  y <- which(tree == ";")
  Ntree <- length(y)
  x <- c(1, y[-Ntree] + 1)
  if (is.na(y[1])) 
    return(NULL)
  STRING <- character(Ntree)
  for (i in 1:Ntree) STRING[i] <- paste(tree[x[i]:y[i]], sep = "",  collapse = "")
  tmp <- unlist(lapply(STRING, unname))
  tmpnames <- tmp[c(TRUE, FALSE)]
  STRING <- tmp[c(FALSE, TRUE)]
  if (is.null(tree.names) && any(nzchar(tmpnames))) 
    tree.names <- tmpnames
  colon <- grep(":", STRING)
  if (!length(colon)) {
    obj <- lapply(STRING, clado.build)
  }
  else if (length(colon) == Ntree) {
    obj <- lapply(STRING, tree.build)
  }
  else {
    obj <- vector("list", Ntree)
    obj[colon] <- lapply(STRING[colon], tree.build)
    nocolon <- (1:Ntree)[!1:Ntree %in% colon]
    obj[nocolon] <- lapply(STRING[nocolon], clado.build)
  }
  for (i in 1:Ntree) {
    ROOT <- length(obj[[i]]$tip.label) + 1
    if (sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1) 
      stop(paste("The tree has apparently singleton node(s): cannot read tree file.\n  Reading Newick file aborted at tree no.",i))
  }
  if (Ntree == 1 && !keep.multi) 
    obj <- obj[[1]]
  else {
    if (!is.null(tree.names)) 
      names(obj) <- tree.names
    class(obj) <- "multiPhylo"
  }
  obj$dollarData <- obj2
  obj$tempchars<- tempchars 
  obj$tempvalues <- tempvalues
  return(obj)
}