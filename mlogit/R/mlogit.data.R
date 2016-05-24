mlogit.data <- function(data, choice, shape = c("wide","long"), varying = NULL,
                        sep = ".", alt.var = NULL, chid.var = NULL, 
                        alt.levels = NULL, id.var = NULL, opposite = NULL,
                        drop.index = FALSE, ranked = FALSE, ...){
  # chid.name, alt.name : the name of the index variables
  # chid, alt : the index variables
  if (shape == "long"){
    if (is.null(chid.var)){
      chid.name <- "chid"
      chid.is.variable <- FALSE
    }
    else{
      chid.name <- chid.var
      chid.is.variable <- ifelse(is.null(data[[chid.var]]), FALSE, TRUE)
    }
    if (!ranked){
      choice.name <- choice
      choice <- data[[choice]]
    }
    if (is.null(alt.var) && is.null(alt.levels))
      stop("at least one of alt.var and alt.levels should be filled")
    
    if (!is.null(alt.levels)){
      J <- length(alt.levels)
      n <- nrow(data)/J
      alt <- factor(rep(alt.levels, n), levels = alt.levels)
      if (!is.null(alt.var) && !is.null(data[[alt.var]])){
        warning(paste("variable",alt.var,"exists and will be replaced"))
        alt.is.variable <- TRUE
      }
      else alt.is.variable <- FALSE
      alt.name <- ifelse(is.null(alt.var), "alt", alt.var)
    }
    else{
      alt.name <- alt.var
      alt.is.variable <- TRUE
      if (!is.factor(data[[alt.name]])) data[[alt.name]] <- factor(data[[alt.name]])
      alt.levels <- levels(data[[alt.name]])
      J <- length(alt.levels)
      alt <- data[[alt.name]]
    }
    n <- nrow(data) / J
    if (!chid.is.variable) chid <- rep(1:n, each = J) else chid <- data[[chid.name]]
    if (!ranked){
      if (!is.logical(data[[choice.name]])){
        if (is.factor(choice) && 'yes' %in% levels(choice))
          data[[choice.name]] <- data[[choice.name]] == 'yes'
        if (is.numeric(choice)) data[[choice.name]] <- data[[choice.name]] != 0
      }
    }
    # remplacer id par chid à gauche
    chid <- as.factor(chid)
    alt <- as.factor(alt)
    row.names(data) <- paste(chid, alt, sep = ".")
  }
  
  if (shape == "wide"){
    if (!ranked){
      choice.name <- choice
      if (is.ordered(data[[choice]])) class(data[[choice]]) <- "factor"
      else data[[choice]] <- as.factor(data[[choice]])
    }
    # this doesn't work for ordered factors which remains ordered
    if (is.null(alt.var)) alt.name <- "alt" else alt.name <- alt.var
    if (is.null(chid.var)) chid.name <- "chid" else chid.name <- chid.var
    if (!is.null(varying)){
      data <- reshape(data, varying = varying, direction = "long", sep = sep,
                      timevar = alt.name, idvar = chid.name,  ...)
    }
    else{
      if (ranked)
        stop("for ranked data in wide format, the varying argument is mandatory and should contain at least the choice variable")
      id.names <- as.numeric(rownames(data))
      nb.id <- length(id.names)
      data[[chid.name]] <- id.names
      lev.ch <- levels(data[[choice]])
      data <- data.frame(lapply(data, rep, length(lev.ch)))
      data[[alt.name]] <- rep(lev.ch, each = nb.id)
      row.names(data) <- paste(data[[chid.name]], data[[alt.name]], sep = ".")
    }
    data <- data[order(data[[chid.name]], data[[alt.name]]), ]
    # remplacer id par chid à gauche
    chid <- as.factor(data[[chid.name]])
    alt <- as.factor(data[[alt.name]])
    if (!is.null(alt.levels)){
      levels(data[[choice]]) <- alt.levels
      levels(alt) <- alt.levels
      row.names(data) <- paste(chid, alt, sep = ".")
    }
    if (!ranked) data[[choice]] <- data[[choice]] == alt
    else{
      if (is.null(data[[choice]])) stop("the choice variable doesn't exist")
    }
  }
  chidpos <- which(names(data) == chid.name)
  altpos <- which(names(data) == alt.name)
  if (!is.null(id.var)){
    idpos <- which(names(data) == id.var)
    id.var <- as.factor(data[[id.var]])
  }
  if (drop.index){
    if (!is.null(id.var)) data <- data[, -c(chidpos, altpos, idpos)]
    else data <- data[, -c(chidpos, altpos)] 
  }
  
  if (!is.null(opposite)){
    for (i in opposite){
      data[[i]] <- - data[[i]]
    }
  }
  index <- data.frame(chid = chid, alt = alt)
  if (!is.null(id.var)) index <- cbind(index, id = id.var)
  rownames(index) <- rownames(data)
  attr(data, "index") <- index
  attr(data, "class") <- c("mlogit.data", "data.frame")
  if (ranked) data <- mlogit2rank(data, choicename = choice)
  if (!ranked) attr(data, "choice") <- choice.name
  data
}


mlogit2rank <- function(x, choicename, ...){
  choicepos <- match(choicename, names(x))
  id <- attr(x, "index")$chid
  lev.id <- levels(id)
  theid <- as.numeric(id)
  oalt <-  attr(x, "index")$alt
  lev.alt <- levels(oalt)
  choice <- x[[choicename]]
  J <- length(unique(choice))
  d <- data.frame()
  chid <- c()
  alt <- c()
  id <- c()
  k <- 0
  for (i in unique(theid)){
    aid <- which(theid == i)
    adata <- x[aid, - choicepos]
    achoice <- choice[aid]
    aalt <- oalt[aid]
    remAlts <- rep(TRUE, J)
    alogchoice <- achoice == 1
    d <- rbind(d, cbind(adata, alogchoice))
    Z <- sum(remAlts)
    k <- k + 1
    chid <- c(chid, rep(k, Z))
    id <- c(id, rep(i, Z))
    alt <- c(alt, aalt)
    for (j in 1:(J - 2)){
      k <- k + 1
      min.index <- achoice == j
      remAlts[min.index] <- FALSE
      Z <- sum(remAlts)
      chid <- c(chid, rep(k, Z))
      alt <- c(alt, aalt[remAlts])
      id <- c(id, rep(i, Z))
      alogchoice <- achoice[remAlts] == j + 1
      d <- rbind(d, cbind(adata[remAlts,], alogchoice))
    }
  }
  colnames(d)[length(d)] <- choicename
  alt <- factor(alt, labels = lev.alt)
  index <- data.frame(chid = chid, alt = alt, id = id)
  rownames(d) <- rownames(index) <- paste(chid, as.character(alt), sep = ".")
  structure(d, index = index, class = c('mlogit.data', 'data.frame'))
}


"[.mlogit.data" <- function(x, i, j, drop = TRUE){
  mydata <- `[.data.frame`(x, i, j, drop = drop)
  index <- "[.data.frame"(attr(x, "index"), i,)
  index <- data.frame(lapply(index, function(x) x[drop = TRUE]), row.names = rownames(mydata))
  
  if (is.null(dim(mydata))){
    structure(mydata,
              index = index,
              class = c("pseries", class(mydata))
              )
  }
  else{
    structure(mydata,
              index = index,
              class = c("mlogit.data", "data.frame"))
  }
}

print.mlogit.data <- function(x, ...){
  attr(x, "index") <- NULL
  class(x) <- "data.frame"
  print(x, ...)
}

"[[.mlogit.data" <- function(x, y){
  index <- attr(x, "index")
  attr(x, "index") <- NULL
  class(x) <- "data.frame"
  result <- x[[y]]
  if (!is.null(result)){
    result <- structure(result,
                        class = c("pseries", class(x[[y]])),
                        index = index,
                        names = row.names(x)
                        )
  }
  result
}  

"$.mlogit.data" <- function(x,y){
  "[["(x, paste(as.name(y)))
}

"$<-.mlogit.data" <- function(object, x, value){
  # object : le data.frame
  # x : la variable
  # value : la nouvelle valeur
  object[[x]] <- value
  object
}

"[[<-.mlogit.data" <- function(object, x, value){
  if (class(value)[1] == "pseries"){
    class(value) <- class(value)[-1]
    attr(value, "index") <- NULL
  }
  object <- "[[<-.data.frame"(object, x, value = value)
  object
}



print.pseries <- function(x, ...){
  attr(x, "index") <- NULL
  attr(x, "class") <- attr(x, "class")[-1]
  if (length(attr(x, "class")) == 1
      && class(x) %in% c("character", "logical", "numeric"))
    attr(x, "class") <- NULL
  print(x, ...)
}

index.mlogit.data <- function(x, ...){
  attr(x, "index")
}
