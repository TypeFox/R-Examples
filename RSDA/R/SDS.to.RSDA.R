SDS.to.RSDA <-
function(file.path, labels=FALSE) {
  
  get.last <- function(X) {
    return (X[[length(X)]])
  }
  
  regex.preprocess <- function (SdsData) {
    SdsData <- paste(SdsData, collapse=" ") # It joins the lines of the file with whitespaces. The processing starts using one string with all the data.
    SdsData <- gsub(pattern='\\\\"', replacement="'", x=SdsData, perl=TRUE) # Replaces escaped double quotes with single quotes.
    SdsData <- gsub(pattern=' ?(=|\\(|\\)|,|\\:)(?=(?:(?:[^"]*"){2})*[^"]*$) ?', replacement=" \\1 ", x=SdsData, perl=TRUE) # It isolates certain special characters like: (, ), : and = from the adjacent strings, so we can later decompose the string in separate tokens (we split the string using whitespaces).
    SdsData <- gsub(pattern='(?:(\\d+)\\s*\\(\\s*((?:1|0)(?:\\.\\d+)?)\\s*\\)\\s*)(?=(?:(?:[^"]*"){2})*[^"]*$)', replacement=' "\\1" = \\2 ', x=SdsData, perl=TRUE) # This converts the data of histograms variables in RECTANGLE_MATRIX to valid JSON objects.
    SdsData <- gsub(pattern='(\\b\\w+\\b) =(?=(?:(?:[^"]*"){2})*[^"]*$)', replacement='"\\1" =', x=SdsData, perl=TRUE) # It quotes the names for properties on objects.
    SdsData <- gsub(pattern='"VAR"\\s*=\\s*(\\d+)\\s*,\\s*(\\d+)(?=(?:(?:[^"]*"){2})*[^"]*$)', replacement='\\1 , \\2 , ', x=SdsData, perl=TRUE) # This convert the data in hierarchies into valid JSON data.
    SdsData <- strsplit(x=SdsData, split='\\s+(?=(?:(?:[^"]*"){2})*[^"]*$)', perl=TRUE) # We split the string into separate tokens.
    SdsData <- sub(pattern='^:$', replacement=",", x=SdsData[[1]], perl=TRUE) # We replace the ":" characters (on intervals data) with commas.
    SdsData <- sub(pattern='^=$', replacement=":", x=SdsData, perl=TRUE) # We replace the "=" characters with ":".
    SdsData <- sub(pattern='^(?:NA|NU)$', replacement="null", x=SdsData, perl=TRUE) # We replace NA and NU with nulls.
    SdsData <- sub(pattern='^\\.(\\d+)$', replacement="0.\\1", x=SdsData, perl=TRUE) # We fill numbers that do not have the integer part (for example we replace .12 with 0.12).
    SdsData <- sub(pattern='^(FILES|HEADER|INDIVIDUALS|VARIABLES|RECTANGLE_MATRIX|DIST_MATRIX|HIERARCHIE|RULES|proba|inter_cont(?:inue)?|nominal|continue?|mult_nominal(?:_Modif)?)$', replacement='"\\1"', x=SdsData, perl=TRUE) # We quote certain special words like: INDIVIDUALS, nominal, continue, among others.
    SdsData <- append(x=SdsData, values="(", after=0) # We add as a first token a "(" character. So it marks the beginning of the data.
    SdsData[length(SdsData)] <- ")" # We replace the final END token with a ")" character, marking the end of the data.
    return (SdsData)
  }
  
  # See step 2 above. The function replaces the () pairs of characters with [] or {} according to the context.
  preprocessed.sds.to.json <- function (SdsData) {
    pStack <- list() # We have a stack that keeps tracking the parenthesis open, but that are not closed.
    for(i in 1:length(SdsData)) { # We iterate over every token in the string, doing the proper replacements according to the context.
      currentToken <- SdsData[i]
      switch (currentToken,
              '(' = {
                pStack[length(pStack) + 1] <- i
              },
              ')' = {
                if (pStack[[length(pStack)]] == "[") {
                  SdsData[i] <- "]"  
                } else if (pStack[[length(pStack)]] == "{") {
                  SdsData[i] <- "}"
                } else {
                  SdsData[pStack[[length(pStack)]]] <- "["
                  SdsData[i] <- "]" 
                }
                pStack[[length(pStack)]] <- NULL
              },
              ':' = {
                if (mode(pStack[[length(pStack)]]) == "numeric") {
                  SdsData[pStack[[length(pStack)]]] <- "{"
                  pStack[[length(pStack)]] <- "{"
                }  
              },
              ',' = {
                if (mode(pStack[[length(pStack)]]) == "numeric") {
                  SdsData[pStack[[length(pStack)]]] <- "["
                  pStack[[length(pStack)]] == "["
                } 
              })
    }
    SdsData <- paste(SdsData, collapse=" ")
    return (fromJSON(SdsData)) # We return the parsed JSON data, ready to be processed. This corresponds to step 3.
  }
  
  process.continue.variable <- function (number.of.rows, data, variable.index, variable.name) {
    aux <- list()
    aux[[1]] <- rep("$C", number.of.rows)
    aux[[2]] <- sapply(X=data$SODAS$RECTANGLE_MATRIX, FUN=function (dat.ind) {
      if (is.null(dat.ind[[variable.index]]))
        return (NA)
      else
        return (round(dat.ind[[variable.index]], 3))
    })
    
    aux <- data.frame(aux)
    colnames(aux) <- c("$C", make.names(names=variable.name))
    return (aux)
  }
  
  process.inter.cont.variable <- function (number.of.rows, data, variable.index, variable.name) {
    aux <- list()
    aux[[1]] <- rep("$I", number.of.rows)
    
    aux[[2]] <- sapply(X=data$SODAS$RECTANGLE_MATRIX, FUN=function (dat.ind) {
      if (is.null(dat.ind[[variable.index]]))
        return (NA)
      else
        return (round(dat.ind[[variable.index]][[1]], 3))
    })
    aux[[3]] <- sapply(X=data$SODAS$RECTANGLE_MATRIX, FUN=function (dat.ind) {
      if (is.null(dat.ind[[variable.index]]))
        return (NA)
      else
        return (round(dat.ind[[variable.index]][[2]], 3))
    })
    
    aux <- data.frame(aux)
    colnames(aux) <- c("$I", make.names(names=variable.name), make.names(names=variable.name))
    return (aux)
  }
  
  process.nominal.variable <- function (labels, number.of.rows, data, variable.index, variable.name) {
    aux <- list()
    aux[[1]] <- rep("$S", number.of.rows)
    
    categories <- sapply(X=get.last(data$SODAS$VARIABLES[[variable.index]]), FUN=function (cat) { cat[[3]] })
    categories.labels <- sapply(X=get.last(data$SODAS$VARIABLES[[variable.index]]), FUN=function (cat) { cat[[2]] })
    
    aux[[2]] <- rep(length(categories), number.of.rows)
    
    categories.data <- sapply(X=data$SODAS$RECTANGLE_MATRIX, FUN=function(dat.ind) {
      if (is.null(dat.ind[[variable.index]]))
        return (rep(x=NA, times=length(categories)))
      else {
        categories.for.individual <- rep(x=0, times=length(categories))
        present.categories <- dat.ind[[variable.index]]
        categories.for.individual[present.categories] <- 1
        return (categories.for.individual)
      }
    })
    
    aux <- data.frame(aux, as.data.frame(t(x=categories.data)))
    if (labels) {
      colnames(aux) <- c("$S", make.names(names=variable.name), make.names(names=categories.labels))
    } else {
      colnames(aux) <- c("$S", make.names(names=variable.name), make.names(names=categories))
    }
    return (aux)
  }
  
  process.mult.nominal.modif.variable <- function (labels, number.of.rows, data, variable.index, variable.name) {
    aux <- list()
    aux[[1]] <- rep("$H", number.of.rows)
    
    categories <- sapply(X=get.last(data$SODAS$VARIABLES[[variable.index]]), FUN=function (cat) { cat[[3]] })
    categories.labels <- sapply(X=get.last(data$SODAS$VARIABLES[[variable.index]]), FUN=function (cat) { cat[[2]] })
    
    aux[[2]] <- rep(length(categories), number.of.rows)
    
    for (i in 1:length(categories)) {
      aux[[i + 2]] <- sapply(X=data$SODAS$RECTANGLE_MATRIX, FUN=function(dat.ind) {
        val <- dat.ind[[variable.index]][as.character(i)]
        if (is.null(val))
          return (NA)
        else
          return (ifelse(test=is.na(val), yes=0.000, no=round(val, 3)))
      })
    }
    
    aux <- data.frame(aux)
    if (labels) {
      colnames(aux) <- c("$H", make.names(names=variable.name), make.names(names=categories.labels))
    } else {
      colnames(aux) <- c("$H", make.names(names=variable.name), make.names(names=categories))
    }
    return (aux)
  }
  
  # --------------------
  # Main function Logic
  # --------------------
  
  data <- readLines(con=file.path, warn=FALSE)
  
  cat("Preprocessing file\n")
  data <- regex.preprocess(data)
  
  cat("Converting data to JSON format\n")
  data <- preprocessed.sds.to.json(data)
  
  if (labels) {
    sym.obj.names <- sapply(X=data$SODAS$INDIVIDUALS, FUN=function (ind) { ind[[2]] })
  } else {
    sym.obj.names <- sapply(X=data$SODAS$INDIVIDUALS, FUN=function (ind) { ind[[3]] })
  }
  
  variables.names <- sapply(X=data$SODAS$VARIABLES, FUN=function (var) { var[[5]] })
  variables.types <- sapply(X=data$SODAS$VARIABLES, FUN=function (var) { var[[2]] })
  result <- data.frame(row.names = make.names(names=sym.obj.names, unique=TRUE))
  number.of.rows <- nrow(result)
  
  for (i in 1:length(variables.types)) {
    
    cat(paste0("Processing variable ", i, ": ", variables.names[[i]],"\n"))
    
    switch (variables.types[[i]],
            'inter_continue' = {
              result <- cbind(result, process.inter.cont.variable(number.of.rows, data, i, variables.names[[i]]))
            },
            'inter_cont' = {
              result <- cbind(result, process.inter.cont.variable(number.of.rows, data, i, variables.names[[i]]))
            },
            'continue' = {
              result <- cbind(result, process.continue.variable(number.of.rows, data, i, variables.names[[i]]))
            },
            'continu' = {
              result <- cbind(result, process.continue.variable(number.of.rows, data, i, variables.names[[i]]))
            },
            'nominal' = {
              result <- cbind(result, process.nominal.variable(labels, number.of.rows, data, i, variables.names[[i]]))
            },
            'mult_nominal' = {
              result <- cbind(result, process.nominal.variable(labels, number.of.rows, data, i, variables.names[[i]]))
            },
            'mult_nominal_Modif' = {
              result <- cbind(result, process.mult.nominal.modif.variable(labels, number.of.rows, data, i, variables.names[[i]]))
            },
            cat(paste0("Variable type not supported:"), variables.types[[i]], "\n"))
  }
  
  return (newSobject(result))
}
