#' Replace symbols in a character vector by other symbols
#' 
#' @param what vector of type character, the symbols to be replaced, e.g. c("A", "B")
#' @param by vector of type character, the replacement, e.g. c("x[0]", "x[1]")
#' @param x vector of type character, the object where the replacement should take place
#' @return vector of type character, conserves the names of x.
#' @examples replaceSymbols(c("A", "B"), c("x[0]", "x[1]"), c("A*B", "A+B+C"))
#' @export
replaceSymbols <- function(what, by, x) {
  
  xOrig <- x
  is.not.zero <- which(x!="0")
  x <- x[is.not.zero]
  
  mynames <- names(x)
  
  x.parsed <- parse(text = x, keep.source = TRUE)
  data <- utils::getParseData(x.parsed)
    
  by <- rep(by, length.out=length(what))
  names(by) <- what
  
  data$text[data$text%in%what] <- by[data$text[data$text%in%what]]
  data <- data[data$token!="expr",]
  
  
  
  breaks <- c(0, which(diff(data$line1) == 1), length(data$line1))
  
  out <- lapply(1:(length(breaks)-1), function(i) {
    
    paste(data$text[(breaks[i]+1):(breaks[i+1])], collapse="")
    
  })
    
  names(out) <- mynames
  out <- unlist(out)
  
  xOrig[is.not.zero] <- out
  
  return(xOrig)
  
  
}



#' Replace a binary operator in a string by a function
#' 
#' @param what character, the operator symbol, e.g. "^"
#' @param by character, the function string, e.g. "pow"
#' @param x vector of type character, the object where the replacement should take place
#' @return vector of type character
#' @examples replaceOperation("^", "pow", "(x^2 + y^2)^.5")
#' @export
replaceOperation <- function(what, by, x) {
  
  mynames <- names(x)
  
  if(is.null(x)) return()
  
  
  names(by) <- what
  if(length(what)>1) {
    cat("Only one operator can be replaced at a time\n")
    return()
  }
  
  containsWhat <- grep(what, x, fixed=TRUE)
  xOrig <- x
  
  x <- paste(x[containsWhat], collapse=";")
  
  exit <- FALSE
  if(x=="") exit <- TRUE
  
  
  
  while(!exit) {
    
    parsed <- parse(text=x, keep.source = TRUE)
    parsData <- utils::getParseData(parsed)
    pres <- parsData[parsData$terminal==TRUE,]
    
    
    #add negative signs to numeric constants
    signs <- pres[1:(nrow(pres)-1), "token"]
    numbers <- pres[2:nrow(pres), "token"]
    index <- which(signs == "'-'" & numbers == "NUM_CONST")
    if(length(index) > 0) {
      pres[index, "text"] <- paste0("-", pres[index+1, "text"])
      pres[index, c("line2", "col2", "parent", "token")] <- pres[index+1, c("line2", "col2", "parent", "token")]
      pres <- pres[-(index+1),]  
    }
    
    positions <- which(pres$text%in%what)
    
    if(1 %in% positions) {
      cat("Binary operator should not be at beginning of the line\n")
      return()
    }
    if(length(positions)==0) {
      x <- as.list(strsplit(x, ";", fixed=TRUE)[[1]])
      return(x)
    }
    
    symbols <- by[pres$text[positions]]
    parents <- cbind(pBefore = pres$parent[positions - 1], pAfter = pres$parent[positions+1])
    
    
    # Get positions of the affiliated two arguments
    getCols <- function(parent) {
      out <- lapply(parent, function(par) {
        entries <- which(pres$parent == par)
        m <- switch(length(entries),
                    "1" = {
                      if(pres$token[entries] == "SYMBOL_FUNCTION_CALL") {
                        cat("Functions in the exponent need to be set in parenteses. Returning NULL.")
                        stop
                      } else c(pres$col1[entries], pres$col2[entries])
                    },
                    "2" = if(pres$token[max(c(entries[1]-1, 1))] == "SYMBOL_FUNCTION_CALL"){c(pres$col1[entries[1]-1], pres$col2[entries[2]])}
                    else{c(pres$col1[entries[1]], pres$col2[entries[2]])})
        
        return(m)
        
      })
      return(do.call(rbind, out))
      
    }
    parentCols <- as.matrix(cbind(getCols(parents[,"pBefore"]), getCols(parents[,"pAfter"])))
    
    
    # Function to determine nesting
    getNesting <- function(M) {
      n <- dim(M)[1]
      mysample <- data.frame(number = 1:n, level = 1)
      l <- 0
      
      for(l in 1:1) {
        
        numbers <- mysample$number[mysample$level==l]
        if(length(numbers) <= 1) break
        
        for(i in numbers) {
          isContained <- any(as.logical((M[i,1] > M[numbers,1])*(M[i,4] < M[numbers,2]) + (M[i,1] > M[numbers,3])*(M[i,4] < M[numbers,4])))
          if(isContained) mysample$level[i] <- mysample$level[i] + 1
        }
        
      }
      
      return(mysample)
      
      
    }
    nesting <- getNesting(parentCols)
    
    # Restrict replacement to lowest level
    lowest <- which(nesting$level == 1)
    if(length(unique(nesting$level)) == 1) exit <- TRUE
    parents <- parents[lowest]
    symbols <- symbols[lowest]
    parentCols <- matrix(parentCols[lowest,], nrow=length(lowest))
    
    
    # Do replacement
    unaffIni <- c(1, parentCols[,4]+1)
    unaffEnd <- c(parentCols[,1]-1, nchar(x))
    nParts <- length(unaffIni)
    parts <- substr(rep(x, nParts), unaffIni, unaffEnd)
    
    replacements <- paste(symbols, "(",
                          substr(rep(x, nParts-1), parentCols[,1], parentCols[,2]),
                          ", ",
                          substr(rep(x, nParts-1), parentCols[,3], parentCols[,4]),
                          ")", sep="")
    
    
    all <- c()
    all[seq(2, 2*length(parts), by=2)-1] <- parts
    all[seq(2, 2*length(replacements), by=2)] <- replacements
    all <- paste(all, collapse="")
    
    x <- all
    
  }
  
  
  #lapply(1:length(x), function(i) system(paste("rm", filenames[i])))
  
  x <- strsplit(x, ";", fixed=TRUE)[[1]]
  
  xOrig[containsWhat] <- x
  
  names(xOrig) <- mynames
  
  return(xOrig)
  
  
}


#' Compute Jacobian of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param variables other variables, e.g. paramters, f depends on. If variables is
#' given, f is derived with respect to variables instead of \code{names(f)}
#' @return named vector of type character with the symbolic derivatives
#' @examples 
#' jacobianSymb(c(A="A*B", B="A+B"))
#' jacobianSymb(c(x="A*B", y="A+B"), c("A", "B"))
#' @export
#' @importFrom stats D splinefun
jacobianSymb <- function(f, variables=NULL) {
  
  if(is.null(variables)) variables <- names(f)
  jacnames <- apply(expand.grid.alt(names(f), variables), 1, paste, collapse=".")
  jacobian <- matrix("0", nrow = length(f), ncol = length(variables))
  
  # Determine possible non-zero elements of the jacobian
  out <- lapply(variables, function(v) which(grepl(v, f)))
  inz <- do.call(rbind, lapply(1:length(variables), function(j) if(length(out[[j]])>0) cbind(i = out[[j]], j = j)))
  
  
  # Commpute derivatives for potential non-zero elements
  if(!is.null(inz)) {
    for(k in 1:dim(inz)[1]) {
      
      i <- inz[k, 1]
      j <- inz[k, 2]
      
      myeq <- parse(text = f[i])
      myvar <-variables[j]
      myderiv <- paste(deparse(D(myeq, myvar)), collapse="")
      
      jacobian[i, j] <- myderiv
      
    }
  }
    
  out <- as.vector(jacobian)
  out <- gsub(" ", "", out, fixed=TRUE)
  out <- gsub("++", "+", out, fixed=TRUE)
  out <- gsub("--", "+", out, fixed=TRUE)
  out <- gsub("-+", "-", out, fixed=TRUE)
  out <- gsub("+-", "-", out, fixed=TRUE)
  
  names(out) <- jacnames
  
  return(out)
  
}

#' Get symbols from a character
#' 
#' @param char Character vector (e.g. equation)
#' @param exclude Character vector, the symbols to be excluded from the return value
#' @return character vector with the symbols
#' @examples getSymbols(c("A*AB+B^2"))
#' @export
getSymbols <- function(char, exclude = NULL) {
  
  char <- char[char!="0"]
  out <- parse(text=char, keep.source = TRUE)
  out <- utils::getParseData(out)
  names <- unique(out$text[out$token == "SYMBOL"])
  if(!is.null(exclude)) names <- names[!names%in%exclude]
  return(names)
  
}


#' Compute matrix product symbolically
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character, the matrix product of M and N
#' @export
prodSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  dimM <- dim(M)
  dimN <- dim(N)
  if(dimM[2] != dimN[1]) {
    cat("Something is wrong with the dimensions of the matrices\n")
    return(NA)
  }
  m <- 1:dimM[1]
  n <- 1:dimN[2]
  grid <- expand.grid(m,n)
  
  MN <- apply(grid, 1, function(ik) {
    
    v <- M[ik[1],]
    w <- N[,ik[2]]
    result <- ""
    for(i in 1:length(v)) {
      if(i==1 & !(0 %in% c(w[i], v[i]))) result <- paste("(", v[i], ") * (", w[i], ")", sep="")
      if(i==1 & (0 %in% c(w[i], v[i]))) result <- "0"
      if(i >1 & !(0 %in% c(w[i], v[i]))) result <- paste(result," + (", v[i], ") * (", w[i], ")", sep="")
      if(i >1 & (0 %in% c(w[i], v[i]))) result <- result
      
    }
    if(substr(result, 1,3)=="0 +") result <- substr(result, 5, nchar(result))
    return(result)})
  
  return(matrix(MN, nrow=dimM[1], ncol=dimN[2]))
  
  
}

#' Compute matrix sumSymbolically
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character, the matrix sum of M and N
#' @export
sumSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  
  if(class(M) == "matrix") dimM <- dim(M) else dimM <- c(length(M),1)
  if(class(N) == "matrix") dimN <- dim(N) else dimN <- c(length(N),1)
  
  M <- as.character(M)
  N <- as.character(N)
  result <- c()
  
  for(i in 1:length(M)) {
    if(M[i] == "0" & N[i] == "0") result[i] <- "0"
    if(M[i] == "0" & N[i] != "0") result[i] <- N[i]
    if(M[i] != "0" & N[i] == "0") result[i] <- M[i]
    if(M[i] != "0" & N[i] != "0") result[i] <- paste(M[i]," + ",N[i], sep="")
  }  
  return(matrix(result, nrow=dimM[1], ncol=dimM[2]))
  
}

# Faster version of expand.grid
# 
# seq1 Vector
# seq2 Vector
expand.grid.alt <- function(seq1,seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}


