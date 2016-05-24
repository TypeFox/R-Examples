sb.extract <- function(x, effect, operation) {
# Arguments:
#   x: output of sb.submit() (html contents)
#   effect: "main", "mplusall", or "mplussome"
#   operation: "construct" or "check"


  calculation.flag <- TRUE


# Reconstruct html contents
  html.content.tree <- htmlTreeParse(x$html)
  html.part         <- html.content.tree$children$html
  body              <- xmlChildren(html.part)$body
  body.value        <- xmlSApply(body, xmlValue)


# Extract inputs
  if (operation == "check") {
    # Choice Sets to be Checked
    nrows <- length(gregexpr("\n", body.value[[10]])[[1]]) + 1
    temp  <- gsub("\n", " ", body.value[[10]]) 
    CSC   <- matrix(as.numeric(unlist(strsplit(temp, " "))), 
                    nrow = nrows, byrow = TRUE)
    TCE <- NULL
    SGE <- NULL

  } else {  # operation == "construct"
    # Choice Sets to be checked
    CSC <- NULL

    # Treatment Combination Entered (TCE)
    temp  <- xmlValue(body[[10]][[1]][[2]]) 
    nrows <- length(gregexpr("\n", temp)[[1]]) + 1
    temp  <- gsub("\n", " ", temp) 
    TCE   <- matrix(as.numeric(unlist(strsplit(temp, " "))), 
                    nrow = nrows, byrow = TRUE)

    # Sets of Generators Entered (SGE)
    temp <- xmlValue(body[[10]][[2]][[2]])
    if (gregexpr("\n", temp)[[1]] > 0) {
      nrows <- length(gregexpr("\n", temp)[[1]]) + 1
    } else {
      nrows <- 1
    }
    temp <- gsub("\n", " ", temp)
    SGE  <- matrix(as.numeric(unlist(strsplit(temp, " "))), 
                   nrow = nrows, byrow = TRUE)
  }

  # Determinant of C and effect (valid for effect == "mplussome")
  if (isTRUE(effect == "mplussome")) {
    DetC  <- xmlValue(body[[11]][[1]][[2]])
    INTER <- xmlValue(body[[11]][[2]][[2]])
  }

  # Determinant of C and effect (valid for effect == "mplus...")
  if (isTRUE(effect == "mplusall")) {
    DetC  <- xmlValue(body[[11]][[1]][[2]])
    INTER <- NULL
  }

  # Determinant of C and effect (valid for effect == "mplussome")
  if (isTRUE(effect == "main")) {
    DetC  <- NULL
    INTER <- NULL
  }


# Extract results of constructing/checking choice sets
  # Message from processing stage
  MSG <- body.value[[14]]

  # Choice sets created
  if (operation == "check") {
    i  <- 16
    CS <- NULL
  } else {
    i <- 18
    nrows <- length(gregexpr("\n", body.value[[16]])[[1]]) + 1
    temp  <- gsub("\n", " ", body.value[[16]])
    temp  <- as.numeric(unlist(strsplit(temp, " ")))
    CS    <- matrix(temp, nrow = nrows, byrow = TRUE)
  }

  # B matrix
  nrows   <- length(gregexpr("\n", body.value[[i]])[[1]]) + 1
  temp    <- gsub("\n", " ", body.value[[i]])
  vec     <- unlist(strsplit(temp, " "))
  vec.chr <- subset(vec, vec != "")
  if (length(vec.chr) == 1) {
    B.chr <- vec.chr
    B.num <- NULL
    calculation.flag <- FALSE
  } else {
    vec.num <- sapply(vec.chr, function(x) eval(parse(text = x)))
    B.chr   <- matrix(vec.chr, nrow = nrows, byrow = TRUE)
    B.num   <- matrix(vec.num, nrow = nrows, byrow = TRUE)
  }

  # Lambda Matrix
  nrows   <- length(gregexpr("\n", body.value[[i + 2]])[[1]]) + 1
  temp    <- gsub("\n", " ", body.value[[i + 2]])
  vec     <- unlist(strsplit(temp, " "))
  vec.chr <- subset(vec, vec != "")
  if (length(vec.chr) == 1) {
    L.chr  <- vec.chr
    L.num  <- NULL
    calculation.flag <- FALSE
  } else {
    vec.num <- sapply(vec.chr, function(x) eval(parse(text = x)))
    L.chr   <- matrix(vec.chr, nrow = nrows, byrow = TRUE)
    L.num   <- matrix(vec.num, nrow = nrows, byrow = TRUE)
  }

  # C Matrix
  nrows   <- length(gregexpr("\n", body.value[[i + 4]])[[1]]) + 1
  temp    <- gsub("\n", " ", body.value[[i + 4]])
  vec     <- unlist(strsplit(temp, " "))
  vec.chr <- subset(vec, vec != "")
  if (length(vec.chr) == 1) {
    C.chr  <- vec.chr
    C.num  <- NULL
    calculation.flag <- FALSE
  } else {
    vec.num <- sapply(vec.chr, function(x) eval(parse(text = x)))
    C.chr   <- matrix(vec.chr, nrow = nrows, byrow = TRUE)
    C.num   <- matrix(vec.num, nrow = nrows, byrow = TRUE)
  }

  # C inverse
  nrows   <- length(gregexpr("\n", body.value[[i + 6]])[[1]]) + 1
  temp    <- gsub("\n", " ", body.value[[i + 6]])
  vec     <- unlist(strsplit(temp, " "))
  vec.chr <- subset(vec, vec != "")
  if (length(vec.chr) == 1) {
    CI.chr  <- vec.chr
    CI.num  <- NULL
    calculation.flag <- FALSE
  } else {
    vec.num <- sapply(vec.chr, function(x) eval(parse(text = x)))
    CI.chr  <- matrix(vec.chr, nrow = nrows, byrow = TRUE)
    CI.num  <- matrix(vec.num, nrow = nrows, byrow = TRUE)
  }

  # Correlation matrix
  nrows   <- length(gregexpr("\n", body.value[[i + 8]])[[1]]) + 1
  temp    <- gsub("\n", " ", body.value[[i + 8]])
  vec     <- unlist(strsplit(temp, " "))
  vec.chr <- subset(vec, vec != "")
  if (length(vec.chr) == 1) {
    COR.chr <- vec.chr
    COR.num <- NULL
    calculation.flag <- FALSE
  } else {
    vec.num <- sapply(vec.chr, function(x) eval(parse(text = x)))
    COR.chr <- matrix(vec.chr, nrow = nrows, byrow = TRUE)
    COR.num <- matrix(vec.num, nrow = nrows, byrow = TRUE)
  }

  
# Return  
  rtn <- list(input  = list(operation         = operation,
                            effect            = effect,
                            choice.sets.check = CSC,
                            treatment         = TCE,
                            generators        = SGE,
                            determinant.c     = DetC,
                            interaction       = INTER),
              messages                        = MSG,
              output = list(choice.sets       = CS,
                            b.chr             = B.chr,
                            b.num             = B.num,
                            l.chr             = L.chr,
                            l.num             = L.num,
                            c.chr             = C.chr,
                            c.num             = C.num,
                            cinv.chr          = CI.chr,
                            cinv.num          = CI.num,
                            correlation.chr   = COR.chr,
                            correlation.num   = COR.num),
              calculation                     = calculation.flag,
              html                            = as.character(x$html),
              version                         = x$version)
  return(rtn)
}
