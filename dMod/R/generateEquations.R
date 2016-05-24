#' Generate equations from data.frame(s)
#' 
#' @param ... one or more \code{data.frame}s with columns "Description" (character), "Rate" (character), 
#' and one column per ODE state with the state names. The state columns correspond to the stoichiometric
#' matrix.
#' @param volumes Named character, volume parameters for species. Names must be a subset of the species.
#' Values can be either characters, e.g. "V1", or numeric values for the volume. If \code{volumes} is not
#' \code{NULL}, missing species are treated like 1.
#' @details This function is supposed to translate a reaction network as being defined in a csv file
#' into the raw equations.
#' 
#' @return An object of class \code{eqnList}, a named vector with the equations. Contains attributes "SMatrix"
#' (the stoichiometric matrix), "species" (the state names), "rates" (the rate expressions) and "description".
#' @examples 
#' #######################################
#' # Example 1
#' #######################################
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
#' f <- generateEquations(reactions)
#' f
#' 
#' #######################################
#' # Example 2
#' #######################################
#' reactions <- data.frame(Description = c("Activation", "Deactivation", "Production", "Degradation"), 
#'                         Rate = c("act*A", "deact*pA", "prod", "degrad*pA"), 
#'                         A=c(-1,1, 1, NA), 
#'                         pA=c(1, -1, NA, -1))
#' volumes <- c(A = "V1", pA = "V2")
#' f <- generateEquations(reactions, volumes = volumes)
#' f[1:length(f)]
#' @export
generateEquations <- function(..., volumes = NULL) {
  
  mylist <- list(...)
  if(length(mylist) > 1) S <- combine(...) else S <- mylist[[1]]
  
  description <- as.character(S$Description)
  rate <- as.character(S$Rate)
  variables <- colnames(S)[-(1:2)]
  SMatrix <- as.matrix(S[,-(1:2)]); colnames(SMatrix) <- variables
  
  volumes.draft <- structure(rep("1", length(variables)), names = variables)
  volumes.draft[names(volumes)] <- volumes
  volumes <- volumes.draft
  
    
  # check for potential errors
  exclmark <- unlist(lapply(1:length(rate), function(i) {
    
    myrate <- rate[i]
    parsedRate <- getParseData(parse(text=myrate))
    symbols <- parsedRate$text[parsedRate$token=="SYMBOL"]
    
    educts <- variables[which(SMatrix[i,]<0)]
    
    !all(unlist(lapply(educts, function(e) any(e==symbols))))
    
  }))
  
  # generate equation expressions
  terme <- lapply(1:length(variables), function(j) {
    v <- SMatrix[,j]
    nonZeros <- which(!is.na(v))
    positives <- which(v > 0)
    negatives <- which(v < 0)
    volumes.destin <- volumes.origin <- rep(volumes[j], length(v))
    if(length(positives) > 0) {
        volumes.origin[positives] <- sapply(positives, function(i) {
        candidates <- which(SMatrix[i,] < 0)
        myvolume <- unique(volumes[candidates])
        if(length(myvolume) > 1) 
          stop("species from different compartments meet in one reaction")
        if(length(myvolume) == 0) myvolume <- volumes[j]
        
        return(myvolume)
      })
    }
    volumes.ratios <- paste0("*(", volumes.origin, "/", volumes.destin, ")")
    #print(volumes.ratios)
    volumes.ratios[volumes.destin == volumes.origin] <- ""
    
    numberchar <- as.character(v)
    if(nonZeros[1] %in% positives){
      numberchar[positives] <- paste(c("", rep("+", length(positives)-1)), numberchar[positives], sep = "") 
    } else {
      numberchar[positives] <- paste("+", numberchar[positives], sep = "")
    }
    paste0(numberchar[nonZeros], "*(",  rate[nonZeros], ")", volumes.ratios[nonZeros])
  })
  
  terme <- lapply(terme, function(t) paste(t, collapse=" "))
  names(terme) <- variables
  
  terme <- do.call(c, terme)
  
  class(terme) <- "eqnList"
  attr(terme, "SMatrix") <- SMatrix
  attr(terme, "species") <- variables
  attr(terme, "volumes") <- volumes
  attr(terme, "rates") <- rate
  attr(terme, "description") <- description
  attr(terme, "exclmarks") <- which(exclmark)
  
  
  if(length(exclmark) == 0) cat("There might be a problem with one or more of the reactions.\n")
  
  return(terme)
  
  
}

#' Print function for eqnList
#' 
#' @param x equation list
#' @param ... Argument not used right now
#' @export print.eqnList
#' @export
print.eqnList <- function(x, ...) {
  
  S <- attr(x, "SMatrix")
  reactions <- apply(S, 1, function(v) {
    
    numbers <- v[which(!is.na(v))]
    educts <- -numbers[numbers < 0]
    educts[which(educts == 1)] <- " "
    products <- numbers[numbers > 0]
    products[which(products == 1)] <- " "
    educts <- paste(paste(educts, names(educts), paste = ""), collapse="+")
    products <- paste(paste(products, names(products), paste = ""), collapse="+")
    
    reaction <- paste(educts, "->", products)
    return(c(educts, products))
    
  })
  
  
  educts <- reactions[1,]
  products <- reactions[2,]

  
  exclMarks <- rep(" ", length(reactions[[1]]))
  exclMarks[attr(x, "exclmarks")] <- "!"
  
  cat("Reaction table:\n")
  out <- data.frame(educts, "->", products, attr(x, "rates"), attr(x, "description"), exclMarks)
  colnames(out) <- c("Educt",  "->",  "Product", "Rate", "Description", "Check")
  print(out)
  
  if(!is.null(attr(x, "observables")))
  cat("\nObservables:\n")
  cat(paste(paste(names(attr(x, "observables")), attr(x, "observables"), sep=" = "), "\n"))
  
  
  S[is.na(S)] <- 0
  v <- MASS::Null(t(S))
  if(ncol(v) > 0) {
    
    v <- round(v/min(abs(v)[round(abs(v), 2) > 0]), 2)
    cq <- sapply(1:ncol(v), function(i) {
      paste(paste(v[,i], names(x), sep = "*"), collapse = "+")
    })
    cat("Conserved quantities:\n")
    cat(paste0("const_", 1:length(cq), " = ", cq, "\n"))
  } 
    
  
}

#' subset of an equation list
#' 
#' @param x the equation list
#' @param ... logical expression for subsetting
#' @details The argument \code{...} can contain "Educt", "Product", "Rate" and "Description".
#' The "%in%" operator is modified to allow searches in Educt and Product (see examples).
#' 
#' @return An object of class \code{eqnList}, a named vector with the equations. Contains attributes "SMatrix"
#' (the stoichiometric matrix), "species" (the state names), "rates" (the rate expressions) and "description".
#' @examples 
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
#' f <- generateEquations(reactions)
#' subset(f, "A"%in%Educt)
#' subset(f, "pA"%in%Product)
#' subset(f, grepl("act", Rate))
#' @export subset.eqnList
#' @export
subset.eqnList <- function(x, ...) {
  
  SMatrix <- attr(x, "SMatrix")
  Educt <- lapply(1:dim(SMatrix)[1], function(row) colnames(SMatrix)[which(SMatrix[row,] < 0)])
  Product <- lapply(1:dim(SMatrix)[1], function(row) colnames(SMatrix)[which(SMatrix[row,] > 0)])
  Rate <- attr(x, "rates")
  Description <- attr(x, "description")
  Volume <- attr(x, "volumes")
  
  "%in%" <- function(x, table) sapply(table, function(mytable) match(x, mytable, nomatch=0) > 0)
  
  ind <- which(eval(substitute(...), list(Description=Description, Educt=Educt, Product=Product, Rate=Rate)))
  
  if(length(ind) == 0) return()
    
  S <- matrix(SMatrix[ind,], nrow = length(ind), dimnames = list(NULL, colnames(SMatrix)))
  ind2 <- which(sapply(1:ncol(S), function(j) !all(is.na(S[, j]))))
  S <- matrix(S[,ind2], ncol = length(ind2), dimnames = list(NULL, colnames(SMatrix)[ind2]))
  rates <- Rate[ind]
  description <- Description[ind]
  volumes <- Volume[colnames(S)]
  
  reactions <- cbind(data.frame(Description = description, Rate = rates), S)

  generateEquations(reactions, volumes = volumes)
   
  
}

#' Combine several data.frames from csv read-in
#' 
#' @param ... data.frames from csv read-in
#' @details This function is useful when separating the model into independent csv model files,
#' e.g.~a receptor model and several downstream paths. Then, the models can be recombined 
#' into one model by \code{combine()}.
#' 
#' @return A \code{data.frame} with columns "Description", "Rate" and one column for each
#' ODE state in the model. See \link{generateEquations}.
#' @export
combine <- function(...) {
  
  mylist <- list(...)
  
  species <- unique(unlist(lapply(mylist, function(S) colnames(S)[-(1:2)])))
  
  
  S <- rep(NA, length(species)+2)
  allnames <- c("Description", "Rate", species)
  names(S) <- allnames
  S <- as.data.frame(as.list(S))
  
  
  for(i in 1:length(mylist)) {
    single <- mylist[[i]]
    snames <- colnames(single)
    tofill <- which(!(allnames%in%snames))
    
    if(length(tofill)>0) {
      fill <- rep(NA, length(tofill))
      names(fill) <- allnames[!(allnames%in%snames)]
      fill <- as.data.frame(as.list(fill))
      single <- cbind(single, fill)
    }
        
    S <- rbind(S, single[,allnames])
  }
  
  S <- S[-1, ]
  return(S)
  
  
}

#' Add observables to ODEs
#' 
#' @param observable named character vector. Names correspond to observable names, the chars 
#' correspond to the observation function
#' @param f equation list
#' @details Observables are translated into an ODE and added to the list of equations
#' @return An object of class \code{eqnList}, a named vector with the equations. Contains attributes "SMatrix"
#' (the stoichiometric matrix), "species" (the state names), "rates" (the rate expressions) and "description".
#' @examples 
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1))
#' f <- generateEquations(reactions)
#' myobs <- c(tA = "s1*(pA + A)", dA = "s2*(pA-A)")
#' f <- addObservable(myobs, f)
#' @export
addObservable <- function(observable, f) {
  
  myattr <- names(attributes(f))[-1]
  allattr <- attributes(f)[myattr]
  
  if(any(names(observable) %in% names(allattr$observables))) {
    
    cat("The following observable(s) are already in the observables list:", paste(names(observable)[names(observable) %in% names(allattr$observables)], collapse=", "), "\nThe original equation list was returned.")
    return(f)
    
  }
  
  # Analyze the observable character expression
  parsedata <- getParseData(parse(text=as.character(observable)))
  symbols <- unique(parsedata[parsedata$token == "SYMBOL","text"])
  species <- symbols[symbols%in%attr(f, "species")]
  derivatives <- lapply(observable, function(obs) {
    out <- lapply(as.list(species), function(x) paste(deparse(D(parse(text=obs), x), width.cutoff = 500),collapse=""))
    names(out) <- species
    return(out)
  })
  
  newodes <- sapply(derivatives, function(der) {
    
    out <- sapply(names(der), function(n) {
      d <- der[n]
      if(d != "0") paste( paste("(", d, ")", sep="") , paste("(", f[names(d)], ")",sep=""), sep="*") else return("0")
    })
    out <- paste(out, collapse = "+")
    
    return(out)
    
  })
  
  
  
  allattr$species <- c(allattr$species, names(observable))
  f <- c(f, newodes)
  for(i in 1:length(myattr)) attr(f, myattr[i]) <- allattr[[i]]
  attr(f, "observables") <- c(attr(f, "observables"), observable)
  
  
  return(f)
  
  
}


#' Add reaction to reaction table
#' 
#' @param from character with the left hand side of the reaction, e.g. "2*A + B"
#' @param to character with the right hand side of the reaction, e.g. "C + 2*D"
#' @param rate named character. The rate associated with the reaction. The name is employed as a description
#' of the reaction.
#' @param f equation list, see \link{generateEquations}
#' @return An object of class \code{eqnList}, a named vector with the equations. Contains attributes "SMatrix"
#' (the stoichiometric matrix), "species" (the state names), "rates" (the rate expressions) and "description".
#' @examples 
#' \dontrun{
#' f <- addReaction("2*A+B", "C + 2*D", "k1*B*A^2", NULL)
#' f <- addReaction("C + A", "B + A", "k2*C*A", f)
#' }
#' @export
addReaction <- function(from, to, rate, f=NULL) {
  
  myattr <- c("class", "SMatrix", "species", "rates", "description", "exclmarks", "observables")
  
  allattr <- attributes(f)[myattr]
  
  if(class(f) != "eqnList") newEqn <- TRUE else newEqn <- FALSE
  
  # Analyze the reaction character expressions
  educts <- getSymbols(from)
  eductCoef <- 0
  if(length(educts)>0) eductCoef <- sapply(educts, function(e) sum(getCoefficients(from, e)))
  products <- getSymbols(to)
  productCoef <- 0
  if(length(products)>0) productCoef <- sapply(products, function(p) sum(getCoefficients(to, p)))

  
  # Species
  species <- unique(c(educts, products))
  
  # Description
  description <- names(rate)
  if(is.null(description)) description <- ""
  
  # Stoichiometric matrix
  SMatrix <- matrix(NA, nrow = 1, ncol=length(species)); colnames(SMatrix) <- species
  if(length(educts)>0) SMatrix[,educts] <- -eductCoef
  if(length(products)>0) {
    filled <- !is.na(SMatrix[,products])
    SMatrix[,products[filled]] <- SMatrix[,products[filled]] + productCoef[filled]
    SMatrix[,products[!filled]] <- productCoef[!filled]  
  }
  
    
  SMatrix[SMatrix == "0"] <- NA
  
    
  # data.frame
  mydata <- cbind(data.frame(Description = description, Rate = as.character(rate)), as.data.frame(SMatrix))
  row.names(mydata) <- NULL
  
  
  if(!is.null(f)) {
    mydata0 <- cbind(data.frame(Description = attr(f, "description"), 
                                Rate = attr(f, "rates")),
                     as.data.frame(attr(f, "SMatrix")))
    
    
    mydata <- combine(mydata0, mydata)
    
  }
  
  
  out <- generateEquations(mydata)
  
  return(out)

  
  
}
