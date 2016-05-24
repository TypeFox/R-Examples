#' Built-in library schemes
#'
#' @name schemes
#' @title Built-in library schemes for peptider
#' @description This data set contains descriptions of amino acid classes several commonly used library schemes: NNN, NNB, NNK, 20/20, and variations of each in which Cysteine is not considered a viable amino acid.

#' @details
#' The schemes are defined as:
#'
#' NNN: All four bases (\"N\" = G/A/T/C) possible at all three positions in the codon.
#' NNB: All four bases in the first two codon positions possible, the third position is restricted to G, T or C (= \"B\")
#' NNK/S: All four bases in the first two codon positions possible, the third position is restricted to G/T (= \"K\") or two C/G (= \"S\").
#' 2020: 20/20 describes the concept that DNA is assembled from prefabricated trimeric building blocks. This allows the generation of libraries from a predefined set of codons and thereby complete exclusion of Stop codons and other unwanted codons.
#' NNN (-C): NNN with Cysteine ignored.
#' NNB (-C): NNB with Cysteine ignored.
#' NNK/SC (-C): NNK/S with Cysteine ignored.
#' 2020 (-C): 20/20 with Cysteine ignored.
#' 
#' The schemes differ in the number of used codons, ranging from 64 (NNN), 48 (NNB), 32 (NNK/S) to 20 or less (20/20). Coding schemes that allow varying ratios of codons/amino acid, result in libraries biased towards amino acids which are encoded more often. Further, the number of Stop codons that can lead to premature termination of the peptide sequence influences the performance of the library.
#' 
#' @docType data
#' @usage data(schemes)
NULL

#' Get the specified library scheme definition
#' 
#' @param name name of the scheme as a character vector
#' @param file CSV file hosting scheme definition, if provided
#' 
#' @return a data frame of peptide classes, amino acids, and size of the classes corresponding to the selected scheme
#' 
#' @export
#' 
#' @importFrom utils data
#' @importFrom utils read.csv
#' 
#' @examples
#' scheme("NNN")
#' scheme("NNK")
scheme <- function(name, file = NULL) {
    if (is.null(file)) {
        schemes <- NULL ## R Check
        data(schemes, envir=environment())
        
        scheme_def <- schemes[[paste(tolower(name), "scheme", sep = "_")]]
        if (is.null(scheme_def)) stop(paste("No library with name", name, "is included in peptider"))
        
        return(scheme_def)
    } else {
        return(read.csv(file))
    }
}

#' Get the specified library scheme
#' 
#' @param schm either a character vector giving the name of a built-in scheme, or a data frame consisting of the scheme definition
#' @param k length of peptide sequences
#' 
#' @return list consisting of a data frame of peptide classes, size of class, and its probabilities, and a list of additional information relating to the library scheme
#' 
#' @export
#' 
#' @examples
#' libscheme("NNN")
#' libscheme("NNK", 2)
#' 
#' # Build a custom 20/20 library
#' custom <- data.frame(class = c("A", "Z"), aacid = c("SLRAGPTVIDEFHKNQYMW", "*"), c = c(1, 0))
#' libscheme(custom)
libscheme <- function(schm, k = 1) {
    if (is.character(schm)) return(libBuild(scheme(schm), k = k))
    else if (is.data.frame(schm)) return(libBuild(k, schm))
    else stop("scheme must be either a character or a data frame")
}

#' Diversity according to peptides paper (Sieber)
#'
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @param N size of the library 
#' @param lib library scheme
#' @param variance return the variance instead of the expected value
#' 
#' @return Expected Diversity of the library
#' 
#' @export
#' 
#' @examples
#' diversity(2, "NNN", 10^3)
#' diversity(2, "NNK", 10^3)
diversity <- function (k, libscheme, N, lib = NULL, variance = FALSE) 
{
    libschm <- as.character(substitute(libscheme))
    if (inherits(try(scheme(libschm), silent = TRUE), "try-error")) 
        libschm <- libscheme
    if (is.null(lib)) 
        lib <- libscheme(libschm, k)
    libdata <- lib$data
    initialloss <- (1 - (lib$info$valid/lib$info$nucleotides)^k)
    libdata$expected <- libdata$probs * N * (1 - initialloss)
    val <- sum(with(libdata, di * choices * (1 - exp(-expected/di))))
    if (variance) {
        cn <- with(libdata, (1 - 2/di)^expected)
        xn <- with(libdata, (1 - 1/di)^expected)
        libdata$var <- with(libdata, (di * (xn - cn) - 
                                          expected * (1 - 2/di)^(expected-1)))
        idx <- which(libdata$di <= 2)
        if (length(idx) > 0)
            libdata$var[idx] <- with(libdata, di*(1-probs))[idx]
        idx <- which(libdata$var < 0) # just a precaution against numerical instabilities
        if (length(idx) > 0) libdata$var[idx] <- 0
        
        val <- with(libdata, sum(choices * var)) + N*(1-initialloss)* initialloss  
    }
    return(as.numeric(val))
}

#' Diversity index according to Makowski
#'
#' The Diversity of a peptide library of length k according to Makowski and colleagues
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' 
#' @details
#' Makowski and colleagues [Makowski, Soares 2003] present another approach by defining functional diversity. They provide the mathematical background to determine the quality of a peptide library based on the probability of individual peptides to appear. In an ideal case, where every peptide has the same frequency the functional diversity is at a maximum of 1. With increasingly skew distributions, this value drops towards a minimum of 0. It is mostly independent of the actual number of sequences in a library but reflects effects caused by the degeneration of the genetic code. In the genetic code the number of codons per amino acid varies from one to six. Therefore random DNA sequences are biased towards encoding peptides enriched in amino acids encoded more frequently, which results in skew distributions of peptide frequencies. 
#' 
#' @return diversity index between 0 and 1
#' @export
#' @examples
#' makowski(2, "NNN")
#' makowski(3, "NNK")
#' makowski(3, "2020")
makowski <- function(k, libscheme) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    scheme_def <- libscheme(libschm, k)
    
    dframe <- scheme_def$data
    info <- scheme_def$info$scheme
    numAA <- sum(info$s[-nrow(info)]) 
    
    with(dframe, 1/(numAA^k*sum(((probs^2)*choices)/di)))
}

#' Coverage as expected number of peptides given all possible peptides
#'
#' Coverage of library of size N given random sampling from the pool of all possible peptides according to probabilities determined according to the library scheme.
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @param N size of the library 
#' @param lib library scheme
#' @param variance return the variance instead of the expected value
#' @return coverage index between 0 and 1
#' @export
#' @examples
#' coverage(2, "NNN", 10^3)
#' coverage(2, "NNK", 10^3)
#' coverage(2, "2020", 10^3) ## 20/20 coverage is not 1 because of random sampling.
coverage <- function(k, libscheme, N, lib=NULL, variance = FALSE) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (is.null(lib)) lib <- libscheme(libschm, k)
    
    s_count <- sum(subset(lib$info$scheme, class != "Z")$s)
    
    val <- min(diversity(k, libscheme, N, lib, variance) / s_count^k, 1)
    if (variance) val <- val / s_count^k
    
    return(val)
}

#' Relative efficiency of a library
#'
#' Relative efficiency of a peptide library, defined as the ratio of expected diversity of a peptide library relative to its overall number of oligonucleotides
#' @param k length of peptide sequences
#' @param libscheme Name (character vector) or definition (data frame) of scheme
#' @param N size of the library 
#' @param lib library, if null, libscheme will be used to create it
#' @param variance return the variance instead of the expected value
#' @return relative efficiency index between 0 and 1
#' @export
#' @examples
#' efficiency(3, "NNN", 10^2)
#' efficiency(3, "NNK", 10^2)
#' efficiency(3, "2020", 10^2) ## 20/20 efficiency is not 1 because of random sampling.
efficiency <- function(k, libscheme, N, lib=NULL, variance = FALSE) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (is.null(lib)) lib <- libscheme(libschm, k)
    libdata <- lib$data
    
    s_count <- sum(subset(lib$info$scheme, class != "Z")$s)
    
    val <- min(diversity(k, libscheme, N, lib, variance), s_count^k) / N
    if (variance) val <- val / N
    
    return(val)
}

#' Build peptide library of k-length sequences according to specified scheme
#' 
#' @import discreteRV
#' 
#' @param k length of peptide sequences
#' @param libscheme library scheme specifying classes of amino acids according to number of encodings
#' last class is reserved for stop tags and other amino acids we are not interested in.
#' @param scale1 Scaling factor for first probs
#' @param scale2 Scaling factor for second probs
#' @return library and library scheme used
#' @examples
#' user_scheme <- data.frame(class=c("A", "B", "C", "Z"),
#'                           aacid=c("SLR", "AGPTV", "CDEFHIKMNQWY", "*"),
#'                           c=c(3,2,1,1))
#' user_library <- libBuild(3, user_scheme)                        
#' @export
libBuild <- function(k, libscheme, scale1 = 1, scale2 = 1) {
    libscheme$class <- as.character(libscheme$class)
    libscheme$s <- nchar(as.character(libscheme$aacid))
    
    seq <- with(libscheme[-nrow(libscheme),], RV(class, scale1 * s*c / sum(s * c), fractions = FALSE, verifyprobs = all(scale1 == 1)))
    d <- with(libscheme[-nrow(libscheme),], RV(class, scale2 * s / sum(s), fractions = FALSE, verifyprobs = all(scale2 == 1)))
    
    d7 <- mult_reduced(d, k)
    seq7 <- mult_reduced(seq, k)
    
    di <- with(libscheme, round(d7$Prob*sum(s[-length(unique(class))])^k,0))
    pi <- seq7$Prob
    mult <- with(libscheme, s*c)
    list(data=data.frame(class = as.vector(d7[,1]), di = di, choices = seq7$Choices, probs = pi),
         info=list(nucleotides=sum(with(libscheme, mult)), 
                   valid=with(libscheme, sum(mult[-length(mult)])),
                   scheme=libscheme))
}

getChoices <- function(str) {
    test <- as.numeric(unlist(strsplit(str, split = ",")))
    
    left <- sum(test)
    total <- 1
    for (i in 1:length(test)) {
        val <- choose(left, test[i])
        left <- left - test[i]
        total <- total * val
    }
    
    return(total)
}

#' @importFrom dplyr filter
mult_reduced <- function(X, n = 2) {
    num_outcomes <- length(X)
    
    str.func <- paste("expand.grid(", paste(rep(paste("0:", n, sep = ""), times = num_outcomes), collapse = ", "), ")")
    
    grid.df <- eval(parse(text = str.func))
    grid.sub <- filter(grid.df, apply(grid.df, 1, sum) == n)
    
    grid.str <- apply(grid.sub, 1, paste, collapse = ",")
    grid.list <- split(grid.sub, 1:nrow(grid.sub))
    grid.prob <- lapply(grid.list, function(x) {
        prod(probs(X)^x)
    })
    
    grid.choices <- lapply(grid.str, getChoices)
    data.frame(Encoding = as.character(grid.str), Prob = as.numeric(grid.prob), Choices = as.numeric(grid.choices))
}

#' Detection probability in a single library of size N
#'
#' The probability that at least one of a number of specific peptide sequences (e. g. the `best' and closely related sequences) is contained in a library
#' @param lib library used in experiment, defaults to NNK with peptide length 7
#' @param size size of the library, defaults to 10^8
#' @return vector of detection probabilities for peptide sequences in each class
#' @export
#' @examples
#' summary(detect())
#'
#' require(ggplot2)
#' lib <- libscheme("NNK", 7)
#' qplot(detect(lib, size=10^8), weight=di, geom="histogram", data=lib$data)
detect <- function(lib = libscheme("NNK", 7), size = 10^8) {
    with(lib$data, 1 - exp(-size*probs/di))
}

#' @importFrom utils data
getNeighborOne <- function(x, blosum=1) {
    ## For CRAN check
    BLOSUM80 <- AA1 <- Blosum <- AA2 <- NULL
    
    data(BLOSUM80, envir=environment())
    
    replacements <- llply(strsplit(x,""), function(y) {
        llply(y, function(z) {
            as.character(subset(BLOSUM80, (AA1 == z) & (Blosum >= blosum)& (AA2 != z))$AA2 )
        })
    })[[1]]
    neighbors <- NULL
    for (i in 1:nchar(x)) {
        neighbors <- c(neighbors, paste(substr(x, 1,i-1), replacements[[i]], substr(x, i+1, nchar(x)), sep=""))
    }
    # check that all neighbors have the correct length
    idx <- which(nchar(neighbors) != nchar(x))
    if (length(idx)>0) neighbors <- neighbors[-idx]
    neighbors <- unique(c(x, neighbors))
    
    return(neighbors)
}

#' Find all neighbors of degree one for a set of peptide sequences
#' 
#' first degree neighbors - a neighbor of a peptide is defined as a peptide sequence that differs in at most one amino acid from a given sequence. 
#' Additionally, we can restrict neighbors to regard only those sequences that have a certain minimal BLOSUM loading. 
#' @import plyr
#' @param x (vector) of character strings of  peptide sequences.
#' @param blosum minimal BLOSUM loading, defaults to 1 for positive loadings only
#' @return list of neighbor sequences
#' @export
#' @examples
#' getNeighbors("APE")
#' getNeighbors(c("HI", "APE"))
#' getNeighbors(c("HI", "EARNEST", "APE"), blosum=3)
#' ## degree 2 neighbors:
#' unique(unlist(getNeighbors(getNeighbors("APE"))))
getNeighbors <- function(x, blosum=1) {
    x <- as.character(x)
    if (length(x) == 1) return(getNeighborOne(x, blosum))
    llply(x, getNeighborOne)
}

#' @importFrom utils data
getNofNeighborsOne <- function(x, blosum = 1, method="peptide", libscheme=NULL) {
    ## For CRAN check
    BLOSUM80 <- AA1 <- Blosum <- AA2 <- NULL
    
    data(BLOSUM80, envir=environment())
    
    if (method == "peptide") {
        replacements <- llply(strsplit(x,""), function(y) {
            llply(y, function(z) {
                as.character(subset(BLOSUM80, (AA1 == z) & (Blosum >= blosum) & (AA2 != z))$AA2)
            })
        })[[1]]
        return(length(unlist(replacements))+1)
    }
    
#     replacements <- llply(strsplit(x,""), function(y) {
#         llply(y, function(z) {
#             as.character(subset(BLOSUM80, (AA1 == z) & (Blosum >= blosum))$AA2)
#         })
#     })[[1]]
    stopifnot(!(is.null(libscheme) & nchar(libscheme) == 0))    
#    lib <- peptider::libscheme(libscheme)$info$scheme
    
    dnas <- sum(codons(getNeighbors(x, blosum=blosum), libscheme=libscheme))
    
#     dnas <- sum(unlist(llply(replacements, function(x) {
#         sum(unlist(llply(strsplit(x, split=""), function(w) {
#             lib[grep(w, lib$aacid, fixed = TRUE), "c"]
#         }))) 
#     })))
 #   dnas <- dnas + sum(unlist(llply(unlist(strsplit(x, split="")), function(w) { 
 #       lib[grep(w, lib$aacid, fixed=TRUE),"c"]
 #   })))
    return(dnas)
}

#' Compute the number of neighbor of degree one for a set of peptide sequences
#' 
#' first degree neighbors - a neighbor of a peptide is defined as a peptide sequence that differs in at most one amino acid from a given sequence. 
#' Additionally, we can restrict neighbors to regard only those sequences that have a certain minimal BLOSUM loading. 
#' Use this function for only a few peptide sequences. Any larger number of peptide sequences will take too much main memory.
#' @param x (vector) of character strings of  peptide sequences.
#' @param blosum minimal BLOSUM loading, defaults to 1 for positive loadings only
#' @param method character string, one of "peptide" or "codon". This specifies the level at which the neighbors are calculated.
#' @param libscheme library scheme under which neighbors are being calculated. this is only of importance, if method="dna"
#' @return vector of numbers of neighbors 
#' @import plyr
#' @export
#' @examples
#' getNofNeighbors("APE")
#' getNofNeighbors(c("NEAREST", "EARNEST"))
#' getNofNeighbors("N")
#' getNofNeighbors("N", method="codon", libscheme="NNK")
getNofNeighbors <- function(x, blosum = 1, method="peptide", libscheme=NULL) {
    data(BLOSUM80, envir=environment())
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    x <- as.character(x)
    if (length(x) == 1) return(getNofNeighborsOne(x, blosum, method, libschm))
    
    return(laply(x, getNofNeighborsOne, blosum, method, libschm))
}

#' Compute the number of codon representations for a (vector of) peptide sequence(s)
#' 
#' use this function for only a few peptide sequences. Any larger number of peptide sequences should be dealt with in the framework of the library scheme and the detect function.
#' @param x (vector) of character strings of  peptide sequences.
#' @param libscheme library scheme under which neighbors are being calculated. this is only of importance, if method="dna"
#' @param flag internal use only: Set to true if calling this from another function
#' @return vector of numbers of codons 
#' @export
#' @import plyr
#' @examples
#' codons("APE", libscheme="NNK")
#' codons("HENNING", libscheme="NNK")
codons <- function(x, libscheme, flag = FALSE) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    if (length(x) == 1) return(codonsOne(x, libschm))
    
    unlist(llply(x, codonsOne, schm=libschm))
}

codonsOne <- function(x, schm) {
    stopifnot(!(is.null(schm) & nchar(schm) == 0))
    lib <- libscheme(schm)$info$scheme
    
    x <- as.character(x)
    prod(unlist(llply(strsplit(x, split="")[[1]], function(w) { 
        lib[grep(w, lib$aacid, fixed=TRUE),"c"]
    })))
}

#' Probability of detection of a peptide sequence 
#' 
#' use this function for only a few peptide sequences. Any larger number of peptide sequences should be dealt with in the framework of the library scheme and the detect function.
#' @param x (vector) of character strings of  peptide sequences.
#' @param libscheme library scheme under which neighbors are being calculated. 
#' @param N number of valid DNA clones investigated
#' @return probability of detection
#' @export
#' @examples
#' ppeptide("APE", libscheme="NNK", N=10^8)
#' ppeptide("HENNING", libscheme="NNK", N=10^8)
ppeptide <- function(x, libscheme, N) {
    libschm <- as.character(substitute(libscheme)) ## Compatibility with old interface
    if (inherits(try(scheme(libschm), silent = TRUE), 'try-error')) libschm <- libscheme
    
    n <- sum(codons(x, libscheme=libschm, flag = TRUE))
    Max <- peptider::libscheme(libschm, 1)$info$valid^nchar(as.character(x[1]))
    1 - exp(-N*n/Max)
}

#' BLOSUM80 matrix
#' 
#' @name BLOSUM80
#' @title BLOSUM80 matrix
#' @description The BLOSUM80 matrix, which stands for Blocks Substitution Matrix, defines log-odds scores for the ratio of the chance of two amino acids appearing in a sequence over the chance that the two amino acids appear in any sequence.  Larger scores indicate a higher probability of substitutions.  This matrix is used in order to compute sequences which are in the neighborhood of other sequences.
#' @docType data
#' @usage data(BLOSUM80)
NULL

#' Calculate neighborhood distribution
#' 
#' Calculate distribution of neighbors under library scheme lib for peptide sequences of length k.
#' @param sch library scheme
#' @param k length of the peptide sequences
#' @return dataset of peptide sequences: AA are amino acid sequences, 
#' c0 are codons for self representation, 
#' cr is the ratio of #neighbors in first degree neighborhood (not counting self representations) and #codons in self representation
#' N1 is the number of neighbors in codon representation (including self representation) 
#' @export
#' @import plyr
#' @importFrom stats xtabs
#' @examples
#' genNeighbors(scheme("NNK"), 2)
#' genNeighbors(scheme("2020"), 2)
genNeighbors <- function(sch, k) {
#    sch <- libscheme(lib,1)$info$scheme
    
    # don't include class Z
    schl <- unlist(strsplit(as.character(sch$aacid[-length(sch$aacid)]), ""))
    schd <- data.frame(AA=schl, C0=codons(schl, libscheme=sch))
    schd$C1 <- getNofNeighbors(schd$AA, method="codon", libscheme=sch) - schd$C0
    schd$C1T0 <- with(schd, C1/C0)
    
    ctabs <- function(values, labels) {
        labels <- as.character(labels)
        
        tv <- xtabs(~values)
        ct <- ldply(names(tv), function(x) {
            data.frame(AA=paste(labels[which(values == x)], collapse=""),
                       c=x)
        })
        ct
    }
    
    user <- with(schd, ctabs(paste(C0, C1, sep=":"), AA))
    user <- data.frame(user, ldply(strsplit(as.character(user$c), ":"), 
                                   function(x) as.numeric(unlist(x))))
    names(user)[3:4] <- c("c0", "c1")
    user$cr <- with(user, c1/c0)
    user$c <- user$c0
    user$s <- nchar(as.character(user$AA))
    
    i <- 1
    C0 <- user$c
    CR <- user$cr
    L <- user$AA
    S <- user$s
    while (i < k) {
        C0 <- outer(C0, user$c, FUN = "*")
        S <- outer(S, user$s, FUN = "*")
        CR <- outer(CR, user$cr, FUN = "+")
        L <- outer(L, user$AA, FUN = "paste", sep = ",")
        
        sprintf("stage %d done",i)
        i <- i + 1
    }
    x <- data.frame(AA=as.vector(L), 
                    c0=as.vector(C0), 
                    cr=as.vector(CR), 
                    s = as.vector(S))
    ## number of neighbors
    x$N1 <- with(x, c0*(cr+1))
    x
}

#' Calculate neighborhood distribution
#' 
#' Calculate distribution of neighbors under library scheme lib for peptide sequences of length k.
#' @param sch library scheme
#' @param k length of the peptide sequences
#' @return dataset of peptide sequences: L are amino acid sequences, 
#' c0 are codons for self representation, 
#' cr is the ratio of #neighbors in first degree neighborhood (not counting self representations) and #codons in self representation
#' N1 is the number of neighbors in codon representation (including self representation) 
#' s is the number of peptide sequences described by the label
#' o is the number of peptide sequences reached by permutations
#' @export
#' @import plyr
#' @importFrom stats xtabs
#' @examples
#' genNeighbors_reduced(scheme("NNK"), 2)
#' genNeighbors_reduced(scheme("2020"), 2)
genNeighbors_reduced <- function(sch, k) {    
    schl <- unlist(strsplit(as.character(sch$aacid[-length(sch$aacid)]), ""))
    schd <- data.frame(AA=schl, C0=codons(schl, libscheme=sch))
    schd$C1 <- getNofNeighbors(schd$AA, method="codon", libscheme=sch) - schd$C0
    schd$C1T0 <- with(schd, C1/C0)
    
    ctabs <- function(values, labels) {
        labels <- as.character(labels)
        
        tv <- xtabs(~values)
        ct <- ldply(names(tv), function(x) {
            data.frame(AA=paste(labels[which(values == x)], collapse=""),
                       c=x)
        })
        ct
    }
    
    user <- with(schd, ctabs(paste(C0, C1, sep=":"), AA))
    user <- data.frame(user, ldply(strsplit(as.character(user$c), ":"), 
                                   function(x) as.numeric(unlist(x))))
    names(user)[3:4] <- c("c0", "c1")
    user$cr <- with(user, c1/c0)
    user$c <- user$c0
    user$s <- nchar(as.character(user$AA))
    user$one <- 1
    
    i <- 1
    dx <- user
    dx$o <- dx$one
    dx$L <- user$AA
    while (i < k) {
        C0 <- as.vector(outer(dx$c0, user$c, FUN = "*"))
        S <- as.vector(outer(dx$s, user$s, FUN = "*"))
        CR <- as.vector(outer(dx$cr, user$cr, FUN = "+"))
        O <- as.vector(outer(dx$o, user$one, FUN = "*"))
        #    L <- as.vector(outer(L, user$AA, FUN = function(x,y) {        
        #        x <- as.character(x)
        #        y <- as.character(y)
        #       paste(pmin(x,y), pmax(x,y), sep=",") 
        #    }))    
        L <- as.vector(outer(dx$L, user$AA, "paste", sep = ","))
        lx <- strsplit(as.character(L), split=",")
        L <- laply(llply(lx, sort), paste, collapse=",")
        
        x <- data.frame(L, C0,CR,S,O)
        dx <- ddply(x, .(L), summarise, 
                    c0=C0[1],
                    cr=CR[1],
                    s=S[1],
                    o=sum(O))
        
        #    save(L,C0,CR,S,O, file=sprintf("%s-%d.RData", lib, i))
        dx$N1 <- with(dx, c0*(cr+1))
        dx$k <- i+1
        #   cat(sprintf("stage %d done\n",i))
        #   write.table(dx, file="neighbors2.csv", sep=",", col.names=!file.exists("neighbors2.csv"), append=TRUE, row.names=FALSE)
        
        i <- i + 1
    }
    dx
}
