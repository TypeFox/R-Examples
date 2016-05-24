"read.nexus.data" <- function (file)
{
    # Simplified NEXUS data parser.
    #
    # Version: 09/13/2006 01:01:59 PM CEST
    #          (modified by EP 2011-06-01)
    #
    # By:      Johan Nylander, nylander @ scs.fsu.edu
    #
    # WARNING: This is parser reads a restricted nexus format,
    #          see README for details.
    #
    # Argument (x) is a nexus formatted data file.
    #
    # Returns  (Value) a list of data sequences each made of a single
    #          vector of mode character where each element is a character.
    #
    # TODO:    Error checking, gap/missing, find.datatype, etc.
    #------------------------------------------------------------------

    "find.ntax" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                                       "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        ntax
    }

    "find.nchar" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                                        "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        nchar
    }

    "find.matrix.line" <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
                matrix.line <- as.numeric(i)
                break
            }
        }
        matrix.line
    }

    "trim.whitespace" <- function (x)
    {
        gsub("\\s+", "", x)
    }

    "trim.semicolon" <- function (x)
    {
        gsub(";", "", x)
    }

    X <- scan(file = file, what = character(), sep = "\n",
              quiet = TRUE, comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0

    for (j in start.reading:NROW(X)) {
        Xj <- trim.semicolon(X[j])
        if(Xj == "") {
            break
        }
        if(any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE))) {
            break
        }
        ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
        if (length(ts) > 2) {
            stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
        }
        if (length(ts) !=2) {
            stop("nexus parser failed to read the sequences (ts!=2)")
        }
        Seq <- trim.whitespace(ts[2])
        Name <- trim.whitespace(ts[1])
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
        if (any(l <- grep(nAME, names(Obj)))) {
            tsp <- strsplit(Seq, NULL)[[1]]
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[l]][p] <- tsp[k]
                chars.done <- k
            }
        }
        else {
            names(Obj)[i] <- Name
            tsp <- strsplit(Seq, NULL)[[1]]
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[i]][p] <- tsp[k]
                chars.done <- k
            }
        }
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax) {
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            if (tot.nchar == nchar*ntax) {
                print("ntot was more than nchar*ntax")
                break
            }
            pos <- tot.nchar
        }
        else {
            i <- i + 1
        }
    }
    if (tot.ntax != 0) {
        cat("ntax:",ntax,"differ from actual number of taxa in file?\n")
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            cat(names(Obj[i]),"has",length(Obj[[i]]),"characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
    Obj <- lapply(Obj, tolower)
    Obj
}

