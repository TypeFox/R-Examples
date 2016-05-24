
### some tools that make life easier

### copy *Rout to *Rout.save
cpRoutsave <- function(Routdir = NULL, Routsavedir = NULL) {

    Routfiles <- list.files(path = Routdir, pattern = "\\.Rout$", 
                            full.names = FALSE)
    srcfiles <- file.path(Routdir, Routfiles)
    destfiles <- file.path(Routsavedir, 
                           paste(Routfiles, ".save", sep = ""))
    file.copy(srcfiles, destfiles, overwrite = TRUE)
}

### attach all data frames in the global environment
gattach <- function() {

    env <- globalenv()
    var <- eval(parse(text = "ls()"), envir = env)
    df <- sapply(var, function(x)
        eval(parse(text = 
            paste("is.data.frame(", x, ")", sep = "", collapse = "")), 
            envir = env))
    if (any(df)) {
      var <- var[df]
      out <- sapply(var, function(x) 
          eval(parse(text = 
              paste("attach(", x, ")", sep = "", collapse = "")),
              envir = env))
    }
}


### extract and check Robject or Rcmd LaTeX markup
extRact <- function(file, what = "Robject") {

    x <- readLines(file)
    indx <- grep(what, x)
    
    out <- sapply(indx, function(i) {
        obj <- NULL
        while (TRUE) {
            where <- regexpr(what, x[i])
            if (where != -1) {
                x[i] <- substring(x[i], where)
                dm <- delimMatch(x[i])
                obj <- c(obj, (substring(x[i], dm + 1, 
                         dm + attr(dm, "match.length") - 2)))
                x[i] <- substring(x[i], dm + attr(dm, "match.length"))
            } else {
                break
            }
        }
        return(obj)
    })
    cmds <- unique(gsub("\\\\", "", out))
    gattach()
    for (cmd in cmds) {
        a <- try(eval(parse(text = cmd)))
        if (class(a) == "try-error") print(a)
    }
    cmds
}

### try to polish S{in,out}put environments, this needs
### manual refinements in some places
prettyS <- function(file, texenvironment = c("Sinput", "Soutput"), 
                    width = 63, split = " ", write = TRUE) {

    ### handle Sinput or Soutput environments
    texenvironment <- match.arg(texenvironment)
    if (texenvironment == "Sinput" && split == " ")
        split <- c(", ", "/", " ")

    ### dirty hack: in `Makefile's I want to call `prettyS'
    ### right after weaving and thus have only `file.Rnw' available
    if (length(grep("Rnw$", file)) > 0) file <- gsub("Rnw$", "tex", file)

    ### read file
    x <- readLines(file)

    ### remove all end-line spaces
    x <- gsub("\\s+$", "", x)

    ### determine begin and end lines of environment
    start <- grep(paste("^\\\\begin\\{", texenvironment, "\\}$", 
                  sep = "", collapse = ""), x)
    end <- grep(paste("^\\\\end\\{", texenvironment, "\\}$", 
                  sep = "", collapse = ""), x)
    if (length(start) == 0) return(NULL)
    if (length(start) != length(end)) 
        stop("unbalanced begin and end statements")
    n <- length(start)

    for (i in 1:n) {
  
        ### iterate over all lines longer than width
        lines <- (start[i]):(end[i])
        lines <- lines[sapply(x[lines], nchar) > width]
        for (line in lines) {
            cat("prettyS: line ", line, " too long: \n", x[line], "\n")
            y <- x[line]
            add <- sapply(split, function(s) 
                ifelse(length(grep(s, y)) > 0, nchar(s), 0))
            if (all(add == 0)) next()
            s <- split[min(which(add > 0))]
            y <- unlist(strsplit(y, split = s))
            nc <- sapply(y, nchar) + add[min(which(add > 0))]
            pos <- cumsum(nc) <= width
            if (!any(pos)) next()
            newline <- cumsum(nc)[max(which(pos))]
            plus <- length(grep("^\\+", x[line])) > 0 && 
                    substr(x[line], newline - 1, newline) != ", "
            x[line] <- paste(substr(x[line], 1, newline), "\n",
                ifelse(texenvironment == "Sinput", options("continue"), ""),
                ifelse(plus, "    ", ""),
                "    ",
                substr(x[line], newline + 1, nchar(x[line])), sep = "", 
                collapse = "")
#            if (length(grep("^\\+", x[line + 1])) > 0 && 
#                (nchar(x[line + 1]) + (nchar(x[line]) - newline) < width)) {
#                y <- x[line + 1]
#                y <- gsub("^\\+   ", "", y)
#                x[line] <- paste(x[line], y, sep = "", collapse = "")
#                x[line + 1] <- ""
#            }
             cat("prettyS: ", x[line], "\n")
        }
    }
    if (write)
        writeLines(x, con = file)
}

### extract all Sinput environments from tex files
chkS <- function(file) {

    texenvironment <- "Sinput"

    ### read file
    x <- readLines(file)

    ### determine begin and end lines of environment
    start <- grep(paste("^\\\\begin\\{", texenvironment, "\\}$", 
                  sep = "", collapse = ""), x)
    end <- grep(paste("^\\\\end\\{", texenvironment, "\\}$", 
                  sep = "", collapse = ""), x)
    if (length(start) == 0) return(NULL)
    if (length(start) != length(end)) 
        stop("unbalanced begin and end statements")
    n <- length(start)

    y <- NULL

    for (i in 1:n) {
  
        ### iterate over all lines longer than width
        lines <- (start[i] + 1):(end[i] - 1)
        x[lines] <- gsub("^R>", "", x[lines])
        x[lines] <- gsub("^\\+", "", x[lines])
        y <- c(y, x[lines])

    }
    y
}


### read in a BibTeX file and return as list
readBibtex <- function(file = NULL) {

    bib <- readLines(file)

    entries <- grep("^@", bib)
    labels <- gsub(",$", "", gsub("^@[A-Za-z].*\\{", "", bib[entries]))

    if (any(duplicated(labels))) {
        print(labels[duplicated(labels)])
        stop("non-unique BibTeX labels in ", file)
    }

    biblist <- vector(mode = "list", length = length(entries))

    for (i in 1:length(entries)) {
        nexte <- ifelse(i == length(entries), length(entries), 
                        entries[i + 1] - 1)
        biblist[[i]] <- bib[entries[i]:nexte]
        empty <- grep("^$", biblist[[i]])
        if (length(empty) > 0)
        biblist[[i]] <- biblist[[i]][-empty]
    }
    names(biblist) <- labels
    class(biblist) <- "txtBibtex"
    return(biblist)
}

### the subset of a BibTeX database actually used in `file'
extractBibtex <- function(file, bibtex) {

    if (class(bibtex) != "txtBibtex")
        bibtex <- readBibtex(bibtex)
    tex <- readLines(file)
    tex <- tex[grep("\\cite", tex)]
    enames <- gsub("\\+", "\\\\+", names(bibtex))
    cited <- sapply(enames, function(name) length(grep(name, tex)) > 0)
    biblist <- bibtex[cited]
    class(biblist) <- "txtBibtex"
    return(biblist)
}

### output to a file
toBibtex.txtBibtex <- function(object, ...) {

    tmp <- lapply(object, function(bib) {
        cat(paste(bib, "\n"))
        cat("\n\n")
    })
}

### set package version in BibTeX (quick'n'dirty hack)
pkgversions <- function(file) {

    x <- readLines(file)
    indx <- grep("VERSION", x)

    for (i in indx) {
        xx <- strsplit(x[i], " ")[[1]]
        xx <- xx[grep("VERSION", xx)]
        pkg <- gsub("[},]", "", gsub("VERSION", "", xx))
        version <- packageDescription(pkg)$Version
        x[i] <- gsub(paste(pkg, "VERSION", sep = "", collapse = ""), version, 
                     x[i])
    }
    class(x) <- "Latex"
    x
}

pkgs <- function()
 c("scatterplot3d", "ape", "coin", "flexmix", "gee", "ipred", "lme4",
   "mclust", "party", "randomForest", "rmeta", "vcd")
