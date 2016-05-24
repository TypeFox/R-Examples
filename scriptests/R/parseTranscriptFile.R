parseTranscriptFile <- function(file, ignoreUpToRegExpr=NULL, ignoreAfterRegExpr=NULL, subst=NULL)
{
    ## This function reads a transcript a file, separates commands and output,
    ## parses the commands, and returns a list with each element containing
    ## a list with some of the following elements:
    ##   input/comment: the exact text that contained the commands (including prompts)
    ##   expr: the parsed expression (can be NULL if source is a comment)
    ##   control: control text
    ##   output: the output that apparently came from the command (can be character(0))
    ##   garbage: lines that couldn't be understood
    ## file : [character] the name of the file to read
    ##
    lines <- readLines(file, -1, warn=FALSE)
    # Ignore intial blank lines
    blank <- grep("^[ \\t]*$", lines, perl=TRUE, invert=TRUE)
    if (length(blank) && blank[1] > 1)
        lines <- lines[-seq(len=blank[1]-1)]
    if (!is.null(ignoreUpToRegExpr)) {
        i <- grep(ignoreUpToRegExpr, lines, perl=TRUE)
        if (length(i)>0)
            lines <- lines[-seq(len=i[1])]
        # Ignore intial blank lines
        blank <- grep("^[ \\t]*$", lines, perl=TRUE, invert=TRUE)
        if (length(blank) && blank[1] > 1)
            lines <- lines[-seq(len=blank[1]-1)]
    }
    if (!is.null(ignoreAfterRegExpr)) {
        i <- grep(ignoreAfterRegExpr, lines, perl=TRUE)
        if (length(i)>0)
            lines <- lines[-seq(len=i[1])]
    }
    if (length(subst)) {
        if (!is.character(subst) || is.null(names(subst)))
            stop("subst must be a named character vector")
        for (i in seq(along=subst)) {
            pattern <- names(subst)[i]
            repl <- subst[i]
            if (pattern=="")
                stop("subst cannot have empty names -- the names are the pattern to replace")
            lines <- gsub(pattern, repl, lines, perl=TRUE)
        }
    }

    ## calculate a code for each line:
    ##   0: output
    ##   1: command
    ##   2: continuation
    ##   3: comment or empty command
    ##   4: control
    ##   5: garbage: not understood
    lineType <- c("output", "command", "continuation", "comment/empty", "control", "garbage")
    code <- ifelse(regexpr("^> ?(#|[[:space:]]*$)", lines)>0, 3,
                   ifelse(regexpr("^> ", lines)>0, 1,
                          ifelse(regexpr("^\\+ ", lines)>0, 2,
                                 ifelse(regexpr("^#@", lines)>0, 4, 0))))
    ## Look for output in the middle of comments -- where these are blank
    ## lines, they can be reclassified as comments.
    ## If we don't do this, it messes up block counting and results in
    ## confusing error messages.
    lastCode <- -1
    for (i in seq(along=lines)) {
        if (lastCode==3) {
            if (lines[i]=="" || regexpr("^ *$", lines)>0)
                code[i] <- 6
        } else {
            lastCode <- code[i]
        }
    }
    if (length(i <- which(code==6))) {
        lines <- lines[-i]
        code <- code[-i]
    }
    ## Identify blocks of contiguous command+continuation, & control-output.
    ## Comments go in their own block.
    ## Insert separators between lines starting with ">"
    ## because these must be separate commands (the first having no output)
    ## Convert command+continuation to all command + separator (-1)
    code2 <- as.vector(rbind(code, ifelse((code==1 | code==2) & c(code[-1],0)==1, -1, -2)))
    code2 <- code2[code2 != -2]
    ## Now that we have continuation blocks separated (by -1), change code 2 to code 1
    code2 <- replace(code2, code2 == 2, 1)
    runs <- rle(code2)
    i <- runs$values != -1
    runs$values <- runs$values[i]
    runs$lengths <- runs$lengths[i]
    nblocks <- sum(is.element(runs$values, c(1, 3)))
    if (! runs$values[length(runs$values)] %in% c(1,3))
        nblocks <- nblocks + 1
    i <- 1 # counter in lines
    j <- 1 # counter in runs
    blocks <- lapply(seq(len=nblocks), function(k) {
        res <- list()
        while (j <= length(runs$values) && !is.element(runs$values[j], c(1,3))) {
            ## warning("not expecting lines of type '", lineType[runs$values[j]+1], "' at line ", i,
            ##         " (", paste("\"", lines[seq(i, len=min(runs$lengths[j], 3))], "\"", collapse=", ", sep=""), ")")
            res$garbage <- c(res$garbage, lines[seq(i,len=runs$lengths[j])])
            i <<- i + runs$lengths[j]
            j <<- j + 1
        }
        if (j <= length(runs$values) && runs$values[j]==3) {
            res$comment <- lines[seq(i,len=runs$lengths[j])]
            i <<- i + runs$lengths[j]
            j <<- j + 1
        } else if (j <= length(runs$values) && runs$values[j]==1) {
            res$input <- lines[seq(i,len=runs$lengths[j])]
            i <<- i + runs$lengths[j]
            j <<- j + 1
            # output and control lines can be mixed up, separate them out
            while (j <= length(runs$values) && !is.element(runs$values[j], c(1,3))) {
                if (runs$values[j]==0) {
                    res$output <- c(res$output, lines[seq(i,len=runs$lengths[j])])
                    i <<- i + runs$lengths[j]
                    j <<- j + 1
                } else if (runs$values[j]==4) {
                    res$control <- c(res$control, lines[seq(i,len=runs$lengths[j])])
                    i <<- i + runs$lengths[j]
                    j <<- j + 1
                } else {
                    warning("not expecting lines of type '", lineType[runs$values[j]+1], "' at line ", i,
                            " (", paste("\"", lines[seq(i, len=min(runs$lengths[j], 3))], "\"", collapse=", ", sep=""), ")")
                    res$garbage <- c(res$garbage, lines[seq(i,len=runs$lengths[j])])
                    i <<- i + runs$lengths[j]
                    j <<- j + 1
                }
            }
            text <- gsub("^[>+] ?", "", res$input, perl=TRUE)
            res$expr <- try(parse(text=text, srcfile=NULL), silent=TRUE)
            ## try to make a syntax error message look like it does at the prompt
            if (is(res$expr, "try-error"))
                res$expr[1] <- gsub("Error in parse\\(text = text, srcfile = NULL\\) :[ \n\t]+", "Error: ", res$expr[1], perl=TRUE)
        }
        while (j <= length(runs$values) && !is.element(runs$values[j], c(1,3))) {
            ## warning("not expecting lines of type '", lineType[runs$values[j]+1], "' at line ", i,
            ##         " (", paste("\"", lines[seq(i, len=min(runs$lengths[j], 3))], "\"", collapse=", ", sep=""), ")")
            res$garbage <- c(res$garbage, lines[seq(i,len=runs$lengths[j])])
            i <<- i + runs$lengths[j]
            j <<- j + 1
        }
        return(res)
    })
    return(blocks)
}

