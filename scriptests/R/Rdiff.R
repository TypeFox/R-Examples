## compares 2 files -- function is the same as in package tools
Rdiff <- function(from, to, useDiff = FALSE)
{
    clean <- function(txt)
    {
        ## remove R header
        if(length(top <- grep("^(R version|R : Copyright)", txt,
                              perl = TRUE, useBytes = TRUE)) &&
           length(bot <- grep("quit R.$", txt, perl = TRUE, useBytes = TRUE)))
            txt <- txt[-(top[1]:bot[1])]
        ## remove BATCH footer
        nl <- length(txt)
        if(grepl("^> proc.time()", txt[nl-2])) txt <- txt[1:(nl-3)]
        ## regularize fancy quotes.
        txt <- gsub("(\xe2\x80\x98|\xe2\x80\x99)", "'", txt,
                      perl = TRUE, useBytes = TRUE)
        if(.Platform$OS.type == "windows") # not entirely safe ...
            txt <- gsub("(\x93|\x94)", "'", txt, perl = TRUE, useBytes = TRUE)
        pat <- '(^Time |^Loading required package|^Package [A-Za-z][A-Za-z0-9]+ loaded|^<(environment|promise|pointer): )'
        txt[!grepl(pat, txt, perl = TRUE, useBytes = TRUE)]
    }

    left <- clean(readLines(from))
    right <- clean(readLines(to))
    if (!useDiff && (length(left) == length(right))) {
        bleft <- gsub("[[:space:]]+", " ", left, perl=TRUE)
        bright <- gsub("[[:space:]]+", " ", right, perl=TRUE)
        if(all(bleft == bright)) return(0L)
        cat("\n")
        diff <- bleft != bright
        ## FIXME do run lengths here
        for(i in which(diff)) {
            cat(i,"c", i, "\n< ", left[i], "\n", "---\n> ", right[i], "\n",
                sep = "")
        }
        return(1L)
    } else {
        ## FIXME: use C code, or something like merge?
        ## The files can be very big.
        if(!useDiff) cat("\nfiles differ in number of lines:\n")
        a <- tempfile()
        b <- tempfile()
        writeLines(left, a)
        writeLines(right, b)
        return(system(paste("diff -bw", shQuote(a), shQuote(b))))
    }
}
