read_stm_CLUTO <-
function(file)
{
    ## Read CLUTO sparse matrix format.
    
    ## Read in the matrix file.
    l <- strsplit(readLines(file, warn = FALSE), "[[:space:]]+")
    l <- lapply(l, as.double)
    l <- lapply(l, na.omit)

    ## Extract header information.
    nRow <- as.integer(l[[1L]][1L])
    nCol <- as.integer(l[[1L]][2L])
    nElem <- l[[1L]][3L]
    ## Remove header
    l <- l[-1L]

    ## Compute i, j, and v slots for a simple_triplet_matrix.
    rowLen <- sapply(l, length)
    l <- unlist(l)
    i <- rep.int(seq_len(nRow), rowLen / 2)
    j <- l[seq.int(1, length(l), by = 2)]
    v <- l[seq.int(2, length(l), by = 2)]

    ## Sanity check
    if(nElem != length(v))
        stop("invalid matrix format")

    ## Generate sparse matrix
    m <- simple_triplet_matrix(i, j, v, nRow, nCol)

    if(is.character(file)) {
        ## Use col labels file if available and valid.
        if(file.exists(f <- sprintf("%s.clabel", file))) {
            lines <- readLines(f)
            if(length(lines) == nCol)
                colnames(m) <- lines
        }
        ## Use row labels file if available and valid.
        if(file.exists(f <- sprintf("%s.rlabel", file))) {
            lines <- readLines(f)
            if(length(lines) == nRow)
                rownames(m) <- lines
        }
        ## Use row class file if available.
        if(file.exists(f <- sprintf("%s.rclass", file))) {
            lines <- readLines(f)
            if(length(lines) == nRow)
                attr(m, "rclass") <- lines
        }
    }
        
    m
}

write_stm_CLUTO <-
function(x, file)
{
    ## Write CLUTO sparse matrix format.
    
    x <- as.simple_triplet_matrix(x)

    ## Generate header.
    header <- paste(x$nrow, x$ncol, length(x$v))

    ## Generate content.
    content <- Map(function(u, v) paste(u, v, collapse = " "),
                   split(x$j, x$i),
                   split(x$v, x$i))

    ## Write out.
    writeLines(c(header, unlist(content)), file)

    if(is.character(file)) {
        if(!is.null(rnms <- rownames(x)))
            writeLines(rnms, sprintf("%s.rlabel", file))
        if(!is.null(cnms <- colnames(x)))
            writeLines(cnms, sprintf("%s.clabel", file))
    }
            
}

read_stm_MC <-
function(file, scalingtype = NULL)
{
    ## Read the CCS format variant employed by MC
    ## (http://www.cs.utexas.edu/users/dml/software/mc/) and related
    ## software projects at cs.utexas.edu such as gmeans.
    ## The main MC web page points to
    ## http://www.cs.utexas.edu/users/jfan/dm/README.html 
    ## which no longer seems to exist: but the MC sources contain a file
    ## README with some information.
    ## The basic CCS format is documented in
    ## http://www.cs.utexas.edu/users/inderjit/Resources/sparse_matrices.
    
    d <- scan(sprintf("%s_dim", file), what = integer(0), quiet = TRUE)
    nr <- d[1L]
    nc <- d[2L]

    i <- scan(sprintf("%s_row_ccs", file), what = integer(0),
              quiet = TRUE)
    p <- scan(sprintf("%s_col_ccs", file), what = integer(0),
              quiet = TRUE)
    if(is.null(scalingtype)) {
        ## The name of the file with the non-zero entries varies with
        ## the t-f-n scaling pattern employed (and possibly an 'i' at
        ## the end indicating that row and columne scaling were
        ## performed independently:
        scalingtype <-
            expand.grid(c("t", "l"),
                        c("x", "f", "e", "1"),
                        c("x", "n", "1"),
                        c("", "i"))
        ## (Not sure whether all combinations really make sense.)
        scalingtype <- 
            apply(scalingtype, 1L, paste, collapse = "")
    }
    files <- sprintf("%s_%s_nz", file, scalingtype)
    pos <- which(file.exists(files))[1L]
    x <- scan(files[pos], what = numeric(0), quiet = TRUE)
    scalingtype <- scalingtype[pos]

    ## Sanity check
    if(d[3L] != length(x))
        stop("invalid matrix format")

    ## In special cases (e.g., when CCS was produced by the MC toolkit,
    ## see http://userweb.cs.utexas.edu/users/jfan/dm/README.html) we
    ## can also infer the row and col names.
    rnms <- if(file.exists(f <- sprintf("%s_words", file))) {
        readLines(f)[seq(from = 2L, length.out = nr)]
    } else NULL
    cnms <- if(file.exists(f <- sprintf("%s_docs", file))) {
        sub("^[^ ]*: ", "", readLines(f))
    } else NULL

    m <- simple_triplet_matrix(i + 1L,
                               rep.int(seq_len(nc), diff(p)),
                               x,
                               nr, nc, list(rnms, cnms))
    attr(m, "scalingtype") <- scalingtype
    m
}

write_stm_MC <-
function(x, file)
{
    ## Write CCS sparse matrix format as used by MC and other software
    ## projects from cs.utexas.edu such as gmeans.

    ## <FIXME>
    ## This said:

    ## Gmeans uses a compressed column storage (CCS)
    ## See http://www.cs.utexas.edu/users/inderjit/Resources/sparse_matrices
    ##
    ## However since Gmeans clusters along columns, and the input for
    ## our skmeans clusters along rows, we would need to transpose the
    ## matrix first, and then write it to CCS. 
    ##
    ## Instead we could directly write to compressed row storage (CRS)
    ## to avoid the transpose
    ## See
    ## http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000

    ## Does this mean we should not transpose in general, but when
    ## writing out for gmeans in skmeans only?
    ## </FIXME>
    
    x <- t(as.simple_triplet_matrix(x))
    # Based on slam/work/Matrix.R
    ind <- order(x$j, x$i)
    write(paste(nrow(x), ncol(x), length(x$v)),
          sprintf("%s_dim", file))
    write(x$i[ind] - 1L,
          sprintf("%s_row_ccs", file), sep = "\n")
    write(c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
          sprintf("%s_col_ccs", file), sep = "\n")
    write(x$v[ind],
          sprintf("%s_tfn_nz", file), sep = "\n")

    ## Could also try to write a _docs file.
    ## But what does the 2nd half of the _words files contain?
}

