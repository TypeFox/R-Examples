readLxb <- function(paths, filter=TRUE, text=FALSE) {
    # Read multiple LXB files and return a list of matrices (one for each LXB).
    #
    # If 'text=TRUE' then each item is a list with a 'text' and 'data' entry.
    # The 'text' is the text segment of the LXB file and the 'data' entry is
    # the data segment (same as the matrix return when 'text=FALSE').
    #
    # If only one LXB is read then the head of the list is returned instead of
    # returning a list with only one element.
    #
    # If 'filter=TRUE' then bad reads are dropped (reads which have no bead ID
    # or which did not pass the doublet discriminator test).  Otherwise all
    # data is included in the output.
    #
    # The name of each LXB file is used to set the 'names' attribute of the
    # returned list.  If the name ends with a letter and a 1-2 digit number
    # then it is assumed that this encodes the row&column of each well on a
    # plate.  In this case the output will be sorted by column and the 'names'
    # attribute is set to the well name instead of the full file name.

    go <- function(filename) {
        x <- .Call("read_lxb", as.character(filename), as.logical(text))
        if (!is.null(x)) {
            if (!is.null(x$data) && filter && 'RID' %in% colnames(x$data)
                    && 'DBL' %in% colnames(x$data))
                x$data <- x$data[x$data[ ,'RID'] != 0 & x$data[ ,'DBL'] != 0, ]

            if (!text)
                x <- x$data
        }
        x
    }

    names <- Sys.glob(paths)
    lxbs  <- lapply(names, go)

    m <- regexec(".*([a-zA-Z])([0-9]+)[.]lxb", names)
    if (all(lapply(m, '[', 1L) != -1)) {
        # All names contain row & column information, extract into data frame
        parts <- do.call(rbind, lapply(regmatches(names, m), '[', c(2L,3L)))
        df <- data.frame(row=parts[ ,1], column=as.integer(parts[ ,2]))

        # Order LXBs by: A1, B1, ..., A2, B2, ...
        idx <- with(df, order(column, row))

        names(lxbs) <- with(df, paste(row, column, sep=""))
        lxbs <- lxbs[idx]
    } else {
        # No row & column information found, use filename
        names(lxbs) <- lapply(names, function(x) sub(".lxb", "", basename(x)))
    }

    if (length(lxbs) == 1)
        lxbs <- lxbs[[1]]
    lxbs
}
