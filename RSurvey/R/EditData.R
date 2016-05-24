# A GUI for viewing and editing table formatted data.

EditData <- function(d, col.names=names(d), row.names=NULL, col.formats=NULL,
                     read.only=FALSE, changelog=NULL, win.title="Data",
                     parent=NULL) {

  ## Additional functions

  # Save table and close
  SaveTable <- function() {
    tclServiceMode(FALSE)
    SaveActiveEdits()
    changelog <<- GetEdits()
    tclServiceMode(TRUE)
    tclvalue(tt.done.var) <- 1
    return()
  }

  # Save edits in active cell
  SaveActiveEdits <- function() {
    cell <- as.character(tkindex(frame3.tbl, "active"))
    ij <- as.integer(strsplit(cell, ",")[[1]])
    i <- ij[1]
    j <- ij[2]
    old.val <- FormatValues(i, j)
    new.val <- paste(as.character(tkget(frame3.tbl, cell)), collapse=" ")
    if (identical(new.val, old.val))
      return()
    SaveEdits(new.val, i, j)
    new.val <- FormatValues(i, j)
    e <- data.frame(timestamp=Sys.time(), cell=cell, old=old.val, new=new.val,
                    stringsAsFactors=FALSE)
    undo.stack <<- rbind(undo.stack, e)
    redo.stack <<- NULL
    return()
  }

  # Get single cell value for table
  GetCellValue <- function(r, c) {
    i <- as.integer(r)
    j <- as.integer(c)
    if (i > 0 && j > 0) {
      val <- FormatValues(i, j, is.fmt=TRUE)
    } else if (i == 0 && j > 0) {
      val <- col.names[j]
    } else if (i > 0 && j == 0) {
      val <- row.names[i]
    } else {
      val <- ""
    }
    return(as.tclObj(val, drop=TRUE))
  }

  # Format values
  FormatValues <- function(i, j, is.fmt=FALSE) {
    fmt.vals <- rep("", length(i))
    for (column in unique(j)) {
      idxs <- j %in% column
      vals <- d[[column]][i[idxs]]
      col.class <- class(vals)
      is.time <- any(c("POSIXt", "Date") %in% col.class)
      fmt <- ifelse(is.fmt || is.time, col.formats[column], "")
      if (is.time) {
        if ("POSIXt" %in% col.class)
          fmt.vals[idxs] <- POSIXct2Character(vals, fmt=fmt)
        else
          fmt.vals[idxs] <- format(vals, format=fmt)
      } else {
        if (fmt == "")
          fmt.vals[idxs] <- format(vals, trim=TRUE, nsmall=ndigits,
                                   scientific=ifelse(is.fmt, NA, FALSE),
                                   drop0trailing=TRUE)
        else
          fmt.vals[idxs] <- sprintf(fmt, vals)
      }
    }
    return(fmt.vals)
  }

  # Save edits to data frame
  SaveEdits <- function(vals, i, j) {
    if (is.character(i)) {
      Fun <- function(x) as.integer(strsplit(x, ",")[[1]])
      ij <- t(vapply(i, Fun, c(0L, 0L)))
      i <- ij[, 1]
      j <- ij[, 2]
    }
    for (column in unique(j)) {
      idxs <- which(j == column)
      col.class <- class(d[[column]])
      new.vals <- vals[idxs]
      new.vals[new.vals == ""] <- NA
      if ("POSIXt" %in% col.class) {
        fmt <- gsub("%OS[[:digit:]]+", "%OS", col.formats[column])
        fmt <- ifelse(fmt == "", "%Y-%m-%d %H:%M:%S", fmt)
        tz <- attr(d[[column]], "tzone")
        new.vals <- strptime(new.vals, format=fmt, tz=tz)
        if ("POSIXct" %in% col.class)
          new.vals <- as.POSIXct(new.vals, tz=tz)
      } else if ("Date" %in% col.class) {
        fmt <- col.formats[column]
        new.vals <- as.Date(new.vals, format=ifelse(fmt == "", "%Y-%m-%d", fmt))
      } else if ("factor" %in% col.class) {
        old.levels <- levels(d[[column]])
        new.levels <- unique(c(old.levels, na.omit(new.vals)))
        if (!all(new.levels %in% old.levels))
          levels(d[[column]]) <<- new.levels
      } else {
        new.vals <- suppressWarnings(as(new.vals, col.class[1]))
      }
      d[[column]][i[idxs]] <<- new.vals
    }
    return()
  }

  # Validate entry value
  ValidateEntryValue <- function(P, S) {
    cell <- as.character(tkindex(frame3.tbl, "active"))
    ij <- as.integer(strsplit(cell, ",")[[1]])
    new.val <- as.character(S)
    if (identical(new.val, CheckEntry(class(d[[ij[2]]][ij[1]]), new.val))) {
      tkset(frame3.tbl, cell, as.character(P))
      is.valid <- TRUE
    } else {
      is.valid <- FALSE
    }
    return(as.tclObj(is.valid))
  }

  # Validate cell value
  ValidateCellValue <- function(s, S) {
    cell <- as.character(tkindex(frame3.tbl, "active"))
    ij <- as.integer(strsplit(cell, ",")[[1]])
    i <- ij[1]
    j <- ij[2]
    new.val <- as.character(S)
    if (identical(new.val, CheckEntry(class(d[[j]][i]), new.val))) {
      tclvalue(value.var) <- as.character(S)
      is.valid <- TRUE
    } else {
      is.valid <- FALSE
    }
    return(as.tclObj(is.valid))
  }

  # Set active cell
  SetActiveCell <- function() {
    cell <- as.character(tkindex(frame3.tbl, "active"))
    ij <- as.integer(strsplit(cell, ",")[[1]])
    fmt.val <- FormatValues(ij[1], ij[2])
    if (!read.only)
      tkset(frame3.tbl, cell, fmt.val)
    tclvalue(value.var) <- as.character(fmt.val)
    return()
  }

  # Change active cell
  ChangeActiveCell <- function(s, S) {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    if (!read.only && s != "") {
      old.cell <- as.integer(strsplit(s, ",")[[1]])
      i <- old.cell[1]
      j <- old.cell[2]
      old.val <- FormatValues(i, j)
      new.val <- paste(as.character(tkget(frame3.tbl, s)), collapse=" ")
      if (!identical(new.val, old.val)) {
        SaveEdits(new.val, i, j)
        e <- data.frame(timestamp=Sys.time(), cell=s, old=old.val,
                        new=FormatValues(i, j), stringsAsFactors=FALSE)
        undo.stack <<- rbind(undo.stack, e)
        redo.stack <<- NULL
      }
      tkset(frame3.tbl, s, FormatValues(i, j, is.fmt=TRUE))
    }
    ij <- as.integer(strsplit(S, ",")[[1]])
    i <- ij[1]
    j <- ij[2]
    if (i == 0 || j == 0) {
      if (i == 0)
        i <- 1L
      if (j == 0)
        j <- 1L
      tkactivate(frame3.tbl, paste(i, j, sep=","))
    }
    tktag.delete(frame3.tbl, "row.idx")
    tktag.delete(frame3.tbl, "col.idx")
    tcl(frame3.tbl, "tag", "cell", "row.idx", paste(i, 0, sep=","))
    tcl(frame3.tbl, "tag", "cell", "col.idx", paste(0, j, sep=","))
    tktag.raise(frame3.tbl, "row.idx")
    tktag.raise(frame3.tbl, "col.idx")
    tktag.configure(frame3.tbl, "row.idx", background="#B3B3B3")
    tktag.configure(frame3.tbl, "col.idx", background="#B3B3B3")
    tktag.raise(frame3.tbl, "active", "sel")
    if (read.only) {
      tktag.delete(frame3.tbl, "active_readonly")
      tcl(frame3.tbl, "tag", "cell", "active_readonly", paste(i, j, sep=","))
      tktag.configure(frame3.tbl, "active_readonly", background="#FBFCD0")
      tktag.raise(frame3.tbl, "active_readonly", "sel")
    }
    SetActiveCell()
    return()
  }

  # Undo edit
  UndoEdit <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    if (is.null(undo.stack) || nrow(undo.stack) == 0)
      return()
    time.stamp <- undo.stack$timestamp
    idxs <- which(time.stamp == time.stamp[nrow(undo.stack)])
    cells <- undo.stack[idxs, "cell"]
    Fun <- function(i) as.integer(strsplit(i, ",")[[1]])
    ij <- t(vapply(cells, Fun, c(0, 0)))
    old.vals <- undo.stack[idxs, "old"]
    SaveEdits(old.vals, cells)
    redo.stack <<- rbind(redo.stack, undo.stack[idxs, , drop=FALSE])
    undo.stack <<- undo.stack[-idxs, , drop=FALSE]
    if (nrow(undo.stack) == 0)
      undo.stack <<- NULL
    fmt.old.vals <- FormatValues(ij[, 1], ij[, 2], is.fmt=TRUE)
    for (i in seq_along(cells))
      tkset(frame3.tbl, cells[i], fmt.old.vals[i])
    SetActiveCell()
    return()
  }

  # Redo edit
  RedoEdit <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    if (is.null(redo.stack) || nrow(redo.stack) == 0)
      return()
    time.stamp <- redo.stack$timestamp
    idxs <- which(time.stamp == time.stamp[nrow(redo.stack)])
    cells <- redo.stack[idxs, "cell"]
    Fun <- function(i) as.integer(strsplit(i, ",")[[1]])
    ij <- t(vapply(cells, Fun, c(0, 0)))
    new.vals <- redo.stack[idxs, "new"]
    SaveEdits(new.vals, cells)
    undo.stack <<- rbind(undo.stack, redo.stack[idxs, , drop=FALSE])
    redo.stack <<- redo.stack[-idxs, , drop=FALSE]
    fmt.new.vals <- FormatValues(ij[, 1], ij[, 2], is.fmt=TRUE)
    for (i in seq_along(cells))
      tkset(frame3.tbl, cells[i], fmt.new.vals[i])
    SetActiveCell()
    return()
  }

  # Bypass copy command
  BypassCopyCmd <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    cells <- as.character(tkcurselection(frame3.tbl))
    Fun <- function(i) as.integer(strsplit(i, ",")[[1]])
    ij <- t(vapply(cells, Fun, c(0, 0)))
    ilim <- range(ij[, 1])
    jlim <- range(ij[, 2])
    ni <- ilim[2] - ilim[1] + 1L
    nj <- jlim[2] - jlim[1] + 1L
    ij.all <- cbind(rep(ilim[1]:ilim[2], nj), rep(jlim[1]:jlim[2], each=ni))
    cells.all <- paste(ij.all[, 1], ij.all[, 2], sep=",")
    if (!all(cells.all %in% cells)) {
      msg <- "The copy command cannot be used on multiple selections."
      tkmessageBox(icon="info", message=msg, title="Copy", type="ok",
                   parent=tt)
      return()
    }
    vals <- FormatValues(ij[, 1], ij[, 2])
    write.table(matrix(vals, nrow=ni, ncol=nj, byrow=TRUE), file="clipboard",
                sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
    return()
  }

  # Bypass cut command
  BypassCutCmd <- function() {
    BypassCopyCmd()
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    cells <- as.character(tkcurselection(frame3.tbl))
    Fun <- function(i) as.integer(strsplit(i, ",")[[1]])
    ij <- t(vapply(cells, Fun, c(0, 0)))
    i <- ij[, 1]
    j <- ij[, 2]
    old.vals <- FormatValues(i, j)
    new.vals <- "NA"
    e <- data.frame(timestamp=Sys.time(), cell=cells, old=old.vals,
                    new=new.vals, stringsAsFactors=FALSE)
    undo.stack <<- rbind(undo.stack, e)
    redo.stack <<- NULL
    SaveEdits(new.vals, i, j)
    for (cell in cells)
      tkset(frame3.tbl, cell, "NA")
    SetActiveCell()
    return()
  }

  # Bypass paste command
  BypassPasteCmd <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    active.cell <- as.character(tkindex(frame3.tbl, "active"))
    args <- list(file="clipboard", colClasses="character", sep="\t",
                 flush=TRUE, fill=TRUE, na.strings=NULL)
    new.vals <- suppressWarnings(try(do.call(read.table, args), silent=TRUE))
    if (inherits(new.vals, "try-error"))
      return()
    match.length <- attr(regexpr("\t", new.vals, fixed=TRUE), "match.length")
    if (!all(match.length == match.length[1])) {
      msg <- paste("Clipboard contains text string copied from multiple",
                   "selections and is unsuitable for pasting.")
      tkmessageBox(icon="info", message=msg, title="Paste", type="ok",
                   parent=tt)
      return()
    }
    active.ij <- as.integer(strsplit(active.cell, ",")[[1]])
    m <- nrow(new.vals)
    n <- ncol(new.vals)
    ij <- cbind(rep(seq_len(m), n), rep(seq_len(n), each=m))
    new.vals <- new.vals[ij]
    ij[, 1] <- ij[, 1] + active.ij[1] - 1L
    ij[, 2] <- ij[, 2] + active.ij[2] - 1L
    i <- ij[, 1]
    j <- ij[, 2]
    cells <- paste(i, j, sep=",")
    old.vals <- FormatValues(i, j)
    e <- data.frame(timestamp=Sys.time(), cell=cells, old=old.vals,
                    new=new.vals, stringsAsFactors=FALSE)
    undo.stack <<- rbind(undo.stack, e)
    redo.stack <<- NULL
    SaveEdits(new.vals, i, j)
    new.vals <- FormatValues(i, j, is.fmt=TRUE)
    for (i in seq_along(cells))
      tkset(frame3.tbl, cells[i], new.vals[i])
    SetActiveCell()
    return()
  }

  # Bypass return command
  BypassReturnCmd <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    active.cell <- as.character(tkindex(frame3.tbl, "active"))
    old.ij <- as.integer(strsplit(active.cell, ",")[[1]])
    i <- old.ij[1]
    if (i > 0 & i < (m + 1L))
      i <- i + 1L
    j <- old.ij[2]
    if (i == 0 | j == 0)
      return()
    new.ij <- paste(i, j, sep=",")
    tkselection.clear(frame3.tbl, "all")
    tkactivate(frame3.tbl, new.ij)
    tkselection.set(frame3.tbl, new.ij)
    tksee(frame3.tbl, new.ij)
    return()
  }

  # Search data table

  CallSearch <- function(is.replace=FALSE) {
    search.defaults$find.what <- as.character(tclvalue(pattern.var))
    ans <- Search(is.replace, defaults=search.defaults, parent=tt)
    if (is.null(ans))
      return()
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"), add=TRUE)

    tclvalue(pattern.var) <- ans$find.what
    if (!identical(ans, search.defaults)) {
      search.defaults <<- ans
      matched.cells <<- NULL
    }

    Find("next", is.replace)

    are.matches <- !is.null(matched.cells) && nrow(matched.cells) > 0
    if (!(is.replace && are.matches))
      return()

    if (ans$is.replace.first)
      matched.cells <<- matched.cells[1, , drop=FALSE]
    cells <- paste0(matched.cells[, 1], ",", matched.cells[, 2])

    old.vals <- FormatValues(matched.cells[, 1], matched.cells[, 2])
    new.vals <- gsub(ans$find.what, ans$replace.with, old.vals,
                     ignore.case=!ans$is.match.case, perl=ans$is.perl,
                     fixed=!ans$is.reg.exps, useBytes=FALSE)

    e <- data.frame(timestamp=Sys.time(), cell=cells, old=old.vals,
                    new=new.vals, stringsAsFactors=FALSE)
    undo.stack <<- rbind(undo.stack, e)
    redo.stack <<- NULL
    matched.cells <<- NULL

    SaveEdits(new.vals, cells)

    for (idx in seq_along(cells))
      tkset(frame3.tbl, cells[idx], new.vals[idx])
    SetActiveCell()
    return()
  }

  # Find value in data table

  Find <- function(direction="next", is.replace=FALSE) {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    pattern <- as.character(tclvalue(pattern.var))
    if (pattern == "")
      return()

    if (is.null(matched.cells)) {
      is.match.word <- search.defaults$is.match.word
      ignore.case <- ifelse(search.defaults$is.match.case, FALSE, TRUE)
      fixed <- ifelse(search.defaults$is.reg.exps, FALSE, TRUE)
      perl <- search.defaults$is.perl
      is.search.sel <- search.defaults$is.search.sel

      if (is.search.sel) {
        cells <- as.character(tcl(frame3.tbl, "tag", "cell", "sel"))
        ncells <- length(cells)
      } else {
        ncells <- m * n
      }

      chunk.size <- 10000L
      nsteps <- ceiling(ncells / chunk.size)
      pb <- ProgressBar(label="Indexing search results\u2026", maximum=ncells,
                        nsteps=nsteps, parent=tt)
      chunks <- split(seq_len(ncells), ceiling(seq_len(ncells) / chunk.size))

      if (is.search.sel) {
        cells <- strsplit(cells, ",", fixed=TRUE)
        cells <- list(i=as.integer(vapply(cells, function(i) i[1], "")),
                      j=as.integer(vapply(cells, function(i) i[2], "")))
        desc.i.order <- order(cells$i)
        cells <- lapply(cells, function(i) i[desc.i.order])
      } else {
        cells <- list(i=rep(seq_len(m), each=n), j=rep(seq_len(n), m))
      }

      if (is.match.word) {
        if (ignore.case)
          Match <- function(i) which(tolower(i) == tolower(pattern))
        else
          Match <- function(i) which(i == pattern)
      } else {
        Match <- function(i) grep(pattern, i, ignore.case=ignore.case,
                                  fixed=fixed, perl=perl)
      }
      Fun <- function(i) {
        ans <- try(SetProgressBar(pb, value=chunks[[i]][1], step=i),
                   silent=TRUE)
        if (inherits(ans, "try-error"))
          return()
        chunks[[i]][Match(FormatValues(cells$i[chunks[[i]]],
                                       cells$j[chunks[[i]]]))]
      }
      matched.idxs <- c(lapply(seq_along(chunks), Fun), recursive=TRUE)

      if (length(matched.idxs) == 0L) {
        msg <- paste0("Search string \'", pattern, "\' not found.")
        tkmessageBox(icon="info", message=msg, title="Find", type="ok",
                     parent=tt)
        return()
      }

      if (is.search.sel) {
        cells <- do.call(cbind, lapply(cells, function(i) i[matched.idxs]))
      } else {
        col.div <- matched.idxs / n
        i <- as.integer(ceiling(col.div))
        j <- as.integer(round(n * (col.div - trunc(col.div))))
        j[j == 0L] <- n
        cells <- cbind(i, j)
      }
      matched.cells <<- cells
    }

    active.i <- as.integer(tcl(frame3.tbl, "tag", "row", "active"))
    active.j <- as.integer(tcl(frame3.tbl, "tag", "col", "active"))
    if (length(active.i) == 0L) {
      active.i <- 1L
      active.j <- 1L
    }

    if (direction == "next") {
      cell.below <- matched.cells[, 1] > active.i |
                   (matched.cells[, 1] == active.i &
                    matched.cells[, 2] > active.j)
      cell.above <- !cell.below
      if (any(cell.below)) {
        cell <- head(matched.cells[cell.below, , drop=FALSE], n=1)
      } else if (any(cell.above)) {
        cell <- head(matched.cells[cell.above, , drop=FALSE], n=1)
      } else {
        return()
      }
    } else {
      cell.above <- matched.cells[, 1] < active.i |
                   (matched.cells[, 1] == active.i &
                    matched.cells[, 2] < active.j)
      cell.below <- !cell.above
      if (any(cell.above)) {
        cell <- tail(matched.cells[cell.above, , drop=FALSE], n=1)
      } else if (any(cell.below)) {
        cell <- tail(matched.cells[cell.below, , drop=FALSE], n=1)
      } else {
        return()
      }
    }

    if (!is.replace) {
      cell.str <- paste(cell[1, 1], cell[1, 2], sep=",")
      tkactivate(frame3.tbl, cell.str)
      tksee(frame3.tbl, cell.str)
    }
  }

  # View data record
  ViewRecord <- function() {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    rec <- as.character(tclvalue(record.var))
    if (is.na(rec))
      return()
    idx <- which(row.names %in% rec)
    if (length(idx) > 0) {
      active.cell <- as.character(tkindex(frame3.tbl, "active"))
      active.col <- as.integer(strsplit(active.cell, ",")[[1]])[2]
      tksee(frame3.tbl, paste(idx[1], active.col, sep=","))
    } else {
      msg <- "Row name (or record number) not found."
      tkmessageBox(icon="info", message=msg, title="View", type="ok",
                   parent=tt)
    }
  }

  # Get edits

  GetEdits <- function() {
    s <- NULL
    if (is.null(changelog) & (is.null(undo.stack) || nrow(undo.stack) == 0))
      return(s)
    if (!is.null(changelog)) {
      changelog$cell <- paste(match(changelog$record, row.names),
                              match(changelog$variable, col.names), sep=",")
      changelog <- changelog[, c("timestamp", "cell", "old", "new")]
      undo.stack <- rbind(changelog, undo.stack)
    }
    for (i in unique(undo.stack$cell)) {
      undo.stack.cell <- undo.stack[undo.stack$cell == i, , drop=FALSE]
      undo.stack.cell <- undo.stack.cell[order(undo.stack.cell$timestamp), ,
                                         drop=FALSE]
      m <- nrow(undo.stack.cell)
      cell <- as.integer(strsplit(undo.stack.cell$cell, ",")[[1]])
      old <- undo.stack.cell$old[1]
      new <- undo.stack.cell$new[m]
      if (identical(old, new))
        next
      lst <- list(timestamp=undo.stack.cell$timestamp[m],
                  record=row.names[cell[1]], variable=col.names[cell[2]],
                  old=old, new=new)
      s <- rbind(s, as.data.frame(lst))
    }
    return(s)
  }

  # View changelog
  ViewChangeLog <- function() {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    ow <- options(width=200)$width
    on.exit(options(width=ow), add=TRUE)
    txt <- paste(capture.output(GetEdits()), collapse="\n")
    EditText(txt, read.only=TRUE, win.title="Change log",
             is.fixed.width.font=TRUE, parent=tt)
    return()
  }

  # View structure and summary information
  ViewInfo <- function(type) {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    ow <- options(width=200)$width
    on.exit(options(width=ow), add=TRUE)
    names(d) <- col.names
    if (type == "Structure")
      txt <- capture.output(str(d))
    else
      txt <- capture.output(lapply(d, summary))
    txt <- paste(c(txt, ""), collapse="\n")
    EditText(txt, read.only=TRUE, win.title=type,
             is.fixed.width.font=TRUE, parent=tt)
    return()
  }

  # Resize column
  ResizeColumn <- function(x, y) {
    tclServiceMode(FALSE)
    on.exit(tclServiceMode(TRUE))
    x <- as.character(x)
    y <- as.character(y)
    border.col <- as.character(tcl(frame3.tbl, "border", "mark", x, y, "col"))
    if (length(border.col) > 0) {
      j <- as.integer(border.col)
      tcl(frame3.tbl, "width", j, col.width[j])
    }
  }

  # Change cell selection
  ChangeCellSelection <- function () {
    is.search.sel <- search.defaults$is.search.sel
    if (is.logical(is.search.sel) && is.search.sel)
      matched.cells <<- NULL
  }

  # Plot histogram
  PlotHistogram <- function() {
    tkconfigure(tt, cursor="watch")
    on.exit(tkconfigure(tt, cursor="arrow"))
    BuildHistogram(d, var.names=col.names, parent=tt)
  }


  ## Main program

  # GUI requires TkTable
  is.tktable <- try(tcl("package", "present", "Tktable"), silent=TRUE)
  if (inherits(is.tktable, "try-error"))
    stop("tkTable is not available")

  # Check validity of data
  if (!inherits(d, c("list", "matrix", "data.frame")))
    stop("invalid class for data table")
  rtn.class <- class(d)
  if (!inherits(d, "list")) {
    d.rownames <- rownames(d)
    d.colnames <- colnames(d)
    if (inherits(d, "matrix")) {
      rownames(d) <- NULL
      colnames(d) <- NULL
      d <- split(d, rep(seq_len(ncol(d)), each=nrow(d)))
    } else {
      d <- as.list(d)
    }
    if (is.null(row.names))
      row.names <- d.rownames
  }

  # Check dimensions of data table
  n <- length(d)
  m <- unique(vapply(d, length, 0))
  if (length(m) != 1)
    stop("components of list have unequal lengths")
  if (m == 0 || n == 0)
    return()

  # Check validity of changelog
  if (!inherits(changelog, c("NULL", "data.frame")))
    stop("invalid changelog")
  changelog.old <- changelog

  # Minimum number of digits to the right of the decimal point
  ndigits <- nchar(format(as.integer(1 / sqrt(.Machine$double.eps))))

  # Number of rows and columns in the viewable table
  nrows <- ifelse(m > 15L, 15L, m)
  ncols <- ifelse(n > 6L, 6L, n)

  # Account for missing arguments

  if (is.character(col.names) && length(col.names) == n) {
    col.names <- col.names[seq_len(n)]
    col.names[is.na(col.names)] <- ""
    col.names <- gsub("(^ +)|( +$)", "", col.names)
  } else {
    Fun <- function(i) paste0(LETTERS, i)
    col.names <- c(LETTERS, as.vector(t(vapply(LETTERS, Fun, rep("", 26)))))
    col.names <- col.names[seq_len(n)]
    names(d) <- col.names
  }

  if (is.character(col.formats) && length(col.formats) == n) {
    col.formats <- as.character(col.formats[seq_len(n)])
    col.formats[is.na(col.formats)] <- ""
  } else {
    col.formats <- rep("", n)
  }

  if (!inherits(row.names, c("integer", "character")) || length(row.names) != m)
    row.names <- seq_len(m)

  # Determine width and height of column 0 and row 0, respectively
  col.0.width  <- max(nchar(row.names)) + 1L
  row.0.height <- max(vapply(strsplit(col.names, "\n"), length, 0L))

  # Determine column widths
  col.width <- NULL
  for (j in seq_len(n)) {
    if (col.names[j] == "")
      nchar.title <- 0
    else
      nchar.title <- max(vapply(strsplit(col.names[j], "\n"), nchar, 0L))
    max.rows <- ifelse(m > 200L, 200L, m)
    nchar.data <- max(nchar(FormatValues(seq_len(max.rows), rep(j, max.rows),
                                         is.fmt=TRUE)))
    col.width[j] <- max(c(nchar.title, nchar.data)) + 1L
    if (col.width[j] < 10L) {
      col.width[j] <- 10L
    } else if (col.width[j] > 50L) {
      col.width[j] <- 50L
    }
  }

  # Assigin global variables
  undo.stack <- NULL
  redo.stack <- NULL
  matched.cells <- NULL
  search.defaults <- list(is.match.word=FALSE, is.match.case=TRUE,
                          is.reg.exps=FALSE, is.search.sel=FALSE,
                          is.perl=FALSE)

  # Assign variables linked to Tk widgets
  table.var   <- tclArray()
  record.var  <- tclVar()
  value.var   <- tclVar()
  pattern.var <- tclVar()
  tt.done.var <- tclVar(0)

  # Open GUI
  tclServiceMode(FALSE)
  tt <- tktoplevel()
  if (!is.null(parent)) {
    tkwm.transient(tt, parent)
    geo <- unlist(strsplit(as.character(tkwm.geometry(parent)), "\\+"))
    tkwm.geometry(tt, paste0("+", as.integer(geo[2]) + 25,
                             "+", as.integer(geo[3]) + 25))
  }
  tktitle(tt) <- win.title

  # Start top menu
  top.menu <- tkmenu(tt, tearoff=0)

  # Edit menu
  menu.edit <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="Edit", menu=menu.edit, underline=0)
  if (!read.only) {
    tkadd(menu.edit, "command", label="Undo", accelerator="Ctrl+z",
          command=UndoEdit)
    tkadd(menu.edit, "command", label="Redo", accelerator="Ctrl+y",
          command=RedoEdit)
    tkadd(menu.edit, "separator")
  }
  tkadd(menu.edit, "command", label="Copy", accelerator="Ctrl+c",
        command=BypassCopyCmd)
  if (!read.only) {
    tkadd(menu.edit, "command", label="Cut", accelerator="Ctrl+x",
          command=BypassCutCmd)
    tkadd(menu.edit, "command", label="Paste", accelerator="Ctrl+v",
          command=BypassPasteCmd)
    tkadd(menu.edit, "separator")
    menu.edit.del <- tkmenu(tt, tearoff=0)
    tkadd(menu.edit.del, "command", label="Character after cursor",
          accelerator="Delete",
          command=function() tkevent.generate(frame3.tbl, "<Delete>"))
    tkadd(menu.edit.del, "command", label="Character before cursor",
          accelerator="Backspace",
          command=function() tkevent.generate(frame3.tbl, "<BackSpace>"))
    tkadd(menu.edit.del, "command", label="All characters after cursor",
          accelerator="Ctrl+k",
          command=function() tkevent.generate(frame3.tbl, "<Control-k>"))
    tkadd(menu.edit, "cascade", label="Inside cell delete", menu=menu.edit.del)
  }
  tkadd(menu.edit, "separator")
  menu.edit.width <- tkmenu(tt, tearoff=0)
  tkadd(menu.edit.width, "command", label="Increase", accelerator="Ctrl+\u003d",
        command=function() tkevent.generate(frame3.tbl, "<Control-equal>"))
  tkadd(menu.edit.width, "command", label="Decrease", accelerator="Ctrl+\u2212",
        command=function() tkevent.generate(frame3.tbl, "<Control-minus>"))
  tkadd(menu.edit, "cascade", label="Column width", menu=menu.edit.width)

# View menu
  menu.view <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="View", menu=menu.view, underline=0)
  tkadd(menu.view, "command", label="Structure",
        command=function() ViewInfo("Structure"))
  tkadd(menu.view, "command", label="Summary",
        command=function() ViewInfo("Summary"))
  if (!read.only) {
    tkadd(menu.view, "separator")
    tkadd(menu.view, "command", label="Change log", command=ViewChangeLog)
  }

  # Search menu
  menu.search <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="Search", menu=menu.search, underline=0)
  tkadd(menu.search, "command", label="Find\u2026", accelerator="Ctrl+f",
        command=function() CallSearch())
  tkadd(menu.search, "command", label="Find next",
        command=function() Find("next"))
  tkadd(menu.search, "command", label="Find previous",
        command=function() Find("prev"))
  if (!read.only)
    tkadd(menu.search, "command", label="Replace\u2026", accelerator="Ctrl+r",
          command=function() CallSearch(is.replace=TRUE))

  # Selection menu
  menu.sel <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="Select", menu=menu.sel, underline=0)
  tkadd(menu.sel, "command", label="Select all cells",
        accelerator="Ctrl+\u2044",
        command=function() tkevent.generate(frame3.tbl, "<Control-slash>"))
  tkadd(menu.sel, "separator")
  menu.sel.extend <- tkmenu(tt, tearoff=0)
  tkadd(menu.sel.extend, "command", label="First cell",
        accelerator="Shift+Ctrl+Home",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Control-Home>"))
  tkadd(menu.sel.extend, "command", label="Last cell",
        accelerator="Shift+Ctrl+End",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Control-End>"))
  tkadd(menu.sel.extend, "separator")
  tkadd(menu.sel.extend, "command", label="Row above",
        accelerator="Shift+\u2191",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Up>"))
  tkadd(menu.sel.extend, "command", label="Row below",
        accelerator="Shift+\u2193",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Down>"))
  tkadd(menu.sel.extend, "command", label="Column left",
        accelerator="Shift+\u2190",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Left>"))
  tkadd(menu.sel.extend, "command", label="Column right",
        accelerator="Shift+\u2192",
        command=function() tkevent.generate(frame3.tbl, "<Shift-Right>"))
  tkadd(menu.sel, "cascade", label="Extend selection to", menu=menu.sel.extend)

  # Navigation menu
  menu.nav <- tkmenu(tt, tearoff=0, relief="flat")
  tkadd(top.menu, "cascade", label="Navigate", menu=menu.nav, underline=0)
  tkadd(menu.nav, "command", label="Move up", accelerator="\u2191",
        command=function() tkevent.generate(frame3.tbl, "<Up>"))
  tkadd(menu.nav, "command", label="Move down", accelerator="\u2193",
        command=function() tkevent.generate(frame3.tbl, "<Down>"))
  tkadd(menu.nav, "command", label="Move left", accelerator="\u2190",
        command=function() tkevent.generate(frame3.tbl, "<Left>"))
  tkadd(menu.nav, "command", label="Move right", accelerator="\u2192",
        command=function() tkevent.generate(frame3.tbl, "<Right>"))
  tkadd(menu.nav, "separator")
  if (read.only) {
    menu.nav.view <- tkmenu(tt, tearoff=0)
    tkadd(menu.nav.view, "command", label="Prior page in view",
          accelerator="PageUp",
          command=function() tkevent.generate(frame3.tbl, "<Prior>"))
    tkadd(menu.nav.view, "command", label="Next page in view",
          accelerator="PageDown",
          command=function() tkevent.generate(frame3.tbl, "<Next>"))
    tkadd(menu.nav, "cascade", label="Move table to have", menu=menu.nav.view)
    tkadd(menu.nav, "separator")
    menu.nav.active <- tkmenu(tt, tearoff=0)
    tkadd(menu.nav.active, "command", label="First cell",
          accelerator="Ctrl+Home",
          command=function() tkevent.generate(frame3.tbl, "<Control-Home>"))
    tkadd(menu.nav.active, "command", label="Last cell", accelerator="Ctrl+End",
          command=function() tkevent.generate(frame3.tbl, "<Control-End>"))
    tkadd(menu.nav, "cascade", label="Move activate cell to",
          menu=menu.nav.active)
  } else {
    menu.nav.view <- tkmenu(tt, tearoff=0)
    tkadd(menu.nav.view, "command", label="First cell in view",
          accelerator="Home",
          command=function() tkevent.generate(frame3.tbl, "<Home>"))
    tkadd(menu.nav.view, "command", label="Last cell in view",
          accelerator="End",
          command=function() tkevent.generate(frame3.tbl, "<End>"))
    tkadd(menu.nav.view, "separator")
    tkadd(menu.nav.view, "command", label="Prior page in view",
          accelerator="Ctrl+PageUp",
          command=function() tkevent.generate(frame3.tbl, "<Control-Prior>"))
    tkadd(menu.nav.view, "command", label="Next page in view",
          accelerator="Ctrl+PageDown",
          command=function() tkevent.generate(frame3.tbl, "<Control-Next>"))
    tkadd(menu.nav, "cascade", label="Move table to have", menu=menu.nav.view)
    tkadd(menu.nav, "separator")
    menu.nav.active <- tkmenu(tt, tearoff=0)
    tkadd(menu.nav.active, "command", label="First cell",
          accelerator="Ctrl+Home",
          command=function() tkevent.generate(frame3.tbl, "<Control-Home>"))
    tkadd(menu.nav.active, "command", label="Last cell", accelerator="Ctrl+End",
          command=function() tkevent.generate(frame3.tbl, "<Control-End>"))
    tkadd(menu.nav.active, "separator")
    tkadd(menu.nav.active, "command", label="Prior page", accelerator="PageUp",
          command=function() tkevent.generate(frame3.tbl, "<Prior>"))
    tkadd(menu.nav.active, "command", label="Next page", accelerator="PageDown",
          command=function() tkevent.generate(frame3.tbl, "<Next>"))
    tkadd(menu.nav, "cascade", label="Move active cell to",
          menu=menu.nav.active)
    tkadd(menu.nav, "separator")
    menu.nav.in <- tkmenu(tt, tearoff=0)
    tkadd(menu.nav.in, "command", label="Left", accelerator="Ctrl+\u2190",
          command=function() tkevent.generate(frame3.tbl, "<Control-Left>"))
    tkadd(menu.nav.in, "command", label="Right", accelerator="Ctrl+\u2192",
          command=function() tkevent.generate(frame3.tbl, "<Control-Right>"))
    tkadd(menu.nav.in, "separator")
    tkadd(menu.nav.in, "command", label="Beggining", accelerator="Ctrl+a",
          command=function() tkevent.generate(frame3.tbl, "<Control-a>"))
    tkadd(menu.nav.in, "command", label="End", accelerator="Ctrl+e",
          command=function() tkevent.generate(frame3.tbl, "<Control-e>"))
    tkadd(menu.nav, "cascade", label="Move inside cell to the",
          menu=menu.nav.in)
  }

  # Graph menu
  menu.graph <- tkmenu(tt, tearoff=0)
  tkadd(top.menu, "cascade", label="Graph", menu=menu.graph, underline=0)
  tkadd(menu.graph, "command", label="Build histogram\u2026",
        command=PlotHistogram)

  # Finish top menu
  tkconfigure(tt, menu=top.menu)

  # Frame 0
  frame0 <- ttkframe(tt, relief="flat")
  frame0.ent.1.1 <- ttkentry(frame0, width=10, font="TkFixedFont",
                             state=if (read.only) "readonly" else "normal",
                             textvariable=value.var, validate="key",
                             validatecommand=function(P, S) {
                                               ValidateEntryValue(P, S)
                                             })
  tkgrid(frame0.ent.1.1, padx=c(10, 25), pady=c(10, 4), sticky="we")
  tkgrid.columnconfigure(frame0, 0, weight=1)
  tkpack(frame0, fill="x", side="top")

  # Frame 1

  frame1 <- ttkframe(tt, relief="flat")

  if (read.only) {
    frame1.but.1.2 <- "x"
    frame1.but.1.3 <- ttkbutton(frame1, width=12, text="Close",
                                command=function() tclvalue(tt.done.var) <- 1)
  } else {
    frame1.but.1.2 <- ttkbutton(frame1, width=12, text="Save",
                                command=SaveTable)
    frame1.but.1.3 <- ttkbutton(frame1, width=12, text="Cancel",
                                command=function() tclvalue(tt.done.var) <- 1)
  }
  frame1.but.1.4 <- ttkbutton(frame1, width=12, text="Help",
                              command=function() {
                                print(help("EditData", package="RSurvey"))
                              })
  frame1.grp.1.5 <- ttksizegrip(frame1)

  tkgrid("x", frame1.but.1.2, frame1.but.1.3, frame1.but.1.4, frame1.grp.1.5)

  tkgrid.columnconfigure(frame1, 0, weight=1)

  if (!read.only)
    tkgrid.configure(frame1.but.1.2, padx=c(10, 0))
  tkgrid.configure(frame1.but.1.3, padx=c(4, 0))
  tkgrid.configure(frame1.but.1.4, pady=10, padx=c(4, 10), columnspan=2)
  tkgrid.configure(frame1.grp.1.5, sticky="se")

  tkraise(frame1.but.1.4, frame1.grp.1.5)

  tkpack(frame1, fill="x", side="bottom", anchor="e")

  # Frame 2

  frame2 <- ttkframe(tt, relief="flat", padding=0, borderwidth=0, height=200)

  frame2.lab.1.1 <- ttklabel(frame2, text="Record")
  frame2.lab.2.1 <- ttklabel(frame2, text="Find")

  frame2.ent.1.2 <- ttkentry(frame2, width=15, font="TkFixedFont",
                             textvariable=record.var)
  frame2.ent.2.2 <- ttkentry(frame2, width=15, font="TkFixedFont",
                             textvariable=pattern.var)

  frame2.but.1.3 <- ttkbutton(frame2, width=4, text="View", command=ViewRecord)
  frame2.but.2.3 <- ttkbutton(frame2, width=2, image=GetBitmapImage("previous"),
                              command=function() Find("prev"))
  frame2.but.2.4 <- ttkbutton(frame2, width=2, image=GetBitmapImage("next"),
                              command=function() Find("next"))

  tkgrid(frame2.lab.1.1, frame2.ent.1.2, frame2.but.1.3, pady=c(0, 4))
  tkgrid(frame2.lab.2.1, frame2.ent.2.2, frame2.but.2.3, frame2.but.2.4)

  tkgrid.configure(frame2.ent.2.2, frame2.ent.1.2, padx=c(0, 2))
  tkgrid.configure(frame2.lab.2.1, frame2.lab.1.1, padx=c(0, 2), sticky="w")
  tkgrid.configure(frame2.but.1.3, columnspan=2, padx=c(0, 70), sticky="we")
  tkgrid.configure(frame2.but.2.4, padx=c(2, 70))

  tkpack(frame2, side="bottom", anchor="nw", padx=c(10, 0))

  # Frame 3

  frame3 <- ttkframe(tt, relief="flat", padding=0, borderwidth=0)

  frame3.tbl <- tkwidget(frame3, "table", rows=m + 1, cols=n + 1,
                         colwidth=-2, rowheight=1,
                         state=ifelse(read.only, "disabled", "normal"),
                         height=nrows + 1, width=ncols + 1, ipadx=3, ipady=1,
                         wrap=0, justify="left", background="#FFFFFF",
                         foreground="#000000", titlerows=1, titlecols=1,
                         multiline=0, resizeborders="col", colorigin=0,
                         bordercursor="sb_h_double_arrow", cursor="plus",
                         colstretchmode="none", rowstretchmode="none",
                         drawmode="single", flashmode=0, rowseparator="\n",
                         colseparator="\t", selectmode="extended",
                         selecttitle=0, insertofftime=0, anchor="nw",
                         highlightthickness=0, cache=1, validate=1,
                         font="TkFixedFont", exportselection=0,
                         browsecommand=function(s, S) ChangeActiveCell(s, S),
                         validatecommand=function(s, S) ValidateCellValue(s, S),
                         command=function(r, c) GetCellValue(r, c),
                         xscrollcommand=function(...) tkset(frame3.xsc, ...),
                         yscrollcommand=function(...) tkset(frame3.ysc, ...))

  frame3.xsc <- ttkscrollbar(frame3, orient="horizontal",
                             command=function(...) tkxview(frame3.tbl, ...))
  frame3.ysc <- ttkscrollbar(frame3, orient="vertical",
                             command=function(...) tkyview(frame3.tbl, ...))

  tcl(frame3.tbl, "width",  0, col.0.width)
  tcl(frame3.tbl, "height", 0, row.0.height)
  for (j in seq_len(n))
    tcl(frame3.tbl, "width", j, col.width[j])

  tkgrid(frame3.tbl, frame3.ysc)
  tkgrid(frame3.xsc, "x")

  tkgrid.configure(frame3.tbl, padx=c(10,  0), pady=c(0, 0), sticky="news")
  tkgrid.configure(frame3.ysc, padx=c( 0, 10), pady=c(0, 0), sticky="ns")
  tkgrid.configure(frame3.xsc, padx=c(10,  0), pady=c(0, 5), sticky="we")

  tktag.configure(frame3.tbl, "active", background="#FBFCD0")
  tktag.configure(frame3.tbl, "sel", background="#EAEEFE",
                  foreground="#000000")
  tktag.configure(frame3.tbl, "title", background="#D9D9D9",
                  foreground="#000000")

  tcl(frame3.tbl, "tag", "row", "coltitles", 0)
  tcl(frame3.tbl, "tag", "col", "rowtitles", 0)

  tktag.configure(frame3.tbl, "coltitles", anchor="nw")
  tktag.configure(frame3.tbl, "rowtitles", anchor="n")

  tktag.raise(frame3.tbl, "title", "sel")

  tkgrid.columnconfigure(frame3, 0, weight=1)
  tkgrid.rowconfigure(frame3, 0, weight=1)

  tkpack(frame3, fill="both", expand=TRUE)

  # Bind events

  tclServiceMode(TRUE)

  tkbind(tt, "<Destroy>", function() tclvalue(tt.done.var) <- 1)

  tkbind(frame0.ent.1.1, "<Return>",
         paste(.Tcl.callback(BypassReturnCmd), "break", sep="; "))
  tkbind(frame0.ent.1.1, "<FocusIn>", function() tksee(frame3.tbl, "active"))

  tkbind(tt, "<Control-f>", function() CallSearch(is.replace=FALSE))
  tkbind(tt, "<Control-r>", function() CallSearch(is.replace=TRUE))

  tkbind(frame3.tbl, "<Control-x>",
         paste(.Tcl.callback(BypassCutCmd), "break", sep="; "))
  tkbind(frame3.tbl, "<Control-v>",
         paste(.Tcl.callback(BypassPasteCmd), "break", sep="; "))
  tkbind(frame3.tbl, "<Return>",
         paste(.Tcl.callback(BypassReturnCmd), "break", sep="; "))
  tkbind(frame3.tbl, "<Control-c>",
         paste(.Tcl.callback(BypassCopyCmd), "break", sep="; "))
  tkbind(frame3.tbl, "<Control-z>", UndoEdit)
  tkbind(frame3.tbl, "<Control-y>", RedoEdit)
  tkbind(frame3.tbl, "<Double-Button-1>", function(x, y) ResizeColumn(x, y))

  tkevent.add("<<TableSelect>>", "<ButtonRelease>",
              "<Control-slash>", "<Shift-Control-Home>", "<Shift-Control-End>",
              "<Shift-Up>", "<Shift-Down>", "<Shift-Right>", "<Shift-Left>",
              "<Up>", "<Down>", "<Right>", "<Left>")
  tkbind(frame3.tbl, "<<TableSelect>>", ChangeCellSelection)

  D <- ""  # hack to force 'D' to be something other than a function
  tkbind(frame3.tbl, "<MouseWheel>",
         function(D) {
           number <- as.integer((-as.integer(D) / 120)^3)
           tkyview(frame3.tbl, "scroll", number, "units")
         })

  tkbind(frame2.ent.2.2, "<KeyRelease>",
         function() {
           matched.cells <<- NULL
         })
  tkbind(frame2.ent.2.2, "<Return>", function() Find("next"))
  tkbind(frame2.ent.2.2, "<Up>",     function() Find("prev"))
  tkbind(frame2.ent.2.2, "<Down>",   function() Find("next"))
  tkbind(frame2.ent.1.2, "<Return>", function() ViewRecord())

  # GUI control

  tkfocus(frame3.tbl)
  tkactivate(frame3.tbl, "origin")
  tkselection.set(frame3.tbl, "active")
  tksee(frame3.tbl, "active")

  tkgrab(tt)
  tkwait.variable(tt.done.var)

  tclServiceMode(FALSE)
  if (identical(changelog.old, changelog)) {
    d <- NULL
  } else {
    if (rtn.class == "data.frame") {
      class(d) <- "data.frame"  # see warning in utils::read.table (R v3.0.2)
      rownames(d) <- d.rownames
      colnames(d) <- d.colnames
    } else if (rtn.class == "matrix") {
      d <- do.call(cbind, d)
      rownames(d) <- d.rownames
      colnames(d) <- d.colnames
    }
  }

  tkgrab.release(tt)
  tkdestroy(tt)
  tclServiceMode(TRUE)

  if (is.null(d))
    invisible(NULL)
  else
    invisible(list(d=d, changelog=changelog))
}
