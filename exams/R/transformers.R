## helper transformator function,
## includes tex2image(), tth(), ttm(), and pandoc_convert() .html conversion
make_exercise_transform_html <- function(converter = c("ttm", "tth", "pandoc", "tex2image"), base64 = TRUE, ...)
{
  ## process converter (plus options for pandoc)
  options <- strsplit(converter, "-", fixed = TRUE)[[1L]]
  converter <- match.arg(options[1L], c("ttm", "tth", "pandoc", "tex2image"))
  options <- options[-1L]
  options <- if(length(options) > 0L) {
    paste0("--", options, collapse = " ")
  } else {
    "--mathml"
  }
  if(converter %in% c("tth", "ttm")) {
    stopifnot(requireNamespace("tth"))
  } else if(converter == "pandoc") {
    stopifnot(requireNamespace("rmarkdown"))
  }

  ## base64 checks
  if(is.null(base64)) base64 <- TRUE
  base64 <- if(isTRUE(base64)) {
    c("bmp", "gif", "jpeg", "jpg", "png")
  } else {
    if(is.logical(base64)) NA_character_  else tolower(base64)
  }
  if(b64 <- !all(is.na(base64))) stopifnot(requireNamespace("base64enc"))

  if(converter == "pandoc") {
    make_exercise_transform_pandoc(to = "html", base64 = base64, options = options, ...)
  } else if(converter == "tex2image") {
    ## transforms the tex parts of exercise to images
    ## of arbitrary format using function tex2image()
    function(x)
    {
      bsname <- if(is.null(x$metainfo$file)) basename(tempfile()) else x$metainfo$file
      sdir <- attr(x$supplements, "dir")
      images <- list(); inames <- NULL
      for(i in c("question", "questionlist", "solution", "solutionlist")) {
        if(!is.null(x[[i]])) {
          if(grepl("list", i)) {
            images <- c(images, as.list(x[[i]]))
            inames <- c(inames, paste(i, 1:length(x[[i]]), sep = "_"))
          } else {
            images <- c(images, list(x[[i]]))
            inames <- c(inames, i)
          }
        }
      }
      names(images) <- inames
      dir <- tex2image(images, idir = sdir, show = FALSE, name = bsname, ...)
      inames <- file_path_sans_ext(basename(dir))
      if(b64) {
        for(i in seq_along(dir))
          dir[i] <- sprintf('<img src="%s" alt="%s" />', base64enc::dataURI(file = dir[i],
            mime = paste('image', format = file_ext(dir[i]), sep = '/')), dir[i])
        for(sf in dir(sdir)) {
          if(length(grep(file_ext(sf), base64, ignore.case = TRUE))) {
            file.remove(file.path(sdir, sf))
            x$supplements <- x$supplements[!grepl(sf, x$supplements)]
            attr(x$supplements, "dir") <- sdir        
          }
        }
      }
      for(i in c("question", "questionlist", "solution", "solutionlist")) {
        if(!is.null(x[[i]])) {
          if(grepl("list", i)) {
            j <- grep(i, inames)
          } else {
            j <- grep(i, inames)
            j <- j[!grepl("list", inames[j])]
          }
          x[[i]] <- if(b64) {
            dir[j]
          } else {
            paste("<img src=\"", file.path(sdir, basename(dir[j])), "\" alt=\"", inames[j], "\" />", sep = "")
          }
          names(x[[i]]) <- inames[j]
        }
      }
      if(!b64) {
        for(i in dir) {
          fp <- file.path(sdir, basename(i))
          file.copy(i, fp)
          if(!(fp %in% x$supplements))
            x$supplements <- c(x$supplements, fp)
        }
        attr(x$supplements, "dir") <- sdir
      }

      x
    }
  } else {
    ## function to apply ttx() on every
    ## element of a list in a fast way
    apply_ttx_on_list <- function(object, converter = "ttm",
      sep = "\\007\\007\\007\\007\\007", ...)
    {
      ## add seperator as last line to each chunk
      object <- lapply(object, c, sep)

      ## call ttx() on collapsed chunks
      rval <- switch(converter,
        "tth" = tth::tth(unlist(object), ...),
	"ttm" = tth::ttm(unlist(object), ...)
      )
      img <- attr(rval, "images")

      ## split chunks again on sep
      ix <- grepl(sep, rval, fixed = TRUE)
      rval <- split(rval, c(0, head(cumsum(ix), -1L)))

      ## FIXME: length of rval may be smaller than the length of object?
      names(rval) <- rep(names(object), length.out = length(rval))

      ## omit sep from last line in each chunk
      cleansep <- function(x) {
        n <- length(x)
	if(n < 1L) return(x)
        if(x[n] == sep) return(x[-n])
        return(c(x[-n], gsub(sep, "", x[n], fixed = TRUE)))
      }
      rval <- lapply(rval, cleansep)

      ## store ttx images
      attr(rval, "images") <- img

      rval
    } 

    ## exercise conversion with ttx()
    function(x)
    {
      owd <- getwd()
      setwd(sdir <- attr(x$supplements, "dir"))

      ## what need to be transormed with ttx()?
      what <- c(
        "question" = list(x$question),
        "questionlist" = as.list(x$questionlist),
        "solution" = list(x$solution),
        "solutionlist" = as.list(x$solutionlist)
      )

      ## transform the .tex chunks
      args <- list("x" = what, ...)
      trex <- apply_ttx_on_list(what, converter, ...)
      namtrex <- names(trex)

      ## base64 image/supplements handling
      if(b64 && length(sfiles <- dir(sdir))) {
        for(sf in sfiles) {
          for(i in seq_along(trex)) {
            if(length(j <- grep(sf, trex[[i]], fixed = TRUE)) && file_ext(sf) %in% base64) {
              base64i <- fileURI(file = sf)
              trex[[i]][j] <- gsub(paste(sf, '"', sep = ''),
                paste(base64i, '"', sep = ""), trex[[i]][j], fixed = TRUE)
              file.remove(file.path(sdir, sf))
              x$supplements <- x$supplements[!grepl(sf, x$supplements)]
            }
          }
        }
        attr(x$supplements, "dir") <- sdir
      }

      ## replace .tex chunks with tth(), ttm() output
      x$question <- trex$question
      x$questionlist <- sapply(trex[grep("questionlist", namtrex)], paste, collapse = "\n")
      x$solution <- trex$solution
      x$solutionlist <- sapply(trex[grep("solutionlist", namtrex)], paste, collapse = "\n")

      for(j in c("question", "questionlist", "solution", "solutionlist")) {
        if(length(x[[j]]) < 1L) x[[j]] <- NULL
      }

      setwd(owd)

      x
    }
  }
}

