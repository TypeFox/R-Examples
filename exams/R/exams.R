## workhorse function for compiling (collections of) exercises
exams <- function(file, n = 1, nsamp = NULL, dir = NULL, template = "plain",
  inputs = NULL, header = list(Date = Sys.Date()), name = NULL,
  quiet = TRUE, edir = NULL, tdir = NULL, control = NULL)
{  
  ## convenience function
  strip_path <- function(file)
    sapply(strsplit(file, .Platform$file.sep), tail, 1)

  ## manage directories: 
  ##   - for producing several files an output directory is required
  if((n > 1 | length(template) > 1) & is.null(dir)) stop("Please specify an output 'dir'.")
  if(!is.null(dir) && !file.exists(dir) && !dir.create(dir))
    stop(gettextf("Cannot create output directory '%s'.", dir))
  ##   - further: dir (output), dir_orig (original), dir_temp (temp), dir_pkg (package)
  if(!is.null(dir)) dir <- file_path_as_absolute(dir)
  dir_orig <- getwd()
  on.exit(setwd(dir_orig))
  dir_temp <- if(is.null(tdir)) tempfile() else tdir
  if(!file.exists(dir_temp) && !dir.create(dir_temp))
    stop(gettextf("Cannot create temporary work directory '%s'.", dir_temp))
  dir_pkg <- find.package("exams")
  
  ## number of available exercises in each element of 'file'
  ## and number of selected samples per element
  nfile <- length(file)
  if(is.null(nsamp)) nsamp <- 1
  if(length(nsamp) < nfile) nsamp <- rep(nsamp, length.out = nfile)
  navail <- sapply(file, length)  
  if(any(navail < nsamp)) {
    ix <- which(navail < nsamp)
    warning(paste("Only", navail[ix], "exercise(s) available in element", ix,
      "of the 'file' argument. Sampling with replacement will be used in order to obtain",
      nsamp[ix], "replications."))
  }
  
  ## file pre-processing:
  ##   - transform to vector (remember grouping IDs)
  ##   - add paths (generate "foo", "foo.Rnw", "foo.tex", and "path/to/foo.Rnw")
  ##   - check existence (use local files if they exist, otherwise take from package)
  ##   - setup sampling (draw random configuration)
  file_id <- rep(seq_along(file), navail)

  file_raw <- unlist(file)
  file_Rnw <- ifelse(
    tolower(substr(file_raw, nchar(file_raw)-3, nchar(file_raw))) != ".rnw",
    paste(file_raw, ".Rnw", sep = ""), file_raw)
  file_base <- file_path_sans_ext(file_Rnw)
  file_tex <- paste(file_base, ".tex", sep = "")
  file_path <- if(is.null(edir)) file_Rnw else file.path(edir, file_Rnw)
  file_path <- ifelse(file.exists(file_path),
    file_path, file.path(dir_pkg, "exercises", file_path))
  if(!all(file.exists(file_path))) stop(paste("The following files cannot be found: ",
    paste(file_raw[!file.exists(file_path)], collapse = ", "), ".", sep = ""))

  sample_id <- function() unlist(lapply(unique(file_id), function(i) {
    wi <- file_id == i
    if(sum(wi) > 1)
      sample(which(wi), nsamp[i], replace = navail[i] < nsamp[i])
    else
      rep(which(wi), length.out = nsamp[i])
  }))
  
  ## similarly: template pre-processing
  template_raw <- template
  template_tex <- template_path <- ifelse(
    tolower(substr(template, nchar(template)-3, nchar(template))) != ".tex",
    paste(template, ".tex", sep = ""), template)
  template_base <- file_path_sans_ext(template_tex)
  template_path <- ifelse(file.exists(template_tex),
    template_tex, file.path(dir_pkg, "tex", template_tex))
  if(!all(file.exists(template_path))) stop(paste("The following files cannot be found: ",
    paste(template_raw[!file.exists(template_path)], collapse = ", "), ".", sep = ""))  

  ## check for using old templates
  if(file.path(dir_pkg, "tex", "exam.tex") %in% template_path) {
    template_path[template_path == file.path(dir_pkg, "tex", "exam.tex")] <- file.path(dir_pkg, "tex", "oexam.tex")
    warning(paste(strwrap(paste(
      "The template exam.tex has been adapted to exams2pdf() and is not fully compatible",
      "with exams() anymore. Template oexam.tex used instead."
      ), exdent = 2), collapse = "\n"))
  }
  if(file.path(dir_pkg, "tex", "solution.tex") %in% template_path) {
    template_path[template_path == file.path(dir_pkg, "tex", "solution.tex")] <- file.path(dir_pkg, "tex", "osolution.tex")
    warning(paste(strwrap(paste(
      "The template solution.tex has been adapted to exams2pdf() and is not fully compatible",
      "with exams() anymore. Template osolution.tex used instead."
      ), exdent = 2), collapse = "\n"))
  }

  ## read template
  template <- lapply(template_path, readLines)
  ## which input types in template?
  input_types <- function(x) {
    x <- x[grep("\\exinput", x, fixed = TRUE)]
    if(length(x) < 1) stop("templates must specify at least one \\exinput{}")
    as.vector(sapply(strsplit(sapply(strsplit(x,
      paste("\\exinput{", sep = ""), fixed = TRUE), tail, 1), "}"), head, 1))
  }
  template_it <- lapply(template, input_types)
  template_has_header <- sapply(template_it, function(x) "header" %in% x)
  template_has_questionnaire <- sapply(template_it, function(x) "questionnaire" %in% x)
  template_has_exercises <- sapply(template_it, function(x) "exercises" %in% x)

  ## output name processing
  if(is.null(name)) name <- strip_path(template_base)
  make_full_name <- function(name, id, type = "")
    paste(name, gsub(" ", "0", format(c(n, id)))[-1], ifelse(type == "", "", "."), type, sep = "")

  ## convenience function for reading metainfo from compiled exercise
  read_metainfo <- function(file) {
    x <- readLines(file)
    get_command <- function(command) {
      cline <- x[grep(command, x, fixed = TRUE)]
      if(length(cline) < 1) NULL else gsub("{", "", strsplit(strsplit(cline[1],
        paste(command, "{", sep = ""), fixed = TRUE)[[1]][2], "}")[[1]], fixed = TRUE)
    }
    type <- match.arg(get_command("\\extype"), c("schoice", "mchoice", "num", "string"))
    sol <- get_command("\\exsolution")
    nam <- get_command("\\exname")
    tol <- get_command("\\extol")
    tol <- rep(if(is.null(tol)) 0 else as.numeric(tol), length.out = 2)
    sol <- switch(type,
      "schoice" = string2mchoice(sol, single = TRUE),
      "mchoice" = string2mchoice(sol),
      "num" = as.numeric(sol),
      "string" = sol
    )
    slength <- length(sol)
    string <- switch(type,
      "schoice" = paste(nam, ": ", mchoice2print(sol), sep = ""),
      "mchoice" = paste(nam, ": ", mchoice2print(sol), sep = ""),
      "num" = if(max(tol) <= 0) {
        paste(nam, ": ", sol, sep = "")
      } else {
        if(slength == 1) {
	  paste(nam, ": ", sol, " (", sol - tol[1], "--", sol + tol[2], ")", sep = "")
	} else {
          paste(nam, ": [", sol[1], ", ", sol[2], "] ([", sol[1] - tol[1], "--", sol[1] + tol[2], ", ",
            sol[2] - tol[1], "--", sol[2] + tol[2], "])", sep = "")
	}
      },
      "string" = paste(nam, ": ", paste(sol, collapse = "\n"), sep = "")
    )

    list(type = type,
         length = slength,
         solution = sol,
	 tolerance = tol,
         string = string)
  }

  control.default <- list(mchoice.print = list(True = letters[1:5], False = rep("", 5)),
                          mchoice.symbol = c(True = "X", False = " "))
  if (is.null(control)) control <- control.default
  else if (is.list(control)) {
    control <- c(control, control.default[!c("mchoice.print", "mchoice.symbol") %in% names(control)])
    if (!all(sapply(control, function(x) identical(c("False", "True"), sort(names(x))))))
      stop("'control' not correctly specified")
    control$mchoice.print <- lapply(control$mchoice.print, function(x) rep(x, 5))
  } else stop("'control' must be NULL or a list")
  
  ## convenience functions for writing LaTeX  
  mchoice2quest <- function(x) paste("  \\item \\exmchoice{",
    paste(ifelse(x, control$mchoice.symbol[["True"]], control$mchoice.symbol[["False"]]), collapse = "}{"), "}", sep = "")
  num2quest <- function(x) {
    rval <-  paste("  \\item \\exnum{", 
      paste(strsplit(format(c(100000.000, x), nsmall = 3, scientific = FALSE)[-1], "")[[1]][-7],
      collapse = "}{"), "}", sep = "")
    if(length(x) > 1) rval <- paste(rval, " \\\\\n        \\exnum{",
      paste(strsplit(format(c(100000.000, x), nsmall = 3, scientific = FALSE)[-1], "")[[2]][-7],
      collapse = "}{"), "}", sep = "")
    rval 
  }
  string2quest <- function(x) paste("  \\item \\exstring{", x, "}", sep = "")
  
  mchoice2print <- function(x) paste(ifelse(x, control$mchoice.print[["True"]], control$mchoice.print[["False"]]), collapse = "")

  ## take everything to temp dir
  file.copy(file_path, dir_temp)
  ## including further inputs (if any)
  if(!is.null(inputs)) {
    inputs_path <- ifelse(file.exists(inputs), inputs, file.path(edir, inputs))
    if(!all(file.exists(inputs_path))) stop(paste("The following inputs cannot be found: ",
      paste(inputs[!file.exists(inputs_path)], collapse = ", "), ".", sep = ""))
    file.copy(inputs_path, dir_temp)
  }
  setwd(dir_temp) 
  on.exit(unlink(dir_temp), add = TRUE)
  
  ## call Sweave and LaTeX, copy and collect results
  metainfo <- list()  
  for(i in 1:n) {
  
    ## select exercise files, run Sweave, collect results
    id <- sample_id()
    for(j in id) Sweave(file_Rnw[j], quiet = quiet) ## FIXME: need envir argument
    metainfo1 <- list()
    for(j in seq_along(id)) metainfo1[[j]] <- read_metainfo(file_tex[id[j]])
    names(metainfo1) <- file_base[id]
    metainfo[[i]] <- metainfo1

    ## assign names
    names(metainfo)[i] <- make_full_name(name[1], i)
    out_tex <- make_full_name(name, i, type = "tex")
    out_pdf <- make_full_name(name, i, type = "pdf")

    ## compile output files for all templates
    for(j in seq_along(template)) {
      tmpl <- template[[j]]

      ## input header
      if(template_has_header[j]) {
        wi <-  grep("\\exinput{header}", tmpl, fixed = TRUE)
        tmpl[wi] <- if(length(header) < 1) "" else paste("\\", names(header), "{",
          sapply(header, function(x) if(is.function(x)) x(i) else as.character(x)), "}", 
	  collapse = "\n", sep = "")
      }

      ## input questionnaire
      if(template_has_questionnaire[j]) {
	typ1 <- sapply(metainfo1, function(x) x[["type"]])
	sol1 <- lapply(metainfo1, function(x) x[["solution"]])
        wi <-  grep("\\exinput{questionnaire}", tmpl, fixed = TRUE)
	tmpl[wi] <- paste(c(
	  "\\begin{enumerate}",
	  sapply(seq_along(typ1), function(i) switch(typ1[i],
	    schoice = mchoice2quest(sol1[[i]]),
	    mchoice = mchoice2quest(sol1[[i]]),
            num =  num2quest(sol1[[i]]),
            string = string2quest(sol1[[i]]))),
          "\\end{enumerate}", ""), collapse = "\n")
      }

      ## input exercise tex
      if(template_has_exercises[j]) {
        wi <-  grep("\\exinput{exercises}", tmpl, fixed = TRUE)
        tmpl[wi] <- paste("\\input{", file_tex[id], "}", sep = "", collapse = "\n")
      }

      ## create and compile output tex
      writeLines(tmpl, out_tex[j])
      texi2dvi(out_tex[j], pdf = TRUE, clean = TRUE, quiet = quiet)
    }

    ## copy to output directory (or show)
    if(!is.null(dir)) {
      file.copy(out_pdf, dir, overwrite = TRUE)
    } else {
      if(.Platform$OS.type == "windows") shell.exec(file.path(dir_temp, out_pdf))
        else system(paste(shQuote(getOption("pdfviewer")), shQuote(out_pdf)), wait = FALSE)
    } 
  }

  ## collect and store meta information
  class(metainfo) <- "exams_metainfo"  
  if(!is.null(dir)) {
    save(metainfo, file = file.path(dir, "metainfo.rda"))
    ## metainfo_df <- as.data.frame(t(sapply(metainfo,
    ##   function(x) as.vector(sapply(x, function(y) y$string)))))
    ## colnames(metainfo_df) <- paste("exercise", gsub(" ", "0", format(1:ncol(metainfo_df))), sep = "")
    ## write.table(metainfo_df, file = file.path(dir, "metainfo.csv"), sep = ",")
  }

  ## return meta information invisibly
  invisible(metainfo)
}

## print exams_metainfo objects
print.exams_metainfo <- function(x, which = NULL, ...) {
  which <- if(is.null(which)) names(x) else {
    if(is.numeric(which)) names(x)[which] else which
  }
  n <- length(x[[1]])
  for(i in which) {
    cat("\n", i, "\n", sep = "")
    for(j in 1:n) {
      cat("    ", format(c(n, j))[-1], ". ", x[[i]][[j]]$string, "\n", sep = "")
    }
  }
  cat("\n")
  invisible(x)
}
