
#' Install a new challenge.
#' @param path         string. install path of the challenge (should be somewhere in your Dropbox).
#' @param out_rmdfile  string. name of the output R Markdown file.
#' @param recursive    logical. should elements of the path other than the last be created? see \code{\link{dir.create}}.
#' @param overwrite    logical. should existing destination files be overwritten? see \code{\link{file.copy}}.
#' @param quiet        logical. deactivate text output.
#' @param showWarnings logical. should the warnings on failure be shown? see \code{\link{dir.create}}.
#' @param data_dir     string. subdirectory of the data.
#' @param submissions_dir string. subdirectory of the submissions. see \code{\link{store_new_submissions}}.
#' @param hist_dir     string. subdirectory of the history. see \code{\link{store_new_submissions}}.
#' @param install_data logical. activate installation of the data files of the template challenge.
#' @param baseline     string. name of the team considered as the baseline.
#' @param add_baseline logical. activate installation of baseline submission files of the template challenge.
#' @param clear_history logical. activate deletion of the existing history folder.
#' @param template     string. name of the template R Markdown script to be installed.
#'   Two choices are available: \code{"en"} (english) and \code{"fr"} (french).
#' @param title        string. title displayed on the webpage.
#' @param author       string. author displayed on the webpage.
#' @param date         string. date displayed on the webpage.
#' @param email        string. email of the challenge administrator.
#' @param date_start   string. start date of the challenge.
#' @param deadline     string. deadline of the challenge.
#' @param data_list    list with members \code{train}, \code{test}, \code{y_test} and
#'   \code{ind_quiz} such as returned by the \code{\link{data_split}} function.
#' @return The path of the created challenge is returned.
#' @export
#' @examples
#' path <- tempdir()
#' wd <- setwd(path)
#' # english version
#' new_challenge()
#' # french version
#' new_challenge(template = "fr")
#' setwd(wd)
#' unlink(path)
new_challenge <- function(path = ".", out_rmdfile = "challenge.rmd", 
                          recursive = FALSE, overwrite = recursive, 
                          quiet = FALSE, showWarnings = FALSE,
                          template = c("en", "fr"),
                          data_dir = "data", 
                          submissions_dir = "submissions",
                          hist_dir = "history", 
                          install_data = TRUE, 
                          baseline = "baseline",
                          add_baseline = install_data,
                          clear_history = overwrite,
                          title = "Challenge",
                          author = "",
                          date = "",
                          email = "EDIT_EMAIL@DOMAIN.com",
                          date_start = format(Sys.Date(), "%d %b %Y"),
                          deadline = paste(Sys.Date()+90, "23:59:59"),
                          data_list = data_split(get_data("german"))) {
  
  dir.create(path, recursive = recursive, showWarnings = showWarnings)
  if (!file.exists(path))
    stop("could not create directory ", path)
  
  stopifnot(is.character(out_rmdfile), length(out_rmdfile)==1, nzchar(out_rmdfile))
  
  stopifnot(is.character(template), nzchar(template))
  template = match.arg(template, c("en", "fr"))
  # currently "challenge_en.rmd" and "challenge_fr.rmd" are available
  
  dir.create(file.path(path, "data"), recursive = recursive, showWarnings = showWarnings)
  
  if (install_data) {
    data_train <- data_list$train
    data_test <- data_list$test
    y_test <- data_list$y_test
    ind_quiz <- data_list$ind_quiz
    
    tmpdir = tempdir()
    
    save('data_train', file = file.path(tmpdir, 'data_train.rda'))
    file.copy(file.path(tmpdir, 'data_train.rda'), file.path(path, data_dir),
              overwrite=overwrite, recursive=recursive)
    
    save('data_test', file = file.path(tmpdir, 'data_test.rda'))
    file.copy(file.path(tmpdir, 'data_test.rda'), file.path(path, data_dir),
              overwrite=overwrite, recursive=recursive)
    
    save('y_test', file = file.path(tmpdir, 'y_test.rda'))
    file.copy(file.path(tmpdir, 'y_test.rda'), file.path(path, data_dir),
              overwrite=overwrite, recursive=recursive)
    
    save('ind_quiz', file = file.path(tmpdir, 'ind_quiz.rda'))
    file.copy(file.path(tmpdir, 'ind_quiz.rda'), file.path(path, data_dir),
              overwrite=overwrite, recursive=recursive)
    
    unlink(tmpdir)
  }
  
  dir.create(file.path(path, submissions_dir), recursive = recursive, showWarnings = showWarnings)
  
  if (install_data && add_baseline) {
    team_dir = new_team(baseline, path=path, submissions_dir = submissions_dir,
                        quiet = TRUE, showWarnings = showWarnings)
    
    # Predict all Good
    y_pred <- rep("Good", nrow(data_test))
    tmpfile = tempfile()
    write(y_pred, file = tmpfile)
    file.copy(tmpfile, file.path(team_dir, 'all_good.csv'), overwrite=overwrite)
    
    # Predict all Bad
    y_pred <- rep("Bad", nrow(data_test))
    write(y_pred, file = tmpfile)
    file.copy(tmpfile, file.path(team_dir, 'all_bad.csv'), overwrite=overwrite)
    
    unlink(tmpfile)
  }
  
  if (clear_history)
    unlink(file.path(path, hist_dir), recursive = TRUE)
  
  dir.create(file.path(path, hist_dir), recursive = recursive, showWarnings = showWarnings)
  
  expr = list(title = title, author = author, date = date, email = email, 
              date_start = date_start, deadline = deadline, baseline = baseline,
              data_dir = data_dir, submissions_dir = submissions_dir,
              hist_dir = hist_dir)
  
  # template files are in "template/challenge_<template>.rmd"
  tpl = system.file('template', paste0("challenge_", template, ".rmd"), package = 'rchallenge')
  if (!nzchar(tpl))
    stop("could not find template ", template)
  
  text = readLines(tpl)
  
  for (n in names(expr))
    text = gsub(paste0("@", toupper(n), "@"), expr[[n]], text)
  
  tmpfile = tempfile()
  
  writeLines(text, tmpfile)
  
  file.copy(tmpfile, file.path(path, out_rmdfile), overwrite=overwrite)
  
  unlink(tmpfile)
  
  if (!quiet) {
    cat('New challenge installed in: "', normalizePath(path), '"\n', sep='')
    cat('Next steps to complete the installation:\n')
    step <- 0
    if (install_data) {
      step <- step + 1
      cat(step, '. Replace the data files in the "data" subdirectory.\n', sep='')
    }
    if (add_baseline) {
      step <- step + 1
      cat(step, '. Replace the baseline predictions in "', file.path(submissions_dir, baseline),'".\n', sep='')
    }
    step <- step + 1
    cat(step, '. Customize the template R Markdown file "', out_rmdfile, '" as needed.\n', sep='')
    step <- step + 1
    cat(step, '. Create and share subdirectories in "', submissions_dir, '" for each team:\n', sep='')
    cat('    rchallenge::new_team("team_foo", "team_bar", path="', path, '", submissions_dir="', submissions_dir, '")\n', sep='')
    step <- step + 1
    cat(step, '. Publish the html page in your "Dropbox/Public" folder:\n', sep='')
    cat('    rchallenge::publish("', file.path(path, out_rmdfile), '")\n', sep='')
    step <- step + 1
    template_html <- paste0(sub("([^.]+)\\.[[:alnum:]]+$", "\\1", basename(out_rmdfile)), ".html")
    cat(step, '. Give the Dropbox public link to "Dropbox/Public/', template_html, '" to the participants.\n', sep='')
    step <- step + 1
    cat(step, '. Automate the updates of the webpage.\n', sep='')
    if (.Platform$OS.type == "unix") {
      cat('   On Unix systems, you can setup the following line to your crontab using "crontab -e":\n', sep='')
      cat('    0 * * * * Rscript -e \'rchallenge::publish("', normalizePath(file.path(path, out_rmdfile)), '")\'\n', sep='')
    }
    if (.Platform$OS.type == "windows") {
      cat('   On Windows systems, you can use the Task Scheduler to create a new task with a "Start a program" action with the settings:')
      cat('   - Program/script: Rscript.exe\n')
      cat('   - options: -e rchallenge::publish(\'', normalizePath(file.path(path, out_rmdfile)), '\')\n', sep='')
    }
  }
  
  invisible(normalizePath(path))
}

#' Create new teams submission folders in your challenge.
#' @param ...          strings. names of the team subdirectories.
#' @param path         string. root path of the challenge. see \code{\link{new_challenge}}.
#' @param submissions_dir string. subdirectory of the submissions. see \code{\link{new_challenge}}.
#' @param quiet        logical. deactivate text output.
#' @param showWarnings logical. should the warnings on failure be shown? see \code{\link{dir.create}}.
#' @return The paths of the created teams are returned.
#' @export
#' @examples
#' path <- tempdir()
#' wd <- setwd(path)
#' new_challenge()
#' new_team("team_foo", "team_bar")
#' setwd(wd)
#' unlink(path)
new_team <- function(..., path = ".", submissions_dir = "submissions", 
                     quiet = FALSE, showWarnings = FALSE) {
  names <- c(...)
  stopifnot(is.character(names))
  if (!file.exists(file.path(path, submissions_dir)))
    stop("could not find submissions directory:", normalizePath(file.path(path, submissions_dir)))
  for (i in seq_along(names)) {
    if (!quiet) cat("Creating team subdirectory:", file.path(submissions_dir, names[i]), "\n")
    dir.create(file.path(path, submissions_dir, names[i]), recursive = FALSE, showWarnings = showWarnings)
  }
  if (!quiet) cat("Next step: share the Dropbox folders with the corresponding teams.\n")
  
  invisible(normalizePath(file.path(path, submissions_dir, names)))
}

#' Publish your challenge R Markdown script to a html page.
#' @param input string. name of the R Markdown input file
#' @param output_file output file. If \code{NULL} then a default based on the name 
#'   of the input file is chosen.
#' @param output_dir string. output directory. default=\code{"~/Dropbox/Public"} 
#'   so that the rendered page can easily be shared on the web with Dropbox.
#' @param quiet      logical. deactivate text output.
#' @param ...        further arguments to pass to \code{\link[rmarkdown]{render}}.
#' @return The compiled document is written into the output file, and the path 
#'   of the output file is returned.
#' @export
#' @seealso \code{\link[rmarkdown]{render}}
#' @importFrom rmarkdown render
#' @examples
#' path <- tempdir()
#' wd <- setwd(path)
#' new_challenge()
#' outdir = tempdir()
#' publish(output_dir = outdir, output_options = list(self_contained = FALSE))
#' unlink(outdir)
#' setwd(wd)
#' unlink(path)
publish <- function(input="challenge.rmd", output_file = NULL, 
                    output_dir = file.path("~/Dropbox/Public"), 
                    quiet = FALSE, ...) {
  wd <- getwd()
  setwd(dirname(input))
  out <- rmarkdown::render(input = basename(input), output_file = output_file,
                           output_dir = output_dir, quiet = quiet, ...)
  setwd(wd)
  
  if (!quiet)
    cat('Next step: give the Dropbox public link to "', file.path(output_dir, basename(out)), '" to the participants.\n', sep='')
  
  invisible(out)
}
