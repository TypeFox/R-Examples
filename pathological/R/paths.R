#' Copy the contents of a directory
#' 
#' Copies the contents of a directory, possibly recursively.
#' @param source_dir String of directory to copy from.
#' @param target_dir String of directory to copy to.
#' @param pattern String regex or \code{NULL}. A filter for filenames, passed  
#' to \code{dir}.
#' @param overwrite Logical value.  Should existing files be overwritten?
#' @param recursive Logical value.  Should subdirectories and their contents 
#' be copied?
#' @param ... Passed from the deprecated \code{dir_copy} to \code{copy_dir}.
#' @note Target directories that don't exist are created, silently (assuming  
#' write permission).
#' @return A logical vector of whether or not each file was successfully  
#' copied is invisibly returned.
#' @seealso \code{\link[base]{basename}}
#' @examples
#' \dontrun{
#' #Copy subdirs by default
#' copy_dir(R.home("etc"), file.path(tempdir(), "etc"))
#' #Just copy the top level
#' copy_dir(R.home("etc"), file.path(tempdir(), "etc2"), recursive = FALSE)
#' #Now copy deeper levels, without overwriting.
#' copy_dir(R.home("etc"), file.path(tempdir(), "etc2"), overwrite = FALSE)
#' #Cleanup
#' unlink(file.path(tempdir(), "etc"), recursive = TRUE)
#' unlink(file.path(tempdir(), "etc2"), recursive = TRUE)
#' }
#' @importFrom assertive.files is_dir
#' @importFrom plyr tryapply
#' @export
copy_dir <- function(source_dir, target_dir, pattern = NULL, overwrite = FALSE, 
                     recursive = TRUE)
{
  #Retrieve all file and directory names
  filenames <- dir(
    source_dir,
    pattern      = pattern, 
    recursive    = recursive, 
    all.files    = TRUE,
    full.names   = FALSE,
    include.dirs = TRUE
  )
  
  #Create missing directories, silently.
  is_directory <- is_dir(file.path(source_dir, filenames))
  directories <- c(target_dir, file.path(target_dir, filenames[is_directory]))
  tryapply(
    directories,
    dir.create,
    showWarnings = FALSE, 
    recursive    = recursive
  )
  
  out_dir <- file.path(target_dir, dirname(filenames[!is_directory]))
  out_dir <- gsub("/\\.$", "", out_dir)   
  
  if(length(out_dir) == 0) return()
  ok <- mapply(
    file.copy,
    from      = file.path(source_dir, filenames[!is_directory]), 
    to        = out_dir,
    overwrite = overwrite, 
    recursive = FALSE
  )
  if(!all(ok))
  {
    warning(
      "The files ", 
      toString(sQuote(filenames[!ok])), 
      " were not copied successfully."
    )
  }
  names(ok) <- filenames[!is_directory]
  invisible(ok)
}

#' @rdname choose_files
#' @importFrom assertive.reflection assert_is_windows
#' @export
choose_dir <- function(default = "", sep = c("/", "\\"))
{
  assert_is_windows()
  sep <- match.arg(sep)
  caption <- gettext("Select a folder", domain = "R-pathological")
  x <- utils::choose.dir(default, caption)
  standardize_path(x, sep, include_names = FALSE)
}

#' Choose files interactively
#' 
#' Choose one or more files or a directory interactively using a pop-up dialog.
#' @param default The default file to be selected.  See the Details section of
#' \code{choose.files} for how to specify.  Only on Windows.
#' @param multi Logical value indicating if multiple files can be selected.
#' Only on Windows.
#' @param sep String separator between directory levels in the output.
#' @return A character vector of standardized file paths that were chosen.
#' @note \code{choose_files} uses \code{choose.files} under Windows 
#' and \code{\link[base]{file.choose}} under other platforms.
#' @importFrom assertive.reflection is_windows
#' @examples 
#' \donttest{
#' if(interactive())
#' {
#'   choose_files()
#'   if(assertive.reflection::is_windows())
#'   {
#'     choose_dir()
#'   }
#' }
#' }
#' @export
choose_files <- function(default = "", multi = FALSE, sep = c("/", "\\"))
{
  if(!interactive())
  {
    stop("You are not running R interactively; use dir() instead.")
  }
  sep <- match.arg(sep)
  x <- if(is_windows())
  {
    filters <- matrix(
      c(
        "All Files", "*",
        "R source files (R, c, cpp, Rnw, Rmd, Rhtml, Rd)", "*.R;*.c;*.cpp;*.Rnw;*.Rmd;*.Rhtml;*.Rd",
        "Delimited text files (csv, dlm, dat)", "*.csv;*.dlm;*.dat",
        "Text files (txt)", "*.txt",
        "Spreadsheets (xlsx, xls, ods)", "*.xlsx;*.xls;*.ods",
        "Archives (zip, tar, tar.gz, tar.xz)", "*.zip;*.tar;*.tar.gz;*.tar.xz",
        "Image files (png, jpeg/jpg, pdf, ps, emf, bmp, svg)", "*.png;*.jpeg;*.jpg;*.pdf;*.ps;*.emf;*.bmp;*.svg",
        "R workspaces and variables (RData, rda, rds)", "*.RData;*.rda;*.rds",
        "Web files (html, xhtml, css, js)", "*.html;*.xhtml;*.css;*.js"
      ),
      ncol = 2,
      byrow = TRUE
    )
    caption <- if(multi)
    {
      gettext("Select at least one file", domain = "R-pathological")
    } else
    {
      gettext("Select a file", domain = "R-pathological")
    }
    # Have to use :: rather than @importFrom because fn only exists on Windows
    utils::choose.files(default, caption, multi, filters = filters, index = 1)
  } else
  {
    if(nzchar(default))
    {
      warning("Setting the default file is only supported under Windows.")
    }
    if(multi)
    {
      warning("Choosing multiple files is only supported under Windows.")
    }
    tryCatch(
      file.choose(),
      error = function(e) character()
    )
  }
  standardize_path(x, sep, include_names = FALSE)
}

#' Create or remove files and directories
#' 
#' A vectorized version of \code{\link[base]{dir.create}}, and
#' \code{\link[base]{file.create}} and \code{\link[base]{unlink}} with more 
#' convenient defaults.
#' @param x A character vector of paths of directories to create/remove. 
#' For \code{create_dirs}, it defaults to a directory inside \code{tempdir()}.
#' @return A logical vector of successes of failures.
#' @note \code{create_dirs} will only attempt to create directories that don't 
#' already exist.
#' @seealso \code{\link[base]{dir.create}}, \code{\link[base]{unlink}}
#' @examples
#' \donttest{
#' dirs <- temp_dir(c("foo", "bar/baz"))
#' create_dirs(dirs)
#' 
#' # Check this worked:
#' assertive.files::assert_all_are_dirs(dirs)
#' 
#' files <- temp_dir("blah/blah/blah", LETTERS)
#' create_files(files)
#' 
#' assertive.files::assert_all_are_existing_files(files)
#' 
#' # Clean up
#' remove_dirs(temp_dir(c("foo", "bar", "blah")))
#' }
#' @importFrom stats setNames
#' @export
create_dirs <- function(x = temp_file(pattern = "dir"))
{
  doesnt_yet_exist <- !file.exists(x)
  yn <- setNames(logical(length(x)), x)
  yn[doesnt_yet_exist] <- vapply(
    x[doesnt_yet_exist],
    dir.create,
    logical(1), 
    recursive = TRUE
  )
  yn
}

#' @rdname create_dirs
#' @export
create_files <- function(x = temp_file())
{
  dirs <- unique(dirname(standardize_path(x)))
  create_dirs(dirs)
  setNames(file.create(x), x)
}

#' Make a path suitable for cygwin
#' 
#' By default, cygwin complains about standard paths.  This function converts 
#' paths to a form that cygwin likes.
#' @param x A character vector of file paths. Defaults to files in the 
#' current directory.
#' @return A character vector of the cygwinified inputs.
#' @seealso \code{standardize_path}
#' @examples
#' \dontrun{
#' cygwinify_path(c("c:/Program Files", "\\\\some/network/drive"))
#' }
#' @importFrom assertive.reflection is_windows
#' @importFrom stringr fixed
#' @importFrom stringr str_detect
#' @importFrom stringr str_split_fixed
#' @export
cygwinify_path <- function(x = dir())
{
  if(!is_windows())
  {
    warning(
      "This function is expecting to be run under Windows, but the OS is ", 
      .Platform$OS.type,
      ".  Returning x untouched."
    )
    return(invisible(x))
  }
  cygwinified_x <- standardize_path(x)
  colon <- fixed(":")
  has_drive <- str_detect(cygwinified_x, colon)
  split_path <- str_split_fixed(cygwinified_x[has_drive], colon, 2L)
  cygwinified_x[has_drive] <- paste0(
    "/cygdrive/",
    split_path[, 1L],
    split_path[, 2L]
  )
  cygwinified_x
}

#' Split a path into its components
#' 
#' \code{decompose_path} splits a path into the directory name, filename 
#' without extension, and extension. \code{strip_extension}, 
#' \code{get_extension} and \code{replace_extension} provide shortcuts to 
#' manipulate the file extension. \code{recompose_path} takes the result of 
#' \code{decompose_path} and returns complete paths.
#' @param x A character vector of file paths. Defaults to files in the 
#' current directory.
#' @param new_extension A new extension to replace the existing ones.
#' @param include_dir Should the directory part of the path be included? If 
#' \code{NA}, the default, keep the directory from the input.  If \code{TRUE},
#' standardize the directory.  If \code{FALSE}, strip the directory.
#' @param ... Not currently used.
#' @return \code{decompose_path} returns a character matrix with three 
#' columns named \code{"dirname"}, \code{"filename"} and \code{"extension"}.
#' \code{strip_extension} returns a character vector of the filename, possibly 
#' with a directory (see \code{include_dir} argument).
#' \code{replace_extension} returns a character vector of the filename with a  
#' newextension, possibly with a directory (see \code{include_dir} argument).
#' \code{get_extension} returns a character vector of the third column.
#' \code{recompose_path} returns a character vector of paths.
#' @examples
#' x <- c(
#'   "somedir/foo.tgz",         # single extension
#'   "another dir\\bar.tar.gz", # double extension
#'   "baz",                     # no extension
#'   "quux. quuux.tbz2",        # single ext, dots in filename
#'   R.home(),                  # a dir
#'   "~",                       # another dir
#'   "~/quuuux.tar.xz",         # a file in a dir
#'   "",                        # empty 
#'   ".",                       # current dir
#'   "..",                      # parent dir
#'   NA_character_              # missing
#' )
#' (decomposed <- decompose_path(x))
#' get_extension(x)
#' strip_extension(x)
#' strip_extension(x, FALSE)
#' recompose_path(decomposed)
#' @importFrom assertive.properties is_empty
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_not_na
#' @importFrom assertive.files is_dir
#' @importFrom assertive.base strip_attributes
#' @importFrom stringr str_detect
#' @importFrom stringr fixed
#' @importFrom stringr str_match
#' @export
decompose_path <- function(x = dir())
{
  if(is_empty(x))
  {
    return(
      structure(
        data.frame(
          dirname          = character(), 
          filename         = character(), 
          extension        = character(),
          stringsAsFactors = FALSE,
          row.names        = character()
        ),
        class = c("decomposed_path", "data.frame")
      )
    )
  }
  original_x <- x <- coerce_to(x, "character")
  x <- standardize_path(x)
  not_missing <- strip_attributes(is_not_na(x))
  is_dir_x <- is_dir(x)
  
  basename_x <- ifelse(
    not_missing,
    ifelse(is_dir_x, "", basename(x)),
    NA_character_
  )
  has_extension <- str_detect(basename_x, fixed("."))
    
  # match one or more letters, numbers and allowed punctuation characters
  # (the filename without extension)
  # then a single period
  # then match one of more letters numbers and periods
  # (the file extension)
  rx <- "^([]\\[[:alnum:] `!@#$%^&()_=+{},.;'-]+?)\\.([[:alnum:].]+)$"
  
  filename_x <- ifelse(not_missing, basename_x, NA_character_)
  extension_x <- ifelse(not_missing, "", NA_character_)
  not_missing_and_has_extension <- not_missing & has_extension
  
  if(any(not_missing_and_has_extension))
  {
    split_name <- str_match(
      basename_x[not_missing_and_has_extension], 
      rx
    )
  
    filename_x[not_missing_and_has_extension] <- split_name[, 2L]
    extension_x[not_missing_and_has_extension] <- split_name[, 3L]
  }
  
  decomposed_x <- data.frame(
    dirname   = ifelse(
      not_missing,
      ifelse(is_dir_x, x, standardize_path(dirname(x))), #restandardisation required
      NA_character_
    ),
    filename  = filename_x, 
    extension = extension_x,
    row.names = ifelse(is.na(original_x), "<NA>", original_x),
    stringsAsFactors = FALSE
  )
  
  structure(decomposed_x, class = c("decomposed_path", "data.frame"))
}

#' @rdname copy_dir
#' @export
dir_copy <- function(...)
{
  .Deprecated("copy_dir")
  copy_dir(...)
}

#' On Windows, return the drive of the path
#' 
#' On a Windows system, this returns the drive letter of the path followed by a 
#' colon.  On other systems, it returns a single forward slash.
#' @param x A character vector of file paths. Defaults to the current directory.
#' @return A character vector of drive paths on Windows systems, or forward 
#' slashes on Unix-based systems.
#' @seealso \code{\link{is_windows_drive}}
#' @examples
#' get_drive(c("~", r_home(), temp_dir()))
#' @importFrom assertive.reflection is_windows
#' @importFrom utils head
#' @export
get_drive <- function(x = getwd())
{
  if(is_windows())
  {
    vapply(strsplit(standardize_path(x), "/"), head, character(1), n = 1)
  } else
  {
    rep.int("/", length(x))
  }
}

#' @rdname decompose_path
#' @importFrom stats setNames
#' @export
get_extension <- function(x = dir())
{
  setNames(decompose_path(x)$extension, x)  
}

#' Get the libraries on your machine
#' 
#' Wrapper to \code{\link[base]{.libPaths}} that gets all the libraries that R 
#' knows about on your machine.
#' @param index A numberic or logical vector specifying the index of the 
#' libraries to return.  By default, all libraries are returned.
#' @param sep String separator between directory levels in the output.
#' @return A character vector of paths to libraries.
#' @seealso \code{\link[base]{.libPaths}}
#' @references \url{http://cran.r-project.org/doc/FAQ/R-FAQ.html#What-is-the-difference-between-package-and-library_003f}
#' @examples
#' get_libraries()
#' get_libraries(1)
#' @export
get_libraries <- function(index = TRUE, sep = c("/", "\\"))
{
  standardize_path(.libPaths()[index], sep = sep, include_names = FALSE)
}

#' Is the path a Windows drive?
#' 
#' Checks to see if the path is a Windows drive.
#' @param x A character vector of file paths. Defaults to files in the 
#' current directory.
#' @return A logical vector, \code{TRUE} when the path is a Windows drive name.
#' On non-Windows machines, the return value is \code{FALSE} everywhere.
#' @note The check is done by regular expression: values are considered to be 
#' Windows drive name if they consist of a letter followed by a colon, 
#' optionally followed by a slash or backslash.
#' Paths are standardardized before checking, so \code{.} and \code{..} are 
#' resolved to their actual locations rather than always returning \code{FALSE}.
#' @seealso \code{\link{get_drive}}
#' @examples
#' x <- c("c:", "c:/", "c:\\", "C:", "C:/", "C:\\", "c:/c", "cc:", NA)
#' # Warnings about OS suppressed so package checks pass on non-Windows systems.
#' suppressWarnings(is_windows_drive(x))
#' @importFrom assertive.reflection is_windows
#' @importFrom assertive.base coerce_to
#' @importFrom stats setNames
#' @importFrom stringr str_detect
#' @export
is_windows_drive <- function(x)
{
  if(!is_windows())
  {
    warning(
      "This function is expecting to be run under Windows, but the OS is ", 
      .Platform$OS.type,
      "."
    )
    return(rep.int(FALSE, length(x)))
  }
  original_x <- x <- coerce_to(x, "character")
  # Want to resolve paths with . or ..
  starts_with_dots <- str_detect(x, "^\\.{1,2}[/\\\\]?")
  # Can't use standardize_path since we want that fn to use this
  x[starts_with_dots] <- normalizePath(x[starts_with_dots]) 
  yn <- str_detect(x, "^[[:alpha:]]:[/\\\\]?$")
  setNames(yn, original_x)
}

#' The OS path 
#' 
#' The locations in the operating system \code{PATH} environment variable.
#' @param sep String separator between directory levels in the output.
#' @param standardize Should the paths be standardized?
#' @param splitter The character to split the PATH environment variable on.
#' Defaults to a semi-colon on Windows systems and a colon elsewhere.
#' @return A character vector of paths.
#' @seealso \code{\link[base]{Sys.getenv}}
#' @examples
#' os_path()
#' @importFrom assertive.reflection is_windows
#' @importFrom assertive.types assert_is_a_bool
#' @importFrom assertive.types assert_is_a_string
#' @export
os_path <- function(sep = c("/", "\\"), standardize = TRUE, 
  splitter = if(is_windows()) ";" else ":")
{
  assert_is_a_bool(standardize)
  assert_is_a_string(splitter)
  
  path <- Sys.getenv("PATH")
  path <- if(!nzchar(path))
  {
    warning("The 'PATH' environment variable is unset or empty.")
    character()
  } else
  {
    strsplit(path, splitter)[[1]]
  }
  if(standardize)
  {
    standardize_path(path, sep = sep, include_names = FALSE)  
  } else path
}

#' Get the parent dir
#' 
#' Gets the parent directory of the input.
#' @param x A character vector of file paths. 
#' @param sep String separator between directory levels in the output.
#' @return A character vector of parent directories of the input.
#' @note Missing values are returned as missing.  On Windows, the parent of a 
#' drive, e.g., \code{"c:/"} is itself.  Likewise, under Unix, the parent of 
#' \code{"/"} is itself.
#' @examples
#' (x <- c(
#'   sys_which("R"),
#'   r_home(),
#'   r_profile_site(),
#'   "c:/",  # different behaviour under Windows/Unix
#'   "~",
#'   "/",
#'   "foo/bar/nonexistent",
#'   NA
#' ))
#' parent_dir(x)
#' @importFrom stats setNames
#' @export
parent_dir <- function(x, sep = c("/", "\\")) 
{
  sep <- match.arg(sep)
  original_x <- x <- coerce_to(x, "character")
  x <- standardize_path(x)
  not_missing <- is_not_na(x)
  win_drive <- suppressWarnings(is_windows_drive(x))
  pdir <- rep.int(NA_character_, length(x))
  pdir[not_missing] <- ifelse(
    win_drive[not_missing],
    x[not_missing],
    ifelse(
      strip_attributes(is_dir(x[not_missing])),
      dirname(x[not_missing]),
      dirname(dirname(x[not_missing]))
    )
  )
  if(sep == "\\")
  {
    pdir[not_missing] <- str_replace_all(pdir[not_missing], "/", "\\")
  }
  setNames(pdir, original_x)
}

#' @rdname r_profile
#' @export
r_environ <- function(sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  # From ?Startup:
  # "The name of the user file can be specified by the R_ENVIRON_USER 
  # environment variable"
  x <- Sys.getenv("R_ENVIRON_USER", NA)
  if(is.na(x))
  {
    # "if this is unset, the files searched for are '.Renviron' in the current"
    x <- if(file.exists(".Renviron"))
    {
      ".Rprofile"
    } else if(file.exists("~/.Renviron"))
    {
      # "or in the user's home directory (in that order)"
      "~/.Rprofile"
    } else 
    {
      NA_character_
    }
  } 
  x <- standardize_path(x, sep = sep, include_names = FALSE)    
  unname(x)
}

#' @rdname r_profile
#' @export
r_environ_site <- function(sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  # From ?Startup:
  # "The name of the site file is the one pointed to by the environment variable 
  # R_ENVIRON"
  x <- Sys.getenv("R_ENVIRON", NA)
  if(is.na(x))
  {
    # "if this is unset, 'R_HOME/etc/Renviron.site' is used"
    x <- r_home("etc", "Renviron.site", sep = sep)
    x <- if(file.exists(x))
    {
      x
    } else
    {
      NA_character_
    }
  } else
  {
    x <- standardize_path(x, sep = sep, include_names = FALSE)
  }
  unname(x)
}

#' The R home directory
#' 
#' Return a path to a file in the R home directory.  A vectorized, standardized
#' version of \code{R.home}.
#' @param component \code{"home"} for the root of the R installation directory,
#' or the name of a subdirectory.
#' @param ... Further subdirectories passed to \code{file.path}.
#' @param sep String separator between directory levels in the output.
#' @return A character vector of paths inside the R installation dir.
#' @seealso \code{\link[base]{R.home}}
#' @examples
#' r_home()
#' r_home("etc", "Rprofile.site")
#' r_home(c("home", "bin", "share"), c("", "i386", "zoneinfo"))
#' @export
r_home <- function(component = "home", ..., sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  roots <- vapply(component, R.home, character(1))
  standardize_path(file.path(roots, ...), sep = sep, include_names = FALSE)
}

#' Get the location of the R profile/environ
#' 
#' Gets the location of the user or site R profile and environ startup files.
#' @param sep String separator between directory levels in the output.
#' @return A string giving the path the \code{".Rprofile"}, \code{".Renviron"}, 
#' \code{"Rprofile.site"}, or \code{".Renviron.site"}.  If the file cannot be 
#' found, NA is returned.
#' @seealso \code{\link[base]{Startup}} for how this is calculated.
#' @examples
#' r_environ()
#' r_environ_site()
#' r_profile()
#' r_profile_site()
#' @aliases startup environ
#' @export
r_profile <- function(sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  # From ?Startup:
  # "The path of this file can be specified by the R_PROFILE_USER environment 
  # variable"
  x <- Sys.getenv("R_PROFILE_USER", NA)
  if(is.na(x))
  {
    # "If this is unset, a file called '.Rprofile' is searched for in the 
    # current directory"
    x <- if(file.exists(".Rprofile"))
    {
      ".Rprofile"
    } else if(file.exists("~/.Rprofile"))
    {
      # "or in the user's home directory (in that order)"
      "~/.Rprofile"
    } else 
    {
      NA_character_
    }
  } 
  x <- standardize_path(x, sep = sep, include_names = FALSE)    
  unname(x)
}

#' @rdname r_profile
#' @export
r_profile_site <- function(sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  # From ?Startup:
  # "The path of this file is taken from the value of the R_PROFILE environment 
  # variable"
  x <- Sys.getenv("R_PROFILE", NA)
  if(is.na(x))
  {
    # "If this variable is unset, the default is 'R_HOME/etc/Rprofile.site'"
    x <- r_home("etc", "Rprofile.site", sep = sep)
    x <- if(file.exists(x))
    {
      x
    } else
    {
      NA_character_
    }
  } else
  {
    x <- standardize_path(x, sep = sep, include_names = FALSE)
  }
  unname(x)
}

#' @rdname decompose_path
#' @export
recompose_path <- function(x, ...)
{
  UseMethod("recompose_path")
}

#' @rdname decompose_path
#' @method recompose_path decomposed_path
#' @importFrom assertive.base is_not_na
#' @export
recompose_path.decomposed_path <- function(x, ...)
{
  not_missing <- is_not_na(x$filename)
  has_an_extension <- nzchar(as.character(x[not_missing, "extension"]))
  path <- rep.int(NA_character_, nrow(x))
  base_x <- ifelse(
    has_an_extension,
    paste(x[not_missing, "filename"], x[not_missing, "extension"], sep = "."),
    x[not_missing, "filename"]
  )
  path[not_missing] <- file.path(x[not_missing, "dirname"], base_x)
  path
}

#' @rdname create_dirs
#' @export
remove_dirs <- function(x)
{
  unlink(x, recursive = TRUE, force = TRUE)
}

#' @rdname decompose_path
#' @importFrom assertive.files is_dir
#' @importFrom assertive.types assert_is_a_bool
#' @importFrom assertive.types assert_is_character
#' @importFrom assertive.base strip_attributes
#' @importFrom stats setNames
#' @export
replace_extension <- function(x = dir(), new_extension, include_dir = NA)
{
  assert_is_character(new_extension)
  assert_is_a_bool(include_dir)
  if(!nzchar(new_extension))
  {
    warning("'new_extension' is empty.  Did you want strip_extension instead?")
  }
  is_dir_x <- strip_attributes(is_dir(x))
  if(any(is_dir_x))
  {
    warning(
      "The directories ", 
      toString(sQuote(x[is_dir_x])), 
      " have no file extensions to replace."
    )
  }
  stripped <- strip_extension(x, include_dir = include_dir)
  setNames(
    ifelse(
      is_dir_x,
      stripped,
      paste(stripped, new_extension, sep = ".")
    ),
    names(stripped)
  )
}

#' Get the RStudio project directory
#' 
#' Gets the current RStudio project directory.
#' @param sep String separator between directory levels in the output.
#' @return A string giving the path to the current RStudio project directory,
#' or character vector with length zero if you are running RStudio without a
#' project open.
#' @note This only works when your IDE is RStudio. Otherwise an error is thrown.
#' @examples 
#' assertive.base::dont_stop(rstudio_project_dir())
#' @importFrom assertive.reflection assert_is_rstudio
#' @export
rstudio_project_dir <- function(sep = c("/", "\\"))
{
  assert_is_rstudio()
  e <- as.environment("tools:rstudio")
  x <- e$.rs.getProjectDirectory()
  standardize_path(x, sep = sep, include_names = FALSE)
}

#' Split a path into directory components
#' 
#' Splits a character vector of paths into directory components.  The opposite  
#' of \code{\link[base]{file.path}}.
#' @param x A character vector of file paths. Defaults to files in the 
#' current directory.
#' @return A named list of character vectors containing the split paths.
#' @note Paths are split on forward and back slashes, except for double forward 
#' or back slashes at the start of (UNC) paths.  These are included in the first 
#' element of that split path.
#' @examples
#' (splits <- split_path(c(getwd(), "~", r_home())))
#' # Reverse the operation
#' sapply(splits, paste, collapse = "/")
#' @importFrom assertive.properties is_empty
#' @importFrom assertive.base coerce_to
#' @importFrom stats setNames
#' @export
split_path <- function(x = dir())
{
  if(is_empty(x))
  {
    return(setNames(list(), character()))
  }
  original_x <- x <- coerce_to(x, "character")
  x <- standardize_path(x)
  split_x <- strsplit(x, "(?<=[^/\\\\])[/\\\\]", perl = TRUE)
  # setting names in a list, as of R3.1.1 processes backslashes cat-style, so 
  # need to duplicate them
  # original_x <- str_replace_all(original_x, fixed("\\"), "\\\\")
  setNames(split_x, original_x)
}

#' Standardize paths
#' 
#' Standardi[sz]e path names so that they can be more easily compared.
#' @param x A character vector of file paths. Defaults to files in the 
#' current directory.
#' @param sep String separator between directory levels in the output.
#' @param include_names A logical value indicating whether the output should be 
#' named with the input file paths.
#' @return A character vector of paths, pointing to the same locations as the
#' input, but in a standardized form.
#' @seealso \code{\link[base]{normalizePath}}, \code{\link[base]{path.expand}},
#' \code{\link[R.utils]{getAbsolutePath}}
#' @examples
#' standardize_path(c(".", "..", "~", R.home(), NA))
#' standardize_path(c(".", "..", "~", R.home(), NA), "\\")
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.properties is_empty
#' @importFrom assertive.reflection is_unix
#' @importFrom assertive.reflection is_windows
#' @importFrom assertive.strings is_non_missing_nor_empty_character
#' @importFrom stats setNames
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_detect
#' @export
standardize_path <- function(x = dir(), sep = c("/", "\\"), include_names = TRUE)
{
  if(is_empty(x))
  {
    return(setNames(character(), character()))
  }
  sep <- match.arg(sep)
  x <- original_x <- coerce_to(x, "character")
  
  ok <- is_non_missing_nor_empty_character(x)
  
  # standardize = expand + normalize
  # normalizePath is uncomfortable with backslashes under Unix.
  x[ok] <- str_replace_all(x[ok], "[/\\\\]", "/")
  x[ok] <- ifelse(
    is.na(x[ok]),
    NA_character_,
    normalizePath(x[ok], "/", mustWork = FALSE)
  )
  
  # again under Unix, normalizePath won't make path absolute
  if(is_unix())
  {
    x[ok] <- ifelse(
      str_detect(x[ok], "^(/|[[:alpha:]]:)"),
      x[ok], 
      file.path(getwd(), x[ok], fsep = "/")
    )
  }
  
  # Under Windows, normalizePath prefixes UNC paths with backslashes rather than 
  # forward slashes
  if(is_windows())
  {
    x[ok] <- str_replace(x[ok], "^\\\\\\\\", "//")
  }
  
  # strip trailing slashes
  x[ok] <- str_replace(x[ok], "/?$", "")  
  
  # Replace / with the chosen slash
  if(sep == "\\")
  {
    x[ok] <- str_replace_all(x[ok], fixed("/"), "\\")
  }
  if(include_names)
  {
    setNames(x, original_x) 
  } else
  {
    unname(x)
  }
}

#' @rdname standardize_path
#' @export
standardise_path <- standardize_path

#' @rdname decompose_path
#' @importFrom assertive.files is_dir
#' @export
strip_extension <- function(x = dir(), include_dir = NA)
{
  # Empty string x gets returned as empty string
  stripped <- character(length(x))
  # Missing x gets returned as NA
  stripped[is.na(x)] <- NA_character_
  
  #Everything else
  ok <- nzchar(x) & !is.na(x)
  decomposed <- decompose_path(x[ok])
  
  stripped[ok] <- if(is.na(include_dir))
  {
    # For include_dir = NA, keep directory same as input
    dirname_x <- ifelse(is_dir(x[ok]), x[ok], dirname(x[ok]))
    ifelse(
      nzchar(decomposed$filename),
      ifelse(
        dirname_x == ".",
        decomposed$filename,                              # no directory
        file.path(dirname_x, decomposed$filename)         # both
      ),
      dirname_x                                           # no filename
    )    
  } else if(include_dir) 
  {
    # For include_dir = TRUE, add standardized directory
    ifelse(
      nzchar(decomposed$filename),
      file.path(decomposed$dirname, decomposed$filename), # both
      decomposed$dirname                                  # no filename
    )
  } else
  {
    # For include_dir = FALSE, strip directory
    decomposed$filename
  }
  setNames(stripped, x)
}

#' Find paths to executables
#' 
#' Wrapper to \code{Sys.which}, that returns standardized paths.
#' @param x A character vector of executables.
#' @param sep String separator between directory levels in the output.
#' @return A character vector of paths to those executables, or \code{NA} if it 
#' doesn't exist. (This behaviour for missing executables differs from 
#' \code{Sys.which}.)
#' @seealso \code{\link[base]{Sys.which}}
#' @examples
#' sys_which("R")              # R executable
#' sys_which(c("make", "gcc")) # tools for running Rcpp
#' @export
sys_which <- function(x, sep = c("/", "\\"))
{
  std_x <- standardize_path(Sys.which(x), sep = sep, include_names = FALSE)
  ifelse(nzchar(std_x), std_x, NA_character_)
}

#' Find a file in a package
#' 
#' Wrapper to \code{system.file} that returns standardized paths.
#' @param ... Character vectors specifying subdirectories and files within some 
#' package. The default, none, returns the root of the package. Wildcards are 
#' not supported.
#' @param package A string with the name of a single package. An error occurs if 
#' more than one package name is given.
#' @param library_location a character vector with path names of R libraries. 
#' See the 'Details section of \code{\link[base]{system.file}} for the meaning 
#' of the default value of NULL.
#' @param must_work If \code{TRUE}, an error is given if there are no matching 
#' files.
#' @param sep String separator between directory levels in the output.
#' @return A character vector of positive length, containing the file paths that 
#' matched \code{...}, or a missing string, \code{NA}, if none matched (unless 
#' \code{mustWork = TRUE}).  (This behaviour for missing paths differs from 
#' \code{system.file}.)
#' If matching the root of a package, there is no trailing separator.
#' system.file() with no arguments gives the root of the base package.
#' @seealso \code{\link[base]{system.file}}
#' @examples
#' # Examples taken from ?system.file
#' system_file()                  # The root of the 'base' package
#' system_file(package = "stats") # The root of package 'stats'
#' system_file("INDEX")
#' system_file("help", "AnIndex", package = "splines")
#' @export
system_file <- function(..., package = "base", library_location = NULL, 
  must_work = FALSE, sep = c("/", "\\"))
{
  paths <- standardize_path(
    system.file(
      ..., 
      package  = package, 
      lib.loc  = library_location, 
      mustWork = must_work
    ), 
    sep = sep, 
    include_names = FALSE
  )
  ifelse(nzchar(paths), paths, NA_character_)
}

#' Return paths to files or dirs within the temp dir
#' 
#' Vectorized wrappers to \code{tempdir} and \code{tempfile} that return  
#' standardized paths.
#' @param ... Character vectors of further directories within the temp 
#' directory. Passed to \code{\link[base]{file.path}}.
#' @param pattern Character vector of prefixes for the temp file name. 
#' Passed to \code{\link[base]{tempfile}}.
#' @param fileext Character vector of file extensions for the temp file. Passed 
#' to \code{\link[base]{tempfile}}.
#' @param sep String separator between directory levels in the output.
#' @return For \code{temp_file} a character vector giving the names of possible 
#' (temporary) files. Note that no files are generated by \code{temp_file}.
#' For \code{temp_dir}, the path of the per-session temporary directory.
#' @seealso \code{\link[base]{tempdir}}
#' @examples
#' temp_dir(c("foo", "bar/baz"))
#' temp_file(c("foo", "bar/baz"), fileext = c(".txt", ".R"))
#' @export
temp_dir <- function(..., sep = c("/", "\\"))
{
  sep <- match.arg(sep)
  tmp <- file.path(tempdir(), ..., fsep = sep)
  standardize_path(tmp, sep = sep, include_names = FALSE)
}

#' @rdname temp_dir
#' @export
temp_file <- function(..., pattern = "file", fileext = "", sep = c("/", "\\"))
{
  x <- temp_dir(...)
  pattern <- rep_len(pattern, length(x))
  fileext <- rep_len(fileext, length(x))
  tmp <- mapply(
    tempfile, 
    x,
    pattern = pattern, 
    fileext = fileext
  )
  standardize_path(tmp, sep = sep, include_names = FALSE)
}
