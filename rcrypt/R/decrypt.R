#' Decrypt a File Using GPG
#'
#' Decrypt a symmetrically encrypted file using GPG.
#'
#' @param input A character string of the file name you wish to decrypt.
#' @param output A character string of the file name that will be created. The default is to create a file with the same name (stripped of the .gpg or .asc file extension) in the same folder.
#' @param passphrase A character string of the passphrase used to decrypt the encrypted file. WARNING: use this to bypass the more secure option of GPG's passphrase popup box. WARNING: the passphrase may be saved in the script as cleartext, saved in the terminal history in cleartext, and/or available in the list of processes in cleartext. The default value is \code{NULL} (Insert the passphrase using GPG's secure pop-up box).
#' @param verbosity An integer \code{0}, \code{1}, \code{2}, or \code{3}. Control GPG's terminal message information in increasing level of detail. A value of \code{0} passes the \code{'--quiet'} flag for a minimum amount of information. A value of \code{1} does not pass any flags. A value of \code{2} passes the \code{'--verbose'} flag. A value of \code{3} passes the \code{'--verbose --verbose'} flag for the most information. The default value is \code{1}.
#' @return A decrypted file.
#' @examples
#' \dontrun{
#' decrypt("path/to/your/file.csv.gpg")
#' decrypt("path/to/your/file.csv.gpg", output = "path/to/your/file.csv")
#' # WARNING: only use the passphrase argument if you understand why it's
#' # not secure.
#' decrypt("path/to/your/file.csv.gpg", passphrase = "your-passphrase")
#' }
#' @export
decrypt <- function(input, output = NULL, passphrase = NULL, verbosity = 1) {
  #-----------------------------------------------------------------------------
  # Check the arguments
  #-----------------------------------------------------------------------------
  if (missing(input)) {
    stop("Check the input argument, it seems to be missing. There's nothing to decrypt.")
  }
  if (!file.exists(input)) {
    stop("Check the input argument, the file name doesn't exist. There's nothing to decrypt.")
  }

  # If output is NULL use input name.
  if (is.null(output)) {
    output <- gsub(".gpg|.asc", "", input)
  }

  # If you try to encrypt to an existing file name.
  if (file.exists(output)) {
    stop("Check the output argument, the file name is already in use! The decrypted file may already exist, or you need to specify a new output file name.")
  }

  # Unix type OS need '--no-tty' for terminal passphrase insertion.
  if (.Platform$OS.type == "unix") {
    tty <- "--no-tty"
  } else{
    tty <- NULL
  }

  # Verbosity
  if (!(verbosity %in% c(0, 1, 2, 3)) ||
      length(verbosity %in% c(0, 1, 2, 3)) == 0) {
    stop("Check the verbosity argument. You've used an invalid value.")
  }
  verbosity <- switch(
    as.character(verbosity),
    "0" = "--quiet",
    "1" = NULL,
    "2" = "--verbose",
    "3" = "--verbose --verbose"
  )

  #-----------------------------------------------------------------------------
  # Decrypt file
  #-----------------------------------------------------------------------------
  if (is.null(passphrase)) {
    # Decrypt with GUI passphrase.
    command <- "gpg"
    system2.args <- c("--output", output, "--decrypt", verbosity, input)
  } else{
    # Decrypt with terminal passphrase insertion.
    # Can't have a space after passphrase when using shell()
    command <- "echo"
    system2.args <- c(
      paste(passphrase, "|", sep = ""),
      "gpg",
      "--passphrase-fd 0",
      "--batch",
      tty,
      "--output",
      output,
      "--decrypt",
      verbosity,
      input
    )
  }

  # Windows doesn't like system2(), need to use shell()
  if (.Platform$OS.type == "unix") {
    what <- "system2"
    args <- list(command = command, args = system2.args)
  } else { # windows
    what <- "shell"
    args <- list(cmd = paste(command,
                             paste(system2.args, collapse = " "),
                             collapse = " "
                             )
                 )
  }

  # Run decryption in terminal
  do.call(what = what, args = args)
}
