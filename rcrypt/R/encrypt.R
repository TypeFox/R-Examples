#' Encrypt a File Using GPG.
#'
#' Symmetric file encryption using GPG. The \code{encrypt} function defaults to
#' the strongest cryptographic flags available for GPG.
#'
#' @param input A character string of the file name you wish to encrypt.
#' @param output A character string of the file name that will be created. The default is to create a file with the same name (with an additional \code{.gpg} or \code{.asc} file extension) in the same folder.
#' @param passphrase A character string of the passphrase used to decrypt the encrypted file. WARNING: use this to bypass the more secure option of GPG's passphrase popup box. WARNING: the passphrase may be saved in the script as cleartext, saved in the terminal history in cleartext, and/or available in the list of processes in cleartext. The default value is \code{NULL} (Insert the passphrase using GPG's secure pop-up box).
#' @param compress A character string of the methods of compression. Possible values are \code{"Uncompressed"}, \code{"ZIP"}, \code{"ZLIB"}, and \code{"BZIP2"}. Values depend on your GPG installation. The default value is \code{"ZLIB"}.
#' @param cipher A character string of the encryption algorithm. Possible values are \code{"AES256"}, \code{"Camellia256"}, \code{"TWOFISH"}, \code{"AES128"}, etc. Values depend on your GPG installation. The default value is \code{"AES256"}.
#' @param armor \code{TRUE} or \code{FALSE}: flag to produce an encrypted ASCII text output file. The default value is \code{FALSE}.
#' @param mdc \code{TRUE} or \code{FALSE}: flag to force the use of modification detection code. It is always used with newer encryption algorithms and recommended to always keep \code{TRUE}. The default value is \code{TRUE}.
#' @param s2k.mode An integer \code{0}, \code{1}, or \code{3}. Sets how passphrases are mangled. A value of \code{0} just uses a plain passphrase (never use). A value of \code{1} will add a salt to the passphrase. A value of \code{3} will salt and iterate the passphrase. It is highly recommended to always use \code{3}. The default value is \code{3}.
#' @param s2k.digest A character string of the digest algorithm used to mangle passphrases. Possible values are \code{"SHA512"}, \code{"SHA384"}, \code{"SHA256"}, etc. The default value is \code{"SHA512"}.
#' @param s2k.count An integer between \code{1024} and \code{65011712}. Specifies how many times the passphrase mangling is repeated. The default value is \code{65011712}.
#' @param verbosity An integer \code{0}, \code{1}, \code{2}, or \code{3}. Control GPG's terminal message information in increasing level of detail. A value of \code{0} passes the \code{'--quiet'} flag for a minimum amount of information. A value of \code{1} does not pass any flags. A value of \code{2} passes the \code{'--verbose'} flag. A value of \code{3} passes the \code{'--verbose --verbose'} flag for the most information. The default value is \code{1}.
#' @return An encrypted file.
#' @examples
#' \dontrun{
#' encrypt("path/to/your/file.csv")
#' encrypt("path/to/your/file.csv", output = "path/to/your/file.csv.gpg")
#' # WARNING: only use the passphrase argument if you understand why it's
#' # not secure.
#' encrypt("path/to/your/file.csv", passphrase = "your-passphrase")
#' }
#' @export
encrypt <- function(input, output = NULL, passphrase = NULL, compress = "ZLIB",
                    cipher = "AES256", armor = FALSE, mdc = TRUE, s2k.mode = 3,
                    s2k.digest = "SHA512", s2k.count = 65011712,
                    verbosity = 1) {
  #-----------------------------------------------------------------------------
  # Check the arguments
  #-----------------------------------------------------------------------------
  if (missing(input)) {
    stop("Check the input argument, it seems to be missing. There's nothing to encrypt.")
  }

  # If output is NULL use input name.
  if (is.null(output)) {
    if (armor) {
      output <- paste(input, ".asc", sep = "")
    } else {
      output <- paste(input, ".gpg", sep = "")
    }
  }

  # If you try to encrypt to an existing file name.
  if (file.exists(output)) {
    stop("Check the output argument, the file name is already in use! The encrypted file may already exist, or you need to specify a new output file name.")
  }

  # Unix OSs need '--no-tty' for terminal passphrase insertion.
  if (.Platform$OS.type == "unix") {
    tty <- "--no-tty"
  } else{
    tty <- NULL
  }

  # If you need an ASCII file, then include the '--armor' flag.
  if (armor == FALSE) {
    armor <- NULL
  } else{
    armor <- "--armor"
  }

  # Force modification detection code.
  if (mdc == TRUE) {
    mdc <- "--force-mdc"
  } else{
    mdc <- NULL
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
  # Encrypt file
  #-----------------------------------------------------------------------------
  if (is.null(passphrase)) {
    # Encrypt with GUI passphrase.
    command <- "gpg"
    system2.args <- c(
      "--output",
      output,
      "--symmetric",
      armor,
      mdc,
      paste("--compress-algo", compress),
      paste("--cipher-algo", cipher),
      paste("--s2k-mode", s2k.mode),
      paste("--s2k-digest-algo", s2k.digest),
      paste("--s2k-count", s2k.count), verbosity, input
    )
  } else{
    # Encrypt with terminal passphrase insertion.
    # Can't have a space after passphrase when using shell()
    command <- "echo"
    system2.args <- c(
      paste(passphrase, "|", sep = ""),
      "gpg",
      "--passphrase-fd 0",
      "--batch",
      tty, "--output", output, "--symmetric", armor, mdc,
      paste("--compress-algo", compress),
      paste("--cipher-algo", cipher),
      paste("--s2k-mode", s2k.mode),
      paste("--s2k-digest-algo", s2k.digest),
      paste("--s2k-count", s2k.count),
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

  # Run encryption in terminal
  do.call(what = what, args = args)
}
