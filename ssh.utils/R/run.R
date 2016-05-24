#-------------------------------------------------------------------------------
#
# Package ssh.utils 
#
# Implementation. 
# 
# Sergei Izrailev, 2011-2014
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------

#' \code{run.withwarn} - Evaluates the expression (e.g. a function call) and returns the
#' result with additional atributes:
#' \itemize{
#' \item num.warnings - number of warnings occured during the evaluation
#' \item last.message - the last warning message
#' }
#' Otherwise, \code{run.withwarn} is similar to \code{base::supressWarnings}
#' 
#' @param expr Expression to be evaluated.
#' @name run.remote
#' @title Functions to run commands remotely via \code{ssh} and capture output.
#' @rdname run.remote
run.withwarn <- function(expr) 
{
   .number_of_warnings <- 0L
   .last_warning_msg <- NULL
   frame_number <- sys.nframe()
   ans <- withCallingHandlers(expr, warning = function(w) 
   {
      assign(".number_of_warnings", .number_of_warnings + 1L, 
            envir = sys.frame(frame_number))
      assign(".last_warning_msg", w$message, 
            envir = sys.frame(frame_number))
      #print(computeRestarts())
      invokeRestart("muffleWarning")
   })
   #message(paste("No. of warnings thrown:", .number_of_warnings, .last_warning_msg))
   attr(ans, "num.warnings") <- .number_of_warnings
   attr(ans, "last.message") <- .last_warning_msg
   ans
}

#------------------------------------------------------------------------------

#' ignore.
#' 
#' \code{run.remote} - Runs the command locally or remotely using ssh. 
#' 
#' In \code{run.remote} the remote commands are enclosed in wrappers that allow to capture output. 
#' By default stderr is redirected to stdout.
#' If there's a genuine error, e.g., the remote command does not exist, the output is not captured. In this case, one can
#' see the output by setting \code{intern} to \code{FALSE}. However, when the command is run but exits with non-zero code, 
#' \code{run.remote} intercepts the generated warning and saves the output. 
#' 
#' The remote command will be put inside double quotes twice, so all quotes in cmd must be escaped twice: \code{\\\"}. 
#' However, if the command is not remote, i.e., \code{remote} is \code{NULL} or empty string, quotes should be escaped 
#' only once.
#' 
#' If the command itself redirects output, the \code{stderr.redirect} flag should be set to \code{FALSE}.
#' @param cmd Command to run. If run locally, quotes should be escaped once. 
#'        If run remotely, quotes should be escaped twice.
#' @param remote Remote machine specification for ssh, in format such as \code{user@@server} that does not 
#'        require interactive password entry. For local execution, pass an empty string "" (default).
#' @param intern Useful for debugging purposes: if there's an error in the command, the output of the remote 
#'        command is lost. Re-running with \code{intern=FALSE} causes the output to be printed to the console.
#'        Normally, we want to capture output and return it.
#' @param stderr.redirect When TRUE appends \code{2>&1} to the command. 
#'        Generally, one should use that to capture STDERR output 
#'        with \code{intern=TRUE}, but this should be set to \code{FALSE} 
#'        if the command manages redirection on its own.
#' @param verbose When \code{TRUE} prints the command.
#' @return \code{run.remote} returns 
#' a list containing the results of the command execution, error codes and messages.
#' \itemize{
#' \item \code{cmd.error} - flag indicating if a warning was issued because command exited with non-zero code
#' \item \code{cmd.out} - the result of the command execution. If there was no error, this contains the output
#'                        as a character array, one value per line, see \code{\link{system}}. If there was 
#'                        an error (as indicated by \code{cmd.error}), this most likely contains the error message
#'                        from the command itself. The \code{elapsed.time} attribute contains the elapsed
#'                        time for the command in seconds.
#' \item \code{warn.msg} - the warning message when \code{cmd.error} is TRUE.
#' }
#' Warnings are really errors here so the error flag is set if there are warnings.
#' 
#' Additionally, \code{cmd.out} has the \code{elapsed.time}, \code{num.warnings} and, if 
#' the number of warnings is greater than zero, \code{last.warning} attributes.
#' @examples
#' \dontrun{
#' ## Error handling:
#' remote = ""
#' command = "ls /abcde"
#' res <- run.remote(cmd=command, remote=remote)
#' if (res$cmd.error)
#' {
#'    stop(paste(paste(res$cmd.out, collapse="\n"), res$warn.msg, sep="\n"))
#' }
#' # Error: ls: /abcde: No such file or directory
#' # running command 'ls /abcde  2>&1 ' had status 1
#' 
#' ## Fetching result of a command on a remote server
#' 
#' # Get the file size in bytes
#' res <- run.remote("ls -la myfile.csv | awk '{print \\$5;}'", remote = "me@@myserver")
#' res
#' # $cmd.error
#' # [1] FALSE
#' # 
#' # $cmd.out
#' # [1] "42"
#' # attr(,"num.warnings")
#' # [1] 0
#' # attr(,"elapsed.time")
#' # elapsed 
#' # 1.063 
#' #
#' # $warn.msg
#' # NULL
#' 
#' file.length <- as.integer(res$cmd.out)
#' }
#' @rdname run.remote
run.remote <- function(cmd, remote = "", intern = T, stderr.redirect = T, verbose = F)
{
   if (stderr.redirect) redir <- " 2>&1 " else redir <- ""
   if (!is.null(remote) && nchar(remote) > 0)
   {
      command <- paste("ssh ", remote, " \"if [ -f  ~/.bash_profile ]; then source ~/.bash_profile; fi; ", cmd, redir, "\"")
   } else
   {
      command <- paste(cmd, redir)      
   }
   if (verbose) print(command)    
   tm <- round(system.time(cmd.out <- run.withwarn(system(command, intern=intern)))[3], 3)
   attr(cmd.out, "elapsed.time") <- tm
   if (attr(cmd.out, "num.warnings") > 0)
   {
      return(list(cmd.error = TRUE, cmd.out = cmd.out, warn.msg = attr(cmd.out, "last.message")))
   }   
   invisible(list(cmd.error = FALSE, cmd.out = cmd.out, warn.msg = NULL))      
}

#-------------------------------------------------------------------------------
