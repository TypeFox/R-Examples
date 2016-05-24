

#' @rdname params
#' @name params
#' @title Setting/loading and extracting various options into the environment
#'
#' @aliases get_opts set_opts print.opts params
#'
#' @description
#' \itemize{
#' \item set_opts(): set options into a custom envir
#' \item get_opts(): extract options
#' \item load_opts(): Read a tab delimted file using \link{read_sheet} and load them as options using \link{set_opts}
#' \item new_opts(): create a options manager to be included in a pacakge
#' \item print.opts(): print pkg options as a pretty table
#'}
#'
#' @param x \itemize{
#' \item get_opts: a character vector of names of options to extract.
#' \item load_opts: path to a configuration file
#' }
#'
#' @param ... set_opts(): a named set of variable/value pairs seperated by comma
#' @param .dots set_opts(): A named list, as a alternative to ...
#' @param envir environ used to store objects. Default is a environ object called opts [params:::opts]
#' @param check load_opts(): in case of a configuration file, whether to check if files defined in parameters exists. [TRUE]
#' @param .parse set_opts(), load_opts(): logical, whether to auto-complete \code{{{myvar}}} using previously defined options. [TRUE]
#' @param verbose load_opts(): Logical variable indicate level of verboseness [TRUE]
#' @param .use.names get_opts(): The resulting vector should be have names (esp. if length(x) is 1). If length(x)>1, this returns a list.
#' @details
#'
#' \strong{Integrating \link{params} in a package:}
#'
#' \emph{create a options manager}:
#'
#' \code{
#' opts_mypkg = new_opts()
#' }
#'
#' The object \code{opts_mypkg} is a list of a few functions, which set, fetch and load
#' options (using a isolated environment). Here are a few examples:
#'
#' \emph{Set some options}:
#'
#' \code{opts_mypkg$set(version = '0.1', name = 'mypkg')}
#'
#' \emph{Fetch ALL options}:
#'
#' \code{opts_mypkg$get()}
#' OR
#' \code{opts_mypkg$get("version")} to fetch a specific option.
#'
#'
#' \strong{Loading configuration files}:
#'
#' \code{load_opts()} OR \code{opts_pkg$load()}:
#'
#'
#' There are cases when options and params are actually paths to scripts or other apps or folders etc.
#' In such cases it might be useful to quickly check if these paths exists on the sytem.
#' As such, \link{load_opts}() automatically checks params ending with \code{path|dir|exe} (if check=TRUE).
#'
#' For example, values for variables like \code{mypath}, \code{my_path}, \code{tool_exe}, etc would be check if they exists
#' and a warning would be shown if they dont exist.
#'
#'
#' Below is a list example options, retrieved via
#'
#' \code{get_opts()}:
#'
#' \preformatted{
#'	|name          |value            |
#	|:-------------|:----------------|
#'	|default_regex |(.*)             |
#'	|my_conf_path  |~/flowr/conf     |
#'	|my_dir        |path/to/a/folder |
#'	|my_path       |~/flowr          |
#'	|my_tool_exe   |/usr/bin/ls      |
#'}
#' @seealso \link{read_sheet} \link{load_opts}
#'
#' @export
#'
#' @examples
#' ## Set options
#' opts = set_opts(flow_run_path = "~/mypath")
#' #OR
#' opts = set_opts(.dots = list(flow_run_path = "~/mypath"))
#'
#' ## printing options, this is internally called by get_opts()
#' print(opts)
#'
#' ## Fetch options
#' get_opts()
#' get_opts("flow_run_path")
#'
#' ## Load options from a file
#' fl = system.file("conf/params.conf", package = "params")
#' load_opts(fl)
#'
#'
#' ## Create a options manager:
#' opts_mypkg = new_opts()
#' ## this provides three functions
#' opts_mypkg$set(version = '0.1', name = 'mypkg')
#' opts_mypkg$load(fl)
#' opts_mypkg$get()
#'
#' ## Additionally, one has the options of using braces ({{}})
#' ## do define nested options:
#'
#' set_opts(first = "John", last = "Doe", full = "{{first}} {{last}}")
#'
new_opts <- function(envir = new.env()){

	get_opts <- function(x){
		params::get_opts(x, envir = envir)
	}

	set_opts <- function(..., .dots){
		params::set_opts(..., .dots = .dots, envir = envir)
	}

	load_opts <- function(...){
		params::load_opts(..., envir = envir)
	}

	list(get=get_opts, set=set_opts, load = load_opts)

}
