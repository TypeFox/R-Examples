#  -----------------------------------------------------------
#  pyImport
#  ========
#' @title Import virtual Python objects to R
#'
#' @description A convenience function to call the Python function 
#'              \bold{import} and creating virtual Python objects for
#'              the imported objects in R.
#' @param import a character giving the names of the objects to import.
#' @param from an optional character string giving the name of the module.
#' @param as an optional string defining an alias for the module name.
#' @param env an optional environment where the virtual Python objects are 
#'            assigned to.
#' @details The function pyImport works like the import function in Python.
#'            
#'          The function pyImport has a special behavior for the packages numpy 
#'          and pandas. For these two packages pyImport does not only import 
#'          numpy but also register their alias in pyOptions. To be found when
#'          pySet is used with the option useNumpy set to TRUE.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' pyImport("os")
#' \dontrun{
#' #NOTE: The following does not only import numpy but also register the
#' #      alias in the options under the name "numpyAlias". 
#' #      The same is done for pandas, the default alias for pandas and numpy 
#' #      are respectively "pandas" and "numpy". The numpyAlias is used 
#' #      when calling pySet with the pyOption useNumpy set to TRUE.
#' pyOptions("numpyAlias")
#' pyImport("numpy", as="np")
#' pyOptions("numpyAlias")
#' pyImport("pandas", as="pd")
#' pyImport(c("getcwd", "sep"), from="os")
#' getcwd()
#' sep
#' sep = "Hello R!"
#' pyExecp("sep")
#' }
# -----------------------------------------------------------
## NOTE: fun <- function(env=parent.env(environment())) would 
##       have a different, in this case not desired behavior!
pyImport <- function(import, from=NULL, as=NULL, env = parent.frame()){
  if ( pyConnectionCheck() ) return(invisible(NULL))
  if (!is.character(import)) stop('"import" has to be of type character')
  if (!is.null(from)) check_string(from)
  if (!is.null(as)) if (!is.character(as)) stop('"as" has to be of type character')
  # import 0 | import + as 1 | import + from 2 | import + as + from 3
  mode <- sum((1:2)[c(!is.null(as), !is.null(from))])

  if (mode >= 2){
    if ( import[1] == "*" ){
      fromLoaded <- from %in% pyDir()
      if ( !fromLoaded ) pyExec(sprintf("import %s", from))
      import <- pyDir(from)
      if ( !fromLoaded ) pyExec(sprintf("del(%s)", from)) ## cleanup
    }
  }

  if (mode == 0){
    pyExec(sprintf("import %s", import))
    pyAttach(import, env=env)
  }else if(mode == 1){
    pyExec(sprintf("import %s as %s", import, as))
    if (import=="numpy") pyOptions("numpyAlias", as)
    if (import=="pandas") pyOptions("pandasAlias", as)
    pyAttach(as, env=env)
  }else if(mode == 2){
    if ( length(import) > 1 ) {
      imp <- paste(import, collapse=", ")
      cmd <- sprintf("from %s import %s", from, imp)
      pyExec(cmd)
      for (imp in import){
        pyAttach(imp, env=env)
      }
    }else{
      pyExec(sprintf("from %s import %s", from, import))
      pyAttach(import, env=env)
    }
  }else if(mode == 3){
    if ( length(as) > 1 ){
      if ( length(as) !=  length(import) ) stop('if length(as) > 1, length(as) has to be equal to length(import)')
        for ( i in seq_along(import) ) {
          pyExec(sprintf("from %s import %s as %s", from, import[i], as[i]))
        }
    }else{
      pyExec(sprintf("import %s as %s", from, as))
      as <- paste(as, import, sep=".")
    }

    for ( i in seq_along(import) ) {
      pyAttach(as[i], env=env)
    }         
  }else{stop(pi)}
  invisible(NULL)
}
