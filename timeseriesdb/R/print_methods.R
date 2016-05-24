#' @export
print.meta_env <- function(x,...){
  env_name <- deparse(substitute(x))
  obj_count <- length(ls(envir = x))
  out <- sprintf('%s contains %s meta information object(s).',
                 env_name,obj_count)
  cat(out)
  
}

