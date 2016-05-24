
check.os <- function(options){
  
  os <- 2
  tmp <- .C('check_os', os = as.integer(os), PACKAGE = "ARTP2")
  if(!(tmp$os %in% 0:1)){
    msg <- 'Failed to indentify the OS. Please contact the authors. '
    stop(msg)
  }
  options$os <- ifelse(tmp$os == 0, 'windows_or_mac', 'linux')
  msg <- paste0('ARTP2 ', ifelse(options$os == 'windows_or_mac', 'does not support ', 'supports '), 'multi-threading on this OS')
  if(options$print) message(msg)
  
  options
  
}
