customStop <-
function(msg, func){
    custom.msg <- paste('In ', func, ', ', msg, sep = '')
    stop(custom.msg, call. = FALSE)
}
