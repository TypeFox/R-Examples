# documentation{{{
#' sets and fixes knitr hooks
#'
#' fixes a well-known bug in the knit_hook 'chunk' and provides a hook entitle
#' 'com'
#'
#' There is a well_known bug in the knit_hook 'chunk' which prevents using
#' results = 'asis' in conjunction with user-defined hooks (including com, to be
#' discussed next). Calling this function allows user-defined hooks to be called
#' with results = 'asis' and get the expected result. This function also
#' provides a knitr hook called 'com', by setting 'com = TRUE' in a knitr chunk,
#' the resulting code is converted to a latex command. For example:
#' '<<mytable, com=TRUE>>=' results in a latex command entitled '\\mytable', which will
#' produce the exact output that would have appeared in the spot of the chunk
#'
#' @examples
#' knitr_sethooks()
#}}}
knitr_sethooks <- function(){

  hook_chunk = knitr::knit_hooks$get('chunk')

  knitr::knit_hooks$set(chunk = function(x, options){
                        x = hook_chunk(x, options)
                        if (options$results == 'asis'){
                          # remove all kframe's
                          gsub('\\\\(begin|end)\\{kframe\\}', '', x)
                        } else x
  })

  knitr::knit_hooks$set(com = function(before, options, envir){
                        if (before){
                          sprintf("\\newcommand\\%s{%%\n", options$label)
                        } else "}"
  })
}
# vim:set foldmethod=marker:
