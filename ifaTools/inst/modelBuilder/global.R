html_list <- function(vars, id) {
  hl <- paste0("<ul id=\'",id,"\' class='stab'>")
  for(i in vars) hl <- paste0(hl, "<li class='ui-state-default stab'><span class='label'>",i,"</span></li>")
  paste0(hl, "</ul>")
}

returnOrder <- function(inputId, vars) {
  tagList(
    singleton(tags$head(tags$script(src = 'js/sort.js'))),
    singleton(tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'sort.css'))),
    HTML(html_list(vars, inputId)),
    tags$script(paste0("$(function() {$( '#",inputId,"' ).sortable({placeholder: 'ui-state-highlight'}); $( '#",inputId,"' ).disableSelection(); });"))
  )
}
