#do not edit, edit noweb/qmrparser.nw
empty   <- function(action = function(s) list(type="empty",value=s),
                    error  = function(p) list(type="empty",pos  =p))
  
  function (stream) return(list(status="ok" ,node=action(""), stream=stream))
