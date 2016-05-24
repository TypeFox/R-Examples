#do not edit, edit noweb/qmrparser.nw
repetition0N    <- function(rpa0, 
                            action = function(s)            list(type="repetition0N",value=s   ), 
                            error  = function(p,h) list(type="repetition0N",pos=p,h=h)) 
 
  option(repetition1N(rpa0),action=action,error=error)

