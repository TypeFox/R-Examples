#do not edit, edit noweb/qmrparser.nw
option   <- function(ap, 
                     action = function(s)   list(type="option",value=s  ),
                     error  = function(p,h) list(type="option",pos=p, h=h)) 
  
  function(stream) {
    cstream <- ap(stream)
    
    if (cstream$status=="ok")                            
      return(list(status="ok",node=action(cstream$node)                ,stream=cstream$stream))
    
    else return(list(status="ok",node=action(list(type="empty" ,value="")),stream=stream))
    
  }

