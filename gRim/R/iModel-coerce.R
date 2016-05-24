## FIXME: We also want to coerce to igraphs; otherwise people get confused...

setOldClass(c("dModel","cModel","mModel"))

setAs("dModel", "graphNEL",
      function(from){
        ugList(from$glist, result="NEL")
      })

setAs("cModel", "graphNEL",
      function(from){
        ugList(from$glist, result="NEL")
      })

setAs("mModel", "graphNEL",
      function(from){
        ugList(from$glist, result="NEL")
      })

setAs("dModel", "matrix",
      function(from){
        ugList(from$glist, result="matrix")
      })

setAs("cModel", "matrix",
      function(from){
        ugList(from$glist, result="matrix")
      })

setAs("mModel", "matrix",
      function(from){
        ugList(from$glist, result="matrix")
      })

setAs("dModel", "igraph",
      function(from){
        ugList(from$glist, result="igraph")
      })

setAs("cModel", "igraph",
      function(from){
        ugList(from$glist, result="igraph")
      })

setAs("mModel", "igraph",
      function(from){
        ugList(from$glist, result="igraph")
      })

