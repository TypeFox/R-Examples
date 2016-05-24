
####  matR is built on the MG-RAST API...
####  ...so some familiarity with the API is helpful.

####  this demo illustrates how to learn about it.

####  first, an overview of API documentation.

doc.MGRAST()

####  next, documentation for a specific resource.

doc.MGRAST (head="matrix")

####  then, documentation for a request within that resource.

doc.MGRAST (depth=2, head=c("matrix", "organism"))

####  finally, a list of available options for that request.

doc.MGRAST (depth=3, head=c("matrix", "organism", "parameters", "options"))

####  help is available in the usual way, too.

?call.MGRAST
