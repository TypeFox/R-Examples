`xport.namestr.header` <-  function( nvar )
  {
    .C("fill_namestr_header",
       nvar = sprintf("%04.4d", as.integer(nvar)),  # Number of variables *AS A STRING*
       PACKAGE="SASxport"
       )
    
    .Call("getRawBuffer", PACKAGE="SASxport")
  }

