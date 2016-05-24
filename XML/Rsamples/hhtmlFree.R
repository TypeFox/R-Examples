library(XML)
invisible(xmlParse(I("<html/>"))) 
.Call("R_isWeakRef", NULL)
gc(); gc()
.Call("R_isWeakRef", NULL)


#
invisible(xmlParse(I("<html/>"), isHTML = TRUE)) 
.Call("R_isWeakRef", NULL)
gc(); gc()
.Call("R_isWeakRef", NULL)

#
invisible(htmlParse(I("<html/>")))
.Call("R_isWeakRef", NULL)
gc(); gc()
.Call("R_isWeakRef", NULL)  # Still in the weak reference list.


#
hParse = function(...) xmlParse(..., isHTML = TRUE)
invisible(hParse(I("<html/>")))
.Call("R_isWeakRef", NULL)
gc(); gc()
.Call("R_isWeakRef", NULL)  # Still in the weak reference list.
