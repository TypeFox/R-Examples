pos.char <-
function(string,char)
{  zz <- stringf(as.character(string))
   if (char%in%zz)
     { z <- which (zz==char)}  else
       { z <- NA }
    # returns all occurrences
   return (z[1])
}
