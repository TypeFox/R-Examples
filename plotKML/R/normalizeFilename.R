# Purpose        : User-friendly characters for "filename"
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Pierre Roudier (pierre.roudier@landcare.nz);
# Status         : Pre-alpha
# Note           : see also: R.utils and utils packages;


## Fix characters following the naming conventions [http://msdn.microsoft.com/en-us/library/windows/desktop/aa365247(v=vs.85).aspx]:
normalizeFilename <- function(x, form = c("default", "8.3")[1], fix.encoding = TRUE, sub.sign = "_"){

   # reserved characters:
   sel = c("<", ">", ":", '\\"', "/", "\\|", "\\?", "\\*", "[[:space:]]", "\\s+$", "\\[", "\\]")
   for(i in sel){ 
      x <- gsub(pattern=i, replacement=sub.sign, x)
   }
   # shorten the path:
   if(form == "8.3"){
      if(.Platform$OS.type == "windows") {
       x <- utils::shortPathName(x)
      }
   }
   if(fix.encoding==TRUE){
      x <- iconv(x, to = "UTF8")
   }

   return(x)
}

# end of script;
