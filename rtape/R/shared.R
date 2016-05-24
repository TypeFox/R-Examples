# Shared, private code

# Decides whether error object was stored or is an end-of-tape marker
.ckErr<-function(v){
 if("try-error"%in%class(v)){
  v<-attr(v,"condition")
  if(identical(v$call,quote(unserialize(con))))
   return(TRUE)
 }
 return(FALSE)
}
