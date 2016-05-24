helloJavaWorld <- function(){
   hjw <- .jnew("HelloJavaWorld") # create instance of HelloJavaWorld class
   out <- .jcall(hjw, "S", "sayHello") # invoke sayHello method
   return(out)
}


