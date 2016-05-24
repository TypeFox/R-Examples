`print.ecespa.kci` <-
function(x, ...)
{
   if(x$type=="Kinhom.log") { 
        if (x$probname=="null" & is.null(x$modtrend)){
            cat ("K function and envelopes of", x$datanamea,
                   "\nsimulated by", x$nsim, "random thinnings of the original pattern,
    \n using ecespa function Kinhom.log\n\n")
        }
   
        if( x$probname=="null" & !is.null(x$modtrend)){
            cat ("inhomogeneous K function and envelopes of", x$datanamea, 
   "\nsimulated by", x$nsim, "random thinnings of the original pattern, 
            \nusing function ecespa Kinhom.log with spatial trend as", deparse(x$modtrend),"\n\n")
        }

        if(x$probname!="null" & !is.null(x$modtrend)){
            cat ("inhomogeneous K function and envelopes of", x$datanamea, 
                   "\nsimulated by", x$nsim, "thinnings of the original pattern,\n", 
           "according to the probability vector",x$probname, 
           ",\n using ecespa function Kinhom.log with spatial trend as", deparse(x$modtrend),"\n\n")
        }

        if( x$probname!="null" & is.null(x$modtrend)){
            cat ("K function and envelopes of", x$datanamea, 
                   "\nsimulated by", x$nsim, "thinnings of the original pattern, 
            \naccording to the probability vector",x$probname, 
           ",\n using ecespa function Kinhom.log\n\n")
        }
   }

   if(x$type=="Ki"){
      cat("K function and envelopes of", x$datanamea,
            "\nbuilt by", x$nsim,"simulations from", x$modnamea,
             "model,\n using ecespa function Ki\n")
  }

  if(x$type=="Kci"){
   cat("Cross K functions and envelopes of", x$datanamea, "and", x$datanameb,
   "\nbuilt by", x$nsim,"simulations from", x$modnamea,"and", x$modnameb,
    "models,\n using ecespa function Ki\n")
  }
}

