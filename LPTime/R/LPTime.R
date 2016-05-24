LPTime <-
function(z, exo=NULL, m = 3, p = 10)
    {
        zt <- matrix(apply(z, 2, LPTrans,m), nrow = nrow(z))
         if(!is.null(exo))
              Xt <- matrix(apply(as.matrix(exo),2, LPTrans,m), nrow = nrow(as.matrix(exo)))
         else
              Xt = NULL
         VAR(zt, exogen = Xt, p)
         }
