#' just supporting function
#'
#' @importFrom matlab flipud
#' @importFrom signal conv

acfpacf.acf <-
function(x,normflg){

  y=matrix(x,length(x),1)
   nr=nrow(y)
   nc=ncol(y)
  if (nc > nr)
   {stop(" 'x' must be a column vector")}

    f=flipud(x)
    r=conv(f,x)
    r=r[nr:length(r)]

   if (normflg==0)
    { r=r/nr
      } else {
        if (normflg==1)
       { r=r/nr
         r=r/r[1]
        }else {
          if(normflg==2)
            { den=t(seq(nr,1,-1))
               r=r/den
               r=t(r)
             } else {
               den=t(seq(nr,1,-1))
               r=r/den
               r=r/r[1]
               r=t(r)
            }
        }
      }
    acf=r
    result = list(acf=acf)
    class(result) = "acfpacf.acf"
    result
}

