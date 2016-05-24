######### R-function dfltLabs  #########
 
# Obtain default axis labels.

# Last changed: 03 AUG 2005

dfltLabs <- function(d,names.x,xlab,ylab,zlab)
{
   if (d==1)
   {
      if (is.null(xlab))
      {
         if (is.null(names.x)) xlab <- ""
         if (!is.null(names.x)) xlab <- names.x 
      }
   }  

   if (d==2)
   {
      if ((is.null(xlab))|(is.null(ylab)))
      {
         if (is.null(names.x)) 
         { 
            if (is.null(xlab)) xlab <- "" 
            if (is.null(ylab)) ylab <- ""
         }
         if (!is.null(names.x))
         { 
            if (is.null(xlab)) xlab <- names.x[1] 
            if (is.null(ylab)) ylab <- names.x[2]
         }
      }
   }

   if (d>=3)
   {
      if ((is.null(xlab))|(is.null(ylab))|(is.null(zlab)))
      {
         if (is.null(names.x)) 
         { 
            if (is.null(xlab)) xlab <- "" 
            if (is.null(ylab)) ylab <- ""
            if (is.null(zlab)) zlab <- ""
         }
         if (!is.null(names.x))
         { 
            if (is.null(xlab)) xlab <- names.x[1] 
            if (is.null(ylab)) ylab <- names.x[2]
            if (is.null(zlab)) zlab <- names.x[3]
         }
      }
   }    

   return(list(xlab=xlab,ylab=ylab,zlab=zlab))  
}

######## End of dfltLabs ########
