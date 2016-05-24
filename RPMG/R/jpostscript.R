`jpostscript` <-
function(file=NULL, P=NULL, w=NULL, h=NULL)
{
########   suffix eps will be added to the file name
  if(missing(file)) { file="JPLOT" }
  if(missing(P)) { P = NULL }
  if(missing(w)) { w = NULL }
  if(missing(h)) { h = NULL }

   if(is.null(P))
    {
      P = par('din')
      P = round(P, digits=2)
      
    }
   if(!is.null(w))
    {
      P[1] = w
    }
   if(!is.null(h))
    {
      P[2] = h
    }
  
  
  psname = local.file(file, "eps")


  ## P = round(par('pin'))
  
  postscript(file=psname , width=P[1], height=P[2], paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)

  return(psname)
  
}

