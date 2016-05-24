# $Id: LaTeX.R 45 2008-05-23 15:50:48Z mentus $ 

as.latex <- function(x,label=NULL,inline=ifelse(is.null(label),TRUE,FALSE),count=ifelse(is.null(label),FALSE,TRUE))
{
  out <- list(alt=x,inline=inline,count=count,label=label)  
  class(out)<-"latex"
  return(out)
}


"HTML.latex" <- function(x,file = HTMLGetFile(),...)
{
  ### Note: no append argument as it could ONLY be happened to work...
  # count: add a (#)
  # label: add before: Equation (#): label
  if (! is(x,"latex")) x <- as.latex(x)
  if (x$inline)
  {
    # inline : directly there in the text (cant be counted or labeled)
    cat(paste("`",x$alt,"`",sep=""),file=file,append=TRUE,sep= " ")
  }
  else
  {
    # not inline: own living space (will be within a table to center it)
     if (x$count) { 
      cat("\n<script>\n nequations=nequations+1;\n document.write(\"<a name='equation\"+nequations+\"'>&nbsp;</a>\")</script>\n",file=file,append=TRUE)

     }
     if (!is.null(x$label)) 
      {
        txt <- "\n<br /><span class='equation'>Equation"
        if (x$count) txt <- paste(txt, "<script>document.write(nequations);</script>")
        txt <- paste(txt, "-",x$label)
        cat(txt,file=file, append=TRUE)     
      }
      if (!x$count) 
          cat(paste("<br /><center><table border=0><td align='center'>`",x$alt,"`</td></table></center><br />",sep=""),file=file,append=TRUE)
      else
          cat(paste("\n<br /><center><table border=0 width='90%'><td width='50'>&nbsp;</td><td align='center'>`",x$alt,"`</td><td align='right' width='50'><script>document.write('('+nequations+')')</script></td></table></center><br />",sep=""),file=file,append=TRUE)
    }
}

