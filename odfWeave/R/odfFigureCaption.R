# Max wants this to be 'stylified'
odfFigureCaption <- function(caption, numformat='1', numlettersync=FALSE,
      formula='Illustration+1', label='Illustration')
{
   # sanity check numformat argument
   if (!any(numformat == c('A', 'a', 'I', 'i', '1')))
      stop('illegal numformat value: ', numformat)

   # allocate a unique refname using information collected during
   # the preprocessing phase
   refname <- uniqueRef('Illustration', 'refIllustration')

   # save all information needed to generate the caption once
   # we are ready to generate the XML for the figure/illustration
   .odfEnv$fig.caption <- list(caption=caption, numformat=numformat,
                               numlettersync=numlettersync, formula=formula,
                               refname=refname, label=label)
   invisible(refname)
}


xRef <- function(name, type = "value", initial = NULL)
  {

    ## todo: check type values
    out <- paste(
                 "<text:sequence-ref text:reference-format=\"",
                 type,
                 "\" text:ref-name=\"",
                 name,
                 "\">",
                 initial,
                 "</text:sequence-ref>",
                 sep = "")
    out

  }
