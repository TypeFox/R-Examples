uniqueRef <- function(seqname, prefix) {
   refnames <- .odfEnv$seqInfo[[seqname]]
   if (is.null(refnames))
      stop('sequence not declared in document: ', seqname)

   i <- 0
   repeat {
      n <- paste(prefix, i, sep='')
      if (!any(n == refnames)) break
      i <- i + 1
   }
   # add this new refname to the list so it isn't used again
   .odfEnv$seqInfo[[seqname]] <- c(refnames, n)
   n
}

# XXX this needs to handle escaping special characters
odfTableCaption <- function(caption, numformat='1', numlettersync=FALSE,
      formula='Table+1', label='Table')
{
   if (!any(numformat == c('A', 'a', 'I', 'i', '1')))
      stop('illegal numformat value: ', numformat)

   refname <- uniqueRef('Table', 'refTable')

   cat(sprintf('<text:p text:style-name="Table">%s ', label))
   cat(sprintf('<text:sequence style:num-format="%s"', numformat))
   if (numlettersync)
      cat(' style:num-letter-sync="true"')
   cat(sprintf(' text:formula="ooow:%s" text:name="Table"', formula))
   cat(sprintf(' text:ref-name="%s">', refname))
   cat(sprintf('999</text:sequence>: %s</text:p>', caption))
   cat('\n')  # see if this removes the warning messages

   # return the reference name so it can be used to create a cross-reference
   invisible(refname)
}
