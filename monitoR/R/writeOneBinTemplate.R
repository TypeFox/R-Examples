# Not meant to be used directly by users
# Modified: 2013 NOV 13


writeOneBinTemplate <-
function(
   template, 
   file
) {
   olddigits <- options(digits=15)
   suppressWarnings(
      {
      write(file=file, 'Original recording')
      write(file=file, template@clip.path, append=TRUE)
      write(file=file, 'Sample rate', append=TRUE)
      write(file=file, template@samp.rate, append=TRUE)
      write(file=file, 'Points on', append=TRUE)
      write.table(file=file, template@pt.on, append=TRUE, row.names=FALSE, quote=FALSE, sep=' ')
      write(file=file, 'Points off', append=TRUE)
      write.table(file=file, template@pt.off, append=TRUE, row.names=FALSE, quote=FALSE, sep=' ')
      write(file=file, 'Time bin step', append=TRUE)
      write(file=file, template@t.step, append=TRUE)
      write(file=file, 'Frequency bin step', append=TRUE)
      write(file=file, template@frq.step, append=TRUE)
      write(file=file, 'Number of time bins spanned', append=TRUE)
      write(file=file, template@n.t.bins, append=TRUE)
      write(file=file, 'First time bin', append=TRUE)
      write(file=file, template@first.t.bin, append=TRUE)
      write(file=file, 'Number of frequency bins spanned', append=TRUE)
      write(file=file, template@n.frq.bins, append=TRUE)
      write(file=file, 'Duration', append=TRUE)
      write(file=file, template@duration, append=TRUE)
      write(file=file, 'Frequency limits', append=TRUE)
      write(file=file, template@frq.lim, append=TRUE)
      write(file=file, 'Window length', append=TRUE)
      write(file=file, template@wl, append=TRUE)
      write(file=file, 'Window overlap', append=TRUE)
      write(file=file, template@ovlp, append=TRUE)
      write(file=file, 'Window type', append=TRUE)
      write(file=file, template@wn, append=TRUE)
      write(file=file, 'Score cutoff', append=TRUE)
      write(file=file, template@score.cutoff, append=TRUE)
      write(file=file, 'Comment', append=TRUE)
      write(file=file, template@comment, append=TRUE)
      }
   )
   options(olddigits)
   message('Template written to file "', file, '"')
}


