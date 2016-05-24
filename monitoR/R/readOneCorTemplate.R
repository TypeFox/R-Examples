# For reading in an individual correlation template
# Should not be called directly
# Modified: 2013 JULY 17

readOneCorTemplate <-
function(file) {
   lns <- readLines(file)
   nr <- length(lns)
   n.pt <- which(lns == "Time bin step") - 7
  
   clip <- lns[2]
   samp.rate <- as.numeric(lns[4])

   pts <- NULL
   for(i in 1:n.pt + 6) {
      pts <- rbind(pts, as.numeric(strsplit(lns[i], ' ')[[1]]))
   }
   colnames(pts) <- c('t', 'frq', 'amp')

   t.step <- as.numeric(lns[n.pt + 8])
   frq.step <- as.numeric(lns[n.pt + 10])
   n.t.bins <- as.numeric(lns[n.pt + 12])
   first.t.bin <- as.numeric(lns[n.pt + 14])
   n.frq.bins <- as.numeric(lns[n.pt + 16])
   duration <- as.numeric(lns[n.pt + 18])
   frq.lim <- as.numeric(strsplit(lns[n.pt + 20], ' ')[[1]])
   wl <- as.numeric(lns[n.pt + 22])
   ovlp <- as.numeric(lns[n.pt + 24])
   wn <- as.character(lns[n.pt + 26])
   score.cutoff <- as.numeric(lns[n.pt + 28])
   comment <- as.character(lns[n.pt + 30])

   template <- new('corTemplate', clip.path=clip, samp.rate=as.integer(samp.rate), pts=pts, t.step=t.step, frq.step=frq.step, n.t.bins=as.integer(n.t.bins), first.t.bin=first.t.bin, n.frq.bins=as.integer(n.frq.bins), duration=duration, frq.lim=frq.lim, wl=as.integer(wl), ovlp=as.integer(ovlp), wn=wn, score.cutoff=score.cutoff, comment=comment)
   return(template)
}
