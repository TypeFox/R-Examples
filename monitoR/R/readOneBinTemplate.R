# For reading in an individual binary template
# Should not be called directly
# Modified: 2013 JULY 17

readOneBinTemplate <-
function(file) {
   lns <- readLines(file)
   nr <- length(lns)
   n.pt.on <- which(lns == "Points off") - 7
   n.pt.off <- which(lns == "Time bin step") - n.pt.on - 9

   clip.path <- lns[2]
   samp.rate <- as.numeric(lns[4])
  
   pt.on <- NULL
   for(i in 1:n.pt.on + 6) {
      pt.on <- rbind(pt.on, as.numeric(strsplit(lns[i], ' ')[[1]]))
   }
   colnames(pt.on) <- c('t', 'frq')
 
   pt.off <- NULL
   for(i in 1:n.pt.off + 8 + n.pt.on) {
      pt.off <- rbind(pt.off, as.numeric(strsplit(lns[i], ' ')[[1]]))
   }
   colnames(pt.off) <- c('t', 'frq')

   t.step <- as.numeric(lns[n.pt.on + n.pt.off + 10])
   frq.step <- as.numeric(lns[n.pt.on + n.pt.off + 12])
   n.t.bins <- as.numeric(lns[n.pt.on + n.pt.off + 14])
   first.t.bin <- as.numeric(lns[n.pt.on + n.pt.off + 16])
   n.frq.bins <- as.numeric(lns[n.pt.on + n.pt.off + 18])
   duration <- as.numeric(lns[n.pt.on + n.pt.off + 20])
   frq.lim <- as.numeric(strsplit(lns[n.pt.on + n.pt.off + 22], ' ')[[1]])
   wl <- as.numeric(lns[n.pt.on + n.pt.off + 24])
   ovlp <- as.numeric(lns[n.pt.on + n.pt.off + 26])
   wn <- as.character(lns[n.pt.on + n.pt.off + 28])
   score.cutoff <- as.numeric(lns[n.pt.on + n.pt.off + 30])
   comment <- as.character(lns[n.pt.on + n.pt.off + 32])

   template <- new('binTemplate', clip.path=clip.path, samp.rate=as.integer(samp.rate), pt.on=pt.on, pt.off=pt.off, t.step=t.step, frq.step=frq.step, n.t.bins=as.integer(n.t.bins), first.t.bin=first.t.bin, n.frq.bins=as.integer(n.frq.bins), duration=duration, frq.lim=frq.lim, wl=as.integer(wl), ovlp=as.integer(ovlp), wn=wn, score.cutoff=score.cutoff, comment=comment)
   return(template)
}
