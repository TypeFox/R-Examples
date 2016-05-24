# Modified: 6 Sept 2015

showPeaks <-
function(
   detection.obj,                               # Complete output from corMatch or binMatch
   which.one=names(detection.obj@templates)[1], # Name or position of list elements that should be used
   fd.rat=4,                                    # Plot frame to template duration frame ratio. Used only if frame is not specified.
   frame=fd.rat*detection.obj@templates[[which.one]]@duration, # Length of time to be plotted.
   id=1:nrow(pks),                 # The row number or name of hits or peaks (within pks) that are to be plotted.
   t.lim,                          # Instead of id numbers, give time limits
   flim=c(0, 20),                   # Frequency limits for the spectrogram.
   point=TRUE,                     # Use an arrow to identify the peak center.
   ask=if(verify) FALSE else TRUE, # If TRUE, pauses between plots.
   scorelim=NULL,                  # Plot limits for scores
   verify=FALSE,                   # Set to true for verification
   what='detections',              # 'detections' for just detections, 'peaks' for all peaks
   box=TRUE,                       # Set to FALSE to suppress box in spectrogram or "template" 
   player='play',                  # Command to call up external wave player, e.g., 'wv_player.exe' (Windows) or 'play' (SoX in Linux)
   spec.col=gray.3(),              # Color palette for spectrogram
   on.col='#FFA50075',             # Color for on points in binary templates
   off.col='#0000FF75',            # Color for off points in binary templates
   pt.col='#FFA50075'              # Color for correlation template points
) {

   # Check arguments
   if(missing(detection.obj)) stop('Required argument detection.obj missing')

   # ask par
   oldask <- par(ask=par('ask'))
   on.exit(par(oldask))

   # Pull out data for plots
   template <- detection.obj@templates[[which.one]]
   pks <- slot(detection.obj, what)[[which.one]]
   survey <- detection.obj@survey
   score <- detection.obj@scores[[which.one]]
   # Pull out spectrogram data from score.obj object
   amp <- detection.obj@survey.data[[which.one]]$amp
   t.bins <- detection.obj@survey.data[[which.one]]$t.bins
   frq.bins <- detection.obj@survey.data[[which.one]]$frq.bins
   which.frq.bins <- which(frq.bins >= flim[1] & frq.bins <= flim[2])
   frq.bins <- frq.bins[which.frq.bins]
   amp <- amp[which.frq.bins, ]
   # Get score limits for plot
   if(is.null(scorelim)) scorelim <- c(0, max(pks$score))

   # For verification
   if(verify) verC <- NULL

   # Length of survey
   dur.survey <- length(survey@left)/survey@samp.rate

   # If t.lim is given, use it to get ids
   if(!missing(t.lim)) id <- which(pks$time>=min(t.lim) & pks$time<=max(t.lim))

   # Exit if there are no hits to show
   if(length(id) == 2) if(all(id[1] == 1, id[2] == 0)) stop('No peaks selected.')
   # THERE IS A FUNCTION NAMED ID, SHOULD CHANGE NAME HERE

   # Loop through all peaks requested
   i <- min(id)
   while(i<=max(id)) {
      x <- 'x'

      t.start <- max(pks$time[i] - frame/2, 0)
      t.end <- min(pks$time[i] + frame/2, dur.survey)

      survey.clip <- cutWave(wave=survey, from=t.start, to=t.end)
      
      par(mfrow=c(2, 1), mar=c(1, 4,4, 1))
      score.clip <- score[score[, 'time']>=t.start & score[, 'time']<=t.end, ]

      # Make spectrogram
      times <- t.bins[t.bins>=t.start & t.bins<=t.end]
      amp.clip <- amp[, t.bins %in% times]
      image(x=times, y=frq.bins, t(amp.clip), col=spec.col, xlab='', ylab='Frequency (kHz)', xaxt='n', las=1, main=paste(if(what == "detections") "Detection" else "Peak", i))
 
      if(box == TRUE & nrow(pks)>0) {
         xleft <- pks$time[i]-template@duration/2
         xright <- pks$time[i]+template@duration/2
         ylwr <- template@frq.lim[1]
         yupr <- template@frq.lim[2]
         polygon(x=c(xleft, xleft, xright, xright), y=c(ylwr, yupr, yupr, ylwr), border='blue')
      } else if(tolower(box) == 'template' & nrow(pks)>0) {
         xleft <- pks$time[i]-template@duration/2
         ylwr <- template@frq.lim[1]
         # Plot template points
         if(class(template) == 'binTemplate') {
            pt.on <- template@pt.on
            pt.off <- template@pt.off
            pt.on[, 't'] <- pt.on[, 't'] + ((xleft-min(times))/template@t.step)
            pt.off[, 't'] <- pt.off[, 't'] + ((xleft-min(times))/template@t.step)
            pt.on[, 'frq'] <- pt.on[, 'frq'] + ylwr - 1 
            pt.off[, 'frq'] <- pt.off[, 'frq'] + ylwr - 1
            
            bin.amp <- 0*amp.clip
            bin.amp[pt.on[, c(2, 1)]] <- 1
            bin.amp[pt.off[, c(2, 1)]] <- 2
      
            image(x=times, y=frq.bins, t(bin.amp), zlim=c(0, 2), col=c('transparent', on.col, off.col), add=TRUE)
      
         } else if (class(template) == 'corTemplate') {
      
            pts <- template@pts
            pts[, 't'] <- pts[, 't'] + ((xleft-min(times))/template@t.step)
            pts[, 'frq'] <- pts[, 'frq'] + ylwr - 1
            bin.amp <- 0*amp.clip
            bin.amp[pts[, c(2, 1)]] <- 1
      
            image(x=times, y=which.frq.bins, t(bin.amp), zlim=c(0, 1), col=c('transparent', pt.col), add=TRUE)
      
         } else stop('Template list class not recognized: ', class(template))
      }

      par(mar=c(4, 4,1, 1))
      plot(score.clip$time, score.clip$score, xlim=c(t.start, t.end), ylim=scorelim, xlab='Time (s)', ylab='Score', type='l', xaxs='i', las=1)
      if(point) {
         t <- pks$time[pks$time>=t.start & pks$time<=t.end]
         score.center <- pks$score[pks$time>=t.start & pks$time<=t.end]
         points(t, score.center, col='black')     
         t <- pks$time[i]
         score.center <- score.clip$score[score.clip$time == pks$time[i]] 
         score.center <- pks$score[i]
         points(t, score.center, pch=21, col='black', bg='red')
         text(x=t, y=score.center, labels=round(score.center, 3), pos=3, offset=1)     
      }
      abline(h=template@score.cutoff, lty=2)

      if(verify) {
         while(length(x) == 0 || !x %in% c('y', 'n', 'r', NA)) {
            cat(paste0('\n', i,'. True detection?\n Enter y for yes, n for no, NA for NA, p to play, r to rewind, or q to exit: '))      
            #x <- tolower(scan(n=1, what='character', quiet=TRUE))
            x <- tolower(readLines(n=1)[1])

            if(length(x) == 0) {
              cat('\nYou didn\'t enter a response.\n')
              next
            }

            if(!is.na(x) && x == 'na') x <- NA
            if(is.na(x)) {
              cat('NA\n')
              break
            }

            cat(switch(x, n=FALSE, y=TRUE, p='Playing clip', r='Previous', q='Exiting', 'Value not recognized. Enter y, n, NA, p, r, or q.'), '\n')

            if(!x %in% c('y', 'n', 'p', 'r', 'q')) next

            if(x == 'q') return()

            if(x == 'p') {
               tuneR::writeWave(object=survey.clip, filename=tempname <- tempfile(fileext='.wav'))
               # Variation on next line may be needed if player is slow
               #Sys.sleep(2) 
               cat("Shell command:", paste(player, tempname), '\n')
               if(tolower(Sys.info()['sysname']) == 'windows') shell(cmd=paste(player, tempname), wait=FALSE)
               else system(command=paste(player, tempname), wait=FALSE)
               # tuneR::writeWave(object=survey.clip, filename=tempname <- paste0('temp', Sys.time(), '.wav'))
               # cat("Shell command:", paste0(player, " \"", tempname, "\""), '\n')
               # system(command=paste0(player, " \"", tempname, "\""), wait=FALSE)
               prev.line <- 0
               for(j in seq(t.start, t.end, length.out=20)) {
                  t1 <- Sys.time()
                  abline(v=prev.line, col='lightgray', lwd=3)
                  abline(v=j, col='black', lwd=2)
                  prev.line <- j
                  delta.t <- Sys.time()-t1
                  Sys.sleep(max(0, (t.end-t.start)/19 - as.numeric(delta.t)))
               }
               # Remove wav file
               file.remove(tempname)
            }
         }
         if(is.na(x) || x != 'r') verC[i] <- x
      }
      par(ask=ask)
      if(!is.na(x) && x == 'r') i <- i-1 else i <- i+1
      if(i<1) i <- 1
   }

   cat("\n")
   if(verify) {
      slot(detection.obj, what)[[which.one]]$true <- NA
      slot(detection.obj, what)[[which.one]]$true[id] <- verC == 'y'
      return(detection.obj)
   }
   return(invisible(NULL))
}
