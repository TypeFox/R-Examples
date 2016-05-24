`find.traj` <-
function(x, y, box, peel.alpha, paste.alpha, mass.min, threshold,
                     d, n, pasting, verbose=FALSE, paste.all=FALSE,
										 peel_crit)
{
 # mass.min <- ceiling(1/peel.alpha)/nrow(x)  #commented out 2009-10-25
                                              #suspect it is incorrect

  peel.traj <- list()

  y.mean <- mean(y)
  mass <- length(y)/n

  if ((y.mean >= threshold) & (mass >= mass.min))
    boxk.peel <- peel.one(x=x, y=y, box=box, peel.alpha=peel.alpha,
                          mass.min=mass.min,threshold=threshold, d=d, n=n,
													peel_crit=peel_crit)
  else
    boxk.peel <- NULL

  boxk.temp <- NULL

  bi <- 0

  while (!is.null(boxk.peel)){
    bi <- bi+1
    boxk.temp <- boxk.peel

    peel.traj[[bi]] <- boxk.temp #$box
		#print("calling peel.one")
		flush.console()
    boxk.peel <- peel.one(x=boxk.temp$x, y=boxk.temp$y, box=boxk.temp$box,
                          peel.alpha=peel.alpha,
                          mass.min=mass.min, threshold=threshold, d=d, n=n,
													peel_crit=peel_crit)
    #cat("mean:",boxk.peel$y.mean, "mass",boxk.peel$mass,"\n") #diagnostic
	
	}

  if (verbose){
		cat("Peeling completed \n")
	}

  if (!pasting){
		paste.traj <- lapply(peel.traj, "[[", "box")
	} else {

    paste.traj <- list()

#    boxk.paste <- boxk.temp

    for (p in 1:length(peel.traj)){
  
      boxk.paste <- peel.traj[[p]]
  
      while (!is.null(boxk.paste))
      {
        boxk.temp <- boxk.paste
        boxk.paste <- paste.one(x=boxk.temp$x, y=boxk.temp$y, box=boxk.temp$box,
                                x.init=x, y.init=y, paste.alpha=paste.alpha,
                                mass.min=mass.min, threshold=threshold, d=d, n=n)
      }
			
      if (verbose) {
				cat("Pasting completed\n")
			}
		
      paste.traj[[p]] <- boxk.temp$box
  
    }

  }

	#browser()
  boxk <- boxk.temp

#Lines added by BB to accomodate peeling trajectory
  paste.seq <- list()
  paste.seq$box <- paste.traj
  paste.seq$num.class <- length(paste.seq$box)  #THIS SEEMS STUPID!

  return(paste.seq)

}

