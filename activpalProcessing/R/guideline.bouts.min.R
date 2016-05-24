guideline.bouts.min <-
function(mets) {
	total.active.time.in.bouts <- 0

	# indices where transitions take place
	mets.length <- length(mets)
	one <- mets[-mets.length]
	two <- mets[-1]

	trans <- c( FALSE, ((one<3)&(two>=3)) | ((one>=3)&(two<3)) )
	trans.inds <- c(1, seq_along(mets)[trans], (mets.length+1))

	# how long are the periods of activity and inactivity
	durations <- trans.inds[-1]-trans.inds[-length(trans.inds)]

	# identify if interval is activity or inactivity (they alternate)
	active.interval <- rep(FALSE,length=length(durations))

	if (mets[1]<3)
		active.interval <- rep(c(FALSE, TRUE),length=length(durations))
	if (mets[1]>=3)
		active.interval <- rep(c(TRUE, FALSE),length=length(durations))


	# Create some empty vectors which will be used to keep track of the
	# start and end points of the bouts in the durations vector.
	bout.starts <- c()
	bout.ends <- c()

	# Create some variables which will be used in constructing the bouts.
	active.inds <- seq_along(durations)[active.interval]
  if(length(active.inds) > 0) {
  	possible.bout.active.inds <- 1
  	possible.bout.active.inds.length <- 1
  	bout.inactivity <- rep(0, durations[active.inds[possible.bout.active.inds]])
    end.possible.bout <- FALSE

  	while(possible.bout.active.inds[possible.bout.active.inds.length] < length(active.inds)) {
  		# Determine if adding the next inactive interval to the current possible bout
  		# would cause the last 10 minutes of the bout to be more than 20% inactive

  		bout.inactivity.with.next.interval <- c(bout.inactivity, rep(1, durations[active.inds[possible.bout.active.inds[possible.bout.active.inds.length]] + 1]))
  		secs.to.extract <- min(length(bout.inactivity.with.next.interval), 10*60)
  		last.10.min.of.possible.bout <- bout.inactivity.with.next.interval[seq_len(secs.to.extract) + (length(bout.inactivity.with.next.interval) - secs.to.extract)]

  		if( sum(last.10.min.of.possible.bout) <= 2*60 ) {
  			# If adding the next inactive interval means the last 10 min has <= 2 min of inactivity,
  			# add it and the next active interval to the possible-bout.
  			bout.inactivity <- c(bout.inactivity.with.next.interval, rep(0, durations[active.inds[possible.bout.active.inds[possible.bout.active.inds.length]] + 2]))
  			possible.bout.active.inds <- c(possible.bout.active.inds, possible.bout.active.inds[possible.bout.active.inds.length] + 1)
  			possible.bout.active.inds.length <- possible.bout.active.inds.length + 1
        if(possible.bout.active.inds[possible.bout.active.inds.length] == length(active.inds)) {
          # If we just added the last active interval to the possible bout, end the possible bout.
          end.possible.bout <- TRUE
        }
  		} else {
    		# If adding the next inactive interval means the last 10 min has > 2 min of inactivity,
  			# stop building this bout.
        end.possible.bout <- TRUE
  		}
      
      if(end.possible.bout) {
  			# If the possible bout is long enough, add it to the list of bouts.
  			# Reset to start building the next possible bout.

  			possible.bout.length <- sum(durations[seq(from=active.inds[possible.bout.active.inds[1]], to=active.inds[possible.bout.active.inds[possible.bout.active.inds.length]])]) 
  			if(possible.bout.length >= 10*60) {
  				# If the total duration of the possible bout is >= 10 min, it is a bout.

  				# The start position of the bout is recorded as an index in the input mets vector.
  				# This is calculated as the sum of the durations of all intervals before the start of this bout, plus 1
  				bout.starts <- c(bout.starts, sum(durations[seq_len(active.inds[possible.bout.active.inds[1]] - 1)]) + 1)

  				# The end position of the bout is recorded as an index in the input mets vector.
  				# This is calculated as the sum of the durations of all intervals up through the end this bout
  				bout.ends <- c(bout.ends, sum(durations[seq_len(active.inds[possible.bout.active.inds[possible.bout.active.inds.length]])]))

				# add the active time in this bout to the total active time in all bouts
				total.active.time.in.bouts <- total.active.time.in.bouts + length(bout.inactivity) - sum(bout.inactivity)
  			} else {
  				# If it is < 10 min, see if we can add on some inactive time to the beginning and/or end
  				# to reach the 10 min threshold and make it into a bout.  Otherwise, discard it.

  				# Calculate the number of seconds needed to fill in to reach the 10 minute threshold.
  				seconds.missing <- 10*60 - possible.bout.length
  				if(seconds.missing + sum(bout.inactivity) <= 2*60) {
  					# If the number of seconds we need to fill in to reach 10 minutes
  					# would not push us over the 10 minute threshold, add in time
  					# from the beginning and then from the end to reach 10 minutes.

  					# Get the current start time and end time of the bout, calculated as above.
  					bout.start <- sum(durations[seq_len(active.inds[possible.bout.active.inds[1]] - 1)]) + 1
  					bout.end <- sum(durations[seq_len(active.inds[possible.bout.active.inds[possible.bout.active.inds.length]])])

  					# Find the number of seconds we can add on to the beginning.
  					# We do not want to overlap with the previous bout.
  					if(length(bout.ends) > 0) {
  						last.bout.end <- bout.ends[length(bout.ends)]
  					} else {
  						last.bout.end <- 0
  					}
  					seconds.to.add <- min(bout.start - last.bout.end + 1, seconds.missing)
  					bout.start <- bout.start - seconds.to.add
  					seconds.missing <- seconds.missing - seconds.to.add

  					# If necessary, add the rest of the time on to the end.
  					# (We know there are "enough" extra seconds at the end because we exceeded the 2 minute threshold
  					# when we tried to add on the whole inactive interval.)
  					if(seconds.missing > 0) {
  						bout.end <- bout.end + seconds.missing
  					}

  					# Add the new bout to the list of bouts.
  					bout.starts <- c(bout.starts, bout.start)
  					bout.ends <- c(bout.ends, bout.end)

					# add the active time in this bout to the total active time in all bouts
					total.active.time.in.bouts <- total.active.time.in.bouts + length(bout.inactivity) - sum(bout.inactivity)
  				}
  			}

  			# Set up to start building the next possible bout.
  			possible.bout.active.inds <- possible.bout.active.inds[possible.bout.active.inds.length] + 1
  			possible.bout.active.inds.length <- 1

        end.possible.bout <- FALSE

        if(possible.bout.active.inds[possible.bout.active.inds.length] < length(active.inds)) {
  				bout.inactivity <- rep(0, durations[active.inds[possible.bout.active.inds]])
  			}
  		}
  	}
  }


		return(sum((bout.ends + 1) - bout.starts)/60)

}

