stimfunction <-
function(totaltime, onsets, durations, accuracy){

  if(max(onsets)>totaltime){
    stop("Mismatch between onsets and totaltime")
  }
  s <- rep(0,totaltime/accuracy)
  os <- onsets/accuracy
  dur <- durations/accuracy
  if(length(durations)==1){
	dur <- dur*rep(1,length(onsets))
  } else if(length(durations) != length(onsets)) {
		stop("Mismatch between number of onsets and number of durations.")
  }
  for(i in (1:length(onsets))){
	if((os[i]+dur[i]) <= totaltime/accuracy){
    		s[c(os[i]:(os[i]+dur[i]))] <- 1
	} else {
		s[c(os[i]:(totaltime/accuracy))] <- 1
	}
  }
  return(s)
}

