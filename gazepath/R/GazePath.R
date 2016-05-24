gazepath <-
function(data, x1, y1, x2 = NULL, y2 = NULL, distance, trial, height_px, height_mm, width_px, width_mm, res_x = 1280, res_y = 1024, samplerate = 500, method = 'Mould', posthoc = FALSE, thres_vel = 35, thres_dur = 100, min_dist = 250){
  ## Check if input is a data frame
  if(!is.data.frame(data)) {
    stop('please insert a data frame and define the column numbers of the variables')
  }
  ## Check distance
  data[,distance] <- ifelse(data[,distance] < min_dist, NA, data[,distance])
  D <- by(data[,distance], data[,trial], data.frame)
  ## Check input 1 or 2 eyes
  if(!is.null(x2) & !is.null(y2)){
    X <- by((data[,x1] + data[,x2]) / 2, data[,trial], data.frame)
    Y <- by((data[,y1] + data[,y2]) / 2, data[,trial], data.frame)
  } else {
    X <- by(data[,x1], data[,trial], data.frame)
    Y <- by(data[,y1], data[,trial], data.frame)
  }
  
  ## Convert single heigths and widths into vectors
  if(length(height_px) == 1) height_px <- rep(height_px, length(unique(data[,trial])))
  if(length(height_mm) == 1) height_mm <- rep(height_mm, length(unique(data[,trial])))
  if(length(width_px) == 1) width_px <- rep(width_px, length(unique(data[,trial])))
  if(length(width_mm) == 1) width_mm <- rep(width_mm, length(unique(data[,trial])))
  
  final <- 'Please insert a correct method'
  s <- NA
  
  if(method == 'Eyelink'){
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      final[[i]] <- Eyelink(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], Hz = samplerate)
    } 
  }
  
  if(method == 'Tobii'){
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      final[[i]] <- Tobii(cbind(X[[i]], Y[[i]]), D[[i]], thres_dur = 100, Hz = samplerate)
    } 
  }
  
  if(method == 'Mould'){
    fix <- thres_vel <- numeric()
    s <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      ## make sure there is at least 1 second of data avaible
      if(length(which(!is.na(X[[i]]))) > samplerate & length(which(!is.na(Y[[i]]))) > samplerate & length(which(!is.na(D[[i]]))) > samplerate){
        ## Calculate speed
        s[[i]] <- Speed_Deg(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], samplerate)
        ## Omit velocities over 1000 deg/s
        s[[i]] <- ifelse(s[[i]] > 1000, NA, s[[i]])
        thres_vel[i] <- Mould_vel(s[[i]], Hz = samplerate)
        fix <- c(fix, possiblefix(s[[i]], thres_vel[i]))
      } else {
        s[[i]] <- NA; thres_vel[i] <- NA
      }
    }
    ## Make sure there are enough data points to make the fit
    fix <- fix[!is.na(fix)]
    if(length(fix) > 100){
      thres_dur <- Mould_dur(fix, Hz = samplerate)
    } else {
      warning('There were not enough data points to derive a duration threshold, the threshold specified is used (default = 100ms)')
    }
    ## Make the final output with the classification of all samples
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      if(!is.na(thres_vel[i])){
        final[[i]] <- fixationANDsaccade(s[[i]], thres_vel[i], thres_dur, Hz = samplerate)
      } else {
        final[[i]] <- NA
      }
    } 
  }
  
  if(method == 'MouldDur'){
    fix <- thres_vel <- numeric()
    s <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      ## make sure there is at least 1 second of data avaible
      if(length(which(!is.na(X[[i]]))) > samplerate & length(which(!is.na(Y[[i]]))) > samplerate & length(which(!is.na(D[[i]]))) > samplerate){
        ## Calculate speed
        s[[i]] <- Speed_Deg(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], samplerate)
        ## Omit velocities over 1000 deg/s
        s[[i]] <- ifelse(s[[i]] > 1000, NA, s[[i]])
        thres_vel[i] <- Mould_vel(s[[i]], Hz = samplerate)
        fix <- c(fix, possiblefix(s[[i]], thres_vel[i]))
      } else {
        s[[i]] <- NA; thres_vel[i] <- NA
      }
    }
    ## Set fixed duration threshold
    thres_dur <- thres_dur
    ## Make the final output with the classification of all samples
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      if(!is.na(thres_vel[i])){
        final[[i]] <- fixationANDsaccade(s[[i]], thres_vel[i], thres_dur, Hz = samplerate)
      } else {
        final[[i]] <- NA
      }
    } 
  }
  
  if(method == 'Mould.all'){
    s <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      ## Calculate speed
      s[[i]] <- Speed_Deg(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], samplerate)
      ## Omit velocities over 1000 deg/s
      s[[i]] <- ifelse(s[[i]] > 1000, NA, s[[i]])
    }
    
    thres_vel <- Mould_vel(unlist(s), plot = F, Hz = samplerate)
    fix <- possiblefix(unlist(s), thres_vel)
    fix <- fix[!is.na(fix)]
    
    thres_dur <- Mould_dur(fix, plot = F, Hz = samplerate)
    
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      final[[i]] <- fixationANDsaccade(s[[i]], thres_vel, thres_dur, Hz = samplerate)
    }
  }
  
  if(method == 'Mould.fix'){
    s <- list()
    for(i in 1:length(unique(data[,trial]))){
      ## Boundary check
      X[[i]] <- Boundary(X[[i]], (res_x - width_px[i]) / 2, res_x - (res_x - width_px[i]) / 2)
      Y[[i]] <- Boundary(Y[[i]], (res_y - height_px[i]) / 2, res_y - (res_y - height_px[i]) / 2)
      ## Calculate speed
      s[[i]] <- Speed_Deg(X[[i]], Y[[i]], D[[i]], height_mm[i], width_mm[i], height_px[i], width_px[i], samplerate)
      ## Omit velocities over 1000 deg/s
      s[[i]] <- ifelse(s[[i]] > 1000, NA, s[[i]])
    }
    
    thres_vel <- Mould_vel(unlist(s), plot = F, Hz = samplerate)
    final <- list()
    for(i in 1:length(unique(data[,trial]))){
      final[[i]] <- fixationANDsaccade(s[[i]], thres_vel, thres_dur, Hz = samplerate)
    }
  }
  
  ## Post-hoc Check
  if(posthoc == TRUE){
    for(i in 1:length(X)) {
      PH <- posthocCheck(final[[i]], X[[i]], Y[[i]])
      final[[i]] <- PH[[1]]
      X[[i]] <- PH[[2]]
      Y[[i]] <- PH[[3]]
    }
  }
  ## Simplify the raw classification
  sim <- list()
  for(i in 1:length(X)){
    sim[[i]] <- simplify(final[[i]], X[[i]], Y[[i]], samplerate, D[[i]], width_px[i], width_mm[i])
  }
  
  ## Determine robustness and precision
  Robustness <- sapply(1:length(X), function(i) robust(X[[i]], samplerate))
  Pre_x <- sapply(1:length(X), function(i) precision(X[[i]], samplerate))
  Pre_y <- sapply(1:length(X), function(i) precision(Y[[i]], samplerate))
  Precision <- (Pre_x + Pre_y) / 2
  
  output <- list(final, X, Y, method, Robustness, Precision, thres_vel, thres_dur, s, samplerate, D, height_px, height_mm, width_px, width_mm, sim)
  class(output) <- 'gazepath'
  return(output)
}
