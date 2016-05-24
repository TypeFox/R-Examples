## functions to extract GPS features

extractFeatsPALMSOneFile = function(inputFile, outputDir, winSize, names=NULL) {
  # splits a single GPS file from PALMS by days
  
  all_GPS = read.csv(inputFile, header=TRUE, stringsAsFactors=FALSE, 
                     colClasses=c(identifier="character"))
  #figure out the sample rate
  t1 = strptime(all_GPS[1, c("dateTime")], "%Y-%m-%d %H:%M:%S")
  t2 = strptime(all_GPS[2, c("dateTime")], "%Y-%m-%d %H:%M:%S")
  sampleRate = as.numeric(t2 - t1)
  ws = winSize / sampleRate
  
  if (is.null(names)) {
    names = unique(all_GPS$identifier)
  }
  for (id in 1:length(names)) {
    outputFile = file.path(outputDir, names[id])
    if (!file.exists(outputFile)) {
      cat(names[id], "...\n")
      cols = c("identifier","dateTime","speed","distance","duration","ele",
              "elevationDelta","lat","lon","nsatUsed","nsatView","snrUsed",
              "snrView","fixType")
      GPS = all_GPS[all_GPS$identifier == names[id], cols]
    
      r = 1
      lastCoordinates = GPS[r, c("lat", "lon")]
      st = strptime(GPS[r, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      newStart = alignStart(winSize, st)
      while (newStart > st) {
        lastCoordinates = GPS[r, c("lat", "lon")]
        r = r + 1
        st = strptime(GPS[r, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      }
      day = st$mday
      out = file.path(outputFile, strftime(st, "%Y-%m-%d"))
      cat(strftime(st, "%Y-%m-%d"), '\n')
      if (!file.exists(outputFile)) {
        dir.create(outputFile, recursive=TRUE)
      }
      cat("timestamp,avgspeed,sdspeed,coefvar,netdistance,totaldistance,ratiodistance,elechange,avgele,nsat,snr\n", file=out, append=TRUE)
    
      # first feature
      while ((r + ws - 1) <= nrow(GPS)) {
        if (st$mday != day) {
          out = file.path(outputFile,strftime(st, "%Y-%m-%d"))
          cat(strftime(st, "%Y-%m-%d"), '\n')
          cat("timestamp,avgspeed,sdspeed,coefvar,netdistance,totaldistance,ratiodistance,elechange,avgele,nsat,snr\n", file=out, append=TRUE)
          day = st$mday
        }
        feat = computeOneGPSFeat(GPS[r:(r + ws - 1), ], lastCoordinates)
        st = strptime(GPS[r, c("dateTime")],"%Y-%m-%d %H:%M:%S")
        cat(strftime(st, "%Y-%m-%d %H:%M:%S,"), file=out, sep = "", append=TRUE)
        cat(feat, file=out, sep=",", append=TRUE)
        cat('\n', file=out, append=TRUE)
        lastCoordinates = GPS[r, c("lat","lon")]
        r = r + ws
      }
    }
  }
}
extractFeatsPALMSDir = function(inputDir, outputDir, winSize, names = NULL) {
  # splits a directory containing GPS files from PALMS by days
  
  if (is.null(names)) {
    names = file_path_sans_ext(list.files(inputDir))
  }
  for (i in 1:length(names)) {
    outputFile = file.path(outputDir, names[i])
    if (!file.exists(outputFile)) {
      cat(names[i], "...\n")
      GPS = read.csv(file.path(inputDir, paste0(names[i], ".csv")), header=TRUE, 
                     stringsAsFactors=FALSE, colClasses=c(identifier="character"))
      #figure out the sample rate
      t1 = strptime(GPS[1, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      t2 = strptime(GPS[2, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      sampleRate = as.numeric(t2 - t1)
      ws = winSize / sampleRate
  
      cols = c("identifier","dateTime","speed","distance","duration","ele",
              "elevationDelta","lat","lon","nsatUsed","nsatView","snrUsed",
              "snrView","fixType")
    
      r = 1
      lastCoordinates = GPS[r, c("lat","lon")]
      st = strptime(GPS[r, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      newStart = alignStart(winSize, st)
      while (newStart > st) {
        lastCoordinates = GPS[r, c("lat", "lon")]
        r = r + 1
        st = strptime(GPS[r, c("dateTime")], "%Y-%m-%d %H:%M:%S")
      }
      day = st$mday
      out = file.path(outputFile, strftime(st, "%Y-%m-%d"))
      cat(strftime(st, "%Y-%m-%d"), '\n')
      if (!file.exists(outputFile)) {
        dir.create(outputFile, recursive=TRUE)
      }
      cat("timestamp,avgspeed,sdspeed,coefvar,netdistance,totaldistance,ratiodistance,elechange,avgele,nsat,snr\n", file=out, append=TRUE)
    
      # first feature
      while ((r + ws - 1) <= nrow(GPS)) {
      
        if (st$mday != day) {
          out = file.path(outputFile, strftime(st, "%Y-%m-%d"))
          cat(strftime(st, "%Y-%m-%d"), '\n')
          cat("timestamp,avgspeed,sdspeed,coefvar,netdistance,totaldistance,ratiodistance,elechange,avgele,nsat,snr\n", file=out, append=TRUE)
          day = st$mday
        }
        feat = computeOneGPSFeat(GPS[r:(r + ws - 1), ], lastCoordinates)
        st = strptime(GPS[r, c("dateTime")], "%Y-%m-%d %H:%M:%S")
        cat(strftime(st, "%Y-%m-%d %H:%M:%S,"), file=out, sep = "", append=TRUE)
        cat(feat, file=out, sep=",", append=TRUE)
        cat('\n', file=out, append=TRUE)
        lastCoordinates = GPS[r, c("lat", "lon")]
        r = r + ws
      }
    }
  }
}
computeOneGPSFeat = function(w, lastCoordinates) {
  # compute one GPS feature vector from a window
  # input: identifier,dateTime,speed,distance,duration,ele,elevationDelta,lat,lon,nsatUsed,nsatView,snrUsed,snrView,fixType
  fAvgSpeed = mean(w[, c("speed")])  # average speed
  fSdSpeed = sd(w[, c("speed")])  # std of speed
  if (fAvgSpeed > 0) {
    fCoefVar = fSdSpeed / fAvgSpeed  # coefficient of variation
  } else {
    fCoefVar = 0
  }
  fNetDistance = distance(lastCoordinates, w[nrow(w), c("lat", "lon")])  # net distance
  d = distance(lastCoordinates, w[1, c("lat", "lon")])
  if (nrow(w) > 1){
    for (t in 1:(nrow(w) - 1)) {
      d = d + distance(w[t, c("lat", "lon")], w[t + 1, c("lat", "lon")])
    }
  }
  fTotalDistance = d # total distance
  if (fTotalDistance > 0) {
    fRatioDistance = fNetDistance / fTotalDistance
  } else {
    fRatioDistance = 0
  }
  fEleChange = sum(abs(w[, c("elevationDelta")]))  # total elevation change
  fAvgEle = mean(w[, c("ele")])  # average elevation
  fNsat = mean(w[, c("nsatView")])  # average nsatView
  fSNR = mean(w[, c("snrView")])  # average snr
  
  return(c(fAvgSpeed, fSdSpeed, fCoefVar, fNetDistance, fTotalDistance, fRatioDistance, fEleChange, fAvgEle, fNsat, fSNR))
}
distance = function(origin, destination) {
  lat1 = origin$lat
  lon1 = origin$lon
  lat2 = destination$lat
  lon2 = destination$lon
  radius = 6371 * 1000 # m
  
  dlat = (lat2 - lat1) * pi / 180
  dlon = (lon2 - lon1) * pi / 180
  a = sin(dlat / 2) * sin(dlat / 2) + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2) * sin(dlon / 2)
  c = 2 * atan2(sqrt(a), sqrt(1 - a))
  d = radius * c
  return(d)
}