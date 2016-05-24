## functions to compute accelerometer features

computeOneAccFeat = function(w, Fs) {
  # function to compute one accelerometer feature vector from a window
  # axis 1: vertical (z)
  # axis 2: horizontal (y)
  # axis 3: perpindicular (x)
  
  # gravity component
#   g = w[1, ]
#   for (n in 1:nrow(w)) {
#     g = 0.9 * g + 0.1 * w[n, ]
#   }
  
  g = matrix(0, nrow(w), 3)
  x = 0.9
  g[1, ] = (1-x) * w[1, ]
  for (n in 2:nrow(w)) {
    g[n, ] = x * g[n-1] + (1-x) * w[n, ]
  }
  g = g[Fs:nrow(g), ] # ignore beggining
  gg = colMeans(g)
  w = w - gg
  
  # v = vector magnitude
  v = sqrt(rowSums(w ^ 2))
  fMean = mean(v)
  fStd = sd(v)
  if (fMean > 0) {
    fCoefVariation = fStd / fMean
  } else {
    fCoefVariation = 0
  }
  fMedian = median(v)
  fMin = min(v)
  fMax = max(v)
  f25thP = quantile(v, 0.25)[[1]]
  f75thP = quantile(v, 0.75)[[1]]
  
  a = acf(v, plot=FALSE)
  fAutocorr = which.max(abs(a$acf[2:length(a$acf)])) / (nrow(w) / Fs)
  if ((sd(w[, 3]) > 0) & (sd(w[, 2]) > 0)) {
    fCorrxy = cor(w[, 3], w[, 2])
  } else {
    fCorrxy = 0
  }
  if ((sd(w[, 3]) > 0) & (sd(w[, 1]) > 0)) {
    fCorrxz = cor(w[, 3], w[, 1])
  } else {
    fCorrxz = 0
  }
  if ((sd(w[, 2]) > 0) & (sd(w[, 1]) > 0)) {
    fCorryz = cor(w[, 2], w[, 1])
  } else {
    fCorryz = 0
  }
  
  if (is.na(fCorrxy)) fCorrxy = 0
  if (is.na(fCorrxz)) fCorrxz = 0
  if (is.na(fCorryz)) fCorryz = 0
  
  fAvgRoll = mean(atan2(w[, 2],w[, 1]))
  fAvgPitch = mean(atan2(w[, 1],w[, 3]))
  fAvgYaw = mean(atan2(w[, 2],w[, 3]))
  
  fSdRoll = sd(atan2(w[, 2],w[, 1]))
  fSdPitch = sd(atan2(w[, 1],w[, 3]))
  fSdYaw = sd(atan2(w[, 2],w[, 3]))
  
  fRollG = atan2(gg[2], gg[1])
  fPitchG = atan2(gg[1], gg[3])
  fYawG = atan2(gg[2], gg[3])
  
  s = specgram(v, n=length(v), Fs=Fs)
  S = abs(s$S)
  f = S / max(S)
  freq = s$f
  f1 = f[freq >= 0.1]
  freq1 = freq[freq >= 0.1]
  fFmax = freq1[which.max(f1)]
  fPmax = max(f1)
  
  band = f[freq > 0.3 & freq < 3]
  fPmaxBand = max(band)
  freqband = freq[freq > 0.3 & freq < 3]
  fFmaxBand = freqband[which.max(band)]
  fEntropy = - sum(f * log(f))
  
  s = specgram(v, n=Fs, Fs=Fs)
  S = abs(s$S)
  f = S / max(S)
  freq = s$f
  f = rowSums(f) / ncol(f)
  FFT0 = f[1]
  FFT1 = f[2]
  FFT2 = f[3]
  FFT3 = f[4]
  FFT4 = f[5]
  FFT5 = f[6]
  FFT6 = f[7]
  FFT7 = f[8]
  FFT8 = f[9]
  FFT9 = f[10]
  FFT10 = f[11]
  FFT11 = f[12]
  FFT12 = f[13]
  FFT13 = f[14]
  FFT14 = f[15]
  
  return(c(fMean, fStd, fCoefVariation, fMedian, fMin, fMax, f25thP, f75thP, fAutocorr, fCorrxy, fCorrxz, fCorryz, fAvgRoll, fAvgPitch, fAvgYaw, fSdRoll, fSdPitch, fSdYaw, fRollG, fPitchG, fYawG, fFmax, fPmax, fFmaxBand, fPmaxBand, fEntropy, FFT0, FFT1, FFT2, FFT3, FFT4, FFT5, FFT6, FFT7, FFT8, FFT9, FFT10, FFT11, FFT12, FFT13, FFT14))
}
extractAccFeatsFile = function(inputFile, outputPath, winSize) {
  # function to extract accelerometer features from raw actigraph file
  
  con = file(inputFile, open = "r")
  line = readLines(con, n = 1)
  Fs = as.numeric(str_match(line, "(\\d+) Hz")[1, 2])
  dateFmt = str_match(line, "date format ([a-z,A-Z,/]*)")[1, 2]
  dateFmt = gsub("yyyy", "%Y", dateFmt)
  dateFmt = gsub("M", "%m", dateFmt)
  dateFmt = gsub("d", "%d", dateFmt)
  line = readLines(con, n = 1)
  line = readLines(con, n = 1)
  StartTime = gsub("Start Time ", "", line)
  line = readLines(con, n = 1)
  StartDate = gsub("Start Date ", "", line)
  line = readLines(con, n = 6)
  st = strptime(paste(StartDate, StartTime), paste(dateFmt, "%H:%M:%S"))
  day = st$mday
  out = file.path(outputPath,strftime(st, "%Y-%m-%d"))
  cat(strftime(st, "%Y-%m-%d"), '\n')
  if (!file.exists(outputPath)) {
    dir.create(outputPath, recursive=TRUE)
  }
  cat("timestamp,mean,sd,coefvariation,median,min,max,25thp,75thp,autocorr,corrxy,corrxz,corryz,avgroll,avgpitch,avgyaw,sdroll,sdpitch,sdyaw,rollg,pitchg,yawg,fmax,pmax,fmaxband,pmaxband,entropy,fft0,fft1,fft2,fft3,fft4,fft5,fft6,fft7,fft8,fft9,fft10,fft11,fft12,fft13,fft14\n", file=out, append=TRUE)
  
  while (length(line <- readLines(con, n = Fs * winSize)) >= Fs * winSize) {
    line = gsub("\"", "", line)
    M = as.matrix(strsplit(line, " "))
    M = sapply(M, strsplit, ",")
    M = sapply(M, as.numeric)
    M = t(M)
    feat = computeOneAccFeat(M, Fs)
    cat(strftime(st, "%Y-%m-%d %H:%M:%S,"), file=out, sep = "", append=TRUE)
    cat(feat, file=out, sep=",", append=TRUE)
    cat('\n', file=out, append=TRUE)
    
    st = as.POSIXlt(st + winSize)
    if (st$mday != day) {
      out = file.path(outputPath,strftime(st, "%Y-%m-%d"))
      cat(strftime(st, "%Y-%m-%d"), '\n')
      cat("timestamp,mean,sd,coefvariation,median,min,max,25thp,75thp,autocorr,corrxy,corrxz,corryz,avgroll,avgpitch,avgyaw,sdroll,sdpitch,sdyaw,rollg,pitchg,yawg,fmax,pmax,fmaxband,pmaxband,entropy,fft0,fft1,fft2,fft3,fft4,fft5,fft6,fft7,fft8,fft9,fft10,fft11,fft12,fft13,fft14\n", file=out, append=TRUE)
      day = st$mday
    }
  }
  close(con)
}
extractAccelerometerFeatures = function(input, output, winSize, names = NULL) {
  # function to compute accelerometer features for all actigraph files in a directory
  if (file.info(input)$isdir){
    if (is.null(names)) {
      names = file_path_sans_ext(list.files(input))
    }
    if (length(names) == 0) {
      stop("couldn't find any accelerometer files\n")
    }
    for (name in names) {
      outputFile = file.path(output, name)
      if (!file.exists(outputFile)) {
        cat(name, "...\n")
        extractAccFeatsFile(file.path(input, paste0(name, ".csv")), outputFile, winSize)
      }
    }
  } else {
    extractAccFeatsFile(input, output, winSize)
  }
}
