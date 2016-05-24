# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

writesound = function (samples, filename = '', fs = 22050){
  if (class(samples) == "sound"){
    if (filename == '') filename = samples$filename
    fs = samples$fs
    samples = samples$sound
  }
  if (class(samples) == "ts") fs = frequency (samples)
  
  if (!is.numeric(samples)) stop("Non-numeric sample values given.")
  if (filename == '') filename = paste (deparse(substitute(samples)), '.wav', sep='')

  
  maxamp = max(abs(samples))
  sound = round((samples/maxamp) * 32767)
  soundfile = file(filename, "wb")
  on.exit(close(soundfile))
  samples = length(sound)
  writeChar("RIFF", soundfile, 4, eos = NULL)
  writeBin(as.integer(samples * 2 + 36), soundfile, size = 4,endian = "little")
  writeChar("WAVEfmt ", soundfile, 8, eos = NULL)
  writeBin(as.integer(16), soundfile, size = 4, endian = "little")
  writeBin(as.integer(1), soundfile, size = 2, endian = "little")
  writeBin(as.integer(1), soundfile, size = 2, endian = "little")
  writeBin(as.integer(fs), soundfile, size = 4, endian = "little")
  writeBin(as.integer(fs * 2), soundfile, size = 4, endian = "little")
  writeBin(as.integer(2), soundfile, size = 2, endian = "little")
  writeBin(as.integer(16), soundfile, size = 2, endian = "little")
  writeChar("data", soundfile, 4, eos = NULL)
  writeBin(as.integer(samples * 2), soundfile, size = 4, endian = "little")
  writeBin(as.integer(sound), soundfile, size = 2, endian = "little")
}

