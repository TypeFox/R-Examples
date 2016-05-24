# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


loadsound = function (filename=0){
  if (mode(filename)!="character") filename = file.choose()

  soundfile = file(filename,"rb")
  readChar(soundfile, nchars = 8)        ## ChunkId and ChunkSize (4,4)
  if(readChar(soundfile, nchars = 4) != 'WAVE'){
    close(soundfile)
    stop ("File provided is not in .wav format.")
  }

  readBin(soundfile, "integer", n = 10, size = 1)  ## Subchunk1ID, Subchunk1Size, AudioFormat (4,4,2)
  
  numChannels = readBin(soundfile, "integer", n = 1, size = 2)
  if (numChannels > 1){
   stop ("This function only loads mono (single-channel) WAV files.")
  }
  sampleRate = readBin(soundfile,"integer", n = 1, size = 4)

  readBin(soundfile,"integer",n= 6,size=1)   ## ByteRate, BlockAlign (4,2)         

  bitsPerSample = readBin(soundfile, "integer", size = 2)
  if (bitsPerSample > 16) stop ("This function only loads 8 and 16 bit WAV files.")

  readBin(soundfile,"integer",n= 4,size=1)  ## Subchunk2ID
  subchunk2Size = readBin(soundfile,"integer", size=4)

  if (bitsPerSample == 8)
      sound = (readBin (soundfile, "integer" , n = subchunk2Size, size = 1, signed = FALSE) - 128) / 255
  if (bitsPerSample == 16)
      sound = readBin(soundfile,"integer", n = subchunk2Size/2, size=2, signed = TRUE) / 32768
 
  close(soundfile)

  numSamples = subchunk2Size / (bitsPerSample/8)
  fs = sampleRate

  output = list (filename = filename, fs = fs, numSamples = numSamples,
  duration = numSamples/fs * 1000, sound = ts(sound, frequency = fs, start=0))
  class(output) = "sound"
  output  
}
