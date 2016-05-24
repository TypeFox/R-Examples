`growth.fnc` <-
function(text = alice, size = 646, nchunks = 40, chunks = 0) {

  
  if (is.numeric(chunks) & (chunks[1] == 0)) {
    if (length(size) > 1) stop("size should be a single number")
    if (size == 1) stop("size should be greater than one")
    if (length(nchunks) > 1) stop("nchunks should be a single number")
    chunks = seq(size, nchunks*size, by = size)
  } else {
    nchunks = length(chunks)
    if (max(chunks) > length(text)) {
      stop("chunk vector exceeds text size\n")
    }
  }

  types = rep(0, nchunks)
  tokens = rep(0, nchunks)
  hapax = rep(0, nchunks)
  tris = rep(0, nchunks)
  dis = rep(0, nchunks)
  yule = rep(0, nchunks)
  zipf = rep(0, nchunks)
  herdan = rep(0, nchunks)
  lognormal = rep(0, nchunks)

  for (i in 1 : nchunks) { 
     spectrum = spectrum.fnc(text[1:chunks[i]])  
		 if (nrow(spectrum) < 3) 
			 stop(paste("too few spectrum elements for chunk ", i, "\n", sep=""))
     types[i] = sum(spectrum$freqOfFreq)
     tokens[i] = chunks[i]
		 if (length(spectrum[spectrum$frequency == 1,]$freqOfFreq) == 0) {
        hapax[i] = 0
		 } else {
        hapax[i] = spectrum[spectrum$frequency == 1,]$freqOfFreq
		 }
		 if (length(spectrum[spectrum$frequency == 2,]$freqOfFreq) == 0) {
        dis[i] = 0
		 } else {
        dis[i] = spectrum[spectrum$frequency == 2,]$freqOfFreq
		 }
		 if (length(spectrum[spectrum$frequency == 3,]$freqOfFreq) == 0) {
        tris[i] = 0
		 } else {
        tris[i] = spectrum[spectrum$frequency == 3,]$freqOfFreq
		 }
     yule[i] = yule.fnc(spectrum) 
     z = zipf.fnc(text[1:chunks[i]])
     herdan[i] = herdan.fnc(text[1:chunks[i]], 
        cumsum(rep(floor(length(text[1:chunks[i]])/40), 40)))$C
     zipf[i] = coef(lm(log(z$frequency) ~ log(z$rank)))[2]
     lognormal[i] = mean(log(table(text[1:chunks[i]])))
     cat(".")
  }
  cat("\n")

  res = data.frame(Chunk = 1:nchunks, 
    Tokens = tokens, Types = types, HapaxLegomena = hapax, 
    DisLegomena = dis, TrisLegomena = tris, Yule = yule, Zipf = zipf)
  res$TypeTokenRatio = types/tokens
  res$Herdan = herdan
  res$Guiraud = types/sqrt(tokens)
  res$Sichel = dis/types
  res$Zipf = zipf
  res$Lognormal = lognormal

  
  return(growthInit(list(data = res)))

}

