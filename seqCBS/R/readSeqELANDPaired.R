readSeqELANDPaired <-
function(filename) {
	seqRaw = read.table(filename, as.is=TRUE)
	seqMatch = seqRaw[seqRaw[,2]==seqRaw[,6],]
	rm(seqRaw)
	seqChr = factor(seqMatch[,2])
	seqFront = seqMatch[,4]
	normalLFflag = (seqFront=="F")
	seqPosF = rep(0, dim(seqMatch)[1])
	seqPosF[normalLFflag] = as.numeric(seqMatch[normalLFflag,3])
	seqPosF[-normalLFflag] = as.numeric(seqMatch[-normalLFflag,7])
	seqPosR = rep(0, dim(seqMatch)[1])
	seqPosR[normalLFflag] = as.numeric(seqMatch[normalLFflag,7])
	seqPosR[-normalLFflag] = as.numeric(seqMatch[-normalLFflag,3])
	rm(seqMatch)
	return(list(seqF=seqPosF, seqR=seqPosR, seqChr=seqChr))
}

