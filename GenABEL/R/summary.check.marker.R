"summary.check.marker" <-
function(object,...) {
	out1 <- rep(NA,25)
	dim(out1) <- c(5,5)
	out1[1,2] <- cross(object$nocall,object$nofreq)	
	out1[1,3] <- cross(object$nocall,object$nohwe)	
	out1[1,4] <- cross(object$nocall,object$redundant)	
	out1[1,5] <- cross(object$nocall,object$Xmrkfail)	
	out1[1,1] <- length(object$nocall) - out1[1,2] - out1[1,3] - out1[1,4] - out1[1,5]
	out1[2,3] <- cross(object$nofreq,object$nohwe)	
	out1[2,4] <- cross(object$nofreq,object$redundant)	
	out1[2,5] <- cross(object$nofreq,object$Xmrkfail)	
	out1[2,2] <- length(object$nofreq) - out1[1,2] - out1[2,3] - out1[2,4] -out1[2,5]
	out1[3,4] <- cross(object$nohwe,object$redundant)	
	out1[3,5] <- cross(object$nohwe,object$Xmrkfail)	
	out1[3,3] <- length(object$nohwe) - out1[1,3] - out1[2,3] - out1[3,4] - out1[3,5]
	out1[4,5] <- cross(object$redundant,object$Xmrkfail)	
	out1[4,4] <- length(object$redundant) - out1[1,4] - out1[2,4] - out1[3,4] - out1[4,5]
	out1[5,5] <- length(object$Xmrkfail) - out1[1,5] - out1[2,5] - out1[3,5] - out1[4,5]
	rownames(out1) <- c("NoCall","NoMAF","NoHWE","Redundant","Xsnpfail")
	colnames(out1) <- rownames(out1)
	out2 <- rep(NA,49)
	dim(out2) <- c(7,7)
	out2[1,2] <- cross(object$idnocall,object$hetfail)
	out2[1,3] <- cross(object$idnocall,object$ibsfail)
	out2[1,4] <- cross(object$idnocall,object$isfemale)
	out2[1,5] <- cross(object$idnocall,object$ismale)
	out2[1,6] <- cross(object$idnocall,object$isXXY)
	out2[1,7] <- cross(object$idnocall,object$otherSexErr)
	out2[1,1] <- length(object$idnocall) - out2[1,2] - out2[1,3] - out2[1,4] - out2[1,4] - out2[1,5] - out2[1,6] - out2[1,7]
	out2[2,3] <- cross(object$hetfail,object$ibsfail)
	out2[2,4] <- cross(object$hetfail,object$isfemale)
	out2[2,5] <- cross(object$hetfail,object$ismale)
	out2[2,6] <- cross(object$hetfail,object$isXXY)
	out2[2,7] <- cross(object$hetfail,object$otherSexErr)
	out2[2,2] <- length(object$hetfail) - out2[1,2] - out2[2,3] - out2[2,4] - out2[2,4] - out2[2,5] - out2[2,6] - out2[2,7]
	out2[3,4] <- cross(object$ibsfail,object$isfemale)
	out2[3,5] <- cross(object$ibsfail,object$ismale)
	out2[3,6] <- cross(object$ibsfail,object$isXXY)
	out2[3,7] <- cross(object$ibsfail,object$otherSexErr)
	out2[3,3] <- length(object$ibsfail) - out2[1,3] - out2[2,3] - out2[3,4] - out2[3,5] - out2[3,6] - out2[3,7]
	out2[4,5] <- cross(object$isfemale,object$ismale)
	out2[4,6] <- cross(object$isfemale,object$isXXY)
	out2[4,7] <- cross(object$isfemale,object$otherSexErr)
	out2[4,4] <- length(object$isfemale) - out2[1,4] - out2[2,4] - out2[3,4] - out2[4,5] - out2[4,6] - out2[4,7]
	out2[5,6] <- cross(object$ismale,object$isXXY)
	out2[5,7] <- cross(object$ismale,object$otherSexErr)
	out2[5,5] <- length(object$ismale) - out2[1,5] - out2[2,5] - out2[3,5] - out2[4,5] - out2[5,6] - out2[5,7]
	out2[6,7] <- cross(object$isXXY,object$otherSexErr)
	out2[6,6] <- length(object$isXXY) - out2[1,6] - out2[2,6] - out2[3,6] - out2[4,6] - out2[5,6] - out2[6,7]
	out2[7,7] <- length(object$otherSexErr) - out2[1,7] - out2[2,7] - out2[3,7] - out2[4,7] - out2[5,7] - out2[6,7]
#	rownames(out2) <- c("IDnoCall","HetFail","IBSFail","Xidfail","isfemale","ismale")
#	rownames(out2) <- c("IDnoCall","HetFail","IBSFail","isfemale","ismale")
	rownames(out2) <- c("IDnoCall","HetFail","IBSFail","isfemale","ismale","isXXY","otherSexErr")
	colnames(out2) <- rownames(out2)
	out <- list()
	out$"Per-SNP fails statistics" <- out1
	out$"Per-person fails statistics" <- out2
	out
}

