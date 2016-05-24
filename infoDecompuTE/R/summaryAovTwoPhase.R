summaryAovTwoPhase <-
		function(design.df, blk.str1, blk.str2, trt.str, var.comp = NA,
      blk.contr = NA, trt.contr = NA, contr.matrix = all(is.na(trt.contr)),  table.legend = FALSE,
      response = NA, latex = FALSE, fixed.names=NA){
    
      #Extract the fixed and random terms

		rT1 = terms(as.formula(paste("~", blk.str1, sep = "")), keep.order = TRUE) #random terms phase 1
		rT2 = terms(as.formula(paste("~", blk.str2, sep = "")), keep.order = TRUE) #random terms phase 2
		fT = terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

		#########################################################################################
#Preparing the block structures
		#write("1. Preparing the block structure.", "")

		blkTerm1 = attr(rT1,"term.labels")
		blkTerm2 = attr(rT2,"term.labels")

		Z1 = makeBlkDesMatrix(design.df, rev(blkTerm1))
		Z2 = makeBlkDesMatrix(design.df, rev(blkTerm2))

		#write("2. Defining the block structures of second Phase.", "")
		Pb <- makeBlockProjectors(Z2)
		if(names(Pb)[1] == "e")
			names(Pb)[1] = "Within"

		#write("3. Defining the block structures of first phase within second Phase.", "")
		trtTerm = attr(rT1,"term.labels")
		effectsMatrix = attr(rT1,"factor")

		T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix, trt.contr = blk.contr, 
      contr.matrix = all(is.na(blk.contr)))
		N =  getIncidenceMatrix(design.df, trtTerm)


		Pb1 <- lapply(Pb,function(z) blkProkMat(z, T, N))

		#########################################################################################
#Preparing the treatment structures
		#write("4. Preparing the treatment structure.", "")

		trtTerm = attr(fT,"term.labels")
		effectsMatrix = attr(fT,"factor")

		T =  makeTreatProjectors(design.df, trtTerm, effectsMatrix, trt.contr, contr.matrix)
		N =  getIncidenceMatrix(design.df, trtTerm)
		Rep = getReplicationList(design.df, trtTerm)
		trt.Coef = getTrtCoef(design.df, trtTerm)

		if(any(grepl("\\.", names(T)))){
			colnames(Rep) = trtTerm
			names(trt.Coef) = trtTerm

			Rep = Rep[,sapply(strsplit(names(T), "\\."), function(x) x[1])]
			trt.Coef = trt.Coef[sapply(strsplit(names(T), "\\."), function(x) x[1])]
		}

		colnames(Rep) = names(T)
		names(trt.Coef) = names(T)

		#########################################################################################
#Start calculating the VCs
#2-phase experiment
		#write("5. Start calculating the variance components.", "")
		#write("6. Pre- and post-multiply NTginvATN by block projection matrices.", "")


		PNTginvATNP <- lapply(Pb1, function(y) lapply(y, function(z) blkProkMat(z, T, N)))

		PNTginvATNP<- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing=TRUE)]


     v.mat = getVMat.twoPhase(Z.Phase1 = Z1, Z.Phase2 = Z2, design.df = design.df, var.comp = var.comp)
     
     ANOVA = getVcCoefMS.twoPhase(PNTginvATNP = PNTginvATNP, design.df = design.df, v.mat = v.mat, response, table.legend)
     
     
     effFactors = lapply(Pb1, function(y) lapply(y, function(z) trtProkMat(z, T, N, Rep)))
		effFactors <- effFactors[sort(1:length(effFactors), decreasing=TRUE)]

     
     EF = getFixEF.twoPhase(effFactors = effFactors, trt.Coef = trt.Coef,  T = T, Rep = Rep, table.legend)
     
      if(latex){
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names) )
     } else{ 
		  return(list(ANOVA = ANOVA, EF = EF))
		}

     
     }
