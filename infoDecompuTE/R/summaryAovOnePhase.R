summaryAovOnePhase <-
		function(design.df, blk.str, trt.str, var.comp = NA,
      blk.contr = NA, trt.contr = NA, contr.matrix = all(is.na(trt.contr)),  table.legend = FALSE,
      response = NA, latex = FALSE, fixed.names=NA){
    
  		#########################################################################################
#Main methods starts here->
#Extract the fixed and random terms

		rT = terms(as.formula(paste("~", blk.str, sep = "")), keep.order = TRUE) #random terms
		fT = terms(as.formula(paste("~", trt.str, sep = "")), keep.order = TRUE)  #fixed terms

		#########################################################################################
#Preparing the block structures

		Z = makeBlkDesMatrix(design.df, rev(attr(rT,"term.labels")))

		Pb <- makeBlockProjectors(Z)
		if(names(Pb)[1] == "e")
			names(Pb)[1] = "Within"

		#########################################################################################
#Prepating the treatment structures


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
#1-phase experiment

		#pre- and post-multiply NTginvATN by block projection matrices
		PNTginvATNP <- lapply(Pb,function(z) blkProkMat(z, T, N))

		#Now construct variance matrices
		PNTginvATNP<- PNTginvATNP[sort(1:length(PNTginvATNP), decreasing=TRUE)]

    v.mat = getVMat.onePhase(Z.Phase1 = Z, design.df = design.df, var.comp = var.comp)
  
    ANOVA = getVcCoefMS.onePhase(PNTginvATNP = PNTginvATNP, design.df = design.df, v.mat = v.mat, response = response, table.legend = table.legend)

		##############################################################################################################


		effFactors = lapply(Pb, function(z) trtProkMat(z, T, N, Rep))
		effFactors <- effFactors[sort(1:length(effFactors), decreasing=TRUE)]

	   EF = getFixEF.onePhase(effFactors = effFactors, trt.Coef = trt.Coef,  T = T, Rep = Rep, table.legend = table.legend)
     if(latex){
        return(toLatexTable(ANOVA = ANOVA, EF = EF, fixed.names = fixed.names) )
     } else{ 
		  return(list(ANOVA = ANOVA, EF = EF))
		}
	}
