parameter.estimate <-
function(sp, method, alpha=NULL, nsnap=1)
  {
    data1 <- sp$nodes.characteristics
    dist1 <- sp$distance.to.neighbours
    A <- as.vector(data1[,"areas"])
    p <- as.vector(data1[,"species"])
    d <- dist(data1[,c("x","y")])
	alpha <- 1/sp$dispersal
if(method == "Rsnap_1")
	{
        d2 <- as.matrix(exp(-alpha*d))
        diag(d2) <- 0
        d2 <- sweep(d2,2,A,"*")
        if(sum(p) == 1) S <- as.vector((d2[,p>0]))
        if(sum(p)>1) S <- as.vector(rowSums(d2[,p>0]))
        model <- glm(p ~ offset(2*log(S+0.001))+log(A), family = binomial)
        beta1 <-coef(model)
        x <- as.numeric(beta1[2])
        A0 <- min(data1[p > 0,3])
        ey <- exp(as.numeric(-beta1[1]))
        e <- A0^x
        y <- ey/e
        par_output <- c(alpha,x,y,e)
        names(par_output) <- c("alpha","x","y","e")
}
if(method == "Rsnap_x")
    {
        P <- rowSums(as.matrix(data1[, 9:(9+(nsnap-1))]))
        d2 <- as.matrix(exp(-alpha*d))
        diag(d2) <- 0
        d2 <- sweep(d2, 2, A, "*")
        S <- as.vector(rowSums(sweep(d2, 2, P/nsnap, "*")))
        mod2 <- glm(cbind(P, nsnap - P) ~ offset(2 * log(S)) + log(A), family = binomial)
        beta1 <- as.numeric(coef(mod2))
        x <- beta1[2]
        ey <- exp(-beta1[1])
        e <- (min(data1[p > 0, 4]))^x
        y <- ey/e
        par_output <- c(alpha, x, y, e)
        names(par_output) <- c("alpha", "x", "y", "e")
}
if(method == "MCsim")
	{
            x_col <- round(data1[, 1])
            y_col <- round(data1[, 2]) 
            area_col <- round(data1[, 3]*10000) 
            dummy_col <- rep(1, nrow(data1)) 
            pnc_col <- rep(1, nrow(data1))
            SIN_col <- rep(1, nrow(data1)) 
            nr_col <- data1[, 8] 
            occ_cols <- data1[, 9:ncol(data1)]
            df_mcsim <- cbind(x_col, y_col, area_col, dummy_col, pnc_col, SIN_col,
                              nr_col, occ_cols)
            write.table(df_mcsim, "inputMCsim.dat", sep='\t', quote=FALSE, row.names=F, col.names=F)
            settings_text <- "Iterations    = 1\nEstimations   = 1\nTime steps    = 1000\nTime to stab. = 100\nInput file    = inputMC.dat\nOutput file   = outfile\nSnapshots     = xx1xx  //max = 10\nIntervals     = 1 1 1 1 1 1 1 1 1\nEstimate xye  = yes\nEstimate alfa = yes\nEstimate A0   = no\nEstimate stoc.= no\nReplicates in est. = 250\nF.evals. in init.  = 300\nF.evals. in est.   = 900\nTurnover max. speed= 0\n\nA0   = xx7xx\nRemote col. prob.  = 0\nEnvironmental var. = 0\nAlpha = xx10xx\nA     = 1\nB     = 0.5\nx     = xx13xx\nY     = xx14xx\ne     = xx15xx \n\nPN class in use = 2\nPNs to iterate  = 1\nPN list         = 43\nMin. PN size    = 5\nMax. PN size    = 300\n\nOpt. method;Mc/Tmc/Nlr/Bnlr = Nlr\nKernel min. match prop.   = 0.75\nShow graph                = No\nShow J / Error            = Err J\nDisplay                   = SVGA\nStep                      = No\nReport time series        = Yes\nReport individual patches = No\nReport p-value distr.     = No\nShow con
 f. 
limits 
        = No\nMinimum big patch area    = 0\n\nUse report area             = No\nArea: left, up, right, down = 0 140 300 0\n"
            settings_text <- sub("xx1xx", as.character(nsnap), settings_text)
			data2 <- data1
			data2[,ncol(data1)+1] <- as.vector(rowSums(data2[,9:ncol(data1)],na.rm=TRUE))
			A0 <- min(data2[data2[,ncol(data1)+1] > 0, ]$area)
            settings_text <- sub("xx7xx", A0, settings_text)
            settings_text <- sub("xx10xx", alpha, settings_text)
            settings_text <- sub("xx13xx", "edit x", settings_text)
            settings_text <- sub("xx14xx", "edit y", settings_text)
            settings_text <- sub("xx15xx", "edit e", settings_text)
            write(settings_text, file = "inputMCsim.set")
            return(cat(paste("You can find data and paremeter files in ", getwd(),
                               "\nPlease, check input files before running mc application (Moilanen, 1999).
							   IMPORTANT: Please follow the recomendations in the readme file available with 
							   the application. Namely the three step (Nlr-Bnlr-Mc) procedure for parameter 
							   estimation. You should change the 'Opt. method;Mc/Tmc/Nlr/Bnlr' in the produced settings file.\n")))
	}
if(method == "rescue")
    {
            x_col <- round(data1[, 1],3)
            y_col <- round(data1[, 2],3)
            area_col <- round(data1[, 3],3)
            parE <- rep(1, nrow(data1))
            parC <- rep(1, nrow(data1))
            occ_cols <- data1[, 9:ncol(data1)]
            df_mcsim <- cbind(x_col, y_col, area_col, parE, parC,occ_cols)
            npatches.line <- paste("npatch ", nrow(df_mcsim), sep="")
            write(npatches.line,file="input_rescue.dat")
            nyear.line <- paste("nyear ", nsnap, sep="")
            write(nyear.line,file="input_rescue.dat",append=TRUE)
            write.table(df_mcsim, "input_rescue.dat", sep='\t', quote=FALSE, row.names=F, col.names=F, append=TRUE)
            parameters_text <- "nhpar 0\nnmpar 6\nnxparE 0\nnxparC 0\nWith_area_effect 0\ne    xx15xx   .1\nx    xx13xx   .1\ny    xx14xx   .1\nz    1.835    2\nalpha  xx10xx .1\nb    0.1      .1\nread_Dij     1\nicov 0\nscale 0\nUse_J_model  1\nL      0 1\nstartyear 1  xx16xx\nseed   123456\nniter  1000 100\nnsubiter 1 1\nlag 1\nbounds_e     1.0e-6 100\nbounds_x     1.0e-6 100\nbounds_y     1.0e-6 100\nbounds_z     1.0e-6 100\nbounds_alpha 1.0e-6 100\nbounds_b     1.0e-6 100\nprior_trans_e     0 100\nprior_trans_x     1 100\nprior_trans_y     0 100\nprior_trans_z     1 100\nprior_trans_alpha  .5 100\nprior_trans_b      .5 1"
            parameters_text <- sub("xx10xx", round(alpha,3), parameters_text)
            parameters_text <- sub("xx13xx", "edit_x", parameters_text)
            parameters_text <- sub("xx14xx", "edit_y", parameters_text)
            parameters_text <- sub("xx15xx", "edit_e", parameters_text)
			parameters_text <- sub("xx16xx", nsnap, parameters_text)
            write(parameters_text, file = "input_rescue.par")
			rl1 <- remove.species(sp)
			distance <- matrix.graph(rl1,"euc_distance")
			ID1 <- rl1$nodes.characteristics$ID
			comb1 <- as.matrix(combn(ID1, 2))
			comp1 <- length(comb1)/2
			distance2 <- as.data.frame(matrix(NA,comp1,2))
			for(i in 1:nrow(distance2))
				{
				distance2[i,1] <- comb1[,i][1]
				distance2[i,2] <- comb1[,i][2]
				}
			distance2 <- cbind(distance2, 0, NA)
			for(i in 1:nrow(distance2))
				{
				a <- as.character(distance2[i,1])
				b <- as.character(distance2[i,2])
				distance2[i,4] <- distance[a,b]
				}
			distance2[,4] <- round(distance2[,4]/1000,3)
			distance2[nrow(distance2),] <- 0
			write.table(distance2, "input_rescue.dis", sep='\t', quote=FALSE, row.names=F, col.names=F, append=TRUE)
			return(cat(paste("You can find data, parameter and distance (no barriers considered) files in ", getwd(),
                             "\nPlease, check input files before running the application (Ter Braak & Etienne, 2003)\n")))
}
if(method == "norescue")
    {
            x_col <- round(data1[, 1],3)
            y_col <- round(data1[, 2],3)
            area_col <- round(data1[, 3],3)
            dummy_col <- rep(1, nrow(data1))
            pnc_col <- rep(1, nrow(data1))
            SIN_col <- rep(1, nrow(data1))
            nr_col <- data1[, 8]
            occ_cols <- data1[, 9:ncol(data1)]
            df_mcsim <- cbind(x_col, y_col, area_col, dummy_col, pnc_col, SIN_col,
                              nr_col, occ_cols)
            npatches.line <- paste("npatch ", nrow(df_mcsim), sep="")
            write(npatches.line,file="input_norescue.dat")
            nyear.line <- paste("nyear ", nsnap, sep="")
            write(nyear.line,file="input_norescue.dat",append=TRUE)
            write.table(df_mcsim, "input_norescue.dat", sep='\t', quote=FALSE, row.names=F, col.names=F, append=TRUE)
            parameters_text <- "nhpar 0\nnmpar 6\nnxparE 0\nnxparC 0\nWith_area_effect 0\ne    xx15xx   .1\nx    xx13xx   .1\ny    xx14xx   .1\nz    1.835    2\nalpha  xx10xx .1\nb    0.5      .1\nread_Dij     1\nicov 0\nscale 0\nUse_J_model  1\nL      0 1\nstartyear 1  xx16xx\nseed   123456\nniter  1000 100\nnsubiter 1 1\nlag 1\nbounds_e     1.0e-6 100\nbounds_x     1.0e-6 100\nbounds_y     1.0e-6 100\nbounds_z     1.0e-6 100\nbounds_alpha 1.0e-6 100\nbounds_b     1.0e-6 100\nprior_trans_e     0 100\nprior_trans_x     1 100\nprior_trans_y     0 100\nprior_trans_z     1 100\nprior_trans_alpha  .5 100\nprior_trans_b      .5 1"
            parameters_text <- sub("xx10xx", round(alpha,3), parameters_text)
            parameters_text <- sub("xx13xx", "edit x", parameters_text)
            parameters_text <- sub("xx14xx", "edit y", parameters_text)
            parameters_text <- sub("xx15xx", "edit e", parameters_text)
			parameters_text <- sub("xx16xx", nsnap, parameters_text)
            write(parameters_text, file = "input_norescue.par")
			rl1 <- remove.species(sp)
			distance <- matrix.graph(rl1,"euc_distance")
			ID1 <- rl1$nodes.characteristics$ID
			comb1 <- as.matrix(combn(ID1, 2))
			comp1 <- length(comb1)/2
			distance2 <- as.data.frame(matrix(NA,comp1,2))
			for(i in 1:nrow(distance2))
				{
				distance2[i,1] <- comb1[,i][1]
				distance2[i,2] <- comb1[,i][2]
				}
			distance2 <- cbind(distance2, 0, NA)
			for(i in 1:nrow(distance2))
				{
				a <- as.character(distance2[i,1])
				b <- as.character(distance2[i,2])
				distance2[i,4] <- distance[a,b]
				}
			distance2[,4] <- round(distance2[,4]/1000,3)
			distance2[nrow(distance2),] <- 0
			write.table(distance2, "input_norescue.dis", sep='\t', quote=FALSE, row.names=F, col.names=F, append=TRUE)
			return(cat(paste("You can find data, parameter and distance (no barriers considered) files in ", getwd(),
                             "\nPlease, check input files before running the application (Ter Braak & Etienne, 2003)\n")))
}
    return(as.data.frame(par_output))
  }
