##
# This function tests if a region located between two stop codons could be a putative CDS
#
# Data used are the factor scores of the CA computed on the windows by recstat function
##
#v.18.08.2011
test.co.recstat <- function(rec, fac = 1, length.min = 150, stop.max = 0.2, win.lim = 0.8, direct = TRUE, level = 0.01) 
{
    if (fac < 0 | 4 < fac) 
	{ # test if factor is between 1 and 4
        print("Factor number is not in 1:4.")
        return()
    }
    seq <- rec[[1]] # recovery of elements of list n
    sizewin <- rec[[2]]
    shift <- rec[[3]]
    seqsize <- rec[[4]]
    vstopd <- rec[[8]]
    vstopr <- rec[[9]]
    resd <- rec[[12]]
    resr <- rec[[13]]
    recd <- rec[[14]]
    recr <- rec[[15]]
	if (seqsize < length.min)
	{
		print("Seqence length is shorter than minimum distance between two Stop codons.")
		return()
	}
    table.recstat <- function(vstop) 
	{
        tabCDS <- numeric() # initialization
        j <- 0
        for (i in 2:length(vstop)) 
		{ # for each stop codons positions vector
            if ((vstop[i] - vstop[i - 1]) > length.min) 
			{ # test if space between codons is above the threshold
                # in each case gets the p-values between each stop codon and range it in vector seg
                seg <- pvalvec[which((vstop[i - 1] - pvalvec[, 2])/(sizewin + 2) <= stop.max &
                    (vstop[i] - pvalvec[, 2])/(sizewin + 2) >= (1 - stop.max)), 1]
                # create a table with calculation on those vectors seg then go to next space inter-codon, each row correspond to a space inter-stop codon
                k <- 0
                for (l in 1:length(seg)) 
				{
                	if (seg[l] < level) 
					{
                		k <- k + 1
                	}     
                	if (k/length(seg) > win.lim) 
					{
                    	result <- 1
                	} 
					else 
					{
                    	result <- 0
                	}
                }
                tabCDS <- c(tabCDS, vstop[i - 1]+3, vstop[i]+2, result)
                j <- j + 1
            }
        }
        tabCDS <- matrix(tabCDS, nrow = j, ncol = 3, byrow = TRUE) # conversion list to table
        return(tabCDS)
    }
    ##
    ##direct strand##
    ##
    if (direct) 
	{ 
    	vstopdindphase <- numeric()
        if (length(vstopd) > 0) 
		{ # test if vector is not empty because problem with modulo
            vstopdindphase <- sapply(1:length(vstopd), function(x)
            { # index vector of reading frame of vector vstopd
                if (vstopd[x]%%3 == 1) 
				{
                    vstopdindphase <- c(vstopdindphase, 1)
                } 
				else 
				{
                    if (vstopd[x]%%3 == 2) 
					{
                        vstopdindphase <- c(vstopdindphase, 2)
                    } 
					else 
					{
                        vstopdindphase <- c(vstopdindphase, 3)
                    }
                }
            })
        }
        vstop1 <- vstopd[vstopdindphase == 1] # vector with only stop codons in reading frame 1
        vstop2 <- vstopd[vstopdindphase == 2] # vector with only stop codons in reading frame 2
        vstop3 <- vstopd[vstopdindphase == 3] # vector with only stop codons in reading frame 3
        vstop1 <- c(vstop1, 1-3, seqsize-(seqsize%%3)-2) # add start and end positions, "-3" and "-2" because of table.recstat()
        vstop2 <- c(vstop2, 2-3, seqsize-((seqsize-1)%%3)-2)
        vstop3 <- c(vstop3, 3-3, seqsize-((seqsize-2)%%3)-2)
        vstop1 <- sort(unique(vstop1)) # sort of the vector
        vstop2 <- sort(unique(vstop2))
        vstop3 <- sort(unique(vstop3))
    	tab <- numeric()
    	for (i in 1:dim(resd)[1]) 
		{
    		if (sum(resd[i, ]) == sizewin/3) 
			{
    			for (j in 1:64) 
				{
    				tab <- c(tab, rep(recd$co[j, 1], resd[i, j]))
    			}
    		} 
    	}
    	tab <- matrix(unlist(tab), byrow = TRUE, ncol = sizewin/3)
    	seqisize <- floor((dim(tab)[1])/3)
    	phase <- c(rep(1, sizewin/3), rep(2, sizewin/3), rep(3, sizewin/3))
    	phase <- as.factor(phase)
    	pvalvec <- numeric()
    	for (i in 1:seqisize) 
		{
    		v1 <- tab[i,]
    		v2 <- tab[seqisize + i,]
    		v3 <- tab[2*seqisize + i,]
    		v <- c(v1, v2, v3)
    		x <- kruskal.test(v ~ phase)$p.value
    		pvalvec <- c(pvalvec, x, (i - 1)*shift + 1)
    	}
    	pvalvec <- matrix(unlist(pvalvec), byrow = TRUE, ncol = 2)
#    	plot((sizewin/2) + (0:(seqisize - 1))*shift + 1, pvalvec[,1], type = "l")
    	tab1 <- table.recstat(vstop1)
        colnames(tab1) <- c("Start", "End", "CDS") # definition of table headings
        tab2 <- table.recstat(vstop2)
        colnames(tab2) <- c("Start", "End", "CDS")
        tab3 <- table.recstat(vstop3)
        colnames(tab3) <- c("Start", "End", "CDS")
        return(list(tab1, tab2, tab3))
	}
	##
    ##reverse strand##
    ##
    if (!direct) 
	{ 
    	vstoprindphase <- numeric()
        if (length(vstopr) > 0) 
		{ 
            vstoprindphase <- sapply(1:length(vstopr), function(x)
            { 
                if (vstopr[x]%%3 == 1) 
				{
                    vstoprindphase <- c(vstoprindphase, 1)
                } 
				else 
				{
                    if (vstopr[x]%%3 == 2) 
					{
                        vstoprindphase <- c(vstoprindphase, 2)
                    } 
					else 
					{
                        vstoprindphase <- c(vstoprindphase, 3)
                    }
                }
            })
        }
        vstop1 <- vstopr[vstoprindphase == 1] 
        vstop2 <- vstopr[vstoprindphase == 2] 
        vstop3 <- vstopr[vstoprindphase == 3] 
        vstop1 <- c(vstop1, 1-3, seqsize-(seqsize%%3)-2) # add start and end positions
        vstop2 <- c(vstop2, 2-3, seqsize-((seqsize-1)%%3)-2)
        vstop3 <- c(vstop3, 3-3, seqsize-((seqsize-2)%%3)-2)
        vstop1 <- sort(unique(vstop1)) 
        vstop2 <- sort(unique(vstop2))
        vstop3 <- sort(unique(vstop3))
        tab <- numeric()
    	for (i in 1:dim(resr)[1]) 
		{
    		if (sum(resr[i, ]) == sizewin/3) 
			{
    			for (j in 1:64) 
				{
    				tab <- c(tab, rep(recr$co[j, 1], resr[i, j]))
    			}
    		} 
    	}
    	tab <- matrix(unlist(tab), byrow = TRUE, ncol = sizewin/3)
    	seqisize <- floor((dim(tab)[1])/3)
    	phase <- c(rep(1, sizewin/3), rep(2, sizewin/3), rep(3, sizewin/3))
    	phase <- as.factor(phase)
    	pvalvec <- numeric()
    	for (i in 1:seqisize) 
		{
    		v1 <- tab[i,]
    		v2 <- tab[seqisize + i,]
    		v3 <- tab[2*seqisize + i,]
    		v <- c(v1, v2, v3)
    		x <- kruskal.test(v ~ phase)$p.value
    		pvalvec <- c(pvalvec, x, (i - 1)*shift + 1)
    	}
    	pvalvec <- matrix(unlist(pvalvec), byrow = TRUE, ncol = 2)
		tab1 <- table.recstat(vstop1)
        colnames(tab1) <- c("Start", "End", "CDS") # definition of table headings
        tab2 <- table.recstat(vstop2)
        colnames(tab2) <- c("Start", "End", "CDS")
        tab3 <- table.recstat(vstop3)
        colnames(tab3) <- c("Start", "End", "CDS")
        return(list(tab1, tab2, tab3))
    }
}