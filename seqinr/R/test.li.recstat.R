##
# This function tests if a region located between two stop codons could be a putative CDS
#
# Data used are the factor scores of the CA computed on the windows by recstat function
##
#v.18.08.2011
test.li.recstat <- function(rec, fac = 1, length.min = 150, stop.max = 0.2, direct = TRUE, level = 0.05)
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
    vdep <- rec[[6]]
    vind <- rec[[7]]
    vstopd <- rec[[8]]
    vstopr <- rec[[9]]
    recd <- rec[[14]]
    recr <- rec[[15]]
    if (seqsize < length.min)
    {
        print("Seqence length is shorter than minimum distance between two Stop codons.")
        return()
    }
    table.recstat <- function(vstop, rec, frame)
    {
        tabCDS <- numeric() # initialization
        j <- 0
        for (i in 2:length(vstop))
        { # for each stop codons positions vector
            if ((vstop[i] - vstop[i - 1]) > length.min)
            { # test if space between codons is above the threshold
                # in each case gets the values between each stop codon for the 3 reading frames and range it in 3 vector seg
                seg1 <- rec$li[which((vstop[i - 1] - vdep[1:seqisize1])/sizewin <= stop.max &
                    (vstop[i] - vdep[1:seqisize1])/sizewin >= (1 - stop.max)), fac]
                seg2 <- rec$li[(which((vstop[i - 1] - vdep[(seqisize1 + 1):(seqisize1 + seqisize2)])/sizewin <= stop.max
                    & (vstop[i] - vdep[(seqisize1 + 1):(seqisize1 + seqisize2)])/sizewin >= (1 - stop.max)) + seqisize1), fac]
                seg3 <- rec$li[(which((vstop[i - 1] - vdep[(seqisize1 + seqisize2 + 1):(length(vdep))])/sizewin <= stop.max
                    & (vstop[i] - vdep[(seqisize1 + seqisize2 + 1):(length(vdep))])/sizewin >= (1 - stop.max)) + seqisize1 + seqisize2), fac]
                # create a table with calculation on those vectors seg then go to next space
                # inter-codon, each row correspond to a space inter-stop codon
                if (frame == 1)
                {
                    test1 <- t.test(seg1, seg2)$p.value
                    test2 <- t.test(seg1, seg3)$p.value
                }
                if (frame == 2)
                {
                    test1 <- t.test(seg2, seg1)$p.value
                    test2 <- t.test(seg2, seg3)$p.value
                }
                if (frame == 3)
                {
                    test1 <- t.test(seg3, seg1)$p.value
                    test2 <- t.test(seg3, seg2)$p.value
                }
                if (test1 < level & test2 < level)
                {
                    result <- 1
                }
                else
                {
                    result <- 0
                }
                tabCDS <- c(tabCDS, vstop[i-1]+3, vstop[i]+2, mean(seg1), mean(seg2), mean(seg3), test1, test2, result)
                j <- j + 1
            }
        }
        tabCDS <- matrix(tabCDS, nrow = j, ncol = 8, byrow = TRUE) # conversion list to table
        return(tabCDS)
    }
    seqisize <- floor((dim(recd$li)[1])/3) # number of window by reading frame, we take the integer part
    if ((dim(recd$li)[1])%%3 == 1) # adaptation of number of window between each reading frame
    {
        seqisize1 <- seqisize + 1 # for fr1
        seqisize2 <- seqisize # for fr2
    }
    if ((dim(recd$li)[1])%%3 == 2)
    {
        seqisize1 <- seqisize + 1
        seqisize2 <- seqisize + 1
    }
    if ((dim(recd$li)[1])%%3 == 0)
    {
        seqisize1 <- seqisize
        seqisize2 <- seqisize
    }
    ##
    ##direct strand##
    ##
    if (direct)
    { 
        vstopdindphase <- numeric()
        if (length(vstopd) > 0)
        { # test if vector is not empty because problem with modulo
            vstopdindphase <- sapply(1:length(vstopd), function(x) { # index vector of reading frame of vector vstopd
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
        tab1 <- table.recstat(vstop1, recd, 1)
        colnames(tab1) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(1,2)", "t(1,3)", "CDS")
        tab2 <- table.recstat(vstop2, recd, 2)
        colnames(tab2) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(2,1)", "t(2,3)", "CDS")
        tab3 <- table.recstat(vstop3, recd, 3)
        colnames(tab3) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(3,1)", "t(3,2)", "CDS")
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
        tab1 <- table.recstat(vstop1, recr, 1)
        colnames(tab1) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(1,2)", "t(1,3)", "CDS")
        tab2 <- table.recstat(vstop2, recr, 2)
        colnames(tab2) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(2,1)", "t(2,3)", "CDS")
        tab3 <- table.recstat(vstop3, recr, 3)
        colnames(tab3) <- c("Start", "End", "Mean 1", "Mean 2", "Mean 3", "t(3,1)", "t(3,2)", "CDS")
        return(list(tab1, tab2, tab3))
    }
}

