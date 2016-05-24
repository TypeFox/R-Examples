##
# This function draws two graphics, one of the CA of a DNA sequence, and one of
# start/stop codons positions in the three reading frames. This for the direct or 
# the reverse strand.
##
#v.18.08.2011
draw.recstat <- function(rec, fac = 1, direct = TRUE, xlim = c(1, seqsize), col = c("red", "blue", "purple"))
{
    if (fac < 0 | 4 < fac)
    { # test if factor is between 1 and 4
        print("Factor number is not in 1:4.")
        return()
    }
    seq <- rec[[1]] # recovery of elements of list rec
    sizewin <- rec[[2]]
    shift <- rec[[3]]
    seqsize <- rec[[4]]
    seqname <- rec[[5]]
    vstopd <- rec[[8]]
    vstopr <- rec[[9]]
    vinitd <- rec[[10]]
    vinitr <- rec[[11]]
    recd <- rec[[14]]
    recr <- rec[[15]]
    if (xlim[1] < 1 | xlim[1] > seqsize)
    {
        xlim <- c(1, xlim[2])
    }
    if (seqsize < xlim[2] | 1 > xlim[2])
    {
        xlim <- c(xlim[1], seqsize)
    }
    par(mfrow = c(2, 1), mar = c(0, 4, 4, 2) + 0.1) # division of the window for a closer between plots
    par(xaxs = "i")
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
        plot((sizewin/2) + (0:(seqisize1 - 1))*shift, recd$li[1:seqisize1, fac], type = "l", lty = 1,
            col = col[1], xlim = xlim, ylim = c(min(recd$li[, fac]), max(recd$li[, fac])),
            main = "Direct strand", xlab = "", ylab = "Factor scores", bty = 'l') # reading frame 1
        lines((sizewin/2) + (0:(seqisize2 - 1))*shift + 1, recd$li[(seqisize1 + 1):(seqisize1 + seqisize2), fac],
            lty = 2, col = col[2], ylab = "2") # reading frame 2
        lines((sizewin/2) + (0:(seqisize - 1))*shift + 2, recd$li[(seqisize1 + seqisize2 + 1):(dim(recd$li)[1]), fac],
            lty = 3, col = col[3], ylab = "3") # reading frame 3
        legend("topleft", legend = c(paste("Sequence name:", seqname), paste("Sequence length:", seqsize, "bp")),
            inset = c(-0.15, -0.2), bty = "n", xpd = TRUE)
        vstopdindphase <- numeric()
        if (length(vstopd) > 0)
        { # test if vector is not empty because problem with modulo
            vstopdindphase <- sapply(1:length(vstopd), function(x)
            { # index vector of reading frame of vector vstopd
                if (vstopd[x]%%3 == 1)
                {
                    vstopdindphase <- c(vstopdindphase, 2.5)
                }
                else
                {
                    if (vstopd[x]%%3 == 2)
                    {
                        vstopdindphase <- c(vstopdindphase, 1.5)
                    }
                    else
                    {
                        vstopdindphase <- c(vstopdindphase, 0.5)
                    }
                }
            })
        }
        vinitdindphase <- numeric()
        if (length(vinitd) > 0)
        { # test if vector is not empty because problem with modulo
            vinitdindphase <- sapply(1:length(vinitd), function(x)
            { # index vector of reading frame of vector vinitd
                if (vinitd[x]%%3 == 1)
                {
                    vinitdindphase <- c(vinitdindphase, 3)
                }
                else
                {
                    if (vinitd[x]%%3 == 2)
                    {
                        vinitdindphase <- c(vinitdindphase, 2)
                    }
                    else
                    {
                        vinitdindphase <- c(vinitdindphase, 1)
                    }
                }
            })
        }
        par(mar = c(5, 4, 3, 2) + 0.1)
        plot(vstopd, vstopdindphase, pch = 25, cex = 0.7, xlim = xlim, ylim = c(0.25, 3), axes = TRUE,
            ann = TRUE, tcl = -0.5, bty = 'l', yaxt = 'n', xlab = "Start/Stop positions (bp)",
            ylab = '', xpd = FALSE) # stop codons positions
        points(vinitd, vinitdindphase, pch = 24, bg = "slategray", cex = 0.7, col = 'slategray') # start codons positions
        abline(h = c(3.1, 2.4, 2.1, 1.4,  1.1, 0.4), col = c(col[1], col[1], col[2], col[2], col[3], col[3]),
            lty = c(1, 1, 2, 2, 3, 3))
        text(x = (xlim[1]-(xlim[2]-xlim[1])*0.75/6), pos = 4, y = c(2.75, 1.75, 0.75), labels = paste("Ph. ", c(0, 1, 2)), xpd = TRUE)
    }
    ##
    ##reverse strand##
    ##
    if (!direct)
    {    
        plot((sizewin/2) + (0:(seqisize1 - 1))*shift, recr$li[1:seqisize1, fac], type = "l", lty = 1,
            col = col[1], xlim = xlim, ylim = c(min(recr$li[, fac]), max(recr$li[, fac])),
            main = "Reverse strand", xlab = "", ylab = "Factor scores", bty = 'l') # reading frame 1
        lines((sizewin/2) + (0:(seqisize2-1))*shift + 1, recr$li[(seqisize1 + 1):(seqisize1 + seqisize2), fac],
            lty = 2, col = col[2], ylab="2") # reading frame 2
        lines((sizewin/2) + (0:(seqisize - 1))*shift + 2, recr$li[(seqisize1 + seqisize2 + 1):(dim(recr$li)[1]), fac],
            lty = 3, col = col[3], ylab = "3") # reading frame 3
        legend("topleft", legend = c(paste("Sequence name:", seqname), paste("Sequence length:", seqsize, "bp")),
            inset = c(-0.15, -0.2), bty = "n", xpd = TRUE)
        vstoprindphase <- numeric()
        if (length(vstopr) > 0)
        { # test if vector is not empty because problem with modulo
            vstoprindphase <- sapply(1:length(vstopr), function(x)
            { # index vector of reading frame of vector vstopr
                if (vstopr[x]%%3 == 1)
                {
                    vstoprindphase <- c(vstoprindphase, 2.5)
                }
                else
                {
                    if (vstopr[x]%%3 == 2)
                    {
                        vstoprindphase <- c(vstoprindphase, 1.5)
                    }
                    else
                    {
                        vstoprindphase <- c(vstoprindphase, 0.5)
                    }
                }
            })
        }
        vinitrindphase <- numeric()
        if (length(vinitr) > 0)
        { # test if vector is not empty because problem with modulo
            vinitrindphase <- sapply(1:length(vinitr), function(x)
            { # index vector of reading frame of vector vinitr
                if (vinitr[x]%%3 == 1)
                {
                    vinitrindphase <- c(vinitrindphase, 3)
                }
                else
                {
                    if (vinitr[x]%%3 == 2)
                    {
                        vinitrindphase <- c(vinitrindphase, 2)
                    }
                    else
                    {
                        vinitrindphase <- c(vinitrindphase, 1)
                    }
                }
            })
        }
        par(mar = c(5, 4, 3, 2) + 0.1)
        plot(vstopr, vstoprindphase, pch = 25, cex = 0.7, xlim = xlim, ylim = c(0.25, 3), axes = TRUE,
            ann = TRUE, tcl = -0.5, bty = 'l', yaxt = 'n', xlab = "Start/Stop positions (bp)",
            ylab = '', xpd = FALSE) # stop codons positions
        points(vinitr, vinitrindphase, pch = 24, bg = "slategray", cex = 0.7, col = 'slategray') # start codons positions
        abline(h = c(3.1, 2.4, 2.1, 1.4,  1.1, 0.4), col = c(col[1], col[1], col[2], col[2], col[3], col[3]),
            lty = c(1, 1, 2, 2, 3, 3))
        text(x = (xlim[1]-(xlim[2]-xlim[1])*0.75/6), pos = 4, y = c(2.75, 1.75, 0.75), labels = paste("Ph. ", c(0, 1, 2)), xpd = TRUE)
    }
}