peakDrawing <-
function(vecIP, vecControl, lineIP, lineControl, lineFC, lineAverage,
                       posInf = 1, posSup = NULL, add = 10, title = ""){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######

    # argument names were changed for user convenience
    vecINPUT  = vecControl
    lineINPUT = lineControl

    #######
    # ----- Parameter details -----------
    #
    # vecIP    : sequencing depth at each nucleotide (from start pos = 1 to end)
    # vecINPUT : sequencing depth at each nucleotide (from start pos = 1 to end)
    # posInf and posSup : first and last positions for the genomic region to be drawn
    # lineIP, lineINPUT, lineFC, lineAverage : threshold values used for peak detection
    # add      : nucleotides can be added before and after the specified genomic region
    # title    : general title for the graphs
    #######

    # by default, the drawing genomic region starts at position 1 and 
    # ends at the final position
    if(is.null(posSup)){
        posSup = length(vecIP)
    }

    # x label for graphs
    xlabels = "chromosomal location"

    # a page is divided into four parts
    par(mfrow = c(2,2))

    if((posInf-add > 0) & (posSup+add <= length(vecIP))){

    plot(vecIP, type = "h", xlab = xlabels, ylab = "T1: IP signal",
            main = paste(title, "\nIP sample (T1)"), col = "red", 
            xlim = c(posInf-add,posSup+add),
            ylim = c(0, max(vecIP, vecINPUT, na.rm = T)))
    abline(v = posInf, lty = "dashed", col = "blue")
    abline(v = posSup, lty = "dashed", col = "blue")
    abline(h = lineIP, col = "red")

    vecLogFC   = log2((vecIP + 1)/(vecINPUT + 1))
 
    plot(vecLogFC, type = "l", xlab = xlabels, ylab = "T3: log2FC",
            main = paste(title, "\nlog2FC (T3)"), col = "orange",
            xlim = c(posInf-add,posSup+add),
            ylim = c(min(vecLogFC, na.rm = T), max(vecLogFC, na.rm = T)))
    abline(v = posInf, lty = "dashed", col = "blue")
    abline(v = posSup, lty = "dashed", col = "blue")
    abline(h = lineFC, col = "red")

    plot(vecINPUT, type = "h", xlab = xlabels, ylab = "T2: control signal",
            main = paste(title, "\ncontrol sample (T2)"), col = "blue",
            xlim = c(posInf-add,posSup+add), 
            ylim = c(0, max(vecIP, vecINPUT, na.rm = T)))
    abline(v = posInf, lty = "dashed", col = "blue")
    abline(v = posSup, lty = "dashed", col = "blue")
    abline(h = lineINPUT, col = "red")

    averageLog = (log2(vecIP + 1) + log2(vecINPUT + 1))/2

    plot(averageLog, type = "l", xlab = xlabels, ylab = "T4: (log2(IP) + log2(control)) / 2",
            main = paste(title, "\naverage log2 signals (T4)"), col = "purple",
            xlim = c(posInf-add,posSup+add),
            ylim = c(min(averageLog, na.rm = T), max(averageLog, na.rm = T)))
    abline(v = posInf, lty = "dashed", col = "blue")
    abline(v = posSup, lty = "dashed", col = "blue")
    abline(h = lineAverage, col = "red")


    }else{
    
        print("!! peakDrawing() error !! check the values of genomic positions")
    }    

# end of function peakDrawing()
}
