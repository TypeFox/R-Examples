plotenz <-
function(sequences, enznames, enzdata, side = TRUE, 
        type = c("RFLP", "TRFLP"), Terminal = c("T5", "T3")) {

    if(length(sequences) >= 2){
	if(!inherits(sequences, "fasta"))
	{stop("The input sequences must be a \"fasta\" object.")}
	}
	match.arg(type)
	 addmark <-
     function(numLane) {
         par(mar = c(0, 0, 0, 0))
         par(bg = "grey30")
         plot(numLane, log(3000), type = "n", ylim = c(log(40), log(3000)), 
             xlim = c(0, numLane), yaxt = "n", xaxt = "n")
         a = c(seq(100, 800, by = 100), 1500, 3000)
         colMar = paste("grey", c(35, 50, 55, 60, 70, 100, 80, 85, 
             90, 100), sep = "")
         text(0.5, log(a), paste(a, "bp"), col = "white")
         segments(1.5 - 0.28, log(a), 1.5 + 0.28, log(a), col = colMar, 
             lwd = 4)
     }

     addband <- 
     function(a, i, tex) {
         cols = ifelse(a > 700, 90, ifelse(a > 600 & a <= 700, 85, 
                 ifelse(a <= 600 & a > 500, 80, ifelse(a <= 500 & a > 
                 400, 70, ifelse(a <= 400 & a > 300, 60, ifelse(a <= 
                 300 & a > 200, 55, ifelse(a <= 200 & a > 100, 50, 
                 ifelse(a <= 100 & a > 60, 35, 32))))))))
         col1 = paste("grey", cols, sep = "")
         segments((1.5 + i) - 0.28, log(a), (1.5 + i) + 0.28, log(a), col = col1, lwd = 4)
         text((1.5 + i), log(1500), tex, col = "white")
     }
	
    numsq = length(sequences)/2
	
    if (side) {
        numLane = 18
        dev = ceiling(numsq/(numLane - 2))
        enzNam = enznames
        nr = 2
        nc = 1
        layout(matrix(1:(nc * nr), nrow = nr, byrow = TRUE))
        for (j in 1:dev) {
            addmark(numLane)
            x = ifelse(j == max(dev), numsq - (j - 1) * (numLane - 
                2), (numLane - 2))
			if(x < 1) x = 1
            for (i in 1:x) {
                SqNum = i + (numLane - 2) * (j - 1)
                DNAsq = sequences[SqNum * 2]
                if (type == "RFLP") 
                  a = enzCut(enznames = enzNam, DNAsq = DNAsq, enzdata = enzdata)$RFLP.frag
                if (type == "TRFLP") 
                  a = enzCut(enznames = enzNam, DNAsq = DNAsq, enzdata = enzdata)$TRFLP[Terminal]
                tex = SqNum
                addband(a, i, tex)
            }
            for (k in 1:length(enzNam)) {
                if (k == 1) 
                  nam1 = enzNam[1]
                if (k > 1) 
                  nam1 = paste(enzNam[k], nam1, sep = "_")
            }
            nam = paste(type, nam1, sep = "  ")
            if (type == "TRFLP") 
                nam = paste(type, Terminal, nam1, sep = "  ")
            text(numLane/2, log(2000), nam, col = "white", font = 2)
            if (j/(nc * nr) - ceiling(j/(nc * nr)) == 0) {
                dev.new()
                layout(matrix(1:(nc * nr), nrow = nr, byrow = TRUE))
            }
        }
    }
	else {
        numLane = length(enznames) + 2
        nc = ifelse(numLane == 4, 6, ifelse(numLane == 5, 5, 
             ifelse(numLane == 6, 4, ifelse(numLane == 7, 3, 2))))
        nr = 3
        layout(matrix(1:(nc * nr), nrow = nr, byrow = TRUE))
        for (r in 1:numsq) {
            addmark(numLane)
            DNAsq = sequences[2 * r]
            nam = gsub(">", "", sequences[2 * r - 1])
            for (i in 1:length(enznames)) {
                enzNam = enznames[i]
                a = enzCut(enznames = enzNam, DNAsq = DNAsq, enzdata = enzdata)$RFLP.frag
                tex = enznames[i]
                addband(a, i, tex)
                text(numLane/2, log(2000), nam, col = "white", 
                  font = 2)
            }
            if (ceiling(r/(nc * nr)) - r/(nc * nr) == 0) {
                dev.new()
                layout(matrix(1:(nc * nr), nrow = nr, byrow = TRUE))
            }
        }
    }
}
