C1C2 <- rbind(c(0,1), c(0,2), c(1,2), c(2,2), c(1,1), c(0,3), c(1,3), c(2,3))
colnames(C1C2) <- c("C1", "C2")
## regDat <- data.frame(C1=C1C2[, "C1"],
##                      C2=C1C2[, "C2"],
##                      chr=c(6, 3, 5, 8, 9, 7, 2, 5),
##                      begin=c(110, 28, 50, 120, 70, 0, 0, 0),
##                      end=c(Inf, 48, Inf, Inf, Inf, 55, 70, 45),
##                      stringsAsFactors=FALSE)

regDat <- data.frame(C1=C1C2[, "C1"],
                     C2=C1C2[, "C2"],
                     chr=c(6, 3, 5, 8, 9, 7, 2, 5),
                     begin=c(110, 55, 120, 50, 70, 20, 10, 0),
                     end=c(145, 90, 150, 95, 100, 55, 40, 45),
                     stringsAsFactors=FALSE)
regDat$type <- sprintf("(%s,%s)", regDat[["C1"]], regDat[["C2"]])
