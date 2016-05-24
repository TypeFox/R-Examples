C1C2 <- rbind(c(0,1), c(1,1), c(1,2), c(0,2), c(2,2), c(0,3), c(0,0))
colnames(C1C2) <- c("C1", "C2")
regDat <- data.frame(C1=C1C2[, "C1"],
                     C2=C1C2[, "C2"],
                     chr=c(18, 21, 2, 4, 9, 3, 13),
                     begin=c(57, 15, 45, 144, 0, 0, 52.14),
                     end=c(72, 35, 64, Inf, 15, 50, 60.4),
                     stringsAsFactors=FALSE)
regDat$type <- sprintf("(%s,%s)", regDat[["C1"]], regDat[["C2"]])

