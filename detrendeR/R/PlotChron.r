PlotChron = 
function (crn, BANDWITH_SPLINE = 10, P_SPLINE = 0.5 ){
   # BANDWITH_SPLINE = 10;
   # P_SPLINE = 0.5;
        nCrn = ncol(crn)
        yr.vec = as.numeric(rownames(crn))
        crn.names <- colnames(crn)
        if (nCrn > 2) {
           # mat = matrix(1:(nCrn-2), (nCrn-2), 1)
           # layout(mat)
                 par(mfrow=c(2,1),mar = c(3, 3, 3, 2), mgp = c(1.25, 0.25, 0), tcl = 0.25)
            
            for (i in 1:(nCrn-2)) {       
                plot(yr.vec, crn[, i], type = "l", xlab = "",
                  ylab = "Tree-ring index", main = crn.names[i], las=1, lab=c(10,2, 0))
                spl = crn[, i]
                tmp = na.omit(spl)
                tmp = SPLINE(tmp, bandwidth = BANDWITH_SPLINE, p = P_SPLINE)
                spl[!is.na(spl)] = tmp
                lines(yr.vec, spl, col = "red", lwd = 2)
                abline(h = 1)
            }
            
                plot(data.frame(yr.vec, crn[,nCrn]), type = "l", xlab = "Years", ylab = "Sample depth",las=1, lab=c(10,5, 0))
                lines(yr.vec, crn[,nCrn-1], col = "red", lwd = 1)
        }
        }
       # PlotChron(crn)