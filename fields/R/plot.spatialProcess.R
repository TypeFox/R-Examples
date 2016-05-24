# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"plot.spatialProcess" <- function(x, digits = 4, which = 1:4, 
    ...) {
    out <- x
    #
    #   don't do plots 2:4 if a fixed lambda
    #
    if (x$fixed.model) {
        which <- 1
    }
    fitted.values <- predict(out)
    std.residuals <- (out$residuals * sqrt(out$weights))/out$shat.GCV
    if (any(which == 1)) {
        temp <- summary(out)
        plot(fitted.values, out$y, ylab = "Y", xlab = " predicted values", 
            bty = "n", ...)
        abline(0, 1)
        hold <- par("usr")
       title("Observations by predicted values")    

    }
    if (any(which == 2)) {
        plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
        hold <- par("usr")
        text(hold[1], hold[4], paste(" RMSE =", format(signif(sqrt(sum(out$residuals^2)/(temp$num.observation - 
            temp$enp)), digits))), cex = 0.8, adj = 0)
    }
    if (any(which == 3)) {
    	mar.old<- par()$mar
    	par( mar= mar.old + c(0,0,0,2) )
        if (nrow(out$gcv.grid) > 1) {
            ind <- out$gcv.grid[, 3] < 1e+19
            out$gcv.grid <- out$gcv.grid[ind, ]
            yr <- range(unlist(out$gcv.grid[, 3:5]), na.rm = TRUE)
            plot(out$gcv.grid[, 2], out$gcv.grid[, 3], xlab = "Eff. number of parameters", 
                ylab = " GCV function", bty = "n", ylim = yr, 
             
                 ...)
            lines(out$gcv.grid[, 2], out$gcv.grid[, 4], lty = 3)
            lines(out$gcv.grid[, 2], out$gcv.grid[, 5], lty = 1)
            xline(out$eff.df, lwd=2, col="grey")
            usr.save<- par()$usr
            usr.save[3:4]<- range( -out$gcv.grid[,7] )
            par( usr= usr.save, ylog=FALSE)
            lines( out$gcv.grid[, 2], -out$gcv.grid[,7] ,
            lty=2, lwd=2, col="blue")
            axis( side=4)
            mtext( side=4, line=2, "log Profile Likelihood(lamdba)",cex=.75)
            title("GCV-points, solid-model, dots- single  \n REML dashed", 
                cex = 0.6)
            box()
            par( mar=mar.old)
        }
    }
    if (any(which == 4)) {
    	plot( out$eval.grid[,c(1,6)], pch=16, xlab="theta (range parameter)", ylab="log Profile Likelihood (theta)")
    	title("Profile likelihood for theta \n (range parameter)")
    	xline( out$theta.MLE, lwd=2, col="grey")
    	xr<- range( out$eval.grid[,1])
    	xg<- seq( xr[1], xr[2],, 100)
    	ind<- !duplicated(out$eval.grid[,1] )
    	yg<- splint(out$eval.grid[ind,1], out$eval.grid[ind,6], xg )
    	lines( xg, yg, col="grey", lty=2, lwd=1)
           }
}
