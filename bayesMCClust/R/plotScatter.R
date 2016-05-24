plotScatter <-
function(outList, thin=1, xi11=c(1,1), xi12=c(2,2), xi21=c(2,2), xi22=c(3,3), xi31=c(1,1), xi32=c(3,3) ) {
    
    existsXim <- is.element("xi.m", names(outList))
    
    minCoord <- min(xi11, xi12, xi21, xi22, xi31, xi32)
    maxCoord <- min(xi11, xi12, xi21, xi22, xi31, xi32)
    
    stopifnot( (minCoord >= 1) & (maxCoord <= (outList$K+1)) )

    dev.new(width=12, height=6)
    par(mfrow=c(1,3))  
    
    plot( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi11[1], xi11[2], 1] } else { outList$xi_h_m[ xi11[1], xi11[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
          if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi12[1], xi12[2], 1] } else { outList$xi_h_m[ xi12[1], xi12[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] },     
          xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.3, col=2, 
          xlab=substitute(xi[aa], list(aa=paste(xi11,collapse=","))), 
          ylab=substitute(xi[aa], list(aa=paste(xi12,collapse=","))), 
          main=substitute(paste("MCMC draws for ", xi[aa], " vs ", xi[bb]), list(aa=paste(xi11,collapse=","), bb=paste(xi12,collapse=",")))
        )        
    if (outList$Prior$H > 1) {
        for (h in 2:outList$Prior$H) {        
            points( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi11[1], xi11[2], h] } else { outList$xi_h_m[ xi11[1], xi11[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
                    if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi12[1], xi12[2], h] } else { outList$xi_h_m[ xi12[1], xi12[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] },                     
                    col=h+1, pch=19, cex=0.3 )
        }
    }

    plot( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi21[1], xi21[2], 1] } else { outList$xi_h_m[ xi21[1], xi21[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
          if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi22[1], xi22[2], 1] } else { outList$xi_h_m[ xi22[1], xi22[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] },           
          xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.3, col=2, 
          xlab=substitute(xi[aa], list(aa=paste(xi21,collapse=","))), 
          ylab=substitute(xi[aa], list(aa=paste(xi22,collapse=","))), 
          main=substitute(paste("MCMC draws for ", xi[aa], " vs ", xi[bb]), list(aa=paste(xi21,collapse=","), bb=paste(xi22,collapse=",")))
        )
    if (outList$Prior$H > 1) {
        for (h in 2:outList$Prior$H) {        
            points( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi21[1], xi21[2], h] } else { outList$xi_h_m[ xi21[1], xi21[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
                    if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi22[1], xi22[2], h] } else { outList$xi_h_m[ xi22[1], xi22[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] },                     
                    col=h+1, pch=19, cex=0.3 )
        }
    }

    plot( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi31[1], xi31[2], 1] } else { outList$xi_h_m[ xi31[1], xi31[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
          if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi32[1], xi32[2], 1] } else { outList$xi_h_m[ xi32[1], xi32[2], 1, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] },           
          xlim=c(0,1), ylim=c(0,1), pch=19, cex=0.3, col=2, 
          xlab=substitute(xi[aa], list(aa=paste(xi31,collapse=","))), 
          ylab=substitute(xi[aa], list(aa=paste(xi32,collapse=","))), 
          main=substitute(paste("MCMC draws for ", xi[aa], " vs ", xi[bb]), list(aa=paste(xi31,collapse=","), bb=paste(xi32,collapse=",")))
        )
    if (outList$Prior$H > 1) {
        for (h in 2:outList$Prior$H) {
            points( if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi31[1], xi31[2], h] } else { outList$xi_h_m[ xi31[1], xi31[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
                    if ( existsXim ) { outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), xi32[1], xi32[2], h] } else { outList$xi_h_m[ xi32[1], xi32[2], h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)] }, 
                    col=h+1, pch=19, cex=0.3)
        }
    }
    
}
