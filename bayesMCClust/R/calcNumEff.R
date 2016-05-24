calcNumEff <-
function( outList, thin=1, printXi=TRUE, printE=TRUE, printBeta=TRUE, grLabels=paste("Group", 1:outList$Prior$H) ) {

    numEffTables <- list()
    
    # outList$xi_h_m
    if ( is.element("xi_h_m", names(outList)) & printXi ) { # 
        
        numEffXihm <- t( mapply(function(x1,x2,x3) c( j=x1, k=x2, h=x3, round( unlist( numEff( outList$xi_h_m[ x1, x2, x3, seq(outList$Mcmc$M0, outList$Mcmc$M, thin) ] ) ), 8 )), 
                   rep(1:(outList$K+1), each=(outList$K+1)), 1:(outList$K+1), rep(1:outList$Prior$H, each=(outList$K+1)*(outList$K+1)) ) )
        for ( h in 1:outList$Prior$H ) {     
            print( xtable( matrix( numEffXihm[ ((outList$K+1)*(outList$K+1)*(h-1)+1):((outList$K+1)*(outList$K+1)*h), outList$K], (outList$K+1), (outList$K+1), byrow=TRUE), 
                   digits=c(1, rep(3, outList$K+1) ),
                   caption=paste("Numerical efficiency: Inefficiency factors of the MCMC draws (based on DMC) obtained for each row $j=1,...,K+1$ 
                                  of the cluster-specific transition matrices $\\boldsymbol{\\xi}_{h,j\\cdot}$ for cluster", h, "(", grLabels[h], ")." ),
                   label=paste("tab:numEffXihm", h, sep="") )) # Group h
            cat("\n") 
        }        
        numEffTables[["numEffXihm"]] <- numEffXihm
    }
    
    # outList$e_h_m
    if ( is.element("e_h_m", names(outList)) & printE ) { # 
    
        numEffEhm <- t( mapply(function(x1,x2,x3) c( j=x1, k=x2, h=x3, round( unlist( numEff( outList$e_h_m[ x1, x2, x3, seq(outList$Mcmc$M0, outList$Mcmc$M, thin) ] ) ), 8 )), 
                   rep(1:(outList$K+1), each=(outList$K+1)), 1:(outList$K+1), rep(1:outList$Prior$H, each=(outList$K+1)*(outList$K+1)) ) )
        for ( h in 1:outList$Prior$H ) {     
            print( xtable( matrix( numEffEhm[ ((outList$K+1)*(outList$K+1)*(h-1)+1):((outList$K+1)*(outList$K+1)*h), outList$K], (outList$K+1), (outList$K+1), byrow=TRUE), 
                   digits=c(1, rep(3, outList$K+1) ),
                   caption=paste("Numerical efficiency: Inefficiency factors of the MCMC draws (based on DMC) obtained for each row $j=1,...,K+1$ 
                                  of the cluster-specific parameter matrices $\\mathbf{\\e}_{h,j\\cdot}$ for cluster", h, "(", grLabels[h], ")." ),
                   label=paste("tab:numEffEhm", h, sep="") )) # Group h
            cat("\n") 
        }
        numEffTables[["numEffEhm"]] <- numEffEhm
    }
    
    # outList$xi.m
    if ( is.element("xi.m", names(outList)) & printXi ) { # 
        
        numEffXim <- t( mapply(function(x1,x2,x3) c( j=x1, k=x2, h=x3, round( unlist( numEff( outList$xi.m[ seq(outList$Mcmc$M0, outList$Mcmc$M, thin), x1, x2, x3  ] ) ), 8 )), 
                   rep(1:(outList$K+1), each=(outList$K+1)), 1:(outList$K+1), rep(1:outList$Prior$H, each=(outList$K+1)*(outList$K+1)) ) )
        for ( h in 1:outList$Prior$H ) {     
            print( xtable( matrix( numEffXim[ ((outList$K+1)*(outList$K+1)*(h-1)+1):((outList$K+1)*(outList$K+1)*h), outList$K], (outList$K+1), (outList$K+1), byrow=TRUE), 
                   digits=c(1, rep(3, outList$K+1) ),
                   caption=paste("Numerical efficiency: Inefficiency factors of the MCMC draws (based on MCC) obtained for each row $j=1,...,K+1$ 
                                  of the cluster-specific transition matrices $\\boldsymbol{\\xi}_{h,j\\cdot}$ for cluster", h, "(", grLabels[h], ")." ),
                   label=paste("tab:numEffXim", h, sep="") )) # Group h
            cat("\n") 
        }
        numEffTables[["numEffXim"]] <- numEffXim
    }

    # outList$Beta.m
    if ( is.element("Beta.m", names(outList)) & printBeta ) { # 
    
        numEffBeta <- matrix(0, dim(outList$Beta.m)[1], outList$Prior$H-1)
        for ( h in 2:outList$Prior$H ) {         
            numEffBeta[,h-1] <- t( sapply( 1:dim(outList$Beta.m)[1], function(x) unlist( numEff( outList$Beta.m[x, h, seq(outList$Mcmc$M0,outList$Mcmc$M,thin) ] ) ) ) )[,2]    
        }
        rownames(numEffBeta) <- dimnames(outList$Beta.m)[[1]]
        colnames(numEffBeta) <- grLabels[-1]

        print( xtable( numEffBeta, digits=c(0, rep(5, outList$Prior$H-1) ), 
                    caption="Numerical efficiency: Inefficiency factors of the MCMC draws obtained for the MNL regression coefficients for each cluster.", 
                    label="tab:numEffBetas" ))
        cat("\n")
        
        numEffTables[["numEffBeta"]] <- numEffBeta
    }
    
    return(invisible(numEffTables))

}
