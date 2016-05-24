   
extractimf <- function(residue, tt=NULL, tol=sd(residue)*0.1^2, max.sift=20, 
                        stoprule="type1", boundary="periodic", sm="none", spar=NULL, alpha=NULL, check=FALSE, weight=NULL) {     
    
    if (boundary == "none")
        minextrema <- 4
    else
        minextrema <- 2
            
    if ((sm == "spline" || sm == "kernel" || sm == "locfit") & (is.null(spar) || spar == 0)) stop("Provide the smoothing parameter.\n") #|| sm == "quantile"
    #if (sm == "quantile" & is.null(alpha)) stop("Provide the quantile for quantile regression.\n")
        
    ndata <- length(residue); ndatam1 <- ndata - 1 
    if(is.ts(residue)) 
        residue <- as.numeric(residue) 
    
    if(is.null(tt)) tt <- 1:length(residue)
    
    emin <- emax <- em <- h <- imf <- NULL
    n2data <- 2*ndata; tt2 <- 1:n2data
    n3data <- 3*ndata; n3datam1 <- n3data-1; tt3 <- 1:n3data

    input <- residue; rangext <- range(residue)
    
    j <- 1
    ### DANIEL BOWMAN'S EDITS JANUARY 15 2013
    s_count <- 0
    prev_excross <- 0
    ### END DANIEL'S EDITS    
    repeat {
        tmp <- extrema(input, ndata, ndatam1)
        
        if (tmp$nextreme <= minextrema) break
        
        if(j == 1 || boundary == "wave") {
            minindex <- unique(c(t(tmp$minindex))); minn <- length(minindex)
            maxindex <- unique(c(t(tmp$maxindex))); maxn <- length(maxindex)          
  
            extminindex <- minindex; extmaxindex <- maxindex
            tmpwavefreq1 <- diff(sort(c(tt[1], tt[minindex[1]], tt[maxindex[1]])))
            tmpwavefreq2 <- diff(sort(c(tt[minindex[minn]], tt[maxindex[maxn]], tt[ndata])))
                                         
            if(input[1] <= input[minindex[1]] && input[1] <= input[maxindex[1]]) {
                extminindex <- c(1, extminindex); wavefreq1 <- 2 * tmpwavefreq1[1] 
            } else if(input[1] >= input[minindex[1]] && input[1] >= input[maxindex[1]]) {
                extmaxindex <- c(1, extmaxindex); wavefreq1 <- 2 * tmpwavefreq1[1] 
            } else if(input[1] >= (input[minindex[1]] + input[maxindex[1]])/2) {
                wavefreq1 <- tmpwavefreq1[2]  + max(tmpwavefreq1[2], 2*tmpwavefreq1[1])        
            } else {
                wavefreq1 <- tmpwavefreq1[2]  + max(tmpwavefreq1[2], round(1.5*tmpwavefreq1[1]))
            }                 
                  
            if(input[ndata] <= input[minindex[minn]] && input[ndata] <= input[maxindex[maxn]]) {
                extminindex <- c(extminindex, ndata); wavefreq2 <- 2 * tmpwavefreq2[2] 
            } else if(input[ndata] >= input[minindex[minn]] && input[ndata] >= input[maxindex[maxn]]) {
                extmaxindex <- c(extmaxindex, ndata); wavefreq2 <- 2 * tmpwavefreq2[2] 
            } else if(input[ndata] >= (input[minindex[minn]] + input[maxindex[maxn]])/2) {
                wavefreq2 <- tmpwavefreq2[1]  + max(tmpwavefreq2[1], 2*tmpwavefreq2[2])        
            } else {
                wavefreq2 <- tmpwavefreq2[1]  + max(tmpwavefreq2[1], round(1.5*tmpwavefreq2[2]))
            } 

            extminn <- length(extminindex); extmaxn <- length(extmaxindex)
            extttminindex <- c(tt[extminindex[1]] - 4:1 * wavefreq1, tt[extminindex], tt[extminindex[extminn]] + 1:4 * wavefreq2)
            extttmaxindex <- c(tt[extmaxindex[1]] - 4:1 * wavefreq1, tt[extmaxindex], tt[extmaxindex[extmaxn]] + 1:4 * wavefreq2)          
            
            if(sm == "none") {        
                f <- splinefun(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)))
                emin <- cbind(emin, f(tt))
                f <- splinefun(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)))
                emax <- cbind(emax, f(tt))
            } else if (sm == "spline") {                
                f <- sreg(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)), lambda = spar)               
                emin <- cbind(emin, predict(f, tt))
                f <- sreg(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)), lambda = spar)               
                emax <- cbind(emax, predict(f, tt)) 
            #} else if (sm == "quantile") {                
            #    f <- qsreg(extttminindex, c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4)), lam = spar, alpha=1 - alpha)               
            #    emin <- cbind(emin, predict(f, tt))
            #    f <- qsreg(extttmaxindex, c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4)), lam = spar, alpha=alpha)
            #    emax <- cbind(emax, predict(f, tt)) 
            } else if (sm == "kernel") {
                f <- loess(yy~xx, data.frame(xx=extttminindex, yy=c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4))), span=spar, family = "gaussian", degree=2, control=loess.control(surface="direct"))               
                emin <- cbind(emin, predict(f, data.frame(xx=tt))) 
                f <- loess(yy~xx, data.frame(xx=extttmaxindex, yy=c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4))), span=spar, family = "gaussian", degree=2, control=loess.control(surface="direct"))  
                emax <- cbind(emax, predict(f, data.frame(xx=tt)))
            } else if (sm == "locfit") {
                f <- locfit(yy~xx, data=data.frame(xx=extttminindex, yy=c(rep(input[extminindex[1]], 4), input[extminindex], rep(input[extminindex[extminn]], 4))), deg=3, kern="gauss", alpha=spar)
                emin <- cbind(emin, predict(f, data.frame(xx=tt)))     
                f <- locfit(yy~xx, data=data.frame(xx=extttmaxindex, yy=c(rep(input[extmaxindex[1]], 4), input[extmaxindex], rep(input[extmaxindex[extmaxn]], 4))), deg=3, kern="gauss", alpha=spar)
                emax <- cbind(emax, predict(f, data.frame(xx=tt)))
            } 
            
            em <- cbind(em, (emin[,j] + emax[,j]) / 2)
            
            if(check){
                plot(tt, input, type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep="")) 
                points(tt[unique(c(t(tmp$minindex)))], input[unique(c(t(tmp$minindex)))], col=4)
                points(tt[unique(c(t(tmp$maxindex)))], input[unique(c(t(tmp$maxindex)))], col=2)
                lines(tt, emin[,j], col=4)
                lines(tt, emax[,j], col=2)
                lines(tt, em[,j]); locator(1)
            }
        } else {
        if(boundary == "none") {
            minindex <- c(1, unique(c(t(tmp$minindex))), ndata)
            maxindex <- c(1, unique(c(t(tmp$maxindex))), ndata)  
        
            if(sm == "none") {  
                fmin <- splinefun(tt[minindex], input[minindex])
                emin <- cbind(emin, fmin(tt))
        
                fmax <- splinefun(tt[maxindex], input[maxindex])
                emax <- cbind(emax, fmax(tt))
            } else if (sm == "spline") {              
                f <- sreg(tt[minindex], input[minindex], lambda = spar)               
                emin <- cbind(emin, predict(f, tt))
                
                f <- sreg(tt[maxindex], input[maxindex], lambda = spar)
                emax <- cbind(emax, predict(f, tt)) 
            #} else if (sm == "quantile") {                
            #    f <- qsreg(tt[minindex], input[minindex], lam = spar, alpha=1 - alpha)               
            #    emin <- cbind(emin, predict(f, tt))
            #    
            #    f <- qsreg(tt[maxindex], input[maxindex], lam = spar, alpha=alpha)
            #    emax <- cbind(emax, predict(f, tt)) 
            } else if (sm == "kernel") {
                f <- loess(yy~xx, data.frame(xx=tt[minindex], yy=input[minindex]), span=spar, family = "gaussian", degree=2, control=loess.control(surface="direct"))               
                emin <- cbind(emin, predict(f, data.frame(xx=tt)))
                
                f <- loess(yy~xx, data.frame(xx=tt[maxindex], yy=input[maxindex]), span=spar, family = "gaussian", degree=2, control=loess.control(surface="direct"))  
                emax <- cbind(emax, predict(f, data.frame(xx=tt)))
            } else if (sm == "locfit") {
                f <- locfit(yy~xx, data=data.frame(xx=tt[minindex], yy=input[minindex]), deg=3, kern="gauss", alpha=spar)
                emin <- cbind(emin, predict(f, data.frame(xx=tt)))
                
                f <- locfit(yy~xx, data=data.frame(xx=tt[maxindex], yy=input[maxindex]), deg=3, kern="gauss", alpha=spar)
                emax <- cbind(emax, predict(f, data.frame(xx=tt)))
            } 
            
            em <- cbind(em, (emin[,j] + emax[,j]) / 2)    
                   
            if(check){
                plot(tt, input, type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep=""))#,  ylim=rangext) 
                points(tt[unique(c(t(tmp$minindex)))], input[unique(c(t(tmp$minindex)))], col=4)
                points(tt[unique(c(t(tmp$maxindex)))], input[unique(c(t(tmp$maxindex)))], col=2)
                lines(tt, emin[,j], col=4)
                lines(tt, emax[,j], col=2)
                lines(tt, em[,j]); locator(1)
            }            
        }   
        if(boundary == "symmetric" || boundary == "periodic") {
            if(boundary == "symmetric") {
                inputext <- c(rev(input[-1]), input, rev(input)[-1])   
                ttext <- c(tt[1] - rev(cumsum(diff(tt))), tt, tt[ndata] + cumsum(rev(diff(tt))))        
                tmp <- extrema(inputext, n3data - 2, n3data - 3) 
            } else if (boundary == "periodic") {
                inputext <- c(input[-ndata], input, input[-1])   
                ttext <- c(rev(tt[1] - cumsum(rev(diff(tt)))), tt, tt[ndata] + cumsum(diff(tt)))        
                tmp <- extrema(inputext, n3data - 2, n3data - 3)            
            }           
            minindex <- unique(c(t(tmp$minindex)))
            maxindex <- unique(c(t(tmp$maxindex)))
                        
            if(sm == "none") {     
                fmin <- splinefun(ttext[minindex], inputext[minindex])
                tmpmin <- fmin(ttext[ndata:(n2data-1)])
                
                fmax <- splinefun(ttext[maxindex], inputext[maxindex])
                tmpmax <- fmax(ttext[ndata:(n2data-1)])
            } else if (sm == "spline") { 
                fmin <- sreg(ttext[minindex], inputext[minindex], lambda = spar)
                tmpmin <- predict(fmin, ttext[ndata:(n2data-1)])          
                
                fmax <- sreg(ttext[maxindex], inputext[maxindex], lambda = spar)
                tmpmax <- predict(fmax, ttext[ndata:(n2data-1)])                        
            #} else if (sm == "quantile") { 
            #    fmin <- qsreg(ttext[minindex], inputext[minindex], lam = spar, alpha=1 - alpha)
            #    tmpmin <- predict(fmin, ttext[ndata:(n2data-1)])          
            #   
            #    fmax <- qsreg(ttext[maxindex], inputext[maxindex], lam = spar, alpha=alpha)
            #    tmpmax <- predict(fmax, ttext[ndata:(n2data-1)])                        
            } else if (sm == "kernel") {                
                fmin <- loess(yy ~ xx, data.frame(xx=ttext[minindex], yy=inputext[minindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct"))               
                tmpmin <- predict(fmin, data.frame(xx=ttext[ndata:(n2data-1)]))
                
                fmax <- loess(yy ~ xx, data.frame(xx=ttext[maxindex], yy=inputext[maxindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct")) 
                tmpmax <- predict(fmax, data.frame(xx=ttext[ndata:(n2data-1)]))
            } else if (sm == "locfit") {
                fmin <- locfit(yy~xx, data=data.frame(xx=ttext[minindex], yy=inputext[minindex]), deg=3, kern="gauss", alpha=spar)
                tmpmin <- predict(fmin, data.frame(xx=ttext[ndata:(n2data-1)]))
                
                fmax <- locfit(yy~xx, data=data.frame(xx=tt[maxindex], yy=input[maxindex]), deg=3, kern="gauss", alpha=spar)
                tmpmax <- predict(fmax, data.frame(xx=ttext[ndata:(n2data-1)]))
            } 
         
            emin <- cbind(emin, tmpmin[1:ndata])
            emax <- cbind(emax, tmpmax[1:ndata])
            
            tmpm <- (tmpmin+tmpmax)/2; em <- cbind(em, tmpm[1:ndata])

            if(check){
                plot(ttext[ndata:(n2data-1)], inputext[ndata:(n2data-1)], type="l", col=3, xlab="", ylab="", main=paste("Boundary = ", boundary, sep=""))
                points(ttext[minindex[minindex > ndata & minindex < n2data]], inputext[minindex[minindex > ndata & minindex < n2data]], col=4)
                points(ttext[maxindex[maxindex > ndata & maxindex < n2data]], inputext[maxindex[maxindex > ndata & maxindex < n2data]], col=2)
                lines(ttext[ndata:(n2data-1)], tmpmin, col=4)
                lines(ttext[ndata:(n2data-1)], tmpmax, col=2)
                lines(ttext[ndata:(n2data-1)], tmpm); locator(1)
            } 
        } else if(boundary == "evenodd") {
            inputeven <- c(input, rev(input), input) 
            ttext <- c(tt, tt[ndata]+tt[ndata]-tt[ndatam1], tt[ndata]+tt[ndata]-tt[ndatam1] + cumsum(rev(diff(tt))))
          
            tmp <- extrema(inputeven, n3data, n3datam1) 

            minindex <- unique(c(t(tmp$minindex)))
            minindex <- minindex[minindex <= n2data + minindex[1]]
            maxindex <- unique(c(t(tmp$maxindex)))
            maxindex <- maxindex[maxindex <= n2data + maxindex[1]]
            
            if(sm == "none") { 
                f <- splinefun(ttext[minindex], inputeven[minindex])
                emineven <- f(ttext[1:ndata]) 
                #lines(tt2, emineven, col=4)

                f <- splinefun(ttext[maxindex], inputeven[maxindex])
                emaxeven <- f(ttext[1:ndata])
                #lines(tt2, emaxeven, col=2)
            } else if (sm == "spline") { 
                fmin <- sreg(ttext[minindex], inputeven[minindex], lambda = spar)
                emineven <- predict(fmin, ttext[1:ndata])          
                
                fmax <- sreg(ttext[maxindex], inputeven[maxindex], lambda = spar)
                emaxeven <- predict(fmax, ttext[1:ndata])                        
            #} else if (sm == "quantile") { 
            #    fmin <- qsreg(ttext[minindex], inputeven[minindex], lam = spar, alpha=1 - alpha)
            #    emineven <- predict(fmin, ttext[1:ndata])          
            #   
            #    fmax <- qsreg(ttext[maxindex], inputeven[maxindex], lam = spar, alpha=alpha)
            #    emaxeven <- predict(fmax, ttext[1:ndata])                        
            }  else if (sm == "kernel") {                
                fmin <- loess(yy ~ xx, data.frame(xx=ttext[minindex], yy=inputeven[minindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct"))               
                emineven <- predict(fmin, data.frame(xx=ttext[1:ndata]))
                
                fmax <- loess(yy ~ xx, data.frame(xx=ttext[maxindex], yy=inputeven[maxindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct")) 
                emaxeven <- predict(fmax, data.frame(xx=ttext[1:ndata]))
            } else if (sm == "locfit") {
                fmin <- locfit(yy~xx, data=data.frame(xx=ttext[minindex], yy=inputeven[minindex]), deg=3, kern="gauss", alpha=spar)
                emineven <- predict(fmin, data.frame(xx=ttext[1:ndata]))
                
                fmax <- locfit(yy~xx, data=data.frame(xx=tt[maxindex], yy=inputeven[maxindex]), deg=3, kern="gauss", alpha=spar)
                emaxeven <- predict(fmax, data.frame(xx=ttext[1:ndata]))             
            }
            
            inputodd <- c(input, -rev(input), input)
            tmp <- extrema(inputodd, n3data, n3datam1) 

            minindex <- unique(c(t(tmp$minindex)))
            minindex <- minindex[minindex <= n2data + minindex[1]]
            maxindex <- unique(c(t(tmp$maxindex)))
            maxindex <- maxindex[maxindex <= n2data + maxindex[1]]            

            if(sm == "none") { 
                f <- splinefun(ttext[minindex], inputodd[minindex])
                eminodd <- f(ttext[1:ndata]) 
                #lines(tt2, eminodd, col=4)

                f <- splinefun(ttext[maxindex], inputodd[maxindex])
                emaxodd <- f(ttext[1:ndata])
                #lines(tt2, emaxodd, col=2)
            } else if (sm == "spline") { 
                fmin <- sreg(ttext[minindex], inputodd[minindex], lambda = spar)
                eminodd <- predict(fmin, ttext[1:ndata])          
                
                fmax <- sreg(ttext[maxindex], inputodd[maxindex], lambda = spar)
                emaxodd <- predict(fmax, ttext[1:ndata])                        
            #} else if (sm == "quantile") { 
            #    fmin <- qsreg(ttext[minindex], inputodd[minindex], lam = spar, alpha=1 - alpha)
            #    eminodd <- predict(fmin, ttext[1:ndata])            
            #   
            #    fmax <- qsreg(ttext[maxindex], inputodd[maxindex], lam = spar, alpha=alpha)
            #    emaxodd <- predict(fmax, ttext[1:ndata])                        
            } else if (sm == "kernel") {                
                fmin <- loess(yy ~ xx, data.frame(xx=ttext[minindex], yy=inputodd[minindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct"))               
                eminodd <- predict(fmin, data.frame(xx=ttext[1:ndata]))
                
                fmax <- loess(yy ~ xx, data.frame(xx=ttext[maxindex], yy=inputodd[maxindex]), span=spar, family = "gaussian", degree=3, control=loess.control(surface="direct")) 
                emaxodd <- predict(fmax, data.frame(xx=ttext[1:ndata]))
            } else if (sm == "locfit") {
                fmin <- locfit(yy~xx, data=data.frame(xx=ttext[minindex], yy=inputodd[minindex]), deg=3, kern="gauss", alpha=spar)
                eminodd <- predict(fmin, data.frame(xx=ttext[1:ndata]))
                
                fmax <- locfit(yy~xx, data=data.frame(xx=tt[maxindex], yy=inputodd[maxindex]), deg=3, kern="gauss", alpha=spar)
                emaxodd <- predict(fmax, data.frame(xx=ttext[1:ndata]))           
            } 
  
            emin <- cbind(emin, (emineven+eminodd)/2)
            emax <- cbind(emax, (emaxeven+emaxodd)/2)
            em <- cbind(em, (emin[,j]+emax[,j])/2)
            
            if(check) {
                plot(tt, (inputeven+inputodd)/2[1:ndata], type="l", col=3, ylab="")
                lines(tt, emin, col=4)
                lines(tt, emax, col=2)
                lines(tt, (emin+emax)/2); locator(1)
            }  
        }
        }
        
        h <- cbind(h, input - em[,j])    

        if (stoprule == "type1" && (all(abs(em[,j]) < tol) || j >= max.sift)) { 
            imf <- h[,j]
            residue <- residue - imf
            break
        } else if (stoprule == "type2" && j >= 2) {
            if (sum((h[2:ndatam1, j-1]-h[2:ndatam1, j])^2/h[2:ndatam1, j-1]^2) < tol || j >= max.sift) { 
                imf <- h[,j]
                residue <- residue - imf
                break
                }
        } else if (stoprule == "type3" && j >= 2) {
            if (sum((h[, j-1]-h[, j])^2)/sum(h[, j-1]^2) < tol || j >= max.sift) { 
                imf <- h[,j]
                residue <- residue - imf
                break
                }
        } else if (stoprule == "type4") {
             if (sum((em[, j])^2)/sum(h[, j]^2) < tol || j >= max.sift) {  
                imf <- h[,j]
                residue <- residue - imf
                break
                }
        } #DANIEL BOWMAN'S EDITS, JANUARY 15 2013
	else if (stoprule == "type5" && j >= 2) { 
             if (abs(tmp$nextreme - tmp$ncross) <= 1 && tmp$nextreme+tmp$ncross == prev_excross) {
                 s_count <- s_count+1
                 if (s_count >= tol || j >= max.sift) {
                     imf <- h[,j]
                     residue <- residue - imf
                     break
                }
            }
            prev_excross <- tmp$nextreme+tmp$ncross
        }
        #END DANIEL'S EDITS
        
        input <- h[,j]
       
        j <- j+1
    }
    
    if(check) list(emin=emin, emax=emax, em=em, h=h, imf=imf, residue=residue, niter=j) else   
    list(imf=imf, residue=residue, niter=j)
}

emd <- function(xt, tt=NULL, tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1", boundary="periodic", 
                sm="none", smlevels=c(1), spar=NULL, alpha=NULL, 
                check=FALSE, max.imf=10, plot.imf=FALSE, interm=NULL, weight=NULL) {

    if(is.ts(xt))
        xt <- as.numeric(xt) 

    if(is.null(tt)) tt <- 1:length(xt)
        
    if(is.null(interm) || all(interm <= 0)) intermtest <- FALSE else intermtest <- TRUE

    if (sm == "spline" || sm == "kernel" || sm == "locfit") { #|| sm == "quantile"
        if (is.null(spar) || spar == 0) stop("Provide the smoothing parameter.\n")
        else if (length(spar) == 1)
            spar <- rep(spar, length(smlevels))  
        else if (length(smlevels) != length(spar)) stop("Provide the smoothing parameter.\n")
    #} else if (sm == "quantile") {
    #    if (is.null(alpha)) stop("Provide the quantile for quantile regression.\n")
    #    else if (length(alpha) == 1)
    #        alpha <- rep(alpha, length(smlevels))  
    #    else if (length(smlevels) != length(alpha)) stop("Provide the quantile for quantile regression.\n")
    }          

    ndata <- length(xt); ndatam1 <- ndata - 1
    residue <- xt; rangext <- range(residue)
    imf <- NULL
    
    j <- 1
    
    repeat {
        if (j > max.imf) break
        if ((any(j == smlevels)) & (sm == "spline" || sm == "kernel" || sm == "locfit")) #|| sm == "quantile" 
            tmp <- extractimf(residue, tt, tol, max.sift, 
                            stoprule=stoprule, boundary=boundary, sm=sm, spar=spar[j], alpha=alpha[j], check=check, weight=NULL) else 
            tmp <- extractimf(residue, tt, tol, max.sift, 
                            stoprule=stoprule, boundary=boundary, sm="none", check=check) 

        if (is.null(tmp$imf))
            break
        
        if(plot.imf) {
            plot(tt, residue, type="l", xlab="", ylab="",
                main=paste(j-1, "-th residue=", j, "-th imf+", j, "-th residue", sep="")); abline(h=0)    
        }
        
        if(intermtest && length(interm) >= j && interm[j] > 0) {  
            tmpimf <- tmp$imf
            tmpresidue <- tmp$residue

            tmpinterm <- extrema(tmpimf, ndata, ndatam1)
            tmpncross <- tmpinterm$ncross
            zerocross <- as.numeric(round(apply(tmpinterm$cross, 1, mean)))

            if(abs(tt[zerocross[3]] - tt[zerocross[1]]) > interm[j]) {
                tmpresidue[1:zerocross[3]] <- tmpresidue[1:zerocross[3]] + tmpimf[1:zerocross[3]]                     
                tmpimf[1:zerocross[3]] <- 0            
            }
            
            for (k in seq(3, tmpncross-3, by=2))         
                if(abs(tt[zerocross[k+2]] - tt[zerocross[k]]) > interm[j]) {
                    tmpresidue[zerocross[k]:zerocross[k+2]] <- tmpresidue[zerocross[k]:zerocross[k+2]] + 
                                tmpimf[zerocross[k]:zerocross[k+2]]            
                    tmpimf[zerocross[k]:zerocross[k+2]] <- 0
                }
            
            if(!(tmpncross %% 2)) {               
                if(abs(tt[zerocross[tmpncross]] -tt[zerocross[tmpncross-1]]) > interm[j]/2) {       
                    tmpresidue[zerocross[tmpncross-1]:ndata] <- 
                            tmpresidue[zerocross[tmpncross-1]:ndata] + tmpimf[zerocross[tmpncross-1]:ndata]                 
                    tmpimf[zerocross[tmpncross-1]:ndata] <- 0            
                } 
            } else {
                if(abs(tt[zerocross[tmpncross]] - tt[zerocross[tmpncross-2]]) > interm[j])       
                    tmpresidue[zerocross[tmpncross-2]:ndata] <- 
                            tmpresidue[zerocross[tmpncross-2]:ndata] + tmpimf[zerocross[tmpncross-2]:ndata]                 
                    tmpimf[zerocross[tmpncross-2]:ndata] <- 0              
            }
            
            tmp$imf <- tmpimf
            tmp$residue <- tmpresidue            
        }
        
        imf <- cbind(imf, tmp$imf)     
        residue <- tmp$residue     
             
        if(plot.imf) {
            plot(tt, imf[,j], type="l", xlab="", ylab="", 
                main=paste(j, "-th imf", sep="")); abline(h=0)
            plot(tt, residue, type="l", xlab="", ylab="",
                main=paste(j, "-th residue", sep="")); abline(h=0); locator(1)
        }
        
        j <- j+1
    }
    list(imf=imf, residue=residue, nimf=j-1)
}  

   
extrema <- function(y, ndata = length(y), ndatam1 = ndata - 1) {

    minindex <- maxindex <- NULL; nextreme <- 0; cross <- NULL; ncross <- 0 
          
    z1 <- sign(diff(y))
    index1 <- seq(1, ndatam1)[z1 != 0]; z1 <- z1[z1 != 0]  

    if (!(is.null(index1) || all(z1==1) || all(z1==-1))) {

        index1 <- index1[c(z1[-length(z1)] != z1[-1], FALSE)] + 1 
        z1 <- z1[c(z1[-length(z1)] != z1[-1], FALSE)]  
        
        nextreme <- length(index1)
             
        if(nextreme >= 2)
            for(i in 1:(nextreme-1)) {
                tmpindex <- index1[i]:(index1[i+1]-1)
                if(z1[i] > 0) {
                    tmpindex <- tmpindex[y[index1[i]] == y[tmpindex]]
                    maxindex <- rbind(maxindex, c(min(tmpindex), max(tmpindex)))
                } else {
                    tmpindex <- tmpindex[y[index1[i]] == y[tmpindex]]
                    minindex <- rbind(minindex, c(min(tmpindex), max(tmpindex)))
                }     
            } 
            
        tmpindex <- index1[nextreme]:ndatam1  
        if(z1[nextreme] > 0) {
            tmpindex <- tmpindex[y[index1[nextreme]] == y[tmpindex]]
            maxindex <- rbind(maxindex, c(min(tmpindex), max(tmpindex)))
        } else {
            tmpindex <- tmpindex[y[index1[nextreme]] == y[tmpindex]]
            minindex <- rbind(minindex, c(min(tmpindex), max(tmpindex)))
        }  
      
        ### Finding the index of zero crossing  
 
        if (!(all(sign(y) >= 0) || all(sign(y) <= 0) || all(sign(y) == 0))) {
            index1 <- c(1, index1)
            for (i in 1:nextreme) {
                if (y[index1[i]] == 0) {
                    tmp <- c(index1[i]:index1[i+1])[y[index1[i]:index1[i+1]] == 0]
                    cross <- rbind(cross, c(min(tmp), max(tmp)))                 
                } else
                if (y[index1[i]] * y[index1[i+1]] < 0) {
                    tmp <- min(c(index1[i]:index1[i+1])[y[index1[i]] * y[index1[i]:index1[i+1]] <= 0])
                    if (y[tmp] == 0) {
                        tmp <- c(tmp:index1[i+1])[y[tmp:index1[i+1]] == 0]
                        cross <- rbind(cross, c(min(tmp), max(tmp))) 
                    } else 
                    cross <- rbind(cross, c(tmp-1, tmp)) 
                }
            }
            #if (y[ndata] == 0) {
            #    tmp <- c(index1[nextreme+1]:ndata)[y[index1[nextreme+1]:ndata] == 0]
            #    cross <- rbind(cross, c(min(tmp), max(tmp)))         
            #} else
            if (any(y[index1[nextreme+1]] * y[index1[nextreme+1]:ndata] <= 0)) {
                tmp <- min(c(index1[nextreme+1]:ndata)[y[index1[nextreme+1]] * y[index1[nextreme+1]:ndata] <= 0])
                if (y[tmp] == 0) {
                    tmp <- c(tmp:ndata)[y[tmp:ndata] == 0]
                    cross <- rbind(cross, c(min(tmp), max(tmp))) 
                } else
                cross <- rbind(cross, c(tmp-1, tmp))
            }
            ncross <- nrow(cross)        
        }
    } 
    
    list(minindex=minindex, maxindex=maxindex, nextreme=nextreme, cross=cross, ncross=ncross)
}    

emddenoise <- function(
xt, tt=NULL, 
cv.index, cv.level, cv.tol=0.1^3, cv.maxiter=20, by.imf=FALSE,
emd.tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1", boundary="periodic", max.imf=10)
{
### Golden section search.

    if(is.ts(xt))
        xt <- as.numeric(xt) 

    if(is.null(tt)) tt <- 1:length(xt)
        
    ndata <- length(xt)
    cv.kfold <- nrow(cv.index)

    tmpemd <- emd(xt, tt, emd.tol, max.sift, stoprule, boundary, sm="none", check=FALSE, max.imf=max.imf, plot.imf=FALSE)  
    tmpnimf <- tmpemd$nimf
        
    cv.ndim <- min(cv.level, tmpnimf)
    
    R <- (sqrt(5)-1)/2 #0.61803399000000003
    C <- 1 - R   
if(by.imf) {
    lambda <- matrix(0, 4, cv.ndim)  
    lambda.range <- 10 * apply(as.matrix(tmpemd$imf[,1:cv.ndim]), 2, function(t) {sqrt(sum(t^2))/sqrt(ndata)})
    
    lambda[1, ] <- rep(0, cv.ndim)
    lambda[4, ] <- lambda.range    
    lambda[2, ] <- lambda[4, ]/2
    lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])    

    perr <- lambdaconv <- NULL 
    optlambda <- lambda[3, ]
    j <- 0
    repeat{

        for (i in 1:cv.ndim) {                                
            optlambda[i] <- lambda[2, i]
            predxt <- tmpxt <- xt
            for (k in 1:cv.kfold) {                
                cvemd <- emd(xt[-cv.index[k,]], tt[-cv.index[k,]], emd.tol, max.sift, 
                                stoprule, boundary, sm="none", check=FALSE, max.imf=max.imf, plot.imf=FALSE)  
                for (m in 1:cv.ndim)
                    cvemd$imf[,m] <- cvemd$imf[,m] * (abs(cvemd$imf[,m]) > optlambda[m])
                tmpxt[-cv.index[k,]] <- apply(cvemd$imf, 1, sum) + cvemd$residue
        
                tmpindex1 <- cv.index[k,] - 1
                if(tmpindex1[1] == 0) tmpindex1 <- c(tmpindex1[2], tmpindex1[-1])  

                tmpindex2 <- cv.index[k,] + 1
                if(tmpindex2[length(tmpindex2)] == ndata+1) 
                    tmpindex2 <- c(tmpindex2[-length(tmpindex2)], tmpindex2[length(tmpindex2)-1])  

                predxt[cv.index[k,]] <- (tmpxt[tmpindex1] + tmpxt[tmpindex2]) / 2
            }
            f2 <- mean((predxt - xt)^2)
                
            optlambda[i] <- lambda[3, i]
            predxt <- tmpxt <- xt
            for (k in 1:cv.kfold) {                
                cvemd <- emd(xt[-cv.index[k,]], tt[-cv.index[k,]], emd.tol, max.sift, 
                                stoprule, boundary, sm="none", check=FALSE, max.imf=max.imf, plot.imf=FALSE)  
                for (m in 1:cv.ndim)
                    cvemd$imf[,m] <- cvemd$imf[,m] * (abs(cvemd$imf[,m]) > optlambda[m])
                tmpxt[-cv.index[k,]] <- apply(cvemd$imf, 1, sum) + cvemd$residue
        
                tmpindex1 <- cv.index[k,] - 1
                if(tmpindex1[1] == 0) tmpindex1 <- c(tmpindex1[2], tmpindex1[-1])  

                tmpindex2 <- cv.index[k,] + 1
                if(tmpindex2[length(tmpindex2)] == ndata+1) 
                    tmpindex2 <- c(tmpindex2[-length(tmpindex2)], tmpindex2[length(tmpindex2)-1])  

                predxt[cv.index[k,]] <- (tmpxt[tmpindex1] + tmpxt[tmpindex2]) / 2
            }
            f3 <- mean((predxt - xt)^2)
                        
            if(f3 < f2) {
                optlambda[i] <- lambda[3, i] 
                optf <- f3
                lambda[1, i] <- lambda[2, i]
                lambda[2, i] <- lambda[3, i]
                lambda[3, i] <- R * lambda[2, i] + C * lambda[4, i]

            }
            else {
                optlambda[i] <- lambda[2, i]
                optf <- f2
                lambda[4, i] <- lambda[3, i]
                lambda[3, i] <- lambda[2, i]
                lambda[2, i] <- R * lambda[3, i] + C * lambda[1, i]
            } 
            perr <- c(perr, optf)  
            lambdaconv <- rbind(lambdaconv, optlambda)
        }
        
        stopping <- NULL
        for (i in 1:cv.ndim) 
            stopping <- c(stopping, abs(lambda[4, i] - lambda[1, i]) / (abs(lambda[2, i]) + abs(lambda[3, i])))
#        if (all(stopping < cv.tol) || abs(perr[cv.ndim*(j+1)+1]-perr[cv.ndim*j+1])/perr[cv.ndim*j+1] < cv.tol || j > cv.maxiter) break 

        if (all(stopping < cv.tol) || j > cv.maxiter) break 

        j <- j + 1 
    }
} else {
    lambda <- matrix(0, 4, 1)  
    lambda.range <- 10 * apply(as.matrix(tmpemd$imf[,1]), 2, function(t) {sqrt(sum(t^2))/sqrt(ndata)})
    
    lambda[1, ] <- rep(0, 1)
    lambda[4, ] <- lambda.range    
    lambda[2, ] <- lambda[4, ]/2
    lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])    

    perr <- lambdaconv <- NULL 
    optlambda <- lambda[3, ]
    j <- 0
    repeat{

#        for (i in 1:cv.ndim) {                                
            optlambda[1] <- lambda[2, 1]
            predxt <- tmpxt <- xt
            for (k in 1:cv.kfold) {                
                cvemd <- emd(xt[-cv.index[k,]], tt[-cv.index[k,]], emd.tol, max.sift, 
                                stoprule, boundary, sm="none", check=FALSE, max.imf=max.imf, plot.imf=FALSE)  
                for (m in 1:cv.ndim)
                    cvemd$imf[,m] <- cvemd$imf[,m] * (abs(cvemd$imf[,m]) > optlambda[1])
                tmpxt[-cv.index[k,]] <- apply(cvemd$imf, 1, sum) + cvemd$residue
        
                tmpindex1 <- cv.index[k,] - 1
                if(tmpindex1[1] == 0) tmpindex1 <- c(tmpindex1[2], tmpindex1[-1])  

                tmpindex2 <- cv.index[k,] + 1
                if(tmpindex2[length(tmpindex2)] == ndata+1) 
                    tmpindex2 <- c(tmpindex2[-length(tmpindex2)], tmpindex2[length(tmpindex2)-1])  

                predxt[cv.index[k,]] <- (tmpxt[tmpindex1] + tmpxt[tmpindex2]) / 2
            }
            f2 <- mean((predxt - xt)^2)
                
            optlambda[1] <- lambda[3, 1]
            predxt <- tmpxt <- xt
            for (k in 1:cv.kfold) {                
                cvemd <- emd(xt[-cv.index[k,]], tt[-cv.index[k,]], emd.tol, max.sift, 
                                stoprule, boundary, sm="none", check=FALSE, max.imf=max.imf, plot.imf=FALSE)  
                for (m in 1:cv.ndim)
                    cvemd$imf[,m] <- cvemd$imf[,m] * (abs(cvemd$imf[,m]) > optlambda[1])
                tmpxt[-cv.index[k,]] <- apply(cvemd$imf, 1, sum) + cvemd$residue
        
                tmpindex1 <- cv.index[k,] - 1
                if(tmpindex1[1] == 0) tmpindex1 <- c(tmpindex1[2], tmpindex1[-1])  

                tmpindex2 <- cv.index[k,] + 1
                if(tmpindex2[length(tmpindex2)] == ndata+1) 
                    tmpindex2 <- c(tmpindex2[-length(tmpindex2)], tmpindex2[length(tmpindex2)-1])  

                predxt[cv.index[k,]] <- (tmpxt[tmpindex1] + tmpxt[tmpindex2]) / 2
            }
            f3 <- mean((predxt - xt)^2)
                        
            if(f3 < f2) {
                optlambda[1] <- lambda[3, 1] 
                optf <- f3
                lambda[1, 1] <- lambda[2, 1]
                lambda[2, 1] <- lambda[3, 1]
                lambda[3, 1] <- R * lambda[2, 1] + C * lambda[4, 1]

            }
            else {
                optlambda[1] <- lambda[2, 1]
                optf <- f2
                lambda[4, 1] <- lambda[3, 1]
                lambda[3, 1] <- lambda[2, 1]
                lambda[2, 1] <- R * lambda[3, 1] + C * lambda[1, 1]
            } 
            perr <- c(perr, optf)  
            lambdaconv <- rbind(lambdaconv, optlambda)
#        }
        
        stopping <- NULL
#        for (i in 1:cv.ndim) 
            stopping <- c(stopping, abs(lambda[4, 1] - lambda[1, 1]) / (abs(lambda[2, 1]) + abs(lambda[3, 1])))
#        if (all(stopping < cv.tol) || abs(perr[cv.ndim*(j+1)+1]-perr[cv.ndim*j+1])/perr[cv.ndim*j+1] < cv.tol || j > cv.maxiter) break 

        if (all(stopping < cv.tol) || j > cv.maxiter) break 

        j <- j + 1 
    }
}   

    for (m in 1:cv.ndim) 
        if(by.imf) 
            tmpemd$imf[,m] <- tmpemd$imf[,m] * (abs(tmpemd$imf[,m]) > optlambda[m])
        else 
            tmpemd$imf[,m] <- tmpemd$imf[,m] * (abs(tmpemd$imf[,m]) > optlambda)
                    
    dxt <- apply(tmpemd$imf, 1, sum) + tmpemd$residue
    tmpemd$imf <- ts(tmpemd$imf, min(tt))    
    tmpemd$residue <- ts(tmpemd$residue, min(tt))
    list(dxt=dxt, optlambda=optlambda, lambdaconv=lambdaconv, perr=perr, demd=tmpemd, niter=j)  
}

cvtype <-
function (n, cv.bsize = 1, cv.kfold, cv.random = FALSE) 
{
    if (n < cv.bsize * cv.kfold) 
        stop("Block size or no. of fold is too large.")
    if (n <= 0) 
        stop("The number of data must be greater than 0.")
    cv.index <- tmp.index <- NULL
    if (cv.random) {
        cv.nblock <- trunc(n/(cv.bsize * cv.kfold))
        for (k in 1:cv.kfold) tmp.index <- rbind(tmp.index, sort(sample(seq(sample(1:cv.bsize, 
            1), n - cv.bsize + 1, by = cv.bsize), cv.nblock)))
    }
    else {
        for (k in 1:cv.kfold) tmp.index <- rbind(tmp.index, seq(1, 
            n - cv.bsize * (cv.kfold - 1), by = cv.bsize * cv.kfold) + 
            cv.bsize * (k - 1))
    }
    if (cv.bsize > 1) {
        cv.index <- tmp.index
        for (i in 1:(cv.bsize - 1)) cv.index <- cbind(cv.index, 
            tmp.index + i)
        cv.index <- t(apply(cv.index, 1, sort))
    }
    else {
        cv.index <- tmp.index
    }
    list(cv.index = cv.index)
}


emd.pred <- function(varpred, trendpred, ci=0.95, figure = TRUE) 
{
    n.ahead <- length(varpred$fcst[[1]][, 1])
    fcst <- lower <- upper <- rep(0, n.ahead)
    for (i in 1:length(varpred$fcst)) {
        fcst <- fcst + varpred$fcst[[i]][, 1]
        lower <- lower + varpred$fcst[[i]][, 2]
        upper <- upper + varpred$fcst[[i]][, 3]
    }
    fcst <- fcst + trendpred$fit
    lower <- lower + trendpred$fit - qt((1+ci)/2, trendpred$df) * trendpred$se.fit
    upper <- upper + trendpred$fit + qt((1+ci)/2, trendpred$df) * trendpred$se.fit

    if (figure) {
        yrange <- range(c(fcst, lower, upper))
        par(mar=0.1+c(4,4,1.0,1))
        plot(1:n.ahead, fcst, type = "l", lwd = 2, col = "green", 
            ylim = yrange, ylab = "Prediction", xlab = "Time horizon")
        lines(1:n.ahead, upper, lty = 2, col = "red")
        lines(1:n.ahead, lower, lty = 2, col = "blue")
    }
    list(fcst = fcst, lower = lower, upper = upper)
}

semd <- function(xt, tt=NULL, cv.kfold, cv.tol=0.1^1, cv.maxiter=20, 
emd.tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1", boundary="periodic", 
smlevels=1, max.imf=10)
{
### Golden section search.
    sm <- "spline"
    
    if(is.ts(xt))
        xt <- as.numeric(xt) 

    if(is.null(tt)) tt <- 1:length(xt)
        
    ndata <- length(xt)

    cv.index <- cvtype(ndata, 1, cv.kfold)$cv.index

    cv.ndim <- min(length(smlevels), max.imf)

    tmp <- sreg(tt, xt)
    spar <- tmp$lambda
   
    if (cv.ndim != 1) {
        tmplambda <- tmp$lambda
        tmpresidue <- xt
        
        for (i in 1:(cv.ndim-1)) {
            tmpemd <- emd(tmpresidue, tt, emd.tol, max.sift, stoprule, boundary, sm=sm, smlevels=1, spar=tmplambda, alpha=NULL, check=FALSE, max.imf=1, plot.imf=FALSE, interm=NULL, weight=NULL) 

            if (is.null(tmpemd$imf)) 
                break
        
            tmpresidue <- tmpemd$residue
            tmplambda <- sreg(tt, tmpresidue)$lambda
            spar <- c(spar, tmplambda)
        }    ### impute by local mean
        cv.ndim <- min(length(spar), max.imf)
    }

    if (cv.ndim == 1) by.imf <- FALSE
    else by.imf <- TRUE

    imputext <- cvimpute.by.mean(y=xt, impute.index=cv.index)$yimpute
    
    R <- (sqrt(5)-1)/2 #0.61803399000000003
    C <- 1 - R
    lambda.range <- c(0.1, 80)
    
    if(by.imf) {
        lambda <- matrix(0, 4, cv.ndim)  
    
        lambda[1, ] <- lambda.range[1] * spar
        lambda[4, ] <- lambda.range[2] * spar   
        lambda[2, ] <- lambda[4, ]/2
        lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])    

        perr <- lambdaconv <- NULL 
        optlambda <- lambda[3, ]

        j <- 0
        repeat{

            for (i in 1:cv.ndim) {                                
                optlambda[i] <- lambda[2, i]
                predxt <- xt
                for (k in 1:cv.kfold) {
                    tmpxt <- xt
                    tmpxt[cv.index[k,]] <- imputext[cv.index[k,]]
                    cvemd <- emd(tmpxt, tt, emd.tol, max.sift, 
                                stoprule, boundary, sm=sm, smlevels=1:cv.ndim, spar=optlambda, alpha=NULL, check=FALSE, max.imf=cv.ndim, plot.imf=FALSE, interm=NULL, weight=NULL)

                    predxt[cv.index[k,]] <- apply(cvemd$imf[,-(1:cv.ndim)], 1, sum)[cv.index[k,]] + cvemd$residue[cv.index[k,]]
                }
                f2 <- mean((predxt - xt)^2)
                
                optlambda[i] <- lambda[3, i]
                predxt <- xt
                for (k in 1:cv.kfold) {
                    tmpxt <- xt
                    tmpxt[cv.index[k,]] <- imputext[cv.index[k,]]
                    cvemd <- emd(tmpxt, tt, emd.tol, max.sift, 
                                stoprule, boundary, sm=sm, smlevels=1:cv.ndim, spar=optlambda, alpha=NULL, check=FALSE, max.imf=cv.ndim, plot.imf=FALSE, interm=NULL, weight=NULL)

                    predxt[cv.index[k,]] <- apply(cvemd$imf[,-(1:cv.ndim)], 1, sum)[cv.index[k,]] + cvemd$residue[cv.index[k,]]
                }
                f3 <- mean((predxt - xt)^2)
                        
                if(f3 < f2) {
                    optlambda[i] <- lambda[3, i] 
                    optf <- f3
                    lambda[1, i] <- lambda[2, i]
                    lambda[2, i] <- lambda[3, i]
                    lambda[3, i] <- R * lambda[2, i] + C * lambda[4, i]
                }
                else {
                    optlambda[i] <- lambda[2, i]
                    optf <- f2
                    lambda[4, i] <- lambda[3, i]
                    lambda[3, i] <- lambda[2, i]
                    lambda[2, i] <- R * lambda[3, i] + C * lambda[1, i]
                } 
                perr <- c(perr, optf)  
                lambdaconv <- rbind(lambdaconv, optlambda)
            }
        
            stopping <- NULL
            for (i in 1:cv.ndim) 
                stopping <- c(stopping, abs(lambda[4, i] - lambda[1, i]) / (abs(lambda[2, i]) + abs(lambda[3, i])))

            #if (all(stopping < cv.tol) || abs(perr[cv.ndim*(j+1)+1]-perr[cv.ndim*j+1])/perr[cv.ndim*j+1] < cv.tol || j > cv.maxiter) break 
            if (all(stopping < cv.tol) || j > cv.maxiter) break 

            j <- j + 1 
        }
    } else {
        lambda <- matrix(0, 4, 1)  
    
        lambda[1, ] <- lambda.range[1] * spar
        lambda[4, ] <- lambda.range[2] * spar
        lambda[2, ] <- lambda[4, ]/2
        lambda[3, ] <- lambda[2, ] + C * (lambda[4, ] - lambda[2, ])    

        perr <- lambdaconv <- NULL 
        optlambda <- lambda[3, ]
        
        j <- 0

        repeat{
            optlambda[1] <- lambda[2, 1]
            predxt <- xt
            for (k in 1:cv.kfold) {
                tmpxt <- xt
                tmpxt[cv.index[k,]] <- imputext[cv.index[k,]]
                cvemd <- emd(tmpxt, tt, emd.tol, max.sift, 
                                stoprule, boundary, sm=sm, smlevels=1, spar=optlambda[1], alpha=NULL, check=FALSE, max.imf=1, plot.imf=FALSE, interm=NULL, weight=NULL)

                predxt[cv.index[k,]] <- cvemd$residue[cv.index[k,]]
            }
            f2 <- mean((predxt - xt)^2)
                
            optlambda[1] <- lambda[3, 1]
            predxt <- xt
            for (k in 1:cv.kfold) {
                tmpxt <- xt
                tmpxt[cv.index[k,]] <- imputext[cv.index[k,]]
                cvemd <- emd(tmpxt, tt, emd.tol, max.sift, 
                                stoprule, boundary, sm=sm, smlevels=1, spar=optlambda[1], alpha=NULL, check=FALSE, max.imf=1, plot.imf=FALSE, interm=NULL, weight=NULL)

                predxt[cv.index[k,]] <- cvemd$residue[cv.index[k,]]
            }
            f3 <- mean((predxt - xt)^2)
                        
            if(f3 < f2) {
                optlambda[1] <- lambda[3, 1] 
                optf <- f3
                lambda[1, 1] <- lambda[2, 1]
                lambda[2, 1] <- lambda[3, 1]
                lambda[3, 1] <- R * lambda[2, 1] + C * lambda[4, 1]
            } else {
                optlambda[1] <- lambda[2, 1]
                optf <- f2
                lambda[4, 1] <- lambda[3, 1]
                lambda[3, 1] <- lambda[2, 1]
                lambda[2, 1] <- R * lambda[3, 1] + C * lambda[1, 1]
            } 
            perr <- c(perr, optf)  
            lambdaconv <- rbind(lambdaconv, optlambda)
        
            stopping <- abs(lambda[4, 1] - lambda[1, 1]) / (abs(lambda[2, 1]) + abs(lambda[3, 1]))
            if (stopping < cv.tol || j > cv.maxiter) break
            #if (j >= 1)
            #if (abs(perr[j+1]-perr[j])/perr[j] < cv.tol || j > cv.maxiter) break 
            
            j <- j + 1 
        }
    }   

    out <- emd(xt, tt, emd.tol, max.sift, stoprule, boundary, sm=sm, smlevels=1:length(optlambda), spar=optlambda, alpha=NULL, check=FALSE, max.imf, plot.imf=FALSE, interm=NULL, weight=NULL)
    
    list(imf=out$imf, residue=out$residue, nimf=out$nimf, optlambda=optlambda, lambdaconv=lambdaconv, perr=perr)  
}


cvimpute.by.mean <- function(y, impute.index)
{
    n <- length(y)
    impute.nrep <- nrow(impute.index)
    yimpute <- y

    for (k in 1:impute.nrep) {
        tmpyimpute <- y
        tmpyimpute[impute.index[k,]] <- NA
        tmpyimpute <- c(y[2], tmpyimpute, y[n-1])
        tmpyimpute <- apply(embed(tmpyimpute, 3), 1, mean, na.rm=TRUE) 
        yimpute[impute.index[k,]] <- tmpyimpute[impute.index[k,]]
    } # End of k-th impute loop

    return(list(yimpute=yimpute))
}
