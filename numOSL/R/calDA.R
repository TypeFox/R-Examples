#####
calDA <- 
function(dose, minGrainSize, maxGrainSize,
         Ucontent, Thcontent, Kcontent, Wct, depth, altitude, 
         latitude, longitude, bulkDensity=2.5, alphaValue=0.03, 
         nsim=10000, rdcf=0.05, rba=0.05, plot=TRUE)   {
    UseMethod("calDA")
} #
### 2015.04.02.
calDA.default <- 
function(dose, minGrainSize, maxGrainSize,
         Ucontent, Thcontent, Kcontent, Wct, depth, altitude, 
         latitude, longitude, bulkDensity=2.5, alphaValue=0.03, 
         nsim=10000, rdcf=0.05, rba=0.05, plot=TRUE)   {
    ### Stop if not.
    stopifnot(length(dose)==2L, is.numeric(dose), all(dose>0.0),
              length(minGrainSize)==1L, is.numeric(minGrainSize), minGrainSize>0.0,
              length(maxGrainSize)==1L, is.numeric(maxGrainSize),
              minGrainSize<maxGrainSize, maxGrainSize <=300.0, 
              length(Ucontent)==2L, is.numeric(Ucontent), all(Ucontent>0.0),
              length(Thcontent)==2L, is.numeric(Thcontent), all(Thcontent>0.0),
              length(Kcontent)==2L, is.numeric(Kcontent), all(Kcontent>0.0),
              length(Wct)==2L, is.numeric(Wct), all(Wct>0.0), 
              length(depth) %in% c(1L,2L), is.numeric(depth), all(depth>0.0),
              length(altitude) %in% c(1L,2L), is.numeric(altitude), all(altitude>0.0),
              length(latitude) %in% c(1L,2L), is.numeric(latitude), all(latitude>0.0),
              length(longitude) %in% c(1L,2L), is.numeric(longitude), all(longitude>0.0),
              length(bulkDensity) %in% c(1L,2L), is.numeric(bulkDensity), all(bulkDensity>0.0),
              length(alphaValue) %in% c(1L,2L), is.numeric(alphaValue), all(alphaValue>=0.0),
              length(nsim)==1L, is.numeric(nsim), nsim>=100L, nsim<=50000L,
              length(rdcf)==1L, is.numeric(rdcf), rdcf>=0.0, rdcf<=0.1,
              length(rba)==1L, is.numeric(rba), rba>=0.0, rba<=0.1,
              length(plot)==1L, is.logical(plot))
    ###
    ### Inner function used for calculating 
    ### attenuation of beta radiation.
    #######################################
    BetaAttenuation <- function(GrainSize1, GrainSize2, rba)  {
            ###
            meanGrainSize <- (GrainSize1+GrainSize2) / 2.0
            ###
            Uabsorption <- 0.0004413559322*meanGrainSize + 0.01877966102
            Thabsorption <- 0.000513220339 *meanGrainSize + 0.03110169492
             Kabsorption <- 0.0003681355932 *meanGrainSize + 0.001677966102
             ###
             if (rba>0)  {
                 Uabsorption <- rnorm(n=1L, mean=Uabsorption, sd=rba*Uabsorption)
                 Thabsorption <- rnorm(n=1L, mean=Thabsorption, sd=rba*Thabsorption)
                 Kabsorption <- rnorm(n=1L, mean=Kabsorption, sd=rba*Kabsorption)
             } # end if.
             out1 <- c("Uabsorption"=Uabsorption, 
                       "Thabsorption"=Thabsorption,
                       "Kabsorption"=Kabsorption)
             return(out1)
     } # end function  BetaAttenuation.
     ### 
     ### Inner function used for calculating 
     ### Beta and Gamma dose rates.
     #######################################
     BetaGammaDoseRate <- function(minGrainSize, maxGrainSize,
                          Ucontent, Thcontent, Kcontent, rdcf, rba)  {
            ### 
            if (rdcf<=0)  {
                UBeta <- 0.146*Ucontent
                UGamma <- 0.113*Ucontent
                ThBeta <- 0.0273*Thcontent
                ThGamma <- 0.0476*Thcontent
                KBeta <- 0.782*Kcontent
                KGamma <- 0.243*Kcontent
            } else {
                UBeta <- rnorm(n=1L, mean=0.146, sd=rdcf*0.146)*Ucontent
                UGamma <- rnorm(n=1L,mean=0.113,sd=rdcf*0.113)*Ucontent
                ThBeta <- rnorm(n=1L, mean=0.0273, sd=rdcf*0.0273)*Thcontent
                ThGamma <- rnorm(n=1L, mean=0.0476, sd=rdcf*0.0476)*Thcontent
                KBeta <- rnorm(n=1L, mean=0.782, sd=rdcf*0.782)*Kcontent
                KGamma <- rnorm(n=1L, mean=0.243, sd=rdcf*0.243)*Kcontent
            } # end if.
            ###
            AttenuationFactors <- 1.0 - BetaAttenuation(minGrainSize, maxGrainSize, rba)
            betaDoseRate <- UBeta*AttenuationFactors[["Uabsorption"]] +
                            ThBeta*AttenuationFactors[["Thabsorption"]] +
                            KBeta*AttenuationFactors[["Kabsorption"]] 
            ###
            gammaDoseRate <- UGamma +ThGamma +KGamma
            ###
            out2 <- c("betaDoseRate"=betaDoseRate,
                      "gammaDoseRate"=gammaDoseRate)
            return(out2)
    } # end function BetaGammaDoseRate.
    ### 
    ### Inner function used for 
    ### calculating Comic dose rate.
    #################################
    ComicDoseRate <- function(depth, altitude, 
                     latitude, longitude, bulkDensity)  {
            ###
            lambda <- asin( 0.203*cos(latitude*pi/180.0)*
                      cos( (longitude-291.0)*pi/180.0)+
                      0.979*sin(latitude*pi/180.0) )*180.0/pi
            ###
            D0 <- (6072.0/(((((depth*bulkDensity)+11.6)^1.68)+75.0)*
                  ((depth*bulkDensity)+212.0)))*
                  exp(-0.00055*(depth*bulkDensity))
            ###
            if (lambda<=35.0)  {
                F <-  0.397889507 + 0.0004800551553*lambda - 
                      0.0001488550036*lambda^2L + 1.925045242e-07*lambda^3L
                ###
                 J <-  0.5192597637 + 0.00182911157*lambda + 
                       4.365380142e-05*lambda^2L + 2.839825082e-06*lambda^3L
                ###
                H <-  4.399084215 - 0.003717607061*lambda + 
                      2.541252939e-05*lambda^2L- 4.807177888e-06*lambda^3L
            } else {
                F <- 0.231
                J <- 0.761
                H <- 4.098
            } # end if. 
            ###
            comicDoseRate <- D0*(F+J*exp(altitude/H/1000.0))
            return(comicDoseRate)
    } # end function ComicDoseRate.
    ### 
    ### Inner function used for calculating
    ### a annual dose rate value.
    ########################################
    calDoseRate <- function(minGrainSize, maxGrainSize,
                   Ucontent, Thcontent, Kcontent, Wct, depth, altitude, 
                   latitude, longitude, bulkDensity, alphaValue, rdcf, rba)  {
            ### 
            bgRate <- BetaGammaDoseRate(minGrainSize, maxGrainSize,
                      Ucontent, Thcontent, Kcontent, rdcf, rba) 
            ###
            comicDoseRate <- ComicDoseRate (depth, altitude, 
                             latitude, longitude, bulkDensity) 
            ###
            betaDoseRate <- bgRate[["betaDoseRate"]] / (1.0+1.25*Wct)
            gammaDoseRate <- bgRate[["gammaDoseRate"]] /(1.0+1.14*Wct)
            ###
            totalDoseRate <- betaDoseRate +
                             gammaDoseRate +
                             comicDoseRate
            ###
            if (minGrainSize>=90.0)  alphaValue <- 0.0
            alphaDoseRate <- alphaValue*
                             (2.18*0.612*Ucontent+0.611*0.557*Thcontent)/
                             (1.0+1.5*Wct)
            totalDoseRate <- totalDoseRate + alphaDoseRate
            return(totalDoseRate)
    } # end function calDoseRate.
    ###
    ###
    DoseRate<- calDoseRate(minGrainSize, maxGrainSize,
               Ucontent[1L], Thcontent[1L], Kcontent[1L], Wct[1L], 
               depth[1L], altitude[1L], latitude[1L], longitude[1L], 
               bulkDensity[1L], alphaValue[1L], rdcf=-100, rba=-100)
    Age <- dose[1L] / DoseRate
    ### 
    Ages <- DoseRates <- vector(length=nsim)
    ###
    simdepth <- depth[1L] 
    simaltitude <- altitude[1L]
    simlatitude <-latitude[1L]
    simlongitude <- longitude[1L]
    simbulkDensity <- bulkDensity[1L]
    simalphaValue <- alphaValue[1L]
    ###
    for (i in seq(nsim))  {
        ###
        simdose <- rnorm(n=1L, mean=dose[1L], sd=dose[2L])
        simUcontent <- rnorm(n=1L, mean=Ucontent[1L], sd=Ucontent[2L])
        simThcontent <- rnorm(n=1L, mean=Thcontent[1L], sd=Thcontent[2L])
        simKcontent <- rnorm(n=1L, mean=Kcontent[1L], sd=Kcontent[2L])
        simWct <- rnorm(n=1L, mean=Wct[1L], sd=Wct[2L])
        ###                            
        if (length(depth)==2L)  {
            simdepth <- rnorm(n=1L, mean=depth[1L], sd=depth[2L])
        } # end if.
        if (length(altitude)==2L)  {
            simaltitude <- rnorm(n=1L, mean=altitude[1L], sd=altitude[2L])
        } # end if.
        if (length(latitude)==2L)  {
            simlatitude <- rnorm(n=1L, mean=latitude[1L], sd=latitude[2L])
        } # end if.
        if (length(longitude)==2L)  {
            simlongitude <- rnorm(n=1L, mean=longitude[1L], sd=longitude[2L])
        } # end if.
        if (length(bulkDensity)==2L)  {
            simbulkDensity <- rnorm(n=1L, mean=bulkDensity[1L], sd=bulkDensity[2L])
        } # end if.
        if (length(alphaValue)==2L)  {
            simalphaValue <- rnorm(n=1L, mean=alphaValue[1L], sd=alphaValue[2L])
        } # end if.
        ###
        DoseRates[i] <- calDoseRate(minGrainSize, maxGrainSize,
                        simUcontent, simThcontent, simKcontent, simWct, 
                        simdepth, simaltitude, simlatitude, simlongitude, 
                        simbulkDensity, simalphaValue, rdcf=rdcf, rba=rba)
        Ages[i] <- simdose / DoseRates[i]   
    } # end for.
    ###
    if (plot==TRUE)  {
        par(mfrow=c(1L,2L))
        dDoseRates <- density(DoseRates)
        plot(dDoseRates$x,dDoseRates$y, xlab="Total dose rate (Gy/ka)",
             ylab="Density", main="Simulated dose rate", type="n")
        xTicks<-axTicks(side=1L)
        maxYx<-dDoseRates$x[which.max(dDoseRates$y)]
        ###
        polygon(dDoseRates$x,dDoseRates$y, col="grey")
        rug(DoseRates, quiet=TRUE)
        ###
        legend(ifelse(maxYx>median(xTicks),"topleft","topright"), 
               legend=c(paste("Simuluation=", nsim, sep=""),
                        paste("DoseRate=",round(DoseRate,2L)," (Gy/ka)", sep=""),
                        paste("s(DoseRate)=",round(sd(DoseRates),2L)," (Gy/ka)",sep="")), 
                        cex=par("cex"), bty="n")
        ###
        ###
        dAges <- density(Ages)
        plot(dAges$x, dAges$y, xlab="Age (ka)", ylab="Density",
             main="Simulated age", type="n")
        xTicks<-axTicks(side=1L)
        maxYx<-dAges$x[which.max(dAges$y)]
        ###
        polygon(dAges$x, dAges$y,  col="grey")
        rug(Ages, quiet=TRUE)
        ###
        legend(ifelse(maxYx>median(xTicks),"topleft","topright"), 
               legend=c(paste("Simuluation=", nsim, sep=""),
                        paste("Age=",round(Age,2L)," (ka)", sep=""),
                        paste("s(Age)=",round(sd(Ages),2L)," (ka)",sep="")), 
                        cex=par("cex"), bty="n")
        ###
        par(mfrow=c(1L,1L))
    } # end if.
    ###
    sdDoseRate <- sd(DoseRates)
    sdAge <- sd(Ages)
    ###
    out <- matrix(nrow=2L, ncol=2L)
    out[1L,] <- c(DoseRate, sdDoseRate)
    out[2L,] <- c(Age, sdAge)
    rownames(out) <- c("DoseRate", "Age")
    colnames(out) <- c("Pars", "stdPars")
    return(out)
} # end function calDA.
#####
