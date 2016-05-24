# source( "C:/_RTOOLS/SWEAVE_WORK/SOIL_TEXTURES/rforge/pkg/soiltexture/R/text.transf.r" ) 


library( "drc" )


TT.check.ps.lim.Xm <- function(# Internal. Check the consistency between 'base.ps.lim' and 'dat.ps.lim'. 
### Check the consistency between 'base.ps.lim' and 'dat.ps.lim'. 
### 4 tests performed.

 base.ps.lim, 

 dat.ps.lim, 

 ps.lim.length=c(4,4)
### vector of 2 integers. Number of particle size classes + 1. c(base,dat)

##author<<Wei Shangguan

){  #
    # if( length( base.ps.lim ) != length( dat.ps.lim ) ) 
    # {   #
    #     stop( paste( 
    #         sep="", 
    #         "The length of the 'base' particle size classes limits must be equal to\n", 
    #         "the length of the 'dat' particle size classes limits.\n", 
    #         "Either check the 'base' particle size classes limits vector,\n", 
    #         "or check number of column in tri.data.\n"  
    #     ) )   #
    # }   #
    #
    if( length( base.ps.lim ) != ps.lim.length[1] ) 
    {   #
        stop( paste( 
            sep="", 
            "The length of the 'base' particle size classes limits must be equal to\n", 
            ps.lim.length[1], " (number of particle size classes+1; from ps min to ps.max)\n", 
            "Actual value: ", length( base.ps.lim ), ".\n", 
            "Either check the 'base' particle size classes limits,\n", 
            "or check number of column in tri.data.\n"  
        ) )   #
    }   #
    #
    if( length( dat.ps.lim ) != ps.lim.length[2] ) 
    {   #
        stop( paste( 
            sep="", 
            "The length of the 'dat' particle size classes limits must be equal to\n", 
            ps.lim.length[2], " (number of particle size classes +1; from ps min to ps.max)\n", 
            "Actual value: ", length( dat.ps.lim ), ".\n", 
            "Either check the 'dat' particle size classes limits,\n", 
            "or check number of column in tri.data.\n"  
        ) )   #
    }   #
    #
    if( base.ps.lim[1] != dat.ps.lim[1] ) 
    {   #
        stop( paste( 
            sep="", 
            "The first value of the 'dat' particle size classes limits must be equal to\n", 
            "the first value of the 'base' particle size classes limits.\n", 
            "Actual value, base: ", base.ps.lim[1], ", dat: ", dat.ps.lim[1]  
        ) )   #
    }   #
    #
    if( base.ps.lim[ps.lim.length[1]] != dat.ps.lim[ps.lim.length[2]] ) 
    {   #
        stop( paste( 
            sep="", 
            "The last value of the 'dat' particle size classes limits must be equal to\n", 
            "the last value of the 'base' particle size classes limits.\n", 
            "Actual value, base: ", base.ps.lim[ps.lim.length[1]], ", dat: ", dat.ps.lim[ps.lim.length[2]]  
        ) )   #
    }   #
    #
#    if( base.ps.lim[1] == 0 ) 
#    {   #
#        if( base.ps.lim[2] < dat.ps.lim[2] )
#        stop( paste( 
#            sep="", 
#            "When the 1st value of 'dat' and 'base' particle size classes limits is 0\n", 
#            "The 2nd value of the 'base' particle size classes limits must higher or equal to\n", 
#            "the 2nd value of the 'dat' particle size classes limits.\n"  
#        ) )   #
#    }   #
}   #






TT.text.transf.Xm <- function(# Transformations of a soil texture data table between 2 particle size systems (X classes), various methods.
### using various Particle Size Distribution (PSD) models including Anderson (AD), Fredlund4P (F4P), Fredlund3P (F3P),
### modified logistic growth (ML), Offset-Nonrenormalized Lognormal (ONL), Offset-Renormalized Lognormal (ORL),
### Skaggs (S), van Genuchten type(VG),van Genuchten modified(VGM), Weibull (W), Logarithm(L),
### Logistic growth (LG), Simple Lognormal (SL),Shiozawa and Compbell (SC).
### The performance of PSD models is influenced by many aspects like soil texture class,
### number and position (or closeness) of observation points, clay content etc.
### The latter four PSD models perform worse than the former ten.
### The AD, F4P, S, and W model is recommended for most of texture classes.
### And it will be even better to compare several different PSD models and using the results of the model
### with the minimum residual sum of squares (or other measures).
### Sometimes, the fitting will failed for the iteration is not converged or some errors happened.
### Transformation of a soil texture data table
### ('tri.data') from one
### particle size system ('dat.ps.lim') into another
### ('base.ps.lim'). No limit in the number of texture classes
### in the input and output texture tables. See TT.text.transf
### for transformation involving only 3 particle classes. 'tri.data'
### can only contain texture data.
### This function requires the "drc" package to be installed.
##author<<Wei Shangguan

tri.data,

base.ps.lim,

dat.ps.lim,

text.sum= NULL,

text.tol= NULL,

tri.sum.tst= NULL,

tri.pos.tst= NULL,

psdmodel= "AD",

omethod= "all",
### see optim for available methods, the default "all" is to run all methods and 
### choose the best results with minimum residual sum of squares (RSS).

tri.sum.norm=FALSE
###Weather the sum of is

){#
    TT.auto.set( set.par = FALSE )
    #
    TT.data.test.X(
        tri.data    = tri.data,
        text.sum    = text.sum,
        text.tol    = text.tol,
        tri.sum.tst = tri.sum.tst,
        tri.pos.tst = tri.pos.tst
    )   #
    #
    tri.data <- t(  apply(
        X       = tri.data,
        MARGIN  = 1,
        FUN     = function(X){
            cumsum(X)/100
        }   #
    )   )   #
    #
    ps.end   <- dim( tri.data )[2] + 1
    #
    TT.check.ps.lim.Xm(
        base.ps.lim     = base.ps.lim,
        dat.ps.lim      = dat.ps.lim,
        ps.lim.length   = c(length(base.ps.lim),ps.end)
    )   #
    #
    if( base.ps.lim[1] != 0 )
    {   #
        tri.data <- cbind(
            "ZERO" = rep(0,dim(tri.data)[1]),
            tri.data
        )   #
        #
        ps.start    <- 1
    }else{
        ps.start    <- 2
    }   #
    #Particle size diameter in micro-meters (will be converted in milli-meters)
    base.ps.lim <- base.ps.lim/1000
    dat.ps.lim <- dat.ps.lim/1000
#    base.ps.lim2 <- TT.dia2phi(base.ps.lim)
#    dat.ps.lim2  <- TT.dia2phi(dat.ps.lim)
    #
#    old.col.nm   <- colnames( tri.data )
    
    # # Added 2010/06/13 Julien Moeys # Removed on 2012/03/07 by Julien Moeys
    # if( !"drc" %in%  as.character( installed.packages()[,1] ) ) 
    # {   #
    #     stop( "The function 'TT.text.transf.Xm' needs the package 'drc'\n Please install it ( install.packages(\"drc\") )" ) 
    # }else{ 
    #     require( "drc" ) 
    # }   #
    
    require( "drc" ) 
    
    fitpsd <- function(
    y,
    xin,
    xout,
    psdmodel,
    method)
    {
        require( "drc" ) # Added 2010/08/11 by JM
        
        #default max and min of initial parameters
        maxspa1 <- 1
        minspa1 <- 0.1
        maxspa2 <- 1
        minspa2 <- 0.1
        #erf error function for ONL and ORL model
        erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
        #r0 for S model
        r0  <- xin[1]
        #dmax, dmin for W model
        dmax <- xin[length(xin)]
        dmin <- r0
        # Particle Size Distribution (PSD) models
        AD  <- function(dose, parm)
            {parm[, 1] + parm[, 2]*atan(parm[, 3]*log10(dose/parm[, 4]))}
        F4P <- function(dose, parm)
            {(1-(log(1+parm[,1]/dose)/log(1+parm[,1]/0.0001))^7)/(log(exp(1)+(parm[,2]/dose)^parm[,3]))^parm[,4]}
        F3P <- function(dose, parm)
            {(1-(log(1+0.001/dose)/log(1+0.001/0.0001))^7)/(log(exp(1)+(parm[,1]/dose)^parm[,2]))^parm[,3]}
        ML  <- function(dose, parm)
            {1/(1+parm[,1]*exp(-parm[,2]*dose^parm[,3]))}
        ONL <- function(dose, parm)
            {
                t   <- (-1)^(log(dose) >= parm[,1]+1)
                (1+t*erf((log(dose)+parm[,1])/parm[,2]*2^0.5))/2+parm[,3]
            }
        ORL <- function(dose, parm)
            {
                t   <- (-1)^(log(dose) >= parm[,1]+1)
                (1-parm[,3])*(1+t*erf((log(dose)+parm[,1])/parm[,2]*2^0.5))/2+parm[,3]
            }
        # no ability to predict the content below r0
        S   <- function(dose, parm)
            {1/(1+(1/y[1]-1)*exp(-parm[,1]*((dose-r0)/r0)^parm[,2]))}
        VG  <- function(dose, parm)
            {(1+(parm[,1]/dose)^parm[,2])^(1/parm[,2]-1)}
        # old form is right
        VGM <- function(dose, parm)
            {y[1]+(1-y[1])*(1+(parm[,1]*dose)^(-parm[,2]))^(-1/parm[,2]-1)}
        #This is the wrong form of VGM
#        VGM <- function(dose, parm)
#            {y[1]+(1-y[1])*(1+(parm[,1]*dose)^(-parm[,2]))^(1/parm[,2]-1)}
        # no ability to predict the content below dmin
        W   <- function(dose, parm)
            {parm[,3]+(1-parm[,3])*(1-exp(-parm[,1]*((dose-dmin)/(dmax-dmin))^parm[,2]))}
        L   <- function(dose, parm)
            {parm[,1]*log(dose)+parm[,2]}
        LG  <- function(dose, parm)
            {1/(1+parm[,1]*exp(-parm[,2]*dose))}
        SC  <- function(dose, parm)
            {
                t   <-(-1)^(log(dose) >= parm[,1]+1)
                (1-parm[,3])*(1+t*erf((log(dose)+parm[,1])/parm[,2]*2^0.5))/2+parm[,3]*(1+t*erf((log(dose)+1.96)/1*2^0.5))/2
            }
        SL  <- function(dose, parm)
            {
                t   <-(-1)^(log(dose) >= parm[,1]+1)
                (1+t*erf((log(dose)+parm[,1])/parm[,2]*2^0.5))/2
            }

        if( psdmodel == "AD" )
        {
            logi    <- AD
            pn      <- 4
            pname   <- c("f0", "b", "c", "r0")
        }
        else if ( psdmodel == "F4P" )
        {
            logi    <- F4P
            pn      <- 4
            pname   <- c("df","a","n","m")
        }
        else if ( psdmodel == "F3P" )
        {
            logi    <- F3P
            pn      <- 3
            pname   <- c("a","n","m")
        }
        else if ( psdmodel == "ML" )
        {
            logi    <- ML
            pn      <- 3
            pname   <- c("a","b","c")
        }
        else if ( psdmodel == "ONL" )
        {
        logi    <- ONL
        pn      <- 3
        pname   <- c("u","o","c")
        }
        else if ( psdmodel == "ORL" )
        {
            logi    <- ORL
            pn      <- 3
            pname   <- c("u","o","e")
        }
        else if ( psdmodel == "S" )
        {
            logi    <- S
            pn      <- 2
            pname   <- c("u","c")
            #S model can not deal with first texture data with zero value
            if(y[1] == 0) y[1] <- 0.0001
        }
        else if ( psdmodel == "VG" )
        {
            logi    <- VG
            pn      <- 2
            pname   <- c("dg","n")
            maxspa2 <- 2
            minspa2 <- 1
        }
        else if ( psdmodel == "VGM" )
        {
            logi    <- VGM
            pn      <- 2
            pname   <- c("dg","n")
            maxspa1 <- 200
            minspa1 <- 4
            maxspa2 <- 2
            minspa2 <- 0.5
        }
        else if ( psdmodel == "W" )
        {
            logi    <- W
            pn      <- 3
            pname   <- c("a","b","c")
        }
        else if ( psdmodel == "L" )
        {
            logi    <- L
            pn      <- 2
            pname   <- c("a","b")
        }
        else if ( psdmodel == "LG" )
        {
            logi    <- LG
            pn      <- 2
            pname   <- c("a","b")
        }
        else if ( psdmodel == "SC" )
        {
            logi    <- SC
            pn      <- 3
            pname   <- c("u","o","e")
        }
        else if ( psdmodel == "SL" )
        {
            logi    <- SL
            pn      <- 2
            pname   <- c("u","o")
        }
        #default lower and upper limit for drc::drm, these values should not set
        #at the beginning of the function for pn is set later
        lowerl <- rep(10e-9,times=pn)
        upperl <- rep(10e+5,times=pn)
        #Initailize spa for drc::drm
        spa <- c(1,1,1,1)
        #methods used in optim() of drc::drm
        meth <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
        
        mdev <- 100
        for( i in 1:5 ) # The nonlinear optimization runs were carried out using at least
                        # five random initial parameter estimates for all soils.
                        #When the final solution for each soil converged to different parameter values,
                        #the parameter values with the best fitting statistics (RSS) were kept.
        {
            if( method == "all" )# using all optim methods
            {
                for( j in 1:5 )
                {
                    countR <- 0  # Added by Julien Moeys on 2011/11/01
                    repeat
                    {   
                        countR <- countR + 1 # Added by Julien Moeys on 2011/11/01
                        
                        spa[1:pn-1]     <- runif(n = pn-1,max = maxspa1,min = minspa1)
                        spa[pn]         <- runif(n = 1,max = maxspa2,min = minspa2)
                        tt<- try( drm(y ~ xin, fct = list( logi, NULL,pname ), # JM:2010/08/11 changed drc::drm to drm alone
                                start   = spa[1:pn],
                                #roust  = "median",
                                lowerl  = lowerl,
                                upperl  = upperl,
                                control = drmc(constr = TRUE,maxIt = 500, # JM:2010/08/11 changed drc::drmc to drmc alone
                                    noMessage   = TRUE,
                                    method      = meth[j],
                                    # method    = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                                    #trace      = TRUE
                                )
                            )
                            , silent = TRUE
                         )
                        if( !inherits(tt, "try-error") )
                        {
                            dev <- sum(residuals(tt)^2)
                            if( mdev > dev )
                            {
                                mdev    <- dev
                                ttbest  <- tt
                            }
                            break
                        } #
                        
                        if( countR >= 100 ){ break } # Added by Julien Moeys on 2011/11/01
                    }
                }
            }
            else
            {
                countR <- 0  # Added by Julien Moeys on 2011/11/01
                repeat
                {
                    countR <- countR + 1 # Added by Julien Moeys on 2011/11/01
                    
                    spa[1:pn-1] <- runif(n=pn-1,max = maxspa1,min = minspa1)
                    spa[pn]     <- runif(n=1,max = maxspa2,min = minspa2)
                    tt  <- try( drm(y ~ xin, fct = list(logi, NULL, pname), # JM:2010/08/11 changed drc::drm to drm alone
                            start   = spa[1:pn],
                            #roust  = "median",
                            lowerl  = lowerl,
                            upperl  = upperl,
                            control = drmc(constr = TRUE,maxIt = 500,# JM:2010/08/11 changed drc::drmc to drmc alone
                                noMessage   = TRUE,
                                method      = method,
                                # method    = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
                                #trace      = TRUE
                            )
                        )
                        , silent = TRUE
                    )
                    if( !inherits(tt, "try-error") )
                    {
                        dev<-sum(residuals(tt)^2)
                        if( mdev>dev )
                        {
                            mdev <- dev
                            ttbest <- tt
                        }
                        break
                    }
                    
                    if( countR >= 100 ){ break } # Added by Julien Moeys on 2011/11/01
                }
            }
            #when the residual sum of error (deviance) is very small, the iteration is stopped to save time
            if(mdev < 0.0001) break
        }
        #predict() has some bug for PSD model to predict the target values
        if( psdmodel == "AD" )
        {
            pre <- coef(ttbest)[1] + coef(ttbest)[2]*atan(coef(ttbest)[3]*log10(xout/coef(ttbest)[4]))
        }
        else if( psdmodel == "F4P" )
        {
            pre <- (1-(log(1+coef(ttbest)[1]/xout)/log(1+coef(ttbest)[1]/0.0001))^7)/(log(exp(1)+(coef(ttbest)[2]/xout)^coef(ttbest)[3]))^coef(ttbest)[4]
        }
        else if( psdmodel == "F3P" )
        {
            pre <- (1-(log(1+0.001/xout)/log(1+0.001/0.0001))^7)/(log(exp(1)+(coef(ttbest)[1]/xout)^coef(ttbest)[2]))^coef(ttbest)[3]
        }
        else if( psdmodel == "ML" )
        {
            pre <- 1/(1+coef(ttbest)[1]*exp(-coef(ttbest)[2]*xout^(coef(ttbest)[3])))
        }
        else if( psdmodel == "ONL" )
        {
           t    <- (-1)^(log(xout) >= coef(ttbest)[1]+1)
           pre  <- (1+t*erf((log(xout)+coef(ttbest)[1])/coef(ttbest)[2]*2^0.5))/2+(coef(ttbest)[3])
        }
        else if( psdmodel == "ORL" )
        {
            t   <- (-1)^(log(xout) >= coef(ttbest)[1]+1)
            pre <- (1-coef(ttbest)[3])*(1+t*erf((log(xout)+coef(ttbest)[1])/coef(ttbest)[2]*2^0.5))/2+coef(ttbest)[3]
        }
        else if( psdmodel == "S" )
        {
            pre <- 1/(1+(1/y[1]-1)*exp(-coef(ttbest)[1]*((xout-r0)/r0)^coef(ttbest)[2]))
        }
        else if( psdmodel == "VG" )
        {
            pre <- (1+(coef(ttbest)[1]/xout)^coef(ttbest)[2])^(1/coef(ttbest)[2]-1)
        }
        else if( psdmodel == "VGM" )
        {
            pre <- y[1]+(1-y[1])*(1+(coef(ttbest)[1]*xout)^(-coef(ttbest)[2]))^(1/coef(ttbest)[2]-1)
        }
        else if( psdmodel == "W" )
        {
            pre <- coef(ttbest)[3]+(1-coef(ttbest)[3])*(1-exp(-coef(ttbest)[1]*((xout-dmin)/(dmax-dmin))^coef(ttbest)[2]))
        }
        else if( psdmodel == "L" )
        {
            pre <- coef(ttbest)[1]*log(xout)+coef(ttbest)[2]
        }
        else if( psdmodel == "LG" )
        {
            pre <- 1/(1+coef(ttbest)[1]*exp(-coef(ttbest)[2]*xout))
        }
        else if( psdmodel == "SC" )
        {
            t   <- (-1)^(log(xout) >= coef(ttbest)[1]+1)
            pre <- (1-coef(ttbest)[3])*(1+t*erf((log(xout)+coef(ttbest)[1])/coef(ttbest)[2]*2^0.5))/2+coef(ttbest)[3]*(1+t*erf((log(xout)+1.96)/1*2^0.5))/2
        }
        else if( psdmodel == "SL" )
        {
            t   <- (-1)^(log(xout) >= coef(ttbest)[1]+1)
            pre <- (1+t*erf((log(xout)+coef(ttbest)[1])/coef(ttbest)[2]*2^0.5))/2
        }
        #pre are the predicted values, coef(ttbest) are the model paremeters,
        out <- c(pre[1],pre[2:length(pre)]-pre[1:length(pre)-1])*100
        #dev is the deviance ( Residual Sum of Squaures)
        c(out,coef(ttbest),dev=mdev*10000)
    }

    results <- t(apply(
        X       = tri.data[1:dim(tri.data)[1],],

        MARGIN  = 1,
        FUN     = fitpsd,
        xin     = dat.ps.lim[ ps.start:ps.end ],
        xout    = base.ps.lim[ ps.start:length(base.ps.lim) ],
        psdmodel= psdmodel,
        method  = omethod)
    )


#    results <- t(apply(
#        X       = tri.data[1:5,],
#        MARGIN  = 1,
#        FUN     = fitpsd,
#        xin     = dat.ps.lim[ ps.start:ps.end ],
#        xout    = base.ps.lim[ ps.start:length(base.ps.lim) ],
#        psdmodel= psdmodel,
#        method  = omethod)
#    )
    colnames(results)[1:(length(base.ps.lim)-ps.start+1)]<-
         paste(sep = "",c(0,base.ps.lim[ ps.start:(length(base.ps.lim)-1)])*1000,"-",base.ps.lim[ ps.start:length(base.ps.lim) ]*1000)
    results
}

#     my.text4 <- data.frame(
#         "CLAY"  = c(05,60,15,05,25,05,25,45,65,75,13,47),
#         "FSILT" = c(02,04,10,15,25,40,35,20,10,05,10,20),
#         "CSILT" = c(03,04,05,10,30,45,30,25,05,10,07,23),
#         "SAND"  = c(90,32,70,70,20,10,10,10,20,10,70,10)
#     )   #
#     TT.text.transf.Xm(
#       tri.data    = my.text4,
#       base.ps.lim = c(0,2,20,50,2000),
#       dat.ps.lim  = c(0,2,20,63,2000),
#       psdmodel    = "S"
#     )   #
#     TT.text.transf.Xm( # JM: does not work on my PC
#       tri.data    = my.text4,
#       base.ps.lim = c(0,1,50,2000),
#       dat.ps.lim  = c(0,2,30,60,2000),
#       psdmodel    = "AD",
#       omethod     = "Nelder-Mead"
#     )


