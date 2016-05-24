F.step.cjs <- function( cap.covars, surv.covars, fit.crit="qaicc", signif.drop=0, signif.increase=0, plot=TRUE, ... ){
#
#   stepwise selection of a CJS model.
#
#   cap.covars = a vector containing covariates to consider for fitting in the capture equation.
#   surv.covars= a vector containing covariates to consider for fitting in the survival equation.
#   fit.crit   = model fit criterion used to rank models in each step.
#   signif.drop = value below which a drop in fit.crit is considered "significant" during forward steps.  This controls stopping. Routine stops
#       when fit.crit on current iteration minus fit.crit on previous iteration is greater than or equal to signif.drop. 
#       This means signif.drop should be less than or equal to 0, unless for some odd reason, 
#       steps that add variables without predictive abilities are desired.
#       If signif.drop = 0 (the default), steps are haulted when fit.crit increases from previous to current step.  If 
#       signif.drop = -2, steps are haulted when fit.crit does not decrease by 2 or more points between the 
#       current and previous steps.  
#   signif.increase = value above which a increase in fit.crit is considered "significant" during backward looks.  This controls elimination 
#       of covariates during backward steps. A variable in the current model is removed if, upon removal, fit.crit 
#       increases by less than or equal to signif.increase.  This means signif.increase should be greater than or equal to 0 to be meaningfull, 
#       unless no variables should be removed, in which case set signif.increase = -Inf (or some other negative number). 
#       If signif.increase = 0 (the default), variables are left in the model if they increase fit.crit when removed.  If 
#       signif.increase = 2, variables in the current model are left in the model if they increase fit.crit by more than 2 units 
#       when removed. 
#   plot = logical scalar for whether or not to plot fit.crit at end of each step.  Occasionally, inspecting the sequences of 
#       drops in fit.crit is informative.
#   ... = additional arguments to F.cjs.estim, like histories= , nhat.v.meth=, c.hat=, etc.

require(mra)

#   Internal function to make a single text string out of vector of string objects
f.myformula <- function( covs ){
    myformula <- "1"   # models always have intercepts.  This routine will not allow a model without an intercept
    for( i in covs ){
        myformula <- paste( myformula, "+", i )
    }
    myformula
}




#   Fit initial model, and initialize things
a1 <- F.cjs.estim( capture= ~1, survival= ~1, ... )   # starting model is null model
cur.min <- a1[[fit.crit]]

cat(paste("========================= Initial Model =========================\n"))
cat(paste( "Surv=", 1, "\tCap=", 1, "\t:", fit.crit, "=", round(cur.min,4), "\n" ))

if( plot ){
    step.val <- 0
    fit.val <- cur.min
    plot( step.val, fit.val, ylab=fit.crit, xlab="Step", type="b" )
}

step <- 1

cur.capmod <- NULL
cur.surmod <- NULL

while( TRUE ){

    cat(paste("========================= Step", step, "================================\n"))

    cat(paste("--- VARIABLE POOL ---\n"))
    cat(paste(" CAPTURE:", f.myformula( cap.covars ), "\n"))
    cat(paste("SURVIVAL:", f.myformula( surv.covars ), "\n"))


    cat(paste("\n------ FORWARD ------\n"))

    #   Make text strings for both current models, to be used in formula objects below
    cur.cap.formula <- f.myformula( cur.capmod )
    

    #   Fix survival, loop through capture covars
    step.advanced <- FALSE
    added.cov <- NULL
    added.mod <- NULL
    cur.sur.formula <- f.myformula( cur.surmod )   
    for( i in cap.covars){
        newmod <- c( cur.capmod, i )
        cur.cap.formula <- f.myformula( newmod )

        cat(paste( "Surv=", cur.sur.formula, "\tCap=", cur.cap.formula, ":", fit.crit, "= " ))

        a <- F.cjs.estim( capture=formula(paste("~", cur.cap.formula )),
                          survival=formula(paste("~", cur.sur.formula )), ... )
        cat(paste( round(a[[fit.crit]],4), "Converged:", as.logical(a$exit.code), "\n"))

        if( ((a[[fit.crit]] - cur.min) <= signif.drop) & (a$exit.code == 1) ){
            #   This model is better than the current minimum
            a1 <- a            
            cur.min <- a[[fit.crit]]
            added.cov <- i
            added.mod     <- "CAPTURE"
            step.advanced <- TRUE
        }
    }

    #   Fix capture, loop through survival covars
    cur.cap.formula <- f.myformula( cur.capmod )   
    for( i in surv.covars){
        newmod <- c( cur.surmod, i )
        cur.sur.formula <- f.myformula( newmod )

        cat(paste( "Cap =", cur.cap.formula, "\tSurv=", cur.sur.formula, ":", fit.crit, "= " ))

        a <- F.cjs.estim( capture=formula(paste("~", cur.cap.formula)),
                          survival=formula(paste("~", cur.sur.formula )), ... )
        cat(paste( round(a[[fit.crit]], 4), "Converged:", as.logical(a$exit.code),  "\n"))
        
        if( ((a[[fit.crit]] - cur.min) <= signif.drop) & (a$exit.code == 1) ){
            #   This model is better than the current minimum
            a1 <- a            
            cur.min <- a[[fit.crit]]
            added.cov <- i
            added.mod     <- "SURVIVAL"
            step.advanced <- TRUE
        }
    }
    
    #   At this point, a1 is the new best model, added.cov is the new covariate added, cur.min is the new minimum fit statistic, mod = model that we added to

    #   Check whether step added a covariate, if not quit here
    if( !step.advanced ){
        cat(paste("\nMinimum", fit.crit, "model found.\n"))
        break
    }
    

    #   Add covariate to approprate model.  Remove it from pool of potential covariates
    if( added.mod == "CAPTURE" ){
        cur.capmod <- c( cur.capmod, added.cov )
        cap.covars <- setdiff( cap.covars, added.cov )   # remove added covariate from the pool
    } else {
        cur.surmod <- c( cur.surmod, added.cov )
        surv.covars <- setdiff( surv.covars, added.cov )
    }

    
    #   Now do the backward look
    removed.cov <- NULL
    step.retreated <- FALSE
    if( (length(cur.capmod) > 1) | (length(cur.surmod) > 1) ){

        cat(paste("\n------ BACKWARD ------\n"))

        #   First, remove any variables in capture model that were not just added    
        for( i in cur.capmod ){
            if( !(i == added.cov & added.mod == "CAPTURE") ){
                newmod <- setdiff( cur.capmod, i )
                cur.cap.formula <- f.myformula( newmod )
                cur.sur.formula <- f.myformula( cur.surmod )   

                cat(paste( "Surv =", cur.sur.formula, "\tCap =", cur.cap.formula, ":", fit.crit, "=" ))
                a <- F.cjs.estim( capture=formula(paste("~", cur.cap.formula)),
                                  survival=formula(paste("~", cur.sur.formula)), ... )
                cat(paste( round(a[[fit.crit]],4), "Converged:", as.logical(a$exit.code),  "\n"))
            
                if( ((a[[fit.crit]] - cur.min) <= signif.increase) & (a$exit.code == 1) ){
                    #   This model is better than the current minimum
                    a1 <- a            
                    cur.min <- a[[fit.crit]]
                    removed.cov <- i
                    removed.mod     <- "CAPTURE"
                    step.retreated <- TRUE
                }
            }
        }
 
        #   Second, remove any variables in survival model that were not just added    
        for( i in cur.surmod ){
            if( !(i == added.cov & added.mod == "SURVIVAL") ){
                newmod <- setdiff( cur.surmod, i )
                cur.cap.formula <- f.myformula( cur.capmod )
                cur.sur.formula <- f.myformula( newmod )   

                cat(paste( "Cap =", cur.cap.formula, "\tSurv =", cur.sur.formula, ":", fit.crit, "=" ))
                a <- F.cjs.estim( capture=formula(paste("~", cur.cap.formula)),
                                  survival=formula(paste("~", cur.sur.formula)), ... )
                cat(paste( round(a[[fit.crit]],4), "Converged:", as.logical(a$exit.code),  "\n"))
            
                if( (a[[fit.crit]] - cur.min <= signif.increase) & (a$exit.code == 1) ){
                    #   This model is better than the current minimum
                    a1 <- a            
                    cur.min <- a[[fit.crit]]
                    removed.cov <- i
                    removed.mod     <- "SURVIVAL"
                    step.retreated <- TRUE
                }
            }
        }
 
 
        if( step.retreated ){
            if( removed.mod == "CAPTURE" ){
                cur.capmod <- setdiff( cur.capmod, removed.cov )
                cap.covars <- c( cap.covars, removed.cov )   #  removed variable goes back into the pool
            } else {
                cur.surmod <- setdiff( cur.surmod, removed.cov )
                surv.covars <- c( surv.covars, removed.cov )
            }
        }
    }

    cat(paste("\n------ STEP",step, "SUMMARY ------\n"))
    cat("        Added: ")
    if( length(added.cov) > 0 ){
        cat(paste( added.cov, "to", added.mod, "\n"))
    } else {
        cat("\n")
    }
    
    cat("      Removed: ")
    if( length(removed.cov) > 0 ){
        cat(paste( removed.cov, "to", removed.mod, "\n"))
    } else {
        cat("\n")
    }
    cat(paste("Current model:\n"))
    cat(paste(format("CAPTURE:", justify="right", width=9), f.myformula( cur.capmod ), "\n"))
    cat(paste(format("SURVIVAL:", justify="right", width=9), f.myformula( cur.surmod ), "\n"))
    cat(paste(format( paste(fit.crit,":",sep=""), justify="right", width=9 ), round(cur.min,4), "\n"))
   
    if( plot ){
        step.val <- c(step.val, step)
        fit.val <- c(fit.val, cur.min) 
        plot( step.val, fit.val, ylab=fit.crit, xlab="Step", type="b" )
    } 
    
    step <- step + 1

}

#   Return the saved best model object
a1

}
