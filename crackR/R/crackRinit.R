crackRinit <-
function(parameters)
{
    ## this function will generate an object which can run a PDTA analysis
    ## the parameter "analysis.type" determines if it's a single type run, multiple type, or cont. damage

    ## single or multiple type; class parametersSing or parametersMult
    ##  1) list of PF run parameters
    ##  2) the initial PF state: a, kc, w, typ (if multiple)
    ##  3) empty results placeholder

    ## continuing damage; class parametersCD
    ##  1) list of PF run parameters
    ##  2) the initial PF state: a.p, a.s, kc, w
    ##  3) empty results placeholder

    Np <- parameters$Np

    ## component of the output list for storing results; this is an empty placeholder
    ## results are in the same form regardless of the analysis type
    results <- list(Np.total=Np,
                    sfpof=data.frame(),
                    pof.int=data.frame(),
                    pcd=data.frame()
                    )
    class(results)    <-  c("crackRresults", "list")

    analysis.type <- parameters$analysis.type
    
    ## single DTA type analysis
    if( analysis.type == "single" )
    {
        a     <- parameters$dta$ifs.rsamp(Np)
        kc    <- parameters$dta$kc.rsamp(Np)
        w     <- parameters$dta$ifs.dactual(a) / parameters$dta$ifs.dsamp(a)
        w     <- w / sum(w)
        state <- list(a  = a,
                      kc = kc,
                      w  = w)
        class(parameters) <- c("parametersSing", "list")
        class(state)      <- c("stateSing", "list")
        obj <- list(parameters=parameters,
                    state=state,
                    results=results)
        class(obj) <-  c("Sing", "list")

    }

    ## multiple DTA type analysis
    if( analysis.type == "multiple" )
    {
        a     <- parameters$dta[[1]]$ifs.rsamp(Np)
        kc    <- parameters$dta[[1]]$kc.rsamp(Np)
        typ   <- rep(1, Np)
        w     <- parameters$dta[[1]]$ifs.dactual(a) / parameters$dta[[1]]$ifs.dsamp(a)
        w     <- w / sum(w)
        state <- list(a   = a,
                      kc  = kc,
                      typ = typ,
                      w   = w)
        class(parameters) <- c("parametersMult", "list")
        class(state)      <- c("stateMult", "list")
        obj <- list(parameters=parameters,
                    state=state,
                    results=results)
        class(obj) <-  c("Mult", "list")
    }
    
    ## continuing damage analysis
    ## currently assuming Kc is common to both locations (not necessary)
    if( analysis.type == "CD" )
    {
        a.p <- parameters$dta.pc$ifs.rsamp(Np)
        a.s <- parameters$dta.sc$ifs.rsamp(Np)
        kc  <- parameters$dta.pc$kc.rsamp(Np)
        w   <- ( parameters$dta.pc$ifs.dactual(a.p) * parameters$dta.sc$ifs.dactual(a.s) ) /
               ( parameters$dta.pc$ifs.dsamp(a.p)   * parameters$dta.sc$ifs.dsamp(a.s) )
        w   <- w / sum(w)
        state <- list(a.p = a.p,
                      a.s = a.s,
                      kc  = kc,
                      w   = w)
        class(parameters) <- c("parametersCD", "list")
        class(state)      <- c("stateCD", "list")
        obj <- list(parameters=parameters,
                    state=state,
                    results=results)
        class(obj) <-  c("CD", "list")

        ## there's currently no bootstrapping with continuing damage
        if( parameters$bootstrap.sfpof ||  parameters$bootstrap.pcd )
            {
                warning("there is no bootstrapping of sfpof or pcd with continuing damage in the current version of crackR")
            }
    }

    return(obj)
}
