
"stampDive" <- function(x, ignoreZ=TRUE)
{
    ## Value: A data frame; stamping each dive with phase number it belongs
    ## to, activity during that phase, and trip start and end time
    ## --------------------------------------------------------------------
    ## Arguments: x=TDRcalibrate object, ignoreZ=ignore Z phases?
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (!is(x, "TDRcalibrate")) stop ("x needs to be a TDRcalibrate object")
    act <- getGAct(x, "activity")
    diveid <- getDAct(x, "dive.id")

    if (ignoreZ) {
        tt <- getTime(getTDR(x))
        interval <- getDtime(getTDR(x))
        act[act == "Z"] <- "L"
        attlist <- .rleActivity(tt, act, interval) # recalculate
        phaseid <- as.numeric(attlist[[1]])  # what phase.id is now
    } else {
        attlist <- getGAct(x)
        phaseid <- getGAct(x, "phase.id")
    }

    beg <- rep(attlist[[3]], table(phaseid))
    end <- rep(attlist[[4]], table(phaseid))
    phase.no <- numeric(length(act))      # vector of 0s
    phaseid[act == "L"] <- 0             # dry phase.id should be 0
    ## make a sequence for phase.id > 0 from 1:number of such phases
    phase.no[act != "L"] <- rep(seq(along=table(phaseid[phaseid > 0])),
                table(phaseid[phaseid > 0]))
    ok <- match(unique(diveid[diveid > 0]), diveid) # required subscripts
    phase.no <-  phase.no[ok]
    activity <- act[ok]

    data.frame(phase.no, activity, beg=beg[ok], end=end[ok])
}
