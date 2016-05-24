AT.SPC.spectrum.at.depth.step <- function(spc, depth.step)
{
    if(depth.step %in% unique(spc$depth.step)){
        ii     <- spc$depth.step == depth.step
        df     <- data.frame( E.MeV.u     = spc$E.MeV.u[ii],
                              particle.no = spc$particle.no[ii],
                              fluence.cm2 = spc$fluence[ii])
        #cat(paste("Returning", sum(ii), "entries in data frame.\n"))
        return(df)
    }else{
        cat("Depth step not found in data. Skipping.")
        return(NULL)
    }
}
