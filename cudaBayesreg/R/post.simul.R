#
# Postprocessing of MCMC simulation
#

post.simul.betadraw <-
function(out, vreg = 2)  
{
    pmeans <- pmeans.hcoef(out$betadraw) 
    px <- regpostsim(pmeans, vreg=vreg, plot=F)
    spma <- px$spma
    spmn <- px$spmn
    #-----------------------------------------
    # Plotting betadraw with emphasis on active and non-active pixels 
    cat("number of thresholded active:", length(spma),"\n")
    cat("number of thresholded non-active:", length(spmn),"\n")
    #  plot arrays of draws of coefs in hier models
    if(length(spma)) # active voxel ccoeffs
        plot.hcoef.post(out$betadraw,spmname="activated",spm=spma) 
    if(length(spmn)) # non-active voxel ccoeffs
        plot.hcoef.post(out$betadraw,spmname="non-activated",spm=spmn) 
}

post.simul.hist <-
function(out, vreg = 2)
{
    pmeans <- pmeans.hcoef(out$betadraw) 
    px <- regpostsim(pmeans, vreg=vreg)
		invisible(px)
}

