
post.ppm <-
function(out, slicedata, ymaskdata, vreg = 2, swap = FALSE, plot = TRUE,
  col=heat.colors(256))
{
    # exitpar <- par(no.readonly=T) 
    # par(mfrow=c(1,2))
    # x11(width=8,height=4.5)
    kin <- ymaskdata$kin # indices of masked out slice times series
    #-----------
    niislicets <- slicedata$niislicets
    niislicets <- niislicets[,,1] # ! radiological convention 
		if(swap)
	      niislicets <- niislicets[dm[1]:1,,1] # ! radiological convention 
    mask <- slicedata$mask
    dm <- dim(mask)
    # cat("mask ",dm,"\n")
#    if(plot) { # initial pre-filtered image
#        # x11(width=4,height=4.5)
#        par(mar=c(0,0,2,0), xaxt="n", yaxt="n", ask=F)
#        image(niislicets, col=heat.colors(256), main="initial pre-filtered image",
#                  axes=F)
#        # image(niislicets[,,1], col=gray((0:255)/256))
#		}
    #-------------------
    pmeans = pmeans.hcoef(out$betadraw) 
    px <- regpostsim(pmeans, vreg=vreg, plot=F)
    pm2 <- pmeans[,vreg]
    spma <- px$spma
    spmn <- px$spmn
		nactive	<- length(spma)
    #---------
#    # mean activation image
#    print("mean activation:")
#    spm <- matrix(0, nrow=dm[1], ncol=dm[2])
#    print(dim(spm))
#    spm[kin] <- pm2
#    if(swap)
#        spm <- spm[dm[1]:1, ] # ! radiological convention 
#    if(plot) {
#        x11(width=4,height=4.5)
#        par(mar=c(0,0,2,0), xaxt="n", yaxt="n", ask=F)
#        image(spm, col=heat.colors(256), main="mean activation", axes=F)
#    }
    #-----------
    # ppm image 
    pm2.act <- pm2[spma]
    ppm <- matrix(0, nrow=dm[1], ncol=dm[2])
    if(length(pm2.act))
        ppm[kin[spma,]] <- pm2.act
    # plot with overlay
    # ppm.m <- ppm+0.3*mask
    ppm.m <- ppm/max(ppm) * 0.65 + niislicets/max(niislicets) * 0.35
    ppm.m <- ppm.m/max(ppm.m) * mask
    if(swap) 
        ppm.m <- ppm.m[dm[1]:1, ] # ! radiological convention 
    if(plot) {
        # x11(width=4,height=4.5)
        par(mar=c(0,0,2,0), xaxt="n", yaxt="n", ask=F)
    		main <- paste("ppm image ; vreg = ",vreg,sep="") 
        image(ppm.m, col=col, main=main, axes=F)
    }
		# par(exitpar)
		invisible(list(ppm=ppm, nactive=nactive))
}

