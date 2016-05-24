
removesusceptibles <- function(epi){
  infected = !(is.na(epi[,3]) & is.na(epi[,4]) & is.na(epi[,5]))
  if (sum(infected) > 1){
    return(epi[infected,])
  } else {
    # make sure epi remains in correct matrix format
    return(t(as.matrix(epi[infected,])))
  }
}

plotepitree <- function(epi, lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75, 
                        zero.at.start = FALSE, main = "Transmission Tree", xlab = "Time", ylab= "",
                        e.col = "black", i.col = "red", lty.transmission = 3, marktransitions = TRUE, 
                        label.trans = "|", cex.trans = 0.5, ...){
  
  epi = removesusceptibles(epi)
  
  # find total number of infecteds -- need at least 2 infecteds in the epidemic
  ninf <- nrow(epi)
  if (ninf < 2) stop("Plot error: need at least two infecteds to plot the epidemic.")
  
  # make sure we have full data
  if(sum(is.na(epi[,3:5])) > 0) stop("Plot error: need full data to plot the epidemic -- run again with all times known.")
  if(sum(is.na(epi[,1])) > 0) stop("Plot error: missing entries in first (Node ID) column.")
  if(sum(is.na(epi[,2])) > 1) stop("Plot error: too many NA values in second (Parent) column -- should only have one initial infected node.")
  
  if (zero.at.start) epi[,3:5] = epi[,3:5] - min(epi[,3:5])
  
  # sort the epidemic by increasing order of exposure
  epi <- epi[order(epi[ ,3]), ]
  
  # count the cumulative number of children each infective has
  nchild = array(0,max(epi[ ,1]))
  for (i in ninf:2) nchild[epi[i,2]] = nchild[epi[i,2]] + nchild[epi[i,1]] + 1
  
  ypos = array(0,max(epi[,1]))
  ypos[epi[1,1]] = 1
  
  # position of the current is position of its parent + gap for as yet unencountered direct offspring of the parent
  for (i in 2:ninf) {  
    ypos[epi[i,1]] = ypos[epi[i,2]] + nchild[epi[i,2]] - nchild[epi[i,1]]
    nchild[epi[i,2]] = nchild[epi[i,2]] - nchild[epi[i,1]] - 1
  }
  
  # set up axes and titles
  plot(x = c(min(epi[ ,3]), max(epi[ ,5])+2), y = c(0.9,ninf ), pch=" ", 
       xlab = xlab, ylab = ylab, main = main, yaxt = "n", bty = "n", ...)      
  
  for (i in 1: ninf)        {
    # plot Exposed portion of infection
    lines(x=c(epi[i,3],epi[i,4]), y=c(ypos[epi[i,1]],ypos[epi[i,1]]), lwd = lwd, col = e.col)
    # plot Infective portion of infection
    lines(x=c(epi[i,4],epi[i,5]), y=c(ypos[epi[i,1]],ypos[epi[i,1]]), lwd = lwd, col = i.col)
    # plot connection between parent and child
    lines(x=c(epi[i,3],epi[i,3]), y=c(ypos[epi[i,2]],ypos[epi[i,1]]), lwd = lwd, lty = lty.transmission)
  }
  
  # plot labels
  if (leaf.labs) text(epi[,5],ypos[epi[,1]],labels = epi[,1],pos = 4, offset = 0.25, cex = leaf.cex)
  
  # Mark transition points
  if (marktransitions)
    for (i in 1:ninf)
      text(epi[i,4],ypos[epi[i,1]], labels=label.trans, cex = cex.trans)
  
}

	
#plotepitreemcmc <- function(mcmcoutput, index = dim(mcmcoutput$transtree)[2], lwd = 1,
#	leaf.labs = TRUE, leaf.cex = 0.75, zero.at.start = FALSE, main = "Transmission Tree",
#	xlab = "Time", ylab= "", e.col = "black", i.col = "red", lty.transmission = 3,
#	marktransitions = TRUE, label.trans = "|", cex.trans = 0.5, ...) {

#   epi <- buildepifromoutput(mcmcoutput, index)
		
#	plotepitree(epi, lwd = lwd, leaf.labs = leaf.labs, leaf.cex = leaf.cex,
#		zero.at.start = zero.at.start, main = main, xlab = xlab, ylab= ylab, e.col = e.col,
#		i.col = i.col, lty.transmission = lty.transmission, marktransitions = marktransitions,
#		label.trans = label.trans, cex.trans = cex.trans, ...)
#}

epi2newickmcmc <- function(mcmcoutput, index = dim(mcmcoutput$transtree)[2]) {
  epi <- buildepifromoutput(mcmcoutput, index)
	return(epi2newick(epi))
}

epi2newick <- function(epi){
  return(epi2newick.int(epi))
}
  
epi2newick.int <- function(epi,node = epi[1,1],lastrow = 1,infectedatrow = 1, rightchild = FALSE){
  nex = "";
  if (lastrow == 1){
    ## at the root of the tree	
    ## make sure there are only infecteds listed
    epi = removesusceptibles(epi)  
    if (nrow(epi) < 2) stop("Need at least two infecteds to convert transmission tree to Newick format")
    ## start recursing
    return(paste('((',epi2newick.int(epi,node=epi[2,"Parent"],lastrow = 2,infectedatrow = 1),',',epi2newick.int(epi,node = epi[2,"Node ID"],lastrow = 2,infectedatrow = 2, rightchild = TRUE),')[&type = "I"]:',epi[2,'Etime'] - epi[1,'Itime'],')[&type = "E"]:',epi[1,'Itime'] - epi[1,'Etime'],sep = ''))
  }
  ## need to determine whether we are at a leaf or coalescence
  atleaf = FALSE
  if (lastrow == nrow(epi)){ 
    ## at end of epidemic so must be at a leaf
    atleaf = TRUE
  } else {
    ## search remaining rows to see if node is involved in another event
    nextrow = match(node,epi[(lastrow+1):nrow(epi),'Parent'])
    atleaf = is.na(nextrow)		
  }
  if (atleaf){
    ## at a leaf
    if (rightchild){
      #  went from E to I along branch -- has form (name[&type = "I"]:time)[&type = "E"]:time
      return(paste('(',node,'[&type = "I"]:',epi[infectedatrow,'Rtime'] - epi[infectedatrow,'Itime'],')[&type = "E"]:',epi[infectedatrow,'Itime'] - epi[lastrow,'Etime'],sep = ''))
    } else {
      # was I along whole branch -- has form name[&type = "I"]:time
      return(paste(node,'[&type = "I"]:',epi[infectedatrow,'Rtime'] - epi[lastrow,'Etime'],sep = ''))
    }
    
  } else { 
    ## at a coalescence 
    # adjust nextrow index since only searched rows below lastrow 
    nextrow = nextrow+lastrow
    if (rightchild){
      #  went from E to I along branch --  has form ((infector,infectee)[&type = "I"]:time):[&type = "E"]:time
      return(paste('((',epi2newick.int(epi,node=node,lastrow = nextrow,infectedatrow = infectedatrow),',',epi2newick.int(epi,node=epi[nextrow,'Node ID'],lastrow = nextrow,infectedatrow = nextrow, rightchild = TRUE),')[&type = "I"]:',epi[lastrow,'Itime'] - epi[lastrow,'Etime'],')[&type = "E"]:',epi[nextrow,'Etime'] - epi[lastrow,'Itime'],sep = ''))
    } else {
      #  was I along whole branch --  has form (infector,infectee)[&type = "I"]:time
      return(paste('(',epi2newick.int(epi,node=node,lastrow = nextrow,infectedatrow = infectedatrow),',',epi2newick.int(epi,node=epi[nextrow,'Node ID'],lastrow = nextrow,infectedatrow = nextrow, rightchild = TRUE),')[&type = "I"]:',epi[nextrow,'Etime'] - epi[lastrow,'Etime'],sep = ''))
    }
  }
}

write.epinet <- function(out, filename){
  
  ## write the parameter file
  paramfile = paste(filename,".log", sep = "")
  
  param = data.frame(Sample = (0:(length(out$llkd)-1))*out$mcmcinfo$thinning, llkd = out$llkd, beta = out$beta, thetai = out$thetai, thetae = out$thetae, ki = out$ki, ke = out$ke)
  cat('# sampled parameters from epinet run. File created: ',  date(),'\n', file = paramfile,append = F)
  cat('# call:',deparse(out$call),'\n',file = paramfile,append = T)
  suppressWarnings(write.table(cbind(param,out$eta), file = paramfile,append = T, row.names = F, quote = F, sep = "\t", eol = "\n"))
    
  ## write the tree file (if transmission tree posteriors exist)
  if (!is.null(out$transtree)) {
    # write headers
    treefile = paste(filename,".trees", sep = "")
    cat('#NEXUS\n\n', file = treefile,append = F)
    #cat('[sampled trees from epinet run. File created: ',  date(),']\n', file = treefile,append = T)
    #cat('[call:',deparse(out$call),']\n\n',file = treefile,append = T)
    cat('Begin trees;\n', file = treefile,append = T)
  
    # write body
    treestepsize =  out$mcmcinfo$thinning*out$mcmcinfo$extrathinning
    for (i in 1:ncol(out$transtree)) cat('tree tree',i*treestepsize - out$mcmcinfo$thinning,' = [&R] ',epi2newick(buildepifromoutput(out,i)),';\n',file = treefile, append = T,sep = '')

    # terminate
    cat('End;\n', file = treefile,append = T)
  }
}

buildepifromoutput <- function(mcmcoutput, index){
  if(is.null(mcmcoutput$transtree) || is.null(index)) stop("Error: Need inferred transmission tree samples for this function")  
  ninf <- dim(mcmcoutput$transtree)[1]	
  lastiteration <- dim(mcmcoutput$transtree)[2]
  if (index > lastiteration) stop("Error: invalid index specified")
  if (lastiteration == 0) stop("Error: need inferred exposure and infectious times for this function")
  epi <- cbind(mcmcoutput$nodeid[1:ninf], mcmcoutput$transtree[ ,index], mcmcoutput$exptimes[ ,index], 
               mcmcoutput$inftimes[ ,index], mcmcoutput$rectimes[1:ninf])
  colnames(epi) <- c("Node ID", "Parent", "Etime", "Itime", "Rtime" )
  return(epi)
}
