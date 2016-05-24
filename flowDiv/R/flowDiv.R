#' @title Cytometric Diversity Indices from 'FlowJo' Workspaces
#'
#' @description Concatenates some 'flowWorkspace','flowCore' 'vegan' and 'gdata' packages functionalities to import 'FlowJo' workspaces and calculates ecological diversity indices for gated populations, based on bidimensional cytograms. Each detection channel is clustered into “n” bins, resulting in n x n categories, similarly to the methodology proposed by Li (1997) for estimating phytoplankton cytometric diversity. Diversity indices are calculated based on the number of cells per category from each sample.
#' @author Bruno M.S. Wanderley, María Victoria Quiroga, André M. Amado, Fernando Unrein
#'
#' @usage flowDiv(myworkspaces=list(), gate.name=NULL, channel1=NULL,
#' channel2=NULL, nbins=16, dil=c(), flowrate=c(), volume=c(), nsamples=NULL,
#' use.beads=FALSE, beads.gate=NULL, doplot=TRUE, file.name=NULL, save.csv=FALSE,
#' ialpha="invsimpson", ibeta="bray", runalpha="TRUE", runbeta="TRUE")
#'
#'
#' @return A list containing alpha index, beta matrix and Pielou's indices for each cytogram.
#'
#'
#' @param myworkspaces  A list containing the paths to FlowJo workspaces which are meant to be analyzed. More than one workspace can be analyzed at the same time. Workspaces should contain .fcs files (versions 2.0 or 3.0) with its original names.
#' @param gate.name Name of the gate to be analyzed. Must be a single-valued string. The gate should be named exactly the same in all samples from the workspaces.
#' @param channel1 Name of channel to be plot as y-axis. Channel name should be exactly the same as that assigned by flow cytometer. Note that this is equipment dependent and can vary.
#' @param channel2 Name of channel to be plot as x-axis. Channel name should be exactly the same as that assigned by flow cytometer. Note that this is equipment dependent and can vary.
#' @param nbins Number of bins to cluster each channel. Default=16.
#' @param dil A vector containing dilution factors for each sample in the workspaces. Its lenght is the same as the total number of samples meant to be analyzed and must follow the exact same order in which samples are presented to the function (the order of samples corresponds both to the order of workspace described in "myworkspaces" parameter and their order in each of these FlowJo workspaces). By default, it is assumed that all samples have no dilutions whatsoever (i.e. all dilutions factors equal 1).
#' @param flowrate A vector containing volumetric flow rate for each sample in the workspaces. Its lenght is the same as the total number of samples meant to be analyzed and must follow the exact same order in which samples are presented to the function (the order of samples corresponds both to the order of workspace described in "myworkspaces" parameter and their order in each of these FlowJo workspaces). Not necessary if "volume" is declared.
#' @param volume A vector containing the total liquid volumes (i.e. flow rate x ellapsed time) of each sample in the workspaces. Its lenght is the same as the total number of samples meant to be analyzed and must follow the exact same order in which samples are presented to the function (the order of samples corresponds both to the order of workspace described in "myworkspaces" parameter and their order in each of these FlowJo workspaces). Not necessary if "flowrate" is declared.
#' @param nsamples Total number of samples to be analyzed. If there is more than one workspace, that number corresponds to the summation of the samples in each workspace altogether.
#' @param use.beads Logical. If “TRUE” , it brings all cytograms to a common point, based on the arithmetic mean of a standard region for all cytograms (usually beads regions), before proceeding to analysis. It is recommended to proceed this way only if samples were analyzed with different settings (i.e. voltages).
#' @param beads.gate Name of the gate describing the standard region to be used (usually bead's regions). Necessary only if use.beads is set to “TRUE”. The beads gate should be named exactly the same in all samples from the workspaces.
#' @param doplot Logical. It plots scatterplots of gated populations and displays grid lines corresponding to the limits of the bins.
#' @param file.name Name of the file in which diversity results should be stored. Required only if save.csv is set to “TRUE”.
#' @param save.csv Logical. Used to export diversity results as .csv file.
#' @param ialpha Method used to calculate alpha diversity index of cytograms (i.e. the degree of similarity between two cytograms). Should be one of "shannon", "simpson" or "invsimpson", as for vegan::diversity function. Default is "invsimpson".
#' @param ibeta Method used to calculate beta diversity index of cytograms. Should be one of  the sixteeen avaiable for vegan::vegdist function. Default is “bray”.
#' @param runalpha Logical. Indicates if alpha index should be calculated. Default is “TRUE”.
#' @param runbeta Logical. Indicates if beta index should to be calculated. Default is “TRUE”.
#'
#' @references Li, W.K.W. (1997). Cytometric diversity in marine ultraphytoplankton. Limnology and Oceanography 42, 874–880.

#' @examples \dontrun{
#'
#' ### Using one workspace ###
#' ## Not run:
#'# Analyzing a .xml FlowJo workspace containg 23 samples using channels FITC-H and SSC-H.
#'
#' indexes.sw <- flowDiv(myworkspaces = list(“my_flowjo_workspace.xml”),
#' gate.name = “my_gate_name”, channel1 = “FITC-H”, channel2 = “SSC-H”, nsamples=23)
#' #assume that the FlowJo workspace is below the current directory
#'
#' ## End (not run)
#'
#'### Using multiple workspaces ###
#'
#'## Not run:
# Analyzing three .wsp FlowJo workspaces (each one containing twelve samples) based on channels FL1-H and SSC-H with and a gate named “Bacteria”.
#'
#'indexes.mw <- flowDiv(myworkspaces = list(”my_flowjo_workspace1.wsp”,
#'“my_flowjo_workspace2.wsp”, “my_flowjo_workspace3.wsp”),
#'gate.name = ”Bacteria”, channel1 = ”FL1-H”, channel2 = ”SSC-H”, nsamples=36)
#'#assume that the FlowJo workspaces are below the current directory
#'
#' ## End (not run)
#'}


flowDiv <- function (myworkspaces=list(), gate.name=NULL, channel1=NULL, channel2=NULL, nbins=16, dil=c(), flowrate=c(), volume=c(), nsamples=NULL, use.beads=FALSE, beads.gate=NULL, doplot=TRUE, file.name=NULL, save.csv=FALSE, ialpha="invsimpson", ibeta="bray", runalpha="TRUE", runbeta="TRUE"){


  samples=list()
  MatList = list()
  id.files=c()
  mean_chan1=c()
  mean_chan2=c()
  timediff.all=c()
  total.flow = c()
  timediff.all=c()
  if(is.null(dil)) dil=rep(1, nsamples)

  for(i in 1:length(myworkspaces)){




    workspace<-flowWorkspace::openWorkspace(myworkspaces[[i]])


    gating.set<-flowWorkspace::parseWorkspace(workspace, name=1)

    flowList = list()
    matrixList = list()
    id = c()
    beadsList = list()
    fmean_chan1=c()
    fmean_chan2=c()
    c1.invTrans = list()
    c2.invTrans = list()
    versions=c()
    invstate="off"


    etim=c()
    btim=c()
    timediff=c()


    for (z in 1:length(gating.set))
    {

      id[z] = sub("\\.fcs", "", tail(unlist(strsplit(flowWorkspace::getData(gating.set[[z]])@description$FILENAME, split="\\/")), n=1))



      if(!is.null(flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)))
      {
        c1.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)
        invstate="on"
      }
      else  c1.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)
      if(!is.null(flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)))
      {
        c2.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)
      }
      else  c2.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)

      versions[z]<-as.numeric(flowWorkspace::getData(gating.set[[z]])@description$FCSversion)
      etim = as.POSIXlt(flowWorkspace::getData(gating.set[[z]])@description$`$ETIM`,format="%H:%M:%S")
      btim = as.POSIXlt(flowWorkspace::getData(gating.set[[z]])@description$`$BTIM`,format="%H:%M:%S")
      timediff[z] = as.numeric(difftime(etim, btim, units="min"))

      nodelist<-flowWorkspace::getNodes(gating.set[[z]], path = 1)
      node<-nodelist[which(nodelist==gate.name)]
      flowList[[z]] = flowWorkspace::getData(gating.set[[z]],node)

      if(use.beads)
      {

        bds<-nodelist[which(nodelist==beads.gate)]
        beadsList[[z]] = flowWorkspace::getData(gating.set[[z]],bds)

      }
    }

    flowWorkspace::closeWorkspace(workspace)


    for (p in 1:length(flowList)){

      if(invstate=="on")
      {

        if(use.beads)
        {


          while(i<2)
          {

            for (e in 1:length(beadsList)) {

              if(mean(versions)==2)
              {

                mean_chan1[e] = mean(log10(c1.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel1])))
                mean_chan2[e] = mean(log10(c2.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel2])))

              }

              else
              {

                mean_chan1[e] = mean(log10(c1.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel1])/26.2144))
                mean_chan2[e] = mean(log10(c2.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel2])/26.2144))

              }

            }

          }







          fmean_chan1[p] = mean_chan1[p]-mean(mean_chan1)
          fmean_chan2[p] = mean_chan2[p]-mean(mean_chan2)


          if(mean(versions)==2)
          {

            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])) + fmean_chan1[p]
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])) + fmean_chan2[p]

          }

          else
          {

            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])/26.144) + fmean_chan1[p]
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])/26.144) + fmean_chan2[p]

          }

        }
        else
        {
          if(mean(versions)==2)
          {
            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1]))
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2]))
          }

          else

          {
            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])/26.2144)
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])/26.2144)

          }

        }

        if(doplot){

          plot(chan2, chan1, pch=".", xlim=c(0,4), ylim=c(0,4), main=id[p], xlab=channel2, ylab=channel1, cex=3, xaxs="i",yaxs="i" )
          grid(nx=nbins, col="blue")

        }



        chan1.cuts <- seq(from = 0, to = 4, length = nbins + 1)
        chan2.cuts <- seq(from = 0, to = 4, length = nbins +1)
        index.chan1 <- cut(chan1, chan1.cuts, include.lowest = TRUE)
        index.chan2 <- cut(chan2, chan2.cuts, include.lowest = TRUE)
        m <- tapply(chan2, list(index.chan1, index.chan2), base::length)
        m[is.na(m)] <- 0
        matrixList[[p]] = m

      }

      if(invstate=="off")
      {

        if(use.beads)
        {


          while(i<2)
          {

            for (e in 1:length(beadsList)) {

              if(mean(versions)==2)
              {

                mean_chan1[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel1]))
                mean_chan2[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel2]))

              }

              else
              {

                mean_chan1[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel1])/26.2144)
                mean_chan2[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel2])/26.2144)

              }

            }

          }







          fmean_chan1[p] = mean_chan1[p]-mean(mean_chan1)
          fmean_chan2[p] = mean_chan2[p]-mean(mean_chan2)


          if(mean(versions)==2)
          {

            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]) + fmean_chan1[p]
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]) + fmean_chan2[p]

          }

          else
          {

            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]/26.144) + fmean_chan1[p]
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]/26.144) + fmean_chan2[p]

          }

        }
        else
        {
          if(mean(versions)==2)
          {
            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1])
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2])
          }

          else

          {
            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]/26.2144)
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]/26.2144)

          }

        }

        if(doplot){

          plot(chan2, chan1, pch=".", xlim=c(0,4), ylim=c(0,4), main=id[p], xlab=channel2, ylab=channel1, cex=3, xaxs="i",yaxs="i" )
          grid(nx=nbins, col="blue")

        }



        chan1.cuts <- seq(from = 0, to = 4, length = nbins + 1)
        chan2.cuts <- seq(from = 0, to = 4, length = nbins +1)
        index.chan1 <- cut(chan1, chan1.cuts, include.lowest = TRUE)
        index.chan2 <- cut(chan2, chan2.cuts, include.lowest = TRUE)
        m <- tapply(chan2, list(index.chan1, index.chan2), base::length)
        m[is.na(m)] <- 0
        matrixList[[p]] = m
      }


    }

    for (k in 1:length(matrixList))
    {
      matrixList[[k]] = as.matrix(as.data.frame(matrixList[[k]]))

    }

    MatList[[i]]=matrixList
    id.files = c(id.files,id)
    timediff.all = c(timediff.all,timediff)

  }



  for(y in 1:length(MatList))
  {

    for(w in 1:length(MatList[[y]]))

    {

      samples[[length(samples)+1]]= MatList[[y]][[w]]

    }

  }

  if(!is.null(flowrate)&&!is.null(volume)) volume=NULL

  if(is.null(flowrate)) vol1=rep(1, nsamples)
  else vol1=max(timediff.all*flowrate)/(timediff.all*flowrate)
  if (is.null(volume)) vol2=rep(1, nsamples)
  else vol2=max(volume)/volume



  for(d in 1:length(samples))

    samples[[d]] = vol1[d]*vol2[d]*dil[d]*samples[[d]]



  unmatrix <-function (x, byrow = FALSE)
  {
    rnames <- rownames(x)
    cnames <- colnames(x)
    if (is.null(rnames))
      rnames <- paste("r", 1:nrow(x), sep = "")
    if (is.null(cnames))
      cnames <- paste("c", 1:ncol(x), sep = "")
    nmat <- outer(rnames, cnames, paste, sep = ":")
    if (byrow) {
      vlist <- c(t(x))
      names(vlist) <- c(t(nmat))
    }
    else {
      vlist <- c(x)
      names(vlist) <- c(nmat)
    }
    return(vlist)
  }



  if (length(samples)==1)runbeta="FALSE"




  if (runalpha){
    alpha <- c()
    pielou <- c()
    for(w in 1:length(samples)){
      samps <- unmatrix(samples[[w]])
      alpha[w]<- vegan::diversity(samps, index=ialpha)
      pielou[w]<- vegan::diversity(samps)/log(vegan::specnumber(samps))
    }
  }


  if(runbeta){

    matriz <- matrix(ncol=length(samples), nrow=length(samples[[1]]))
    colnames(matriz) <- id.files

    for(g in 1:length(samples)) matriz[,g] <- unmatrix(samples[[g]], byrow="FALSE")


    beta = vegan::vegdist(t(matriz), method=ibeta)
    beta = as.matrix(beta)
  }


  indices <- list()

  if(runalpha){
    indices$alpha=alpha
    names(indices$alpha)=id.files
    indices$pielou=pielou
    names(indices$pielou)=id.files

  }


  if(runbeta){
    indices$beta=beta

  }

  if(save.csv) write.csv(indices, file=file.name, row.names=TRUE)

  return(indices)

}
