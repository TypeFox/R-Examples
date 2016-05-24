# CHNOSZ/anim.R
# functions to make various animations

anim.TCA <- function(redox=list(O2=c(-95,-60)),high.T=FALSE,
  nframes=100,pHlim=c(0,10),width = 420, height = 320) {
  # tricarboxylic acids 20100701, 20110204 jmd
  # we depend on an empty png directory
  if(!"png" %in% dir()) stop("directory 'png' not present")
  else if(length(dir("png")) > 0) stop("directory 'png' not empty")
  # add supplementary data (from default location of data/OBIGT-2.csv)
  # which includes properties for the metabolites
  add.obigt()
  # expand default logfO2 range if we're at high temperature
  if(high.T & missing(redox)) redox <- list(O2=c(-100,-40))
  # the name of 'redox' should be either O2 or H2
  basis(c("CO2","H2O",names(redox),"H+"))
  # load the species in order of increasing degree of ionization
  species(c("pyruvic acid","oxaloacetic acid","malic acid", # Z = 0
    "fumaric acid","a-ketoglutaric acid","citric acid",
    "pyruvate","H-oxaloacetate","H-malate","H-fumarate",    # Z = -1
    "H-a-ketoglutarate","H2-citrate",
    "oxaloacetate-2","malate-2","fumarate-2",               # Z = -2
    "a-ketoglutarate-2","H-citrate-2",
    "citrate-3"),rep("aq",18))                              # Z = -3
  # start the plot device - multiple png figures
  png(filename="png/Rplot%04d.png",width=width,height=height)
  par(mar=c(3,3.5,1.5,1),mgp=c(1.9,1,0))
  # leadin/out frames
  # the variable(s) that change with each frame
  pH <- seq(pHlim[1],pHlim[2],length.out=nframes)
  pHres <- length(pH)
  pHlim <- range(pH)
  baseres <- 128
  # calculate the affinity
  redox[[1]] <- c(redox[[1]],baseres)
  affargs <- c(redox,list(H2O=c(-10,10,baseres),pH=c(pHlim,pHres)))
  a.cold <- do.call(affinity,affargs)
  if(high.T) {
    affargs <- c(affargs,list(T=100))
    a.hot <- do.call(affinity,affargs)
  }
  # put the movie together
  nframes.in <- 15
  nframes.out <- 20
  iframes <- c(rep(1,nframes.in),1:length(pH),rep(pHres,nframes.out))
  for(i in 1:length(iframes)) {
    print(paste("frame",i,"of",length(iframes),"at pH",pH[iframes[i]]))
    # don't fill fields with heat colors if we're doing a high-T overlay
    if(high.T) color <- "lightgrey" else color <- "heat"
    # take a single O2-H2O slice along pH
    diagram(slice.affinity(a.cold,3,iframes[i]),fill=color,cex=1.5)
    ltext <- substitute(paste(italic(T)==x~degree*C),list(x=25))
    lcol <- "black"
    # overlay the high-T diagram
    if(high.T) {
      ltext <- c(ltext,substitute(paste(italic(T)==x~degree*C),list(x=100)))
      lcol <- c(lcol,"red")
      diagram(slice.affinity(a.hot,3,iframes[i]),add=TRUE,col="red",cex=1.5)
    }
    if(high.T) legend("topleft",legend=as.expression(ltext),lty=1,col=lcol)
    title(main=paste("TCA cycle reactants; pH =",
      format(round(pH[iframes[i]],2))),cex.main=1)
  }
  # close the plot device - convert to movie
  dev.off()
  cat("anim.TCA: converting to animated GIF...\n")
  outfile <- "TCA.gif"
  # convert tool from imagemagick.org; -loop 0 for infinite loop
  syscmd <- paste("convert -loop 0 -delay 10 png/*.png png/",outfile,sep="")
  cat(paste(syscmd,"\n"))
  # we need shell() on windows for the paths to include ImageMagick
  if(.Platform$OS.type=="unix") sres <- system(syscmd)
  else sres <- shell(syscmd)
  if(sres==0) cat(paste("anim.TCA: animation is at png/",outfile,"\n",sep=""))
  else {
    cat("anim.TCA: error converting to animated GIF\n")
    cat("anim.TCA: check that 'convert' tool from ImageMagick is in your PATH\n")
  }
}

anim.plasma <- function(width=480, height=480) {
  ## animate relative metastabilities of plasma proteins,
  ## as a function of chemical activities of H2 and O2
  ## start with 71 proteins; one protein is removed in 
  ## each successive frame
  ## 20110804 jmd
  # make sure we have an empty png directory
  if(!"png" %in% dir()) stop("directory 'png' not present")
  else if(length(dir("png")) > 0) stop("directory 'png' not empty")
  # get list of proteins
  f <- system.file("extdata/abundance/AA03.csv", package="CHNOSZ")
  pdata <- read.csv(f, as.is=TRUE)
  notna <- !is.na(pdata$name)
  pname <- pdata$name[notna]
  # set up the system; use O2 aq instead of gas
  basis(c("CO2","NH3","H2S","H2","O2","H+"))
  basis("O2","aq")
  basis(c("CO2","NH3","H2S","H+"),c(-3,-3,-10,-7))
  species(pname,"HUMAN")
  a <- affinity(H2=c(-20,0), O2=c(-80,-60))
  # start with all species
  ispecies <- 1:length(pname)
  # open plot file
  png("png/plasma%03d.png",width=width, height=height)
  # numbers of proteins for each frame
  np <- length(pname):2
  # add some lead-in and lead-out frames
  np <- c(rep(head(np,1),8),np,rep(tail(np,1),8))
  # now loop until we get to two proteins
  for(i in 1:length(np)) {
    d <- diagram(a, groups=as.list(ispecies), names=pname[ispecies], normalize=TRUE, cex=1.5)
    # note that the darker colors go with higher abundances
    # as reported by Anderson and Anderson, 2003
    # how many show up on the diagram?
    nshow <- length(unique(as.numeric(d$predominant)))
    title(main=paste(np[i]," human plasma proteins (",nshow,
      " showing)",sep=""),cex.main=1)
    if(c(0,diff(np))[i] != 0) {
      # identify the protein with the greatest
      # area on the diagram
      imost <- which.max(tabulate(as.numeric(d$predominant)))
      # take out that proteins
      ispecies <- ispecies[-imost]
    }
  }
  # close PNG plot device
  dev.off()
  # make animated GIF using ImageMagick
  cat("anim.plasma: converting to animated GIF...\n")
  outfile <- "plasma.gif"
  syscmd <- paste("convert -loop 0 -delay 25 png/*.png png/", outfile, sep = "")
  cat(paste(syscmd,"\n"))
  if(.Platform$OS.type=="unix") sres <- system(syscmd)
  else sres <- shell(syscmd)
  if(sres==0) cat(paste("anim.plasma: animation is at png/",outfile,"\n",sep=""))
  else {
    cat("anim.plasma: error converting to animated GIF\n")
    cat("anim.plasma: check that 'convert' tool from ImageMagick is in your PATH\n")
  }
}

anim.carboxylase <- function(T=25:125,ntop=5,lcex=0.8,width=420,height=320) {
  # animate rank-activity diagrams over a temperature
  # and logaH2 gradient, or plot a single one for a single temperature
  # plot rank-activity diagrams for 24 carboxylases;
  # 12 ribulose phosphate carboxylase
  # 12 acetyl-coenzyme A carboxylase
  # 6 of each type are nominally from mesophilic organisms
  # and 6 from thermophilic organisms
  # arranged here in order of increasing growth temperature
  rubisco <- c("RBL_BRAJA","A6YF84_9PROT","A1E8R4_9CHLO","A8C9T6_9MYCO","A3EQE1_9BACT","A5CKC7_9CHRO",
    "RBL_SYNJA","Q6JAI0_9RHOD","RBL_METJA","A3DND9_STAMF","A1RZJ5_THEPD","RBL_PYRHO")
  rubisco.organisms <- c("a-proteobacterium-R","b-proteobacterium","Bracteacoccus","Mycobacterium",
    "Leptospirillum","Cyanobium","Synechococcus","Cyanidiales",
    "Methanococcus-R","Desulfurococcus","Thermofilum","Pyrococcus")
  accoaco <- c("Q9F7M8_PRB01","ACCA_DEIRA","A6CDM2_9PLAN","A4AGS7_9ACTN","ACCA_CAUCR","A1VC70_DESVV",
    "A6VIX9_METM7","Q2JSS7_SYNJA","A0GZU2_9CHLR","A7WGI1_9AQUI","Q05KD0_HYDTH","ACCA_AQUAE")
  accoaco.organisms <- c("g-proteobacterium","Deinococcus","Planctomyces","Actinobacterium",
    "a-proteobacterium-A","d-proteobacterium","Methanococcus-A","Synechococcus",
    "Chloroflexus","Hydrogenobaculum","Hydrogenobacter","Aquifex")
  # assemble them all
  organisms <- c(rubisco.organisms,accoaco.organisms)
  # new scheme 20090611: red for hot, blue for cold
  # open for rubisco, filled for accoaco
  col <- rep(c(rep("blue",6),rep("red",6)),2)
  pch <- c(rep(c(0:2,5:7),2),rep(c(15:20),2))
  # how many frames do we want?
  res <- length(T)
  if(res==1) ido <- 1
  else {
    # check for png directory
    if(!"png" %in% dir()) stop("directory 'png' not present")
    else if(length(dir("png")) > 0) stop("directory 'png' not empty")
    # start the plot device - multiple png figures
    png(filename="png/Rplot%04d.png",width=width,height=height)
    # add counters for lead-in and lead-out frames
    ido <- c(rep(1,15),1:res,rep(res,20))
  }
  # set up system
  basis(c("CO2","H2O","NH3","H2","H2S","H+"),
    c("aq","liq","aq","aq","aq","aq"),c(-3,0,-4,-6,-7,-7))
  species(c(rubisco,accoaco))
  # equation for logaH2 as a function of temperature
  # from Dick and Shock, 2011
  # http://dx.plos.org/10.1371/journal.pone.0022782
  get.logaH2 <- function(T) return(-11+T*3/40)
  H2 <- get.logaH2(T)
  # calculate affinities
  if(res==1) {
    basis("H2",H2)
    a <- affinity(T=T)
  } else a <- affinity(T=T,H2=H2)
  # calculate activities
  e <- equilibrate(a, normalize=TRUE)
  # for each point make a rank plot
  rank <- 1:length(e$loga.equil)
  for(i in 1:length(ido)) {
    # print some progress
    if(i%%20 == 0) cat("\n") else cat(".")
    # keep track of positions of previous points
    loga <- numeric()
    for(j in 1:length(e$loga.equil)) loga <- c(loga, e$loga.equil[[j]][ido[i]])
    if(i > 4) myrank4 <- myrank3
    if(i > 3) myrank3 <- myrank2
    if(i > 2) myrank2 <- myrank1
    if(i > 1) myrank1 <- myrank
    order <- order(loga,decreasing=TRUE)
    myrank <- rank(loga)
    cex <- rep(1.2,24)
    # show changes by increasing point size
    # any points that changed on the step before the step 
    # before the step before?
    if(i > 4) {
      ichanged <- myrank3 != myrank4
      cex[ichanged[order]] <- cex[ichanged[order]] + 0.1
    }
    # any points that changed on the step before the step before?
    if(i > 3) {
      ichanged <- myrank2 != myrank3
      cex[ichanged[order]] <- cex[ichanged[order]] + 0.2
    }
   # any points that changed on the step before?
    if(i > 2) {
      ichanged <- myrank1 != myrank2
      cex[ichanged[order]] <- cex[ichanged[order]] + 0.3
    }
    # any points that changed on this step?
    if(i > 1) {
      ichanged <- myrank != myrank1
      cex[ichanged[order]] <- cex[ichanged[order]] + 0.4
    }
    plot(rank,loga[order],col=col[order],pch=pch[order],
      ylab=expression(log~italic(a)),cex=cex,cex.main=1,cex.lab=1,cex.axis=1)
    myT <- format(round(T,1))[ido[i]]
    myH2 <- format(round(H2,2))[ido[i]]
    title(main=substitute(list(X~degree*C, log*italic(a)[paste(H2)]==Y),
      list(X=myT,Y=myH2)))
    # legends showing highest and lowest few
    legend("topright",legend=c(paste("top",ntop),organisms[order[1:ntop]]),
      pch=c(NA,pch[order[1:ntop]]),col=c(NA,col[order[1:ntop]]),
      pt.cex=c(NA,cex[1:ntop]),cex=lcex)
    order <- order(loga)
    legend("bottomleft",legend=c(paste("low",ntop),organisms[order[ntop:1]]),
      pch=c(NA,pch[order[ntop:1]]),col=c(NA,col[order[ntop:1]]),
      pt.cex=c(NA,cex[24:(24-ntop+1)]),cex=lcex)
  }
  # finish up animation stuff
  if(res > 1) {
    # finish progress report
    cat("\n")
    # close PNG plot device
    dev.off()
    # make animated GIF using ImageMagick
    cat("anim.carboxylase: converting to animated GIF...\n")
    outfile <- "carboxylase.gif"
    syscmd <- paste("convert -loop 0 -delay 10 png/*.png png/", outfile, sep = "")
    cat(paste(syscmd,"\n"))
    if(.Platform$OS.type=="unix") sres <- system(syscmd)
    else sres <- shell(syscmd)
    if(sres==0) cat(paste("anim.carboxylase: animation is at png/",outfile,"\n",sep=""))
    else {
      cat("anim.carboxylase: error converting to animated GIF\n")
      cat("anim.carboxylase: check that 'convert' tool from ImageMagick is in your PATH\n")
    }
  }
}


