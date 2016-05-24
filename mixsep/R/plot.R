addBpDye <- function(x,kit=c("ID","Mini","Plus","COfiler","SGM","SE","Profiler","ESI17")){
  ## ABI = Applied Biosystems kit; Promega = Promega PowerPlex kit
  ## locus: Locus designation; offset: bp = offset+4*allele; dye: fluorescent dye/colour [B(lue), G(reen), Y(ellow), R(ed)]
  ## Offsets were computed based on precision data reported by ABI (online user manuals, www.appliedbiosystems.com)
  ## Loci implemented: Deselected below by loci[TRUE/FALSE]
  loci <- c("CSF1PO","D5S818","D7S820","D13S317","TPOX","D3S1358","D8S1179", 
            "D16S539","D18S51","D21S11","FGA","TH01","VWA","D2S1338","D19S433", ## Note: vWA is designated VWA
            "D1S1656","D2S441","D10S1248","D12S391","D22S1045","SE33","AMEL")
  kit <- kit[1]
  kits <- list("ID"=data.frame(locus=loci[c(rep(T,15),rep(F,6),T)], ## Identifiler, Plus & Direct (ABI)
                 offset=c(280,105,231,184,198,63,90,232,234,88,146,146,110,246,65,106), ## Amel-bp is for X, Y=X+6
                 dye=c("B","R","B","G","Y","G","B","G","Y","B","R","G","Y","G","Y","R")), 
               "Mini"=data.frame(locus=loci[c(T,F,T,T,F,F,F,T,T,T,T,F,F,T,rep(F,7),T)], ## MiniFiler (ABI)
                 offset=c(62,125,71,55,96,91,82,60,101),
                 dye=c("R","B","B","Y","Y","G","R","G","G")),
               "NGM"=data.frame(locus=loci[c(rep(F,5),rep(T,17))], ## NGM and NGMSElect (ABI)
                 offset=c(85,90,207,233,87,164,163,107,228,90,137,42,43,172,58,293,100),
                 dye=c("R","G","B","G","G","Y","Y","B","B","Y","R","R","B","R","Y","R","G")),
               "Plus"=data.frame(locus=loci[c(F,T,T,T,F,T,T,F,T,T,T,F,T,rep(F,8),T)], ## Profiler Plus & ID (ABI)
                 offset=c(103,232,173,63,91,234,89,144,110,103),
                 dye=c("Y","Y","Y","B","G","G","G","B","B","G")),
               "COfiler"=data.frame(locus=loci[c(T,F,T,F,T,T,F,T,F,F,F,T,rep(F,9),T)], ## Cofiler (ABI)
                 offset=c(256,232,191,63,209,147,103),
                 dye=c("G","Y","G","B","B","G","G")),
               "SGM"=data.frame(locus=loci[c(rep(F,5),rep(T,10),rep(F,6),T)], ## SGM Plus (ABI)
                 offset=c(63,91,209,234,89,144,147,110,230,66,103),
                 dye=c("B","G","B","G","G","Y","Y","B","B","Y","G")),
               "SE"=data.frame(locus=loci[c(rep(F,5),rep(T,10),rep(F,5),T,T)], ## SEfiler (ABI)
                 offset=c(61,90,207,233,89,142,143,107,228,64,183,102),
                 dye=c("B","G","B","R","R","Y","Y","B","B","Y","G","G")),
               "Profiler"=data.frame(locus=loci[c(rep(T,6),rep(F,4),rep(T,3),rep(F,8),T)], ## Profiler (ABI)
                 offset=c(256,103,232,173,191,63,144,147,110,103),
                 dye=c("G","Y","Y","Y","G","B","B","G","B","G")),
               "ESI17"=data.frame(locus=loci[c(rep(F,5),rep(T,17))], ## ESI17 (Promega)
                 offset=NA,
                 dye=c("B","R","G","G","Y","R","Y","Y","B","B","G","G","G","Y","B","R","B")))
  kits <- kits[[kit]]
  nx <- names(x)
  x <- merge(x,kits,by="locus")
  if(any(ame <- grepl("AME",toupper(x$locus)))){
    xame <- x[ame,] ## Amelogenin data
    xame <- merge(xame,data.frame(allele=c("X","Y"),bp=rep(xame$offset[1],2)+c(0,6)),by="allele")
    xame$offset <- NULL
    x <- x[!ame,] ## Autosomal data
    x$allele <- as.numeric(paste(x$allele))
  }
  else x$allele <- as.numeric(paste(x$allele))
  x$bp <- x$offset+4*( round(x$allele) + (x$allele-round(x$allele))*2.5 )
  x$offset <- NULL
  if(any(ame)) x <- rbind(x,xame)
  if(is.element("R",x$dye)) return(x[order(factor(x$dye,c("B","G","Y","R")),x$bp),c("locus","allele","height","area","bp","dye")])
  else return(x[order(factor(x$dye,c("B","G","Y")),x$bp),c("locus","allele","height","area","bp","dye")])
}

plotEPG <- function(x,color=TRUE,justdata=FALSE,addProfile=FALSE,profiles=NULL,contributor=NULL,...){
  kitcolor <- data.frame("dye"=c("B","G","Y","R"),baseline=seq(from=620,by=-200,len=4),
                         color=c("#1e7ba6","#2ea61e","#ffe51d","#ef2f2f")) ## yellow was: #ffef41, blue was: #0076c9, green was: #1eb04e
  kitcolor <- kitcolor[is.element(kitcolor$dye,x$dye),]
  while(min(kitcolor$baseline>20)) kitcolor$baseline <- kitcolor$baseline-200
  DD <- kitcolor$baseline
  ## Plot in colors
  if(color) x <- merge(x,kitcolor,by="dye",all.x=TRUE)
  else x$color <- "#999999"
  if(justdata) x$expHeight <- x$height ## No expected heights given for 'justdata'
  x$expHeight <- x$exp*(x$height/x$area) ## gives the same proportionality between hat(h) and hat(A) as for observed h and A
  ## was: mh <- max(c(x$area,x$exp),na.rm=TRUE)
  mh <- max(c(x$height,x$expHeight),na.rm=TRUE)
  ## was: x$plotObs <- x$area*170/mh
  x$plotObs <- x$height*170/mh
  ## was: x$plotExp <- x$exp*170/mh
  x$plotExp <- x$expHeight*170/mh
  rx <- range(x$bp*4,na.rm=TRUE)
  par(mar=c(0,0,0,0))
  topskip <- rep(c(0,40,80,120),each=7) ## to make room for profiles in top of plot
  plot(c(rx[1]-40,rx[2]+40),c(min(x$baseline,.na.rm=TRUE)-20,max(x$baseline,na.rm=TRUE)+170+topskip[length(profiles)]*addProfile),type="n",xlab="",ylab="",axes=FALSE,...)
  N <- nrow(x)
  yy <- rep(x$baseline,each=4)
  yy[seq(from=4,by=4,len=N)] <- NA
  obsPeak <- yy
  expPeak <- yy
  obsPeak[seq(from=2,by=4,len=N)] <- obsPeak[seq(from=2,by=4,len=N)]+x$plotObs
  expPeak[seq(from=2,by=4,len=N)] <- expPeak[seq(from=2,by=4,len=N)]+x$plotExp
  xx <- rep(x$bp,each=4)*4+rep(c(-8,0,8,NA),N)
  ## Adding 'grid' lines to the plot indicating the actual intensities of the peaks
  gl <- seq(from=250,to=mh,by=250)
  if(length(gl)<4) gl <- seq(from=100,to=mh,by=100)
  if(length(gl)>6) gl <- gl[floor(seq(from=1,to=length(gl),len=6))]
  gl <- c(50,gl)
  gridLines <- gl*170/mh
  abline(h=(hh <- rep(gridLines,length(DD))+rep(DD,each=length(gridLines))),col="#efefef")
  text(rep(rx[1]-30,length(hh)),hh,gl,cex=0.60,adj=c(1,0.5))
  polygon(xx,obsPeak,col=paste(x$color),border=NA)
  if(!justdata) polygon(xx,expPeak,border=1,lty=1,lwd=1)
  abline(h=DD)
  ## Adding allele and locus designations to plot
  loci <- aggregate(cbind(bp,baseline)~locus,data=x,FUN=mean)
  text(loci$bp*4,loci$baseline-25,paste(loci$locus),cex=0.75,adj=c(0.5,1)) ## was: -15
  text(x$bp*4,x$baseline-17,paste(x$allele),cex=0.75) ## was: -8
  #box()
  if(addProfile){ ## add profile table to plot
    addprofiles2plot(x=mean(rx),y=par("usr")[4],profiles=profiles,contributor=contributor,just=c(0.5,-0.6))
    ## Jakob plot
    addHooks2plot(profiles,data=x[,c("locus","allele","bp","baseline")])
    ## 
  }
}

computeExpArea <- function(x,y,tauhat){
  alts <- apply(x$result$profiles,2,paste,collapse="/")
  alts <- rbind(alts,x$result$alternatives)
  comb <- diag(alts[as.numeric(y)+1,1:length(y)])
  names(comb) <- dimnames(x$result$profiles)[[2]]
  comb <- lapply(comb,function(z){
    pp <- lapply(strsplit(z,"/"),strsplit,split=",")[[1]]
    alleles <- sort(unique(unlist(pp)))
    loc <- as.data.frame(cbind(alleles,matrix(0,length(alleles),length(pp))))
    names(loc) <- c("allele",paste("P",1:length(pp),sep=""))
    for(i in 1:length(pp)) loc[,i+1] <- (loc$allele==pp[[i]][1])+(loc$allele==pp[[i]][2])
    loc
  }) 
  comb <- as.data.frame(cbind(locus=rep(names(comb),unlist(lapply(comb,nrow))),do.call("rbind",comb)))
  data <- x$result$data #[,datacols]
  names(data) <- c("locus","allele","height","area")
  d <- merge(data,comb,by=c("locus","allele"))
  ds <- split(d,d$locus)
  m <- ncol(comb)-2
  N <- nrow(data)-length(ds)-(m-1)
  aa <- ahat(ds,m=m)
  tt <- that(ds,alpha=aa,m=m)
  R2 <- tauhat/tt
  r2 <- (tauhat/tt)^N
  expArea <- expectedAreas(split(d,d$locus),alpha=aa,m=m)
  list(data=expArea[,c("locus","allele",paste("P",1:m,sep=""),"exp")],alpha=aa,tau=tt,R2=R2,r2=r2,N=N)
}

## Modified from function addtable2plot in package 'plotrix'
addprofiles2plot <- function (x, y, profiles, contributor, cex = c(0.85,1), just = c(0,1)) {
  tabdim <- dim(profiles)
  if(is.null(tabdim)){ ## profiles is a vector
    loci <- names(profiles)
    profiles <- do.call("cbind",strsplit(profiles,"/"))
    colnames(profiles) <- loci
    rownames(profiles) <- ifelse(rep(is.null(contributor),nrow(profiles)),"",contributor)
    tabdim <- dim(profiles)
  }
  else{ ## profiles is a table
    loci <- colnames(profiles)
    contributor <- rownames(profiles)
  }
  if(is.null(tabdim)) return(NULL) ## something went wrong (profiles is neither table nor vector)
  mwidth <- strwidth("M", cex = cex[1])
  cellwidth <- max(strwidth(c(loci, contributor, as.vector(unlist(profiles))), cex = cex[1])) + mwidth
  nvcells <- tabdim[1] + 1
  nhcells <- tabdim[2] + 1
  cellheight <- max(strheight(c(loci, contributor, as.vector(unlist(profiles))), cex = cex[1])) * 2
  xleft <- x - just[1] * nhcells * cellwidth
  ytop <- y + just[2] * nvcells * cellheight
  for (row in 1:tabdim[1]) {
    if (row <= nvcells - 1){ 
      segments(xleft, ytop - row * cellheight, xleft + nhcells * cellwidth, ytop - row * cellheight, lwd = 1, col = 1)
      text(xleft + 0.5 * cellwidth, ytop - (row + 0.5) * cellheight, contributor[row], cex = cex[1], col = row)
      for (col in 1:tabdim[2]){
        text(xleft + (col + 0.5) * cellwidth, ytop - (row + 0.5) * cellheight, profiles[row, col], cex = cex[1], col = 1)
        text(xleft + (col + 0.5) * cellwidth, ytop - 0.5 * cellheight, loci[col],  cex = cex[1], col = 1, font = 2)
        segments(xleft, ytop - cellheight, xleft + nhcells * cellwidth, ytop - cellheight, lwd = 1, col = 1)
      }
    }
  }
  text(xleft + (nhcells * cellwidth)/2, ytop + cellheight/2, "Plotted profiles", cex = cex[2], col = 1, font = 2)
  segments(xleft, ytop, xleft + nhcells * cellwidth, ytop, lwd = 2, col = 1)
  segments(xleft, ytop - nvcells * cellheight, xleft + nhcells * cellwidth, ytop - nvcells * cellheight, lwd = 2, col = 1)
}

profileHook <- function(x,data){ ## x, kit=
  allele <- as.data.frame(cbind(locus=names(x),do.call("rbind",strsplit(x,","))))
  allele1 <- allele[,c(1,2)] ## locus and first allele
  allele2 <- allele[,c(1,3)] ## locus and second allele
  names(allele1) <- names(allele2) <- c("locus","allele")
  allele1 <- merge(allele1,data,by=c("locus","allele"),all.x=TRUE) ## allele2bp(allele1,kit=kit)
  allele2 <- merge(allele2,data,by=c("locus","allele"),all.x=TRUE) ## allele2bp(allele2,kit=kit)
  names(allele1) <- paste(names(allele1),"1",sep="")
  names(allele2) <- paste(names(allele2),"2",sep="")
  m.var <- c("locus","baseline")
  allele <- merge(allele1,allele2,by.x=paste(m.var,"1",sep=""),by.y=paste(m.var,"2",sep=""))
  names(allele) <- c("locus","baseline","allele1","bp1","allele2","bp2")
  convertTab(convertTab(allele),mode="num",subset=c("bp","baseline"))
}

addHooks2plot <- function(x,data){ 
  x <- do.call("cbind",strsplit(x,"/"))
  x <- apply(x,1,profileHook,data=data)
  m <- length(x)
  bpskip <- ifelse(m==rep(2,m),c(-2,2),c(-2.5,0,2.5))
  for(j in 1:m){
    for(i in 1:nrow(x[[j]])){
      with(x[[j]][i,],lines(c(bp1,bp1,bp2,bp2)*4+bpskip[j],baseline-c(0,4,4,0)*j,col=j,lwd=2))
      if(x[[j]]$allele1[i]==x[[j]]$allele2[i]) with(x[[j]][i,],points(bp1*4+bpskip[j],baseline-4*j,pch=16,cex=0.8,col=j)) ## hom
    }
  }
}
