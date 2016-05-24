#' Function for reading Bruker OPUS files from MPA, HTS-xt and Alpha
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}
#' Note: Tested using several thousands raw spectra from three instruments. More testing is needed however.

if(!isGeneric("read.opus")){
  setGeneric("read.opus", function(file.name, ...){standardGeneric("read.opus")})
}

.read.opus <- function(file.name, sp=NULL, codes=c("ZFF","RES","SNM","DAT","LWN","FXV","LXV","NPT","MXY","MNY","END","TIM"), plot.spectra=FALSE, print.progress=FALSE, speclib="ICRAF", signif.digit=get("signif.digit", spec.opts), MID, st.wavenumbers=wavenumbers)
{
  if(!(speclib=="ICRAF"|speclib=="New")){ stop("'speclib' must be one of the following: 'ICRAF' or 'New'") }
  if(file.exists(file.name)){
    ## Read metainfo
	  try( pa <- hexView::readRaw(file.name, offset = 0, nbytes = file.info(file.name)$size, human = "char", size = 1, endian = "little"), silent=TRUE )
      if(!class(.Last.value)[1]=="try-error"){
        pr <- pa$fileRaw
	   		## Get source of instrument
		  	ins <- grepRaw("INS", pr, all=TRUE)
		 	ins <- readRaw(file.name, offset = ins[length(ins)]+7, nbytes = 3, human = "char", size = 1, endian = "little")
		 	ins <- blockString(ins)
		 	## Get source of infrared to know if NIR or MIR
		    src <- grepRaw("SRC", pr, all=TRUE)
		    src <- readRaw(file.name, offset = src[length(src)]+4, nbytes = 3, human = "char", size = 1, endian = "little")
		    src <- blockString(src)
  		  	instr.range <- tolower(paste(ins, src, sep="-"))
  		  ## Get Beam Splitter
  		 	bms <- grepRaw("BMS", pr, all=TRUE)
		    bms <- readRaw(file.name, offset = bms[length(bms)]+4, nbytes = 4, human = "char", size = 1, endian = "little")
		    bms <- blockString(bms)
        		## Wavenumbers for MIR spectra from Tensor are assigned prefix "m", MIR spectra from Alpha prefixed "a"; for NIR MPA "n"
        		pref <- ifelse(instr.range=="ten-mir", "m", ifelse(instr.range=="alp-mir", "a", ifelse(instr.range=="mpa-nir", "n", "X")))
  		if(instr.range=="ten-off"){ instr.range="ten-mir"} ## AS: Old ten-mir written as tensor-27
			  if(speclib=="ICRAF"){ 
          		if(instr.range=="ten-mir"){
					  wb <- rev(seq(get("icraf.htsxt", spec.opts)[3], get("icraf.htsxt", spec.opts)[2], (get("icraf.htsxt", spec.opts)[2]-get("icraf.htsxt", spec.opts)[3])/(get("icraf.htsxt", spec.opts)[1]-1)))
				  										}
				  if(instr.range=="alp-mir"){
				  wb <- rev(seq(get("icraf.alpha", spec.opts)[3], get("icraf.alpha", spec.opts)[2], (get("icraf.alpha", spec.opts)[2]-get("icraf.alpha", spec.opts)[3])/(get("icraf.alpha", spec.opts)[1]-1)))
          												}
			    if(instr.range=="mpa-nir"){
				wb <- rev(seq(get("icraf.mpa", spec.opts)[3], get("icraf.mpa", spec.opts)[2], (get("icraf.mpa", spec.opts)[2]-get("icraf.mpa", spec.opts)[3])/(get("icraf.mpa", spec.opts)[1]-1)))
 				  										}
    				if(bms=="ZnSe"){
			  	pref="a"
       		 	wb <- rev(seq(499.8151, 3996.4810, (3996.4810-499.8151)/(1715-1)))
	 		    							}
	 	    								}
	 	     if(!(instr.range=="ten-mir"|instr.range=="alp-mir"|instr.range=="mpa-nir"|instr.range=="ten-off"|bms=="ZnSe")){ stop("Unknown file format. See '?read.opus' for more info.") }  

		 	  ## speclib="New"
			   if(speclib=="New"){
      	  			if(instr.range=="ten-mir"){ pref="m" }
			 		if(instr.range=="alp-mir"){ pref="a"}
					if(instr.range=="mpa-nir"){ pref="n"}
					if(bms=="ZnSe"){ pref="a"}
		    								}
	      ## Get positions where the following parameters are found in the file  
     		z <- grepRaw(codes[1],pr,all=TRUE)[1]+5
  			re <- grepRaw(codes[2],pr,all=TRUE)[1]+5
  	 		snm <- grepRaw(codes[3],pr,all=TRUE)[1]+7
  	 		dat <- grepRaw(codes[4],pr,all=TRUE)[1]+7
  	 		lwn <- grepRaw(codes[5],pr,all=TRUE)[1]+7
  	 		fx <- grepRaw(codes[6],pr,all=TRUE)[3]+7
     		lx <- grepRaw(codes[7],pr,all=TRUE)[3]+7
     		npt0 <- grepRaw(codes[8],pr,all=TRUE)[2]+3
     		npt1 <- grepRaw(codes[8],pr,all=TRUE)[3]+7
     		mxy <- grepRaw(codes[9],pr,all=TRUE)[1]+7 
     		mny <- grepRaw(codes[10],pr,all=TRUE)[3]+7 
     		end <- grepRaw(codes[11],pr,all=TRUE)+11
     		tim <- grepRaw(codes[12],pr,all=TRUE)+11
        		## calculate end and start of each block:
 			offs <- sapply(5:10, function(x){end[x]})
     		byts <- diff(offs)
  	 		ZFF <- readRaw(file.name, offset=z, nbytes=4, human="int", size=2)[[5]][1]
  	 		RES <- readRaw(file.name, offset=re, nbytes=4, human="int", size=2)[[5]][1]
  			snm.lab.material <- blockString(readRaw(file.name, offset = snm, nbytes = 22, human = "char", size = 1, endian = "little"))
  			if(!nzchar(snm.lab.material)){
  		 	  SSN <- ""
  		  	  Material <- ""
  		  	  warning("Product name not found inside OPUS file...")
  		       } else {
            if(!length(grep(snm.lab.material, pattern=";"))==0){
            	snm.lab.material <- as.vector(strsplit(snm.lab.material,";"))[[1]]
           		SSN <- paste0(snm.lab.material[2], snm.lab.material[1])
           		Material <- snm.lab.material[3]
   			  				} 
   			  		else {
       		 				 if(!length(grep(snm.lab.material, pattern="_"))==0){ ## AS: ICR_02182
          	  				SSN <- sub("_", "", snm.lab.material)
           					Material <- ""          
    								}
    						else {
    							if(!length(snm.lab.material)==0){ 
         		 	  		SSN <- snm.lab.material
         		 			Material <- ""          
    				  					}
       							}
    		  					}   
   						}
        ## Set three SSN first three characters to lower
     	SSN <- paste0(tolower(substr(SSN,1,3)), substr(SSN,4,20))
  	 	Scandate <- blockString(readRaw(file.name, offset = dat, nbytes = 10, human = "char", size = 1, endian = "little"))
  	 	Scantime <- blockString(readRaw(file.name, offset = tim[2]-4, nbytes = 8, human = "char", size = 1, endian = "little"))
  	 	Scandate <- paste(Scandate,Scantime)
  		LWN <- readRaw(file.name, offset=lwn, nbytes=8, human="real", size=8)[[5]][1]
  	 	## Combine the above parameters
		spectrum.meta <- c(SSN, Material, Scandate, ZFF, RES, LWN)
		## Get number of data points for each spectra data block   		
  	 	NPT0 <- readRaw(file.name, offset=npt0, nbytes=12, human="int", size=4)[[5]][2]
  	 	NPT1 <- readRaw(file.name, offset=npt1, nbytes=4, human="int", size=4)[[5]][1]
  	 	fxv <- readRaw(file.name, offset=fx, nbytes=16, human="real", size=8)[[5]][1] ## fxv:	Frequency of first point
  	 	lxv <- readRaw(file.name, offset=lx, nbytes=16, human="real", size=8)[[5]][1] ## lxv:	Frequency of last point
  	 	Wavenumbers <- rev(seq(lxv, fxv, (fxv-lxv)/(NPT1-1)))
  	 	## Read all through all the data blocks inside the OPUS file:
  	 	nbytes1 <- NPT0*4 ## initial parameters
  	 	smxa <- c()
  	 	smna <- c()
  	 	nbytes.f <- NPT1*4
              if(offs[1]<2000){
       		offs.f<-offs[3]
       		}
       		
       		 if(offs[1]>20000){
       		offs.f<-offs[2]
       		}
  	   ## Selected spectra block
    	   opus.p <- readRaw(file.name,width=NULL,offset=offs.f-4,nbytes=nbytes.f,human="real",size=4,endian="little") 
      spectra <- opus.p[[5]]
       
  
    		## Make compatible to ICRAF spectra:
		  if(speclib=="ICRAF"){
		## TH: Is spline fitting necessary? AS: Yes, this standardizes all spectral data points to conform to ICRAF spectral library
        spectra <- spline(Wavenumbers, spectra, xout=wb, method="natural")$y
		Wavenumbers <- wb
      }
      ## Specify if graphics showing spectra being converted is displayed
     	if(plot.spectra==TRUE){
        plot(Wavenumbers, spectra, ylab="Absorabance", xlab=expression("Wavenumbers cm"^-1), type="l")   
        mtext(paste("File source: ", getwd(),file.name,sep="/"), side=3,line=2,cex=1.2)
 	    }
  	  ## Print progress of conversion
      if(print.progress==TRUE){
        message(paste("Converting ", file.name, " file", sep=""))
      }
	    ## Add meta ID
  	  if(missing(MID)){ MID <- paste(Sys.getenv(c("USERDNSDOMAIN"))[[1]], instr.range, sep="_") } 
  	  ## create data.frames:
      samples <- data.frame(SAMPLEID=spectrum.meta[1], Material=spectrum.meta[2], Zero.Filing=spectrum.meta[4], Resolution=spectrum.meta[5], LWN=spectrum.meta[6], DateTime=as.POSIXct(spectrum.meta[3], format="%d/%m/%Y %H:%M:%S "), MID=MID)
      ab <- data.frame(as.list(signif(spectra, signif.digit)))
      ## standardize the BAND 'names' so they correspond to the names in the reference models:
      rn <- rank(Wavenumbers, ties.method="first")
      if(missing(st.wavenumbers)){
        st.wavenumbers <- wavenumbers[wavenumbers$TYPE==instr.range,]
     	}
      w.s <- cut(Wavenumbers[rn], breaks=c(st.wavenumbers$LOWER, st.wavenumbers$UPPER[nrow(st.wavenumbers)]), labels=st.wavenumbers$BAND, include.lowest=TRUE)
      w.s <- as.numeric(sapply(w.s, function(x){gsub("[^0-9.]", "", x)}))
      names(ab)[rn] <- paste0(pref, ifelse(is.na(w.s), round(Wavenumbers[rn], 1), w.s))
      ab <- cbind(data.frame(SAMPLEID=spectrum.meta[1]), ab)
      message('Creating object of type "SpectraPoints"...')
      out <- new("Spectra", samples=samples, ab=ab)
      if(is.null(sp)){ 
        sp <- SpatialPoints(data.frame(lat=0,lon=0), proj4string=CRS(as.character(NA)))
      }
      out <- SpectraPoints(Spectra=out, sp=sp)
      out@metadata$MID <- MID 
      out@metadata$Instrument_name <- instr.range
      out@metadata$Material_class <- spectrum.meta[2]
      out@metadata$Laboratory_URL <- Sys.getenv(c("USERDNSDOMAIN"))[[1]] 
      return(out)
    }
  } else {
	    warning(paste("File",file.name,"does not exist"))
	}
}
	
setMethod("read.opus", signature(file.name = "character"), .read.opus)

## read in a loop:
setMethod("read.opus", signature(file.name = "list"), function(file.name, ...){
     xx <- lapply(file.name, .read.opus, ...)
     Spectra <- new("Spectra", samples=do.call(rbind, lapply(xx, function(x){slot(slot(x, "data"), "samples")})), ab=do.call(rbind, lapply(xx, function(x){slot(slot(x, "data"), "ab")})))
     metadata <- do.call(rbind, lapply(xx, function(x){slot(x, "metadata")}))
     metadata <- metadata[!duplicated(metadata$MID),]
     sp <- do.call(rbind.SpatialPoints, lapply(xx, function(x){slot(x, "sp")}))
     out <- new("SpectraPoints",data=Spectra, sp=sp, metadata=metadata)
})

## End of script;







