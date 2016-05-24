# A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra (fwhm).
# read LC/GC-MS spectra(.netCDF, .mzXML, mzData)
# use functions in XCMS package
read_ms_spec<-function(folder.name, metadata = NULL, profmethod='bin', fwhm=30, bw=30){
  files <- list.files(folder.name, recursive=T, full.names=TRUE);
  if (is.null(metadata)) {
	xset <- xcms::xcmsSet(files, profmethod = "bin", fwhm=fwhm)
  } else {
	xset = xcms::xcmsSet(files, profmethod = "bin", fwhm = fwhm, sclass = metadata[,1])
  }
  xset<-xcms::group(xset, bw=bw);
  xset
}



# retention time correction for LC/GC-MS spectra
ms_rt_correction<-function(xset, bw=30){
  xset2<-xcms::retcor(xset)
  # re-group peaks after retention time correction
  xset2<-xcms::group(xset2, bw=30)
  xset2
}



# fill in missing peaks
ms_fill_peaks<-function(xset){
  xset2<-xcms::fillPeaks(xset);
  xset2
}


# into:  integrated area of original (raw) peak
# intf:  integrated area of filtered peak
# maxo:  maximum intensity of original (raw) peak
# maxf:  maximum intensity of filtered peak

ms_create_matrix<-function(xset, intvalue = "into"){
  values <- xcms::groupval(xset, "medret", value = intvalue);
  values
}

# A commonly used value is 30 (seconds) for LC-MS and 4 (seconds) for GC-MS spectra (fwhm).
read_ms_spectra = function(folder.name, type = "undefined", filename.meta= NULL, description = "", prof.method='bin', fwhm=30, bw=30, intvalue = "into", header.col.meta = TRUE, header.row.meta = TRUE, sep.meta = ","){
	if (!is.null(filename.meta))
		metadata = read_metadata(filename.meta, header.col = header.col.meta, header.row = header.row.meta, sep = sep.meta)
	else metadata = NULL
	xset = read_ms_spec(folder.name, metadata = metadata, prof.method, fwhm, bw)
	xset = ms_rt_correction(xset, bw)
	xset = ms_fill_peaks(xset)
	mat = ms_create_matrix(xset, intvalue)
	dataset = create_dataset(mat, type = type, metadata = metadata, description = description, 
							 label.x = "mz/rt", label.values = "intensity", xSet = xset)
	dataset
  
}


