setClass("theta", representation(kinpar = "vector", specpar = "list", 
        irfpar = "vector", parmu = "list", partau = "vector", dummy="vector",
        clpequ = "vector", specdisppar = "list", kinscal = "vector", 
	prel = "vector", oscpar="vector", coh="vector", drel = "vector", cohirf="vector", 
	kinpar2 = "vector", kin2scal="vector", jvec = "vector", 
	anipar = "vector", peakpar = "list", amplitudes = "vector", 
	shift = "vector", amps = "list", kinscalspecial = "list",
	eigenvaluesK ="vector"), 
	prototype = list(kinpar = vector(), specpar = list(), dummy=vector(), 
	irfpar = vector(), parmu = list(), partau = vector(), 
	clpequ = vector(), specdisppar = list(), drel = vector(), oscpar=vector(),
	coh = vector(), kinscal = vector(), prel = vector(),
	jvec = vector(), cohirf= vector(), kinpar2 = vector(), 
	kin2scal = vector(), anipar = vector(), peakpar = list(),
	amplitudes=vector(), shift = vector(), amps = list(),
          kinscalspecial = list(), eigenvaluesK = vector() ))



