setClass("dat", representation(psi.df = "matrix", psi.weight = "matrix", 
        x = "vector", nt = "integer", x2 = "vector", nl = "integer", 
        C2 = "matrix", E2 = "matrix", sigma = "numeric", 
	mod_type = "character", parnames = "vector",dummy="vector",
        simdata = "logical", weightpar = "list", 
        weight = "logical", weightM = "matrix", weightsmooth = "list", 
        fixed = "list", clp0 = "list", makeps = "character", 
        clpequspec = "list", lclp0 = "logical", lclpequ = "logical", 
        title = "character", mhist = "list", datCall = "list", 
        dscalspec = "list", mvecind = "vector", 
        drel = "vector", clpequ = "vector", scalx = "numeric", 
        prel = "vector", prelspec = "list", fvecind = "vector", 
        pvecind = "vector", nvecind = "vector",
	iter = "numeric", free = "list",
        clpCon = "list", ncomp = "numeric", clpdep = "logical", 
	inten = "matrix", positivepar="vector", constrained ="list", 
	clinde = "list", chinde = "list", highcon = "vector", 
	lowcon = "vector", datafile = "character", getX = "logical",
	clpType = "character", clpequspecBD = "list", compnames = "vector", 
	getXsuper = "logical", usecompnames0 = "logical", 
	usecompnamesequ = "logical", autoclp0 = "list", cohcol = "numeric",
        weightList = "list", outMat = "matrix", satMat = "matrix"), 
        prototype = list(psi.df = matrix(), psi.weight = matrix(), 
            x = vector(), nt = integer(), x2 = vector(), nl = integer(), 
            C2 = matrix(), E2 = matrix(), sigma = numeric(), 
            mod_type = character(), simdata = logical(), 
            weightpar = list(), weight = FALSE, weightM = matrix(), 
            weightsmooth = list(), fixed = list(), clp0 = list(), 
            clpequspec = list(), clpCon = list(), lclp0 = logical(), 
            lclpequ = logical(), makeps = character(), title = character(), 
            mhist = list(), datCall = list(), nvecind = vector(), 
            dscal = list(), drel = vector(),  mvecind = vector(),
            scalx = 1, prel = vector(), prelspec = list(), fvecind = vector(), 
            pvecind = vector(), clpequ = vector(), free = list(),
            iter = 1, ncomp = numeric(), clpdep = logical(),inten = matrix(), 
	    parnames=vector(), dummy=vector(), positivepar=vector(), constrained =list(), 
	    clinde = list(), chinde = list(), highcon = vector(), 
	    lowcon = vector(), datafile = "", getX = FALSE, 
	    clpType = "", clpequspecBD = list(), compnames = vector(),
	    getXsuper = FALSE, usecompnames0 = FALSE,  
	    usecompnamesequ = FALSE, autoclp0 = list(), cohcol = 0,
            weightList = list(), outMat = matrix(), satMat = matrix() ))

setClass("kin", representation("dat", kinpar = "vector", specpar =
"list", seqmod = "logical", irf = "logical", mirf = "logical", reftau
= "numeric", measured_irf = "vector", convalg =
"numeric", irffun = "character", irfpar = "vector", cohirf = "vector",
dispmu = "logical", dispmufun = "character", anipar = "vector", parmu
= "list", disptau = "logical", disptaufun = "character", partau =
"vector", fullk = "logical", kmat = "array", jvec = "vector", 
anispec = "list", ncolc = "vector", 
kinscal = "vector", kmatfit = "array", cohspec = "list", coh = "vector",
oscspec = "list", oscpar = "vector",
wavedep = "logical", lambdac = "numeric", speckin2 = "list", 
usekin2 = "logical", kinpar2 = "vector",
kin2scal = "vector", amplitudes = "vector", streakT = "numeric", 
streak="logical", doublegaus = "logical", multiplegaus = "logical", fixedkmat="logical",
kinscalspecialspec ="list", kinscalspecial = "list",
lightregimespec = "list",
numericalintegration = "logical", initialvals = "vector",
reactantstoichiometrymatrix = "vector",
stoichiometrymatrix = "vector"),
         prototype = list( kinpar = vector(), seqmod =
TRUE, irf = FALSE, mirf = FALSE, measured_irf = vector(), convalg = 1,
cohirf = vector(), irffun = "gaus", anispec = list(), irfpar =
vector(), dispmu = FALSE, dispmufun = "poly", parmu = list(), anipar =
vector(), disptaufun = "poly", reftau = 0, specpar = list(), 
partau = vector(), posk = FALSE, disptau = FALSE, fullk = FALSE,
kmat = array(), jvec = vector(), 
ncolc = vector(), kinscal = vector(), kmatfit = array(), 
cohspec = list(), coh = vector(), oscspec = list(), oscpar = vector(), wavedep = logical(), 
lambdac = numeric(), speckin2 = list(), usekin2 = FALSE,
kinpar2 = vector(), kin2scal = vector(), amplitudes = vector(), 
streakT = 0, streak = FALSE, doublegaus = FALSE, multiplegaus = FALSE, fixedkmat=FALSE,
           kinscalspecialspec = list(),  kinscalspecial = list(),
           lightregimespec = list(), numericalintegration = FALSE,
           initialvals = vector(),
           reactantstoichiometrymatrix = vector(), 
           stoichiometrymatrix = vector()
           ))

setClass("mass", 
representation("kin", 
peakpar = "list", 
amplitudes = "vector", 
peakfunct = "character",
lzerofile = "character", 
extracomp = "logical", 
shift = "vector"),
prototype = list( 
peakpar = list(), 
peakfunct = "expmodgaus", 
lzerofile = "", 
amplitudes = vector(), 
getX=TRUE, 
extracomp = TRUE,
shift = vector() )
)

setClass("spec", representation("dat", clpequ = "vector", 
        specpar = "list", specfun = "character", specref = "numeric", 
        specCon = "list", ncole = "vector", specdisp = "logical", 
        specdisppar = "list", specdispindex = "vector", nupow = "numeric", 
        timedep = "logical", parmufunc = "character"), 
        prototype = list(specpar = list(), specfun = "gaus", 
          specCon = list(), ncole = vector(), specdisp = FALSE, 
          specdisppar = list(), specdispindex = vector(), specref =
          numeric(), nupow = 5, clpequ = vector(), 
          timedep = FALSE, parmufunc = "poly")) 

setClass("amp", representation("dat", amps = "list"),
         prototype = list(amps = list(), clpdep = FALSE))


