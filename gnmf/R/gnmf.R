# These usage instructions are from the stand-alone binary.
# They are kept here for informational purposes only,
# but the R implementation should be similar.
#
# Usage: hpcNMF -i input_data_matrix
#
# General suggestion: For an input matrix containing no zero entries, it's recommended that idealization be set to 0 (default),
#  scaling to 'T' and normalizing to 'T'. For an input matrix containing zero entries, idealization should be set to 0.1,
#   scaling to 'F' and normalizing to 'F'. If idealization is set to 0.1, setting both scaling and normalizing to 'T' 
#   is recommended for faster convergence
#
# Other options:
#   -s: v, 'KL', 'Renyi', 'ED', 'GammaJD', 'GammaKL', 'DivComb', 'Div2', 'Pareto', 'InverseLink_GammaKL', 'BD', 'NBD', 'ODP'
#   -n: update steps, default 2000
#   -c: repeats, default 20
#   -k: rank range, default 2-2
#   -t: clustering target, default 'PATTERN (H matrix)'
#   -cs: clustering scheme, default 'Binary', could be 'PearsonHC'
#   -r: reference file, string
#   -scaling: initial scaling for faster convergence, default 'F'
#   -normalizing: H matrix normalization, default 'F'
#   -alphas: csv string, valid only for scheme Renyi and Pareto, default 1.0
#   -RunType: default 'whole', could be 'simulation' or 'evaluation'
#   -cts: set convergence test step size, default 20
#   -idealization, default 0, can be set to 0.1

gnmf <-
function(V, scheme, nsteps=2000, repeats=20, ranks=2,
         cltarget="PATTERN", clscheme="Binary", reffile="", scaling="F",
         normalizing="F", alphas=1, runtype="simulation", cstepsize=20, idealization=1) {

    ###################################################################
    # DATA CHECKS.
    ###################################################################
    # Make sure V is numeric.
    if ( ! is.numeric(V) ) {
        print("hpcnmf(): V must be numeric!")
        return(NaN);
    }

    # Make sure V is a matrix.
    if ( ! is.matrix(V) ) {
        print("hpcnmf(): V must be a matrix!")
        return(NaN);
    }

    # Check 'scheme'.
    if ( ( scheme != "KL"    )
       & ( scheme != "Renyi" )
       & ( scheme != "ED"    ) ) {
        print("hpcnmf(): scheme must be KL, Renyi, or ED!")
        return(NaN);
    }

    # Check 'nsteps'.
    if ( nsteps < 1 ) {
        print("hpcnmf(): nsteps must be at least 1!")
        return(NaN);
    }

    # Check 'repeats'.
    if ( repeats < 1 ) {
        print("hpcnmf(): repeats must be at least 1!")
        return(NaN);
    }

    # Check 'ranks'.
    matrixRank <- qr(V)$rank
    if ( max(ranks) > matrixRank ) {
        print("hpcnmf(): values in ranks cannot be greater than the rank of V!")
        return(NaN);
    }
    if ( min(ranks) < 2 ) {
        print("hpcnmf(): values in ranks must be at least 2!")
        return(NaN);
    }
    numRanks = length(ranks)
    if ( numRanks != 1 ) {
        print("gnmf(): ranks must be a scalar!")
        return(NaN);
    }

    # Check 'cltarget'.
    if ( ( cltarget != "PATTERN"   )
       & ( cltarget != "ALTERNATE" ) ) {
        print("hpcnmf(): cltarget must be PATTERN or ALTERNATE!")
        return(NaN);
    }

    # Check 'clscheme'.
    if ( ( clscheme != "Binary"  )
       & ( clscheme != "PearsonHC" ) ) {
        print("hpcnmf(): clscheme must be Binary or PearsonHC!")
        return(NaN);
    }

    # Check 'reffile'.
    if ( ! is.character(reffile) ) {
        print("hpcnmf(): reffile must be the name of a file!")
        return(NaN);
    }

    # Check 'scaling'.
    if ( ( scaling != "F"  )
       & ( scaling != "T" ) ) {
        print("hpcnmf(): scaling must be either T or F!")
        return(NaN);
    }

    # Check 'normalizing'.
    if ( ( normalizing != "F"  )
       & ( normalizing != "T" ) ) {
        print("hpcnmf(): normalizing must be either T or F!")
        return(NaN);
    }

    # Check 'alphas'.
    # If SCHEME is NOT "Renyi", then 'alphas' does not apply.
    # In that case, force alpha to scalar 1.0.
    if ( scheme != "Renyi" ) {
        alphas = 1
    } else {
        if ( ! is.numeric(alphas) ) {
            print("hpcnmf(): alphas must be numeric!")
            return(NaN);
        }
    }
    numAlphas = length(alphas)
    if ( numAlphas != 1 ) {
        print("gnmf(): alphas must be a scalar!")
        return(NaN);
    }

    # Check 'runtype'.
    if ( ( runtype != "whole"      )
       & ( runtype != "simulation" )
       & ( runtype != "evaluation" ) ) {
        print("hpcnmf(): runtype must be whole, simulation, or evaluation!")
        return(NaN);
    }

    # Check 'cstepsize'.
    if ( cstepsize < 1 ) {
        print("hpcnmf(): cstepsize must be at least 1!")
        return(NaN);
    }

    # Make sure idealization is numeric.
    if ( ! is.numeric(idealization) ) {
        print("hpcnmf(): idealization must be numeric!")
        return(NaN);
    }

    ###################################################################
    # ALLOCATE MEMORY FOR RESULTS, THEN INVOKE C++ CODE.
    ###################################################################

    # Create two matrices of zeros which will hold the NMF result.
    # This step allocates memory that the C++ code will use.
    # V ~ W * H
    nrowV          <- nrow(V)
    ncolV          <- ncol(V)
    sizeAlphaBlock <- repeats * sum(min(ranks):max(ranks)) * length(alphas)
    H              <- matrix(rep(0,ncolV*sizeAlphaBlock),nrow=sizeAlphaBlock)
    W              <- matrix(rep(0,nrowV*sizeAlphaBlock),ncol=sizeAlphaBlock)

    # Invoke C function 'MainOrig'.
    # The first argument of .C() is the name of the C wrapper function.
    # The rest of the arguments provide pointers to various R variables,
    # including pointers to matrix data.
    output =.C("CWrapper",
        matrixV      = as.double(V),
        nRows        = as.integer(nrowV),
        nCols        = as.integer(ncolV),
        scheme       = as.character(scheme),
        nsteps       = as.integer(nsteps),
        repeats      = as.integer(repeats),
        rankrange    = as.integer(ranks),
        numRankRange = as.integer(numRanks),
        cltarget     = as.character(cltarget),
        clscheme     = as.character(clscheme),
        reffile      = as.character(reffile),
        scaling      = as.character(scaling),
        normalizing  = as.character(normalizing),
        alphas       = as.double(alphas),
        nalphas      = as.integer(numAlphas),
        runtype      = as.character(runtype),
        cstepsize    = as.integer(cstepsize),
        idealization = as.double(idealization),
        matrixH      = as.double(H),
        matrixW      = as.double(W)
    )

    ###################################################################
    # REFORMAT RESULTS AND RETURN.
    ###################################################################

    # Reshape 'H' and 'W' back into matrices.
    # Then return the matrices in a named list.
    Htmp             <- matrix(output$matrixH,nrow=sizeAlphaBlock)	
    Wtmp             <- matrix(output$matrixW,ncol=sizeAlphaBlock)
    H                <- gnmfLevel1Parse(Htmp,1,nrowV,ncolV,alphas,ranks,repeats)
    W                <- gnmfLevel1Parse(Wtmp,2,nrowV,ncolV,alphas,ranks,repeats)
    returnVal        <- list(H,W)
    names(returnVal) <- c("H","W")
    return(returnVal)
}
