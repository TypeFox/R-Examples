### --- Test setup ---
#
# regression tests
#

if(FALSE) {
    ## Not really needed, but can be handy when writing tests
    library(RUnit)
    library(GenABEL)
}

### do not run
#stop("SKIP THIS TEST")
###

### ---- common functions and data -----

#source(paste("../inst/unitTests/shared_functions.R"))
#source(paste(path,"/shared_functions.R",sep=""))

### --- Test functions ---

test.descriptives.trait <- function()
{
  	require(GenABEL.data)
    data(ge03d2ex)
# this works
    descriptives.trait(ge03d2ex,by.var=phdata(ge03d2ex)$sex)
# bug [#1184]
# and this does not!
    checkException(descriptives.trait(ge03d2ex,by.var=sex))
    attach(phdata(ge03d2ex))
    descriptives.trait(ge03d2ex,by.var=sex)
    detach(phdata(ge03d2ex))
    descriptives.trait(ge03d2ex,by.var="sex")
    phdata(ge03d2ex)$sex[2] <- NA
# bug [#1185]
# and this does not!
    descriptives.trait(ge03d2ex,by.var=phdata(ge03d2ex)$sex)
    descriptives.trait(ge03d2ex,by.var="sex")
# bug [#1259]
    convert.snp.illumina(infile="test_markers", outfile="test.raw")
    test = load.gwaa.data(pheno="test_phenos", geno="test.raw", force=T)
    attach(phdata(test))
    descriptives.trait(data=test, by = bt)
    detach(phdata(test))
}
