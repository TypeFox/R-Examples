# Automatically generated from all.nw using noweb
alignped2 <- function(x, dad, mom, level, horder, packed,
                      spouselist) {
    x <- x[order(horder[x])]  # Use the hints to order the sibs
    rval <- alignped1(x[1],  dad, mom, level, horder, packed, 
                      spouselist)
    spouselist <- rval$spouselist

    if (length(x) >1) {
        mylev <- level[x[1]]
        for (i in 2:length(x)) {
            rval2 <-  alignped1(x[i], dad, mom, level,
                                horder, packed, spouselist)
            spouselist <- rval2$spouselist
            
            # Deal with the unusual special case:
            if ((rval2$n[mylev] > 1) || 
                          (is.na(match(x[i], floor(rval$nid[mylev,])))))
                rval <- alignped3(rval, rval2, packed)
            }
        rval$spouselist <- spouselist
        }
    rval
    }
