setMethod("checkIC", signature(IC = "IC", L2Fam = "L2RegTypeFamily"), 
    function(IC, L2Fam){ 
        if(dimension(Domain(IC@Curve[[1]])) != dimension(img(L2Fam@RegDistr)) + dimension(img(L2Fam@ErrorDistr)))
            stop("dimension of 'Domain' of 'Curve' !=\n", 
                 "dimension of 'img' of 'RegDistr' + dimension of 'img' of 'ErrorDistr'")

        trafo <- L2Fam@param@trafo
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        cent <- E(L2Fam, IC1)
        cat("precision of centering:\t", cent, "\n")

        dims <- length(L2Fam@param)
        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")
        consist <- E(L2Fam, IC1 %*% t(L2deriv))
        cat("precision of Fisher consistency:\n")
        print(consist - trafo)
        
        return(invisible())
    })
