modelFunction <- function(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, retFct, 
                          doseScaling, respScaling, isFinite, pshifts = NULL)
{
    if (!is.null(retFct))
    {
        drcFct <- retFct(doseScaling, respScaling)
    }
    drcFct1 <- function(dose, parm)
    {
        parmVal <- parm2mat(parm)
#        print(c(dim(pshifts), dim(parmVal)))
        if ((!is.null(pshifts)) & all(dim(pshifts) == dim(parmVal))) 
        {
            parmVal <- parmVal + pshifts
        }     
#        drcFct(dose, (parm2mat(parm))[isFinite, , drop = FALSE])
        drcFct(dose, parmVal[isFinite, , drop = FALSE])
    }

    if (is.null(cm))
    {
        multCurves <- function(dose, parm)
        {
           drcFct1(dose, parm)
        }
    } else {  # not adapting to scaling (not using drcFct1)!!!
        iv <- isFinite & (assayNoOld == cm)
        niv <- !iv
        fctEval <- rep(0, length(dose))

        multCurves <- function(dose, parm)
        {
            parmVal <- (parm2mat(parm))[isFinite, , drop = FALSE]
#            print(c(dim(pweights), dim(parmVal)))
            if ((!is.null(pshifts)) & all(dim(pshifts) == dim(parmVal))) 
            {
                parmVal <- parmVal + pshifts
            }
            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , drop = FALSE])

            fctEval
        }
    }
    
    multCurves
}

