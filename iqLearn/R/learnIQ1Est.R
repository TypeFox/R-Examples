learnIQ1Est <- function (mainObj, cmObj, sigObj, dens){
    A1 = cmObj$A1;
    
    if (dens=="norm"){
        normal = T;
    }
    else if (dens=="nonpar"){
        normal = F;
    }
    else{
        stop ("dens must be one of {norm, nonpar}");
    }
    ## NEXT: estimate optimal treatment for each patient in training
    ## set
    
    ## contrast function mean estimate
    muPos = cmObj$cmPos;
    muNeg = cmObj$cmNeg;
    
    ## lhat is the estimate of the expectation of the main effect term
    lhatPos = mainObj$mainPos; 
    lhatNeg = mainObj$mainNeg;
    n = length (lhatPos);
    ## standard deviation estimates
    if (sigObj$homo){
        sigPos = rep (sigObj$sigPos, n);
        sigNeg = rep (sigObj$sigNeg, n);
    } else {
        sigPos = sigObj$sigPos;
        sigNeg = sigObj$sigNeg;
    }
    ## Estimate Q1 for both A1=1 and A1=-1 using either normal density
    ## or empirical estimate
    if (normal){
        q1HatPos = lhatPos + muPos*(1-2*pnorm (-muPos/sigPos)) +
            sqrt(2/pi)*sigPos*exp (-muPos^2/(2*sigPos^2)); 
        q1HatNeg = lhatNeg + muNeg*(1-2*pnorm (-muNeg/sigNeg)) +
            sqrt(2/pi)*sigNeg*exp (-muNeg^2/(2*sigNeg^2)); 
    }
    else{
        n = length (lhatPos)
        q1HatPos = rep (NA, n);
        q1HatNeg = rep (NA, n);
        for (i in 1:n){
            opPos = sum (abs (muPos[i] +
                sigPos[i]*sigObj$stdResids));  
            opNeg = sum (abs (muNeg[i] +
                sigNeg[i]*sigObj$stdResids));  
            q1HatPos[i] = lhatPos[i] + opPos;
            q1HatNeg[i] = lhatNeg[i] + opNeg;
        }
    }
    
    ## Vector of optimal first-stage txts
    optA1 = sign (q1HatPos - q1HatNeg);
    
    list ("optA1"=optA1)
}
