"baseIRF" <-
function (irfvec, indexlow, indexhigh, removeNeg = FALSE)
{
        irfvec <- irfvec - mean(irfvec[indexlow:indexhigh])
        if(removeNeg)
                      irfvec[which(irfvec < 0 )] <- 0
        irfvec
}
