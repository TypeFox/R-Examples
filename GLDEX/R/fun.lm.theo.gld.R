fun.lm.theo.gld <-
function (L1, L2, L3, L4, param) 
{
    if (length(L1) > 1) {
        L4 <- L1[4]
        L3 <- L1[3]
        L2 <- L1[2]
        L1 <- L1[1]
    }
    result <- rep(NA, 4)
    if (tolower(param) == "rs") {
        result[1] <- L1 + 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)
        result[2] <- fun.Lm.gt.2.rs(L3, L4, 2)/L2
        result[3] <- fun.Lm.gt.2.rs(L3, L4, 3)/L2
        result[4] <- fun.Lm.gt.2.rs(L3, L4, 4)/L2
    }
    if (tolower(param) == "fmkl" | tolower(param)=="fkml") {
        if (L3 != 0 & L4 != 0) {
            result[1] <- L1 - 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)
        }
        if (L3 == 0 & L4 == 0) {
            result[1] <- fun.fmkl0(1)/L2 + L1
        }
        if (L3 != 0 & L4 == 0) {
            result[1] <- fun.fmkl.L40(1, L3)/L2 + L1
        }
        if (L3 == 0 & L4 != 0) {
            result[1] <- fun.fmkl.L30(1, L4)/L2 + L1
        }
        result[2] <- fun.Lm.gt.2.fmkl(L3, L4, 2)/L2
        result[3] <- fun.Lm.gt.2.fmkl(L3, L4, 3)/L2
        result[4] <- fun.Lm.gt.2.fmkl(L3, L4, 4)/L2
    }
    return(result)
}
