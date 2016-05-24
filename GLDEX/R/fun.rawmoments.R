
"fun.rawmoments"<-
function (L1, L2, L3, L4, param = "fmkl") 
{
    if (param == "rs") {
        v1 <- fun.rsb(L3, L4, 1)
        v2 <- fun.rsb(L3, L4, 2)
        v3 <- fun.rsb(L3, L4, 3)
        v4 <- fun.rsb(L3, L4, 4)
        m1 <- L1 + 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)
        central.m2 <- (v2 - v1^2)/L2^2
        central.m3 <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3))
        central.m4 <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 3 * 
            v1^4)/((L2^4))
    }
    if (param == "fmkl") {

        if(L3==0 & L4==0){
       
        v1 <- 0
        v2 <- fun.fmkl0(2)
        v3 <- 0
        v4 <- fun.fmkl0(4)
        m1 <- v1/L2+L1}


        if(L3==0 & L4!=0){
        
        v1 <- fun.fmkl.L30(1,L4)
        v2 <- fun.fmkl.L30(2,L4)
        v3 <- fun.fmkl.L30(3,L4)
        v4 <- fun.fmkl.L30(4,L4)
        m1 <- v1/L2+L1}


        if(L3!=0 & L4==0){
        
        v1 <- fun.fmkl.L40(1,L3)
        v2 <- fun.fmkl.L40(2,L3)
        v3 <- fun.fmkl.L40(3,L3)
        v4 <- fun.fmkl.L40(4,L3)
        m1 <- v1/L2+L1}


        if(L3!=0 & L4!=0){
        
        v1 <- fun.fmklb(L3, L4, 1)
        v2 <- fun.fmklb(L3, L4, 2)
        v3 <- fun.fmklb(L3, L4, 3)
        v4 <- fun.fmklb(L3, L4, 4)
        m1 <- L1 - 1/L2 * ((L3 + 1)^-1 - (L4 + 1)^-1)}

        central.m2 <- (v2 - v1^2)/L2^2
        central.m3 <- (v3 - 3 * v1 * v2 + 2 * (v1)^3)/((L2^3))
        central.m4 <- (v4 - 4 * v1 * v3 + 6 * (v1^2) * v2 - 3 * 
            v1^4)/((L2^4))
    }
    ex1 <- m1
    ex2 <- central.m2 + m1^2
    ex3 <- central.m3 + 3 * central.m2 * m1 + (m1)^3
    ex4 <- central.m4 + 4 * central.m3 * m1 + 6 * central.m2 * 
        m1^2 + m1^4
    return(cbind(ex1, ex2, ex3, ex4))
}
