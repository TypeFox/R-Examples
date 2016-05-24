chao1 <-
function (x) {
    Sobs <- sum(x>0)
    n <- sum(x)
    F1 <- as.numeric(table(x)[names(table(x))=="1"])
    F2 <- as.numeric(table(x)[names(table(x))=="2"])
    F1 <- ifelse(length(F1)==0, 0, F1)
    F2 <- ifelse(length(F2)==0, 0, F2)

    if (F1>0 && F2>0) {

        Chao1 <- Sobs + (((n-1)/n)*(F1^2/(2*F2)))
        var <- F2*((((n-1)/n)*((F1/F2)^2)) + ((((n-1)/n)^2)*((F1/F2)^3)) + (0.25*(((n-1)/n)^2)*((F1/F2)^4)))
        T <- Chao1 - Sobs # This equation equals 0 when F1 = 1 and F2 = 1
        K <- exp(1.96*(log(1+(var/T^2)))^0.5) # This equation outputs Na when F1 = 1 and F2 = 1
        low95 <- Sobs + (T/K)
        upp95 <- Sobs + (T*K)

    } else if (F1>1 && F2==0) {

        Chao1 <- Sobs + (((n-1)/n)*(F1*(F1-1)/(2*(F2+1))))
        var <- (((n-1)/n)*(F1*(F1-1)/2) + ((((n-1)/n)^2)*((F1*(((2*F1)-1)^2))/4)) - (((n-1)/n)^2)*(F1^4/(4*Chao1)))
        T <- Chao1 - Sobs
        K <- exp(1.96*(log(1+(var/T^2)))^0.5)
        low95 <- Sobs + (T/K)
        upp95 <- Sobs + (T*K)

    } else {

        Chao1 <- Sobs
        z <- table(x)
        z <- z[!names(z)=="0"]
        F <- as.numeric(z)
        i <- as.numeric(names(z))
        var <- (sum(F*(exp(-i)-exp(-2*i))) - ((1/n)*sum(i*exp(-i)*F)^2))
        P <- sum(F*exp(-i))/Sobs
        low95 <- max(Sobs, ((Sobs/(1-P)) - ((1.96*(var)^0.5)/(1-P))))
        upp95 <- (Sobs/(1-P)) + ((1.96*(var)^0.5)/(1-P))

    }
    res <- data.frame(S.obs = Sobs, S.chao1 = Chao1, var, low95, upp95)
}
