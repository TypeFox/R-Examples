ntaste <- data.frame(score = tastetest$score,
                     scr   = as.numeric(tastetest$scr) - 1,
                     liq   = as.numeric(tastetest$liq) - 1,
                     scrliq = ( as.numeric(tastetest$scr) -1 ) * 
                              ( as.numeric(tastetest$liq) -1 )
                     );  ntaste

Omega <- lm(score~scr*liq,data=tastetest)
M <- model.matrix(Omega)
M2 <- cbind(M[,3], M[,2] - 2 * M[,4])
M3 <- cbind(M[,2], M[,3] - 2 * M[,4])

omega1 <- lm(score~scr+liq,data=tastetest)
omega2 <- lm(score~M2,tastetest)
omega2a <- lm(score~ liq + I(scr - 2 * scrliq),data=ntaste)
omega3 <- lm(score~M3,tastetest)
omega3a <- lm(score~ scr + I(liq - 2 * scrliq),data=ntaste)

anova(omega1,Omega)   # test for interaction
# test main effect for scr
# anova(omega2a,Omega)  # this gives the same result as line below
anova(omega2,Omega)   
# test main effect for liq
# anova(omega3a,Omega)  # this gives the same result as line below
anova(omega3,Omega)   
