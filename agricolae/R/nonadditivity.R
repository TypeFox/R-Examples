`nonadditivity` <-
function (y,factor1, factor2, df, MSerror)
{
name.y <- paste(deparse(substitute(y)))
conjunto <- subset(data.frame(y, factor1,factor2), is.na(y) == FALSE)
zz<-tapply(conjunto[,1],conjunto[,c(2,3)],mean)
nn<-tapply(conjunto[,1],conjunto[,c(2,3)],length)
r<-1/mean(1/nn) # Repeticiones
sc<-df*MSerror
x1 <- rownames(zz)
y1 <- colnames(zz)
fila <- length(x1)
col <- length(y1)
total <- fila * col
x <- character(length = total)
y <- character(length = total)
z <- numeric(length = total)
k <- 0
for (i in 1:fila) {
    for (j in 1:col) {
        k <- k + 1
        x[k] <- x1[i]
        y[k] <- y1[j]
        z[k] <- zz[i, j]
    }
}
yy <- data.frame(x, y, z)
yy[,1] <- as.factor(yy[,1])
yy[,2] <- as.factor(yy[,2])
modelo.y<-lm(yy[,3] ~ yy[,1]+ yy[,2])
y.est<-predict(modelo.y)
residual  <-residuals(modelo.y)
q<-y.est^2
modelo.q <- lm(q ~ yy[,1]+ yy[,2])
Nonadditivity  <-residuals(modelo.q)
modelo.e <- lm( residual ~ Nonadditivity )
anva <- anova(modelo.e)
beta <-coef(modelo.e)[2]
P <-  r*anva[1,2] / beta
Q <- P/beta
# 1 grado de libertad para la no-aditividad
Sna <-as.numeric(P^2/Q)
if( r == 1) {
sc<-sc-Sna
df<-df-1
}
# Construccion del cuadro de no-additivity
anva<-anova(modelo.e)
anva["Residuals",1]<-df
anva["Residuals",2]<-sc
anva["Residuals",3]<-sc/df
anva["Nonadditivity",2:3]<-r*anva["Nonadditivity",2:3]
anva["Nonadditivity",4]<-anva["Nonadditivity",3]/anva["Residuals",3]
anva["Nonadditivity",5]<- 1- pf( anva["Nonadditivity",4],1, df)
cat("\nTukey's test of nonadditivity\n")
cat(name.y,"\n")
cat("\nP :", P)
cat("\nQ :", Q,"\n\n")
print(anva)
output <- list(P=P, Q=Q, anva=anva)
return(output)
}

