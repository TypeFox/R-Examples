descomponer <- function (y,frequency,type) {
  # Author: Francisco Parra Rodriguez
  # http://rpubs.com/PacoParra/24432
  # date:"y", frequency:"frequency". 
  # Use 7 for frequency when the data are sampled daily, and the natural time period is a week, 
  # or 4 and 12 when the data are sampled quarterly and monthly and the natural time period is a year.
  n <- length(y)
  y <- matrix(y,ncol=1)
  f1 <- NULL
  if(n%%2==0) {f2 <- n/(2*frequency)} else {
    f2 <- (n-1)/(2*frequency)}
  #Modelo para obtener serie con tendencia
  c <- seq(from=2, to=(2+(n/frequency) ))
  #Use the "sort.data.frame" function, Kevin Wright. Package taRifx
  i <- seq(1:n)
  i2 <- i*i  
  if (type==1)
  {eq <- lm(y~i)  
   z <- eq$fitted} else {
     if (type==2) eq <- lm(y~i+i2) 
     z <- eq$fitted} 
  cx <- cdf(z)
  id <- seq(1,n)
  S1 <- data.frame(cx)
  S2 <- S1[1:(2+(n/frequency)),]
  X <- as.matrix(S2)
  cy <- gdf(y)
  B <- solve(X%*%t(X))%*%(X%*%cy)
  Y <- t(X)%*%B
  BTD <- B
  XTD <- t(MW(n))%*%t(X)
  TD <- gdt(Y)
  # Genero la serie residual
  IRST <- y-TD
  # Realizo la regresion dependiente de la frecuenca utilizando como explicativa IRST.
  # modelo para obtener serie con  estacionalidad con trunc.
  frecuencia <- seq(0:(n/2)) 
  frecuencia <- frecuencia-1
  S <- data.frame(f1=frecuencia)
  sel <- subset(S,f1==trunc(2*f2))
  c <- seq(from=2,to=(n/f2))
  for (i in c) {sel1 <- subset(S,f1==i*trunc(2*f2))
                sel <- rbind(sel,sel1)}
  m1 <- c(sel$f1 * 2)
  m2 <- c(m1+1)
  c <- c(m1,m2)
  n3 <- length(c)
  d <- rep(1,n3)
  s <- data.frame(c,d)
  S <- sort.data.frame (s,formula=~c)
  #Use the "sort.data.frame" function, Kevin Wright. Package taRifx
  l <- frequency*trunc(n/frequency)
  i <- seq(1:l)
  i2 <- i*i  
  if (type==1)
  {eq <- lm(y[1:l]~i)  
   z <- eq$fitted} else {
     if (type==2) eq <- lm(y[1:l]~i+i2) 
     z <- eq$fitted} 
  cx <- cdf(z)
  id <- seq(1,l)
  S1 <- data.frame(cx,c=id)
  S2 <- merge(S,S1,by.x="c",by.y="c")
  S3 <- rbind(c(1,1,cx[1,]),S2) 
  m <- l+2
  X1 <- S3[,3:m]
  # matriz de regresores a l
  X1 <- as.matrix(X1)
  # la paso al dominio del tiempo
  X2 <- data.frame(t(MW(l))%*%t(X1))
  if (n==l) X3 <- X2 else
    X3 <- rbind(X2,X2[1:(n-l),])
  # la paso al dominio de la frecuencia
  X4 <-MW(n)%*%as.matrix(X3)
  cy <- gdf(IRST)
  B1 <- solve(t(X4)%*%X4)%*%(t(X4)%*%cy)
  Y <- X4%*%B1
  BST <- B1
  XST <- t(MW(n))%*%X4
  ST <- gdt(Y)
  TDST <- TD+ST
  IR <- IRST-ST  
  data <- data.frame(y,TDST,TD,ST,IR)
  regresoresTD <- data.frame(XTD)
  regresoresST <- data.frame(XST)
  list(datos=data,regresoresTD=regresoresTD,regresoresST=regresoresST,coeficientesTD=BTD,coeficientesST=BST)}
