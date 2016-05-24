## ----,eval=FALSE---------------------------------------------------------
#     #exists in basic packages: you can use their names directly or by a new name
#     qt
#     qnorm
#     b <- qpois
#     #or from other packages  : Here you must give it a new name
#     #you can replace the package::functionname on your needs.
#     f <- stats::qweibull
#     invcdfnames <- c("qt","qnorm","qpois","f")
#    #or invcdfnames <- c("qt","qnorm","b","f") but never
#    #invcdfnames <- c("qt","qnorm","qpois","stats::qweibull")

## ----,eval=FALSE---------------------------------------------------------
#    #always you can use the following way, the inner lists' names should match the
#    #above functions' arguments names.
#    paramslists <- list(
#               m1 = list(df = 5 ),
#               m2 = list(mean = 0, sd = 1),
#               m3 = list(lambda = 3),
#               m4 = list(shape = 1, scale = 1)
#               )
#    #if you are lazy,e.g. qnorm using the default values, then you can use the following way:
#   paramslists2 <- list(
#               m1 = list(df = 5 ),
#               m3 = list(lambda = 3),
#               m4 = list(shape = 1, scale = 1)
#               )
#   defaultindex <- c(2)

## ----,eval=FALSE---------------------------------------------------------
#   #If you are familiar with the bounding RA algorithm, you can set the functions' arguments
#   #on your needs. e.g. let m1 = 80, sigma0 = 0.001 will be ok if you know the smaller
#   #sigma0  the more time will be costed. But if you don't familiar with it, you'd better
#   #use the default values

## ------------------------------------------------------------------------
cor_matrix <- matrix(c(1.0,-0.4,0.1,-0.2,-0.4,
                       1.0,0.8,0.6,0.1,0.8,1.0,
                       0.5, -0.2,0.6,0.5,1.0
                       ),4,4)

## ----,eval=TRUE----------------------------------------------------------
  f <- stats::qweibull
  invcdfnames <- c("qt","qnorm","qpois","f")
  paramslists <- list(
             m1 = list(df = 5 ),
             m2 = list(mean = 0, sd = 1),
             m3 = list(lambda = 3),
             m4 = list(shape = 1, scale = 1)                 
             )
  cor_matrix <- matrix(c(1.0,-0.4,0.1,-0.2,-0.4,
                       1.0,0.8,0.6,0.1,0.8,1.0,
                       0.5, -0.2,0.6,0.5,1.0
                       ),4,4) 
  cor_matrix
  res <- NORTARA::genNORTARA(10000,cor_matrix,invcdfnames,paramslists)
  head(res,5)
  cor(res)
  paramslists2 <- list(
             m1 = list(df = 5 ),             
             m3 = list(lambda = 3),
             m4 = list(shape = 1, scale = 1)                 
             )
  defaultindex <- c(2)
  res2 <- NORTARA::genNORTARA(10000,cor_matrix,invcdfnames,paramslists2,defaultindex)
  head(res2,5)
  cor(res2)

