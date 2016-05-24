#' Runs the Wang-Wei-Zhu decomposition
#' 
#' This function runs the Wang-Wei-Zhu decomposition.
#' 
#' @param x an object of the class decompr
#' @return the decomposed table
#' @author Bastiaan Quast
#' @details Adapted from code by Fei Wang.
#' @references Wang, Zhi, Shang-Jin Wei, and Kunfu Zhu.
#' Quantifying international production sharing at the bilateral and sector levels. 
#' No. w19677. National Bureau of Economic Research, 2013.
#' @export
#' @examples
#' # load example data
#' data(leather)
#' 
#' # create intermediate object (class decompr)
#' decompr_object <- load_tables_vectors(inter,
#'                                       final,
#'                                       countries,
#'                                       industries,
#'                                       out        )
#' 
#' # run the WWZ decomposition on the decompr object
#' wwz(decompr_object)


wwz <- function( x ) {
  
  # Part 1: Decomposing Export into VA( 16 items )
  # defining ALL to contain all decomposed results
  ALL  <- array( 0, dim=c( x$GN, x$G, 19 ) )
  
  # the order of 16 items and exp, expint, expfd
  decomp19 <- c(  "DVA_FIN",
                  "DVA_INT",
                  "DVA_INTrexI1",
                  "DVA_INTrexF",
                  "DVA_INTrexI2",
                  "RDV_INT",
                  "RDV_FIN",
                  "RDV_FIN2",
                  "OVA_FIN",
                  "MVA_FIN",
                  "OVA_INT",
                  "MVA_INT",
                  "DDC_FIN",
                  "DDC_INT",
                  "ODC",
                  "MDC",
                  "texp",
                  "texpint",
                  "texpfd"
  )
  
  ALL[ ,,17 ] <- x$ESR
  ALL[ ,,18 ] <- x$Eint
  ALL[ ,,19 ] <- x$Efd
  x$Eint      <- NULL
  x$Efd       <- NULL
  
  
  ####### THIS ONE TAKES A LONG TIME ##########
  # Part 2-1 == H10-(1): DVA_FIN=[ VsBss#Ysr ]
  # OK
  for ( r in 1: x$G ) {
    z1           <- matrix( x$Ym[ , r], nrow = x$GN, ncol = x$GN )
    ALL[ ,r,1 ]  <- colSums( x$Vhat %*% x$Bd * t(z1) )
  }
  # View( ALL[ ,,1 ] )   # DVA_FIN
  
  
  # Part 2-2: H10-(2): DVA_INT
  # OK
  VsLss <- x$Vhat %*% x$L
  z1    <- x$Am %*% x$Bd %*% x$Yd
  for ( r in 1: x$G ) {
    # r=2
    ALL[ ,r,2 ] <- colSums(VsLss)*t(z1[ ,r ])
  }
  # View( ALL[ ,,2 ] )   #  DVA_INT
  
  
  # Part 2-3: H10-(3): DVA_INTrexI1
  # OK
  z1 <- matrix(rowSums( x$Yd), nrow = x$GN, ncol = x$GN )
  for ( tt in 1:x$G ) {
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z1[ m:n,m:n ] <- 0
  }
  
  z2 <- x$Bm %*% z1
  for ( tt in 1:x$G ) {
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z2[ m:n,m:n ] <- 0
  }
  
  z3 <- x$Am * t( z2 )
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,3 ] <- colSums(VsLss) * (rowSums(z3[ ,m:n ]) )
  }
  # View( ALL[ ,,3 ] )  #  DVA_INTrexI1
  
  
  # Part 2-4: H10-(4): DVA_INTrexF
  # OK
  z <- matrix( 0, nrow = x$GN, ncol = x$GN )
  z1 <- rowSums( x$Ym )
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z[,m:n]      <- z1 - x$Ym[ ,tt ]
    z[ m:n,m:n ] <- 0
  }
  
  z2 <- x$Am * t(x$Bd %*% z)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,4 ] <- colSums(VsLss) * (rowSums(z2[ ,m:n ]) )
  }

  
  # Part 2-5: H10-(5): DVA_INTrexI2
  # OK !
  z1 <- t( x$Bm %*% z )
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z1[ m:n,m:n ] <- 0
  }
  
  z2 <- x$Am * z1
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,5 ] <- colSums(VsLss) * (rowSums(z2[ ,m:n ]) )
  }
  #  View( ALL[ ,,5 ] )  #  DVA_INTrexI2
  
  
  # Part 2-6: H10-(6): RDV_FIN
  # OK !
  z <- matrix( 0, nrow = x$GN, ncol = x$GN )
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z[,m:n] <-  x$Ym[ ,tt ]
  }
  
  z1 <- x$Am * t( x$Bd %*% z )

  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,7 ] <- colSums(VsLss) * (rowSums(z1[ ,m:n ]) )
  }
  # View( ALL[ ,,7 ] )  #  RDV_FIN
  
  
  # Part 2-7 == H10-(7): RDV_FIN2
  # OK !
  z1 <- x$Bm %*% z
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z1[ m:n,m:n ] <- 0
  }
  
  z2 <- x$Am * t(z1)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,8 ] <- colSums(VsLss)*(rowSums(z2[ ,m:n ]) )
  }
  # rm(z2)
  # View( ALL[ ,,8 ] )  #  RDV_FIN2

  
  # Part 2-8 == H10-(8): RDV_INT
  # OK !
  z <- matrix( 0, nrow = x$GN, ncol = x$GN )
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z[,m:n] <-  x$Yd[ ,tt ]
  }
  
  z1 <- x$Am * t( x$Bm %*% z )
  for ( r in 1: x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,6 ] <- colSums(VsLss)*(rowSums(z1[ ,m:n ]) )
  }
  # View( ALL[ ,,6 ] )  #  RDV_INT
  
  
  # Part 2-9 == H10-(9): DDC_FIN
  # OK !
  z <- matrix( 0, nrow = x$GN, ncol = x$GN)
  for ( tt in 1:x$G ){
    m <- 1   + (tt-1) * x$N
    n <- x$N + (tt-1) * x$N
    z[m:n,m:n] <-  rowSums( x$Ym[ m:n, ] )
  }
  
  z1 <- x$Am * t(x$Bm %*% z)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,13 ] <- colSums(VsLss)*(rowSums(z1[ ,m:n ]))
  }
  # View( ALL[ ,,13 ] )  #   DDC_FIN
  

  ####### THIS ONE TAKES A LONG TIME ##########
  # Part 2-10 == H10-(10): DDC_INT
  # OK !
  z <- x$Am %*% diag(x$X)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    ALL[ ,r,14 ] <- colSums(diag(x$Vc) %*% x$Bd - VsLss) * (rowSums( z[ ,m:n ]) )
  }
  # rm(VsLss)
  # View( ALL[ ,,14 ] )  # DDC_INT
  
  # Part 2-11 == H10-(11): MVA_FIN  =[ VrBrs#Ysr ]
  #    H10-(14): OVA_FIN  =[ Sum(VtBts)#rYsr ]
  # OK !
  VrBrs <- x$Vhat %*% x$Bm
  YYsr  <- matrix(0, nrow = x$GN, ncol = x$GN)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    YYsr[ ,1:x$GN ] <- x$Ym[ ,r ]
    z <- VrBrs * t( YYsr )
    ALL[ ,r,10 ] <- colSums( z[ m:n, ]       ) # MVA_FIN[ ,r ]
    ALL[ ,r,9 ]  <- colSums( z[ -c( m:n ), ] ) # OVA_FIN[ ,r ]
  }
  # rm( YYsr )
  # View( ALL[ ,,9 ]  ) # OVA_FIN
  # View( ALL[ ,,10 ] ) # MVA_FIN
  
  ####### THIS ONE TAKES A VERY LONG TIME ##########
  # Part 2-12 == H10-(12): MVA_INT  =[ VrBrs#AsrLrrYrr ]
  #     H10-(15): OVA_INT  =[ Sum(VtBts)#AsrLrrYrr ]
  # OK !
  YYrr <- matrix(0, nrow = x$GN, ncol = x$GN)
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    YYrr[ , 1:x$GN ] <- x$Yd[ , r ]
    z <- VrBrs * t(x$Am %*% x$L %*% YYrr)
    ALL[ ,r,12 ] <- colSums( z[ m:n, ]        )  #  MVA_INT[ ,r ]
    ALL[ ,r,11 ] <- colSums( z[ -c( m:n ) , ] )  #   OVA_INT[ ,r ]
  }
  # rm( YYrr )
  # View( ALL[ ,,11 ] )  #  OVA_INT
  # View( ALL[ ,,12 ] )  #  MVA_INT
  
  ####### THIS ONE TAKES A VERY LONG TIME ##########
  # Part 2-13 == H10-(13): MDC  =[ VrBrs#AsrLrrEr* ]
  # ==  H10-(16): ODC  =[ Sum(VtBts)#AsrLrrEr* ]
  # OK !
  for ( r in 1:x$G ) {
    m <- 1   + (r-1) * x$N
    n <- x$N + (r-1) * x$N
    EEr <- matrix(0, nrow = x$GN, ncol = x$GN)
    EEr[ m:n, 1:x$GN ] <- x$E[ m:n, 1 ]
    z <- VrBrs * t(x$Am %*% x$L %*% EEr)
    ALL[ ,r,16 ] <- colSums( z[ m:n, ]        )   # MDC[ ,r ]
    ALL[ ,r,15 ] <- colSums( z[ -c( m:n ) , ] )   # ODC[ ,r ]
  }
  # rm( EEr, z, VrBrs )
  
  
  ###### try to calculate DViX_Fsr
  DViX_Fsr <- VsLss %*% x$ESR
  DViX_Fsr <- t(DViX_Fsr)
  dim(DViX_Fsr) <- c(x$GN*x$G, 1)
  
  
  dimnames( ALL )  <-  list( x$rownam, x$k, decomp19)  
  
  
  # Part 3 Putting all results in one sheet
  
  for ( u in 1:x$GN ){
    if ( u==1 ){
      #ALLandTotal <- rbind( ALL[ u,, ], colSums( ALL[ u,, ] ))
      ALLandTotal <- rbind( ALL[ u,, ])
    }
    else {
      #ALLandTotal <- rbind( ALLandTotal, ALL[ u,, ], colSums( ALL[ u,, ] ))
      ALLandTotal <- rbind( ALLandTotal, ALL[ u,, ])
    }
  }
  # rm(ALL)
  rownames( ALLandTotal ) <- NULL #x$bigrownam
  
  
  
  # Part 4  checking the differences resulted in texp, texpfd, texpintdiff
  
  texpdiff        <- rowSums( ALLandTotal[ ,1:16 ]) - ALLandTotal[ ,17 ]
  
  texpfddiff      <- ALLandTotal[ ,1] + ALLandTotal[ ,9 ] + ALLandTotal[ ,10 ] - ALLandTotal[ ,19 ]

  texpintdiff     <- rowSums( ALLandTotal[ ,2:8 ] ) + rowSums( ALLandTotal[ ,11:16 ] ) - ALLandTotal[ ,18 ]
  
  texpdiffpercent <-  texpdiff / ALLandTotal[ ,17 ] * 100
  texpdiffpercent[ is.na(texpdiffpercent) ] <- 0
  texpfddiffpercent <-  texpfddiff/ALLandTotal[ ,19 ] * 100
  texpfddiffpercent[ is.na(texpfddiffpercent) ] <- 0
  texpintdiffpercent <-  texpintdiff/ALLandTotal[ ,18 ] * 100
  texpintdiffpercent[ is.na(texpintdiffpercent) ] <- 0
  texpdiff <- round( texpdiff, 4)
  texpfddiff <- round( texpfddiff,4)
  texpintdiff <-  round( texpintdiff, 4)
  texpdiffpercent <- round( texpdiffpercent, 4)
  texpfddiffpercent <- round( texpfddiffpercent,4)
  texpintdiffpercent <-  round( texpintdiffpercent, 4)
  
  ALLandTotal <-   data.frame( rep(x$k,                       each=length(x$k)*length(x$i) ),
                               rep(x$i,  times = length(x$k), each=length(x$k) ),
                               rep(x$k,  times = length(x$k)*length(x$i) ),
                               ALLandTotal,
                               texpdiff,
                               texpfddiff,
                               texpintdiff,
                               texpdiffpercent,
                               texpfddiffpercent,
                               texpintdiffpercent,
                               DViX_Fsr)
  
  names(ALLandTotal)[1:3] <- c("Exporting_Country", "Exporting_Industry", "Importing_Country")
                          
  # dim( ALLandTotal )
  
  attr(ALLandTotal, "decomposition") <- "wwz"
  
  return(ALLandTotal)
  
}