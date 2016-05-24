"lmomtri" <-
function(para, paracheck=TRUE, nmom=c("3", "5")) {
    nmom <- match.arg(nmom)
    nmom <- as.numeric(nmom)

    if(paracheck) {
       if(! are.partri.valid(para)) return()
    }
    attributes(para$para) <- NULL

    MIN  <- para$para[1]
    MODE <- para$para[2]
    MAX  <- para$para[3]
    A <- (MAX-MIN); B <- (MODE-MIN); C <- (MAX-MODE)

    L1 <- (MIN + MODE + MAX)/3
    L2 <- ((MIN-MODE)^2/(MAX-MIN) - (MIN+MODE) + 2*MAX)/15

    tmp1 <- 2*MIN*(B/A)^3 + (12/7)*A^(-3)*B^4 + 2*MAX-2*MAX*(B/A)^3 - (12/7)*A^(-3)*C^4
    tmp2 <- (12/5)*2*A^(-2)*C^3 - (12/3)*A^(-1)*C^2
    tmp3 <- -3*((1/A^2) * (MIN*B^2 + (4/5)*B^3 ) + MAX - (1/A^2) * (MAX*B^2 - (4/5)*C^3 + (4/3)*A*C^2))
    L3   <- tmp1 + tmp2 + tmp3 + L1
    LCV  <- L2/L1
    TAU3 <- L3/L2

    E11  <- L1
    E22  <- L2 + E11
    E33  <- (L3 + 3*E22 - E11)/2
    if(nmom == 3) {
       z <- list(lambdas=c(L1, L2, L3, NA, NA),
                 ratios=c(NA, LCV, TAU3, NA, NA),
                 trim=0, rightrim=0, leftrim=0, E33err=NA,
                 source="lmomtri")
       return(z)
    } else {
       E33p <- E44 <- E55 <- NULL
       try(E33p <- expect.max.ostat(3, para=para, qua=quatri), silent=TRUE)
       try(E44  <- expect.max.ostat(4, para=para, qua=quatri), silent=TRUE)
       try(E55  <- expect.max.ostat(5, para=para, qua=quatri), silent=TRUE)
       if(is.null(E33p)) {
          warning("E33p could not be computed by expect.max.ostat() numerical integration ",
                  "so E33err will be NA")
          E33p <- NA
       }
       if(is.null(E44)) {
          warning("E44 could not be computed by expect.max.ostat() numerical integration ",
                  "so Tau4 will be NA")
          E44 <- NA
       }
       if(is.null(E55)) {
          warning("E55 could not be computed by expect.max.ostat() numerical integration ",
                  "so Tau5 will be NA")
          E55 <- NA
       }
       L4 <-           5*E44 - 10*E33 +  6*E22 - 1*E11
       L5 <- 14*E55 - 35*E44 + 30*E33 - 10*E22 + 1*E11
       TAU4 <- L4/L2; TAU5 <- L5/L2

       E33percent.error <- 100*((E33 - E33p)/E33)

       z <- list(lambdas=c(L1, L2, L3, L4, L5),
                 ratios=c(NA, LCV, TAU3, TAU4, TAU5),
                 trim=0, rightrim=0, leftrim=0, E33err=E33percent.error,
                 source="lmomtri")
       return(z)
    }
}

