Sim_Data <-
function(Obs, VG=1, VP1=0.05, VP2=0.1, VP3=0.3){

  Obs_orig<- Obs

  if (Obs/4 != round(Obs/4)) (Obs<- ceiling(Obs/4)*4)

  # Scaling variable X1: X1 \in [0,100]
    X1_ja  <- rnorm(n=1*Obs/4, mean=80, sd=5)
    MIN1   <- min(X1_ja)
    MAX1   <- MIN1*0.9
    X1_nein<- runif(n=3*Obs/4, min=0, max=MAX1)
    X1     <- c(X1_ja, X1_nein)

  # Scaling variable X2: X2 \in [0,100]
    X2_ja1 <- rnorm(n=1*Obs/4, mean=80, sd=5)
    X2_ja2 <- rnorm(n=1*Obs/4, mean=80, sd=5)
    MIN2   <- min(X2_ja1, X2_ja2)
    MAX2   <- MIN2*0.9
    X2_nein<- runif(n=2*Obs/4, min=0, max=MAX2)
    X2     <- c(X2_ja1, X2_ja2, X2_nein)

  # Scaling variable X3: X3 \in [0,100]
    X3_ja   <- rnorm(n=1*Obs/4, mean=80, sd=5)
    MIN3    <- min(X3_ja)
    MAX3    <- MIN3*0.9
    X3_nein1<- runif(n=2*Obs/4, min=0, max=MAX3)
    X3_nein2<- runif(n=1*Obs/4, min=0, max=MAX3)
    X3     <- c(X3_nein1, X3_ja, X3_nein2)


 ####################
 # Class assignment #
 ####################
    Klasse<- c()
    # Class 1
      T1<- X1 >= MIN1
      K1<- which(T1==1)
      Klasse[K1]<- 1

    # Class 2
      T1<- X1 <  MIN1
      T2<- X2 >= MIN2
      K2<- which(T1*T2==1)
      Klasse[K2]<- 2

    # Class 3
      T1<- X1 <  MIN1
      T3<- X3 >= MIN3
      K3<- which(T1*T3==1)
      Klasse[K3]<- 3

    # Class 4
      T1<- X1 < MIN1
      T2<- X2 < MIN2
      T3<- X3 < MIN3
      K4<- which(T1*T2*T3==1)
      Klasse[K4]<- 4

  # Class assignment
    Y <- rbinom(n=Obs, size=1, prob=VG)*Klasse

  # Observations without assigned class
    m1<- which(Y==0)
    m2<- sample(c(rep(1,Obs/4), rep(2,Obs/4), rep(3,Obs/4), rep(4,Obs/4)))
    for (i in 1:length(m1)) Y[m1[i]]<- m2[i]
    Y2<- c()
    Y2[which(Y==1)]<- "A"
    Y2[which(Y==2)]<- "B"
    Y2[which(Y==3)]<- "C"
    Y2[which(Y==4)]<- "D"

  # Class assignments and addition of noise to the variables X3, X4, and X5
    X4<- X3; X5<- X3; X6<- X3
    for(i in 1:Obs){ 
        if (Y[i]==3){
            b1   <- rbinom(n=1, size=1, prob=VP1)
            if (b1==1) X4[i]<- runif(n=1, min = 0, max = MAX3)  # X3 <-> X4

            b2   <- rbinom(n=1, size=1, prob=VP2)
            if (b2==1) X5[i]<- runif(n=1, min = 0, max = MAX3)  # X3 <-> X5

            b3   <- rbinom(n=1, size=1, prob=VP3)
            if (b3==1) X6[i]<- runif(n=1, min = 0, max = MAX3)  # X3 <-> X6
        }
    }

  # More variables
    X7 <- rnorm(n=Obs, mean = 5, sd = 1)
    X8 <- rnorm(n=Obs, mean = 5, sd = 1)
    X9 <- rnorm(n=Obs, mean = 5, sd = 1)
    X10<- rnorm(n=Obs, mean = 5, sd = 1)
#    X11<- rnorm(n=Obs, mean = 5, sd = 1)
#    X12<- rnorm(n=Obs, mean = 5, sd = 1)
#    X13<- rnorm(n=Obs, mean = 5, sd = 1)
#    X14<- rnorm(n=Obs, mean = 5, sd = 1)
#    X15<- rnorm(n=Obs, mean = 5, sd = 1)
#    X16<- rnorm(n=Obs, mean = 5, sd = 1)
#    X17<- rnorm(n=Obs, mean = 5, sd = 1)
#    X18<- rnorm(n=Obs, mean = 5, sd = 1)
#    X19<- rnorm(n=Obs, mean = 5, sd = 1)
#    X20<- rnorm(n=Obs, mean = 5, sd = 1)
    X  <- as.data.frame(cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10))
#                             X11,X12,X13,X14,X15,X16,X17,X18,X19,X20))

   DATEN<- cbind(Y2, X)
   
   if (nrow(DATEN) > Obs_orig) DATEN<- DATEN[-sample(1:Obs, 2),]

  # Output
    return(DATEN)
}
