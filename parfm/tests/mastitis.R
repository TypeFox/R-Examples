library(parfm)
data(mastitis)
head(mastitis)

mastitis$timeto <- as.numeric((mastitis$Midpoint * 4 / 365.25))

set.seed(1)
mastitis <- mastitis[sample(1:nrow(mastitis), 200), ]


################################################################################
#Example 4.4: The gamma frailty model for the udder quarter infection data     #
#Duchateau and Janssen (2008, page 136)                                        #
################################################################################
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster="Cowid",
                  dist="weibull",
                  frailty="gamma",
                  data=mastitis)
modParfm


################################################################################
#Example 4.7 The inverse Gaussian frailty model for the udder quarter          #
#infection data                                                                #
#Duchateau and Janssen (2008, page 156)                                        #
################################################################################
mastitis <- data.frame(mastitis, 
                       timeto=as.numeric((mastitis$Midpoint *4 / 365.25)))
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster="Cowid",
                  dist="weibull",
                  frailty="ingau",
                  data=mastitis)
modParfm


################################################################################
#Example 4.10 The positive stable frailty model for the udder quarter          #
#infection data                                                                #
#Duchateau and Janssen (2008, page 169)                                        #
################################################################################
mastitis <- data.frame(mastitis, 
                       timeto=as.numeric((mastitis$Midpoint *4 / 365.25)))
modParfm <- parfm(Surv(timeto, Status) ~ Heifer,
                  cluster="Cowid",
                  dist="weibull",
                  frailty="possta",
                  data=mastitis)
modParfm
