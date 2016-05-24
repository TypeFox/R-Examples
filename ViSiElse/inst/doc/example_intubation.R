## ---- fig.show='hold',fig.width=7----------------------------------------
vars <- c("time_in_intub","time_insert_probe","time_out_intub","delay_intub_prob") 
label <- c( "Blade in the mouth", "Tube insert in the mouth",
           "Blade out of the mouth","Blade in - Tube inserted")
typeA <- c( "p","p","p","l")
showorder <- c( 1,3,4,2)
deb   <- c( NA, NA,NA,"time_in_intub")
fin <- c( NA,NA,NA,"time_insert_probe")
bookdataframe <- data.frame(vars,label,typeA,showorder,deb,fin)
library(ViSiElse)
bookintubation <- ConvertoViSibook(bookdataframe)
plot(bookintubation)

## ---- fig.show='hold',fig.width=7,fig.height=7---------------------------
x <- buildViSiGrid( X = intubation, 
                              book = bookintubation ,
                              pixel = 2 )
plot(x, main = "Intubation, t_0 : Blade in the mouth")

## ---- fig.show='hold',fig.width=7, fig.height=7--------------------------
#bookintubation <- changeShoworder( bookintubation , c(2,3,4))
#Vintubation <- buildViSiGrid( X = intubation, 
#                              book = bookintubation ,
#                              t_0 = "time_in_intub",pixel = 1 )
#plot(Vintubation, main = "Intubation, t_0 : Blade in the mouth")

