t_mal <-
function(u, g, ts, tq, d)# ts:age at starting smoking, tq:age at quitting smoking
 {
 if (!is.na(ts) & !is.na(tq) & ts>tq) stop("Error: ts>tq")
 if ( is.na(ts) & !is.na(tq) ) stop("Error: missing ts")
 if(is.na(ts)) {x <- "never smoker"} else if (is.na(tq)) {x <- "current smoker"} else {x <- "former smoker"} 
 ### Here we determine the smoking status x. ###
 ### CAUTION   : the three possible statuses are:
 ### "Never"   : ts==NA and tq==NA
 ### "Current" : ts     and tq==NA
 ### "Former"  : ts     and tq
 Tmal<-99999
 #-----------------------#
 #-----------------------#
 if (x=="never smoker")
   {
    st <- function(t,a){exp(-ht_mal_int(0, t, g, d, "N"))-a}
   }   else if     
 #-----------------------#
   (x=="current smoker") {

      i1  <- ht_mal_int(0, ts, g, d, "N")

      st1 <- function(t,a) {exp( - ht_mal_int(0, t, g, d, "N"))-a}
      st2 <- function(t,a) {exp( - i1 - ht_mal_int(ts, t, g, d, "Y"))-a}

   if (u<exp(-i1)) st<-st2 else st<-st1
   }    else   {
 #-----------------------#
 #  (x=="former smoker")

      i1  <- ht_mal_int(0, ts, g, d, "N")
      i2  <- ht_mal_int(ts, tq, g, d, "Y")

      st1 <- function(t,a) {exp( - ht_mal_int(0, t, g, d, "N"))-a}
      st2 <- function(t,a) {exp( - i1 - ht_mal_int(ts, t, g, d, "Y"))-a}
      st3 <- function(t,a) {exp( - i1 - i2 - ht_mal_int(tq, t, g, d, "N"))-a}

   if (u<exp(-i1-i2)) st<-st3 else if (u<exp(-i1)) st<-st2 else st<-st1
   }
 #-----------------------#
 #-----------------------#

 st_try  <- try(uniroot (st, c(0,150), u), T)
 if(class(st_try)!="try-error")   Tmal<-st_try$root   else   Tmal=150

return(Tmal)
 }
