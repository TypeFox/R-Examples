weigh.SD <-
function(w, data, vari, na.rm)
{
  dacc <- paste0("data$", vari)  
  
  DA  <- eval(parse(text=dacc))
  nas <- is.na(DA)
  
  # NIMMT NA's raus wenn welche in den PVs drin sind!!
  if(any(nas) & na.rm)
    {#NA's rausnehmen
      wna <- which(nas)  
      w1 <- w[-wna]
      
      wm  <- sum(w1*DA[-wna]/sum(w1))
      #erg <- 1/(sum(w1) - mean(w1)) * sum(w1 * (DA[-wna] - wm)^2) 
      erg <- sum(w1)/(sum(w1)^2 - mean(w1^2)) * sum(w1 * (DA[-wna] - wm)^2)
                 
      return(sqrt(erg))
    } else {
            #wm  <- sum(w*DA[-wna]/sum(w))
            wm  <- sum(w*DA/sum(w))
            #erg <- 1/(sum(w) - mean(w)) * sum(w * (DA - wm)^2)
            #erg <- sum(w1)/(sum(w1)^2 - mean(w1^2)) * sum(w1 * (DA[-wna] - wm)^2)
            erg <- sum(w)/(sum(w)^2 - mean(w^2)) * sum(w * (DA - wm)^2)
                       
            return(sqrt(erg))
           }

}
