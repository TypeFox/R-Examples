.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius
arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}

TGammaStar <- function(Tleaf, Patm,
                       Egamma=37830.0, 
                       value25=42.75){  

  value25*arrh(Tleaf,Egamma)*Patm/100
}

TKm <- function(Tleaf, Patm,
                Oi = 210,      # O2 concentration (mmol mol-1)
                Ec = 79430.0,  # activation energy for Kc 
                Eo = 36380.0,  # activation energy for Ko
                Kc25 = 404.9,  # Kc at 25C
                Ko25 = 278.4  # Ko at 25C
                ){
  
  Oi <- Oi * Patm / 100
  
  Ko <- Ko25*arrh(Tleaf, Eo)
  Kc <- Kc25*arrh(Tleaf, Ec)
  Km <- Kc * (1.0 + Oi / Ko)
  
return(Km)
}

# Vcmax temperature response (Arrhenius)
TVcmax <- function(Tleaf, EaV, delsC, EdVC){
  
  if(EdVC > 0){
    V1 <- 1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15)))
    V2 <- 1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tk(Tleaf))))
    f <- V1/V2
  } else f <- 1
  
  exp((Tleaf-25)*EaV/(.Rgas()*Tk(Tleaf)*Tk(25))) * f
}

# Jmax temperature response (Arrhenius)
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/.Rgas()/298.15)
  J2 <- 1+exp((Tk(Tleaf)*delsJ-EdVJ)/.Rgas()/Tk(Tleaf))
  exp(EaJ/.Rgas()*(1/298.15 - 1/Tk(Tleaf)))*J1/J2
}





