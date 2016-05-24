# Funktionen zur Erzeugung von Zykloiden v1.01
# Peter Biber 9.11.2008


#-------------------------------------------------------------------------------

# Groesster gemeinsamer Teiler zweier ganzer Zahlen (Bronstein S. 333)
# Euklidischer Algorithmus
ggT <- function(a, b) {
       if (((a == round(a)) & (b = round(b))) &
           ((a > 0) & (b > 0))) {

           a1 <- max(a, b)
           zahl <- min(a, b)

           repeat {
                 rest <- a1 %% zahl
                 if (rest > 0) {
                    a1 <- zahl
                    zahl <- rest
                 }
                 else break
           } # repeat
           return(zahl)
       } # if: Wenn Zahlen integer und positiv
       else  return(NA)
} # function ggT

#-------------------------------------------------------------------------------

# Kleinstes gemeinsames Vielfaches (Bronstein S. 334)
kgV <- function(a, b) {
       if (((a == round(a)) & (b = round(b))) &
           ((a > 0) & (b > 0))) {

           zahl <- a * b / ggT(a, b)
           
           return(zahl)

       } # if: a,b, integer und positiv
       else return(NA)
} # function kgV

#-------------------------------------------------------------------------------

# npeaks: Gibt die Anzahl der Spitzen bzw. Schlaufen des Zykloiden zurück.
# A, a: natürliche Zahlen, kennzeichnen Verhältnis der Radien des fixen (A)
# und des bewegten Kreises (a)

npeaks <- function(A, a) {
          if (((A == round(A)) & (a = round(a))) &
              ((A > 0) & (a > 0))) {

             n <- kgV(A, a) / a
             return(n)
          }
          else return(NA)

} # function npeaks

#-------------------------------------------------------------------------------

# zykloid: Liefert einen data frame zurück, der die x-y Werte eines
# Hypo- oder Epizykloiden enthält.
# A: Radius des festen Kreises    A, a: Natürliche Zahlen
# a: Radius des rollenden Kreises
# lambda: Relative Position des Spurpunktes auf a
# hypo: TRUE - rollender Kreis läuft im Inneren des festen Kreises, 
#       FALSE - rollender Kreis läuft auf der Außenseite des festen Kreises
# steps: Zahl der Winkelschritte pro Umlauf für die Bewegung des Mittelpunkts
# des rollenden Kreises.
# start: Startwinkel im festen Kreis gegen den Uhrzeigersinn

# TODO: Teste Startpunkt für negatives lambda.

zykloid <- function(A, a, lambda, hypo = TRUE, steps = 360, start = pi/2) {

           if (((A == round(A)) & (a == round(a))) &
               ((A > 0) & (a > 0))) {

               # Zahl der Umläufe des rollenden Kreises bis zum Schluss
               # der Figur.
               n <- kgV(A, a) / A
               
               # Leeren, aber richtig dimensionierten data.frame erzeugen
               # + 1 Zeile ist notwendig wegen Schluss
               nr <- rep(NA, n * steps + 1)
               M  <- data.frame(x = nr, y = nr)
               colnames(M) <- c("x", "y")
               
               if (hypo) { a <- -1 * a }
               
               i   <- c(1:(n * steps + 1))      # index
               phi <- (i - 1) * 2 * pi / steps  # phi: Winkel der Position des Mittelpunktes des rollenden Kreises
                                                # i - 1, damit mit 0 begonnen wird 
               
               M[, "x"]   <- (A + a) * cos(phi + start) - lambda * a * cos((A + a)/a * phi + start)
               M[, "y"]   <- (A + a) * sin(phi + start) - lambda * a * sin((A + a)/a * phi + start)
               
               return(M)
           }
           else return(NA)
} # function zykloid

#-------------------------------------------------------------------------------

# zykloid.scaleA: Erzeugt einen Zykloiden wobei der Radius des fixen Kreises
# vorgegeben werden kann: RadiusA (default: Einheitskreis).
# Ebenso dessen Mittelpunkt mit den Koordinaten Cx und Cy (default Ursprung).
# Weitere Einstellungen wie funktion zykloid.

zykloid.scaleA <- function(A, a, lambda, hypo = TRUE, Cx = 0, Cy = 0, RadiusA = 1, steps = 360, start = pi/2) {

           # Rohzykloid in M
           M <- zykloid(A, a, lambda, hypo, steps, start)

           # Eigentliche Umrechnung
           # Ein data.frame kommt nur heraus, wenn es geklappt hat (sonst NA)
           if (is.data.frame(M)) { 

               M[, "x"] <- M[, "x"] * RadiusA / A + Cx
               M[, "y"] <- M[, "y"] * RadiusA / A + Cy               

               return(M)
               
           } # if 
           else return(NA)

} # function zykloid.scaleA

#-------------------------------------------------------------------------------

# zykloid.scaleAa: Erzeugt einen Zykloiden wobei der Radius des Umkreises
# der Figur vorgegeben werden kann: RadiusAa (default: Einheitskreis).
# Ebenso dessen Mittelpunkt mit den Koordinaten Cx und Cy (default Ursprung).
# Alle weiteren Einstellungen wie funktion zykloid.scaleA

zykloid.scaleAa <- function(A, a, lambda, hypo = TRUE, Cx = 0, Cy = 0, RadiusAa = 1, steps = 360, start = pi/2) {

           # Rohzykloid in M
           M <- zykloid(A, a, lambda, hypo, steps, start)

           # Eigentliche Umrechnung
           # Ein data.frame kommt nur heraus, wenn es geklappt hat (sonst NA)
           if (is.data.frame(M)) { 
           
               if (hypo) Norm <- abs(A - a) + abs(lambda) * a # Gilt auch für a > A
               else      Norm <-     A + a  + abs(lambda) * a

               M[, "x"] <- M[, "x"] * RadiusAa / Norm + Cx
               M[, "y"] <- M[, "y"] * RadiusAa / Norm + Cy               

               return(M)
               
           } # if 
           else return(NA)

} # function zykloid.scaleAa

#-------------------------------------------------------------------------------

# zykloid.scaleP: Erzeugt einen Zykloiden wobei der Radius des Kreises, auf dem
# die Spitzen der Figur liegen, vorgegeben werden kann: RadiusP (default: 
# Einheitskreis).
# Ebenso dessen Mittelpunkt mit den Koordinaten Cx und Cy (default Ursprung).
# Alle weiteren Einstellungen wie funktion zykloid.scaleA

zykloid.scaleP <- function(A, a, lambda, hypo = TRUE, Cx = 0, Cy = 0, RadiusP = 1, steps = 360, start = pi/2) {

           # Rohzykloid in M
           M <- zykloid(A, a, lambda, hypo, steps, start)

           # Eigentliche Umrechnung
           # Ein data.frame kommt nur heraus, wenn es geklappt hat (sonst NA)
           if (is.data.frame(M)) { 
           
               if (hypo) Norm <-   abs(A + a * (abs(lambda) - 1))
               else      Norm <-   abs(A + a * (1 - abs(lambda)))

               M[, "x"] <- M[, "x"] * RadiusP / Norm + Cx
               M[, "y"] <- M[, "y"] * RadiusP / Norm + Cy               

               return(M)
               
           } # if 
           else return(NA)

} # function zykloid.scaleP


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------