`Abschlussnote` <-
structure(function(x,y,z){     
   x.note <- (x/100)*30
   y.note <- (y/100)*30
   z.note <- (z/100)*40
   abschluss <- x.note + y.note + z.note
   cat (c("Abschlussnote: "), abschluss)
 }
, comment = "Funktion zur Errechnung einer gedachten Abschlussnote")
