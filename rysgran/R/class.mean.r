class.mean<-function(mean, output, lang){
  if (output == "phi"){
    if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
      if (mean <= -8) {return("Boulder")}
      if ((mean > -8) & (mean <= -6)) {return("Cobble")}
      if ((mean > -6) & (mean <= -2)) {return("Pebble")}
      if ((mean > -2) & (mean <= -1)) {return("Granules")}
      if ((mean > -1) & (mean <= 0)) {return("Very coarse sand")}
      if ((mean > 0) & (mean <= 1)) {return("Coarse sand")}
      if ((mean > 1) & (mean <= 2)) {return("Medium sand")}
      if ((mean > 2) & (mean <= 3)) {return("Fine sand")}
      if ((mean > 3) & (mean <= 4)) {return("Very fine sand")}
      if ((mean > 4) & (mean <= 5)) {return("Coarse silt")}
      if ((mean > 5) & (mean <= 6)) {return("Medium silt")}
      if ((mean > 6) & (mean <= 7)) {return("Fine silt")}
      if ((mean > 7) & (mean <= 8)) {return("Very fine silt")}
      if (mean > 8) {return("Clay")}
      }
    if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
      if (mean <= -8) {return("Matac\u00E3o")}
      if ((mean > -8) & (mean <= -6)) {return("Calhau")}
      if ((mean > -6) & (mean <= -2)) {return("Seixo")}
      if ((mean > -2) & (mean <= -1)) {return("Gr\u00E2nulo")}
      if ((mean > -1) & (mean <= 0)) {return("Areia muito grossa")}
      if ((mean > 0) & (mean <= 1)) {return("Areia grossa")}
      if ((mean > 1) & (mean <= 2)) {return("Areia m\u00E9dia")}
      if ((mean > 2) & (mean <= 3)) {return("Areia fina")}
      if ((mean > 3) & (mean <= 4)) {return("Areia muito fina")}
      if ((mean > 4) & (mean <= 5)) {return("Silte grosso")}
      if ((mean > 5) & (mean <= 6)) {return("Silte m\u00E9dio")}
      if ((mean > 6) & (mean <= 7)) {return("Silte fino")}
      if ((mean > 7) & (mean <= 8)) {return("Silte muito fino")}
      if (mean > 8) {return("Argila")}
      }
    }
  if (output == "metric"){
    if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
      if (mean >= 256000) {return("Boulder")}
      if ((mean > 64000) & (mean <= 256000)) {return("Cobble")}
      if ((mean > 4000) & (mean <= 64000)) {return("Pebble")}
      if ((mean > 2000) & (mean <= 4000)) {return("Granules")}
      if ((mean > 1000) & (mean <= 2000)) {return("Very coarse sand")}
      if ((mean > 500) & (mean <= 1000)) {return("Coarse sand")}
      if ((mean > 250) & (mean <= 500)) {return("Medium sand")}
      if ((mean > 125) & (mean <= 250)) {return("Fine sand")}
      if ((mean > 63) & (mean <= 125)) {return("Very fine sand")}
      if ((mean > 31) & (mean <= 63)) {return("Coarse silt")}
      if ((mean > 16) & (mean <= 31)) {return("Medium silt")}
      if ((mean > 8) & (mean <= 16)) {return("Fine silt")}
      if ((mean > 4) & (mean <= 8)) {return("Very fine silt")}
      if (mean <= 4) {return("Clay")}
      }
    if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
      if (mean >= 256000) {return("Matac\u00E3o")}
      if ((mean > 64000) & (mean <= 256000)) {return("Calhau")}
      if ((mean > 4000) & (mean <= 64000)) {return("Seixo")}
      if ((mean > 2000) & (mean <= 4000)) {return("Gr\u00E2nulo")}
      if ((mean > 1000) & (mean <= 2000)) {return("Areia muito grossa")}
      if ((mean > 500) & (mean <= 1000)) {return("Areia grossa")}
      if ((mean > 250) & (mean <= 500)) {return("Areia m\u00E9dia")}
      if ((mean > 125) & (mean <= 250)) {return("Areia fina")}
      if ((mean > 63) & (mean <= 125)) {return("Areia muito fina")}
      if ((mean > 31) & (mean <= 63)) {return("Silte grosso")}
      if ((mean > 16) & (mean <= 31)) {return("Silte m\u00E9dio")}
      if ((mean > 8) & (mean <= 16)) {return("Silte fino")}
      if ((mean > 4) & (mean <= 8)) {return("Silte muito fino")}
      if (mean <= 4) {return("Argila")}
      }
    }
  }
