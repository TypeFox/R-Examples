class.skew<-function(skew, method, output, lang){
  if (method == "moment"){
    if (output == "phi"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (skew <= -1.3) {return("Very negative")}
        if ((skew > -1.3) & (skew <= -0.43)) {return("Negative")}
        if ((skew > -0.43) & (skew <= 0.43)) {return("Approximately symmetrical")}
        if ((skew > 0.43) & (skew <= 1.3)) {return("Positive")}
        if (skew > 1.3) {return("Very positive")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (skew <= -1.3) {return("Muito negativa")}
        if ((skew > -1.3) & (skew <= -0.43)) {return("Negativa")}
        if ((skew > -0.43) & (skew <= 0.43)) {return("Aproximadamente sim\u00E9trica")}
        if ((skew > 0.43) & (skew <= 1.3)) {return("Positiva")}
        if (skew > 1.3) {return("Muito positiva")}
        }
      }
    if (output == "metric"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (skew <= -1.3) {return("Very positive")}
        if ((skew > -1.3) & (skew <= -0.43)) {return("Positive")}
        if ((skew > -0.43) & (skew <= 0.43)) {return("Approximately symmetrical")}
        if ((skew > 0.43) & (skew <= 1.3)) {return("Negative")}
        if (skew > 1.3) {return("Very negative")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (skew <= -1.3) {return("Muito positiva")}
        if ((skew > -1.3) & (skew <= -0.43)) {return("Positiva")}
        if ((skew > -0.43) & (skew <= 0.43)) {return("Aproximadamente sim\u00E9trica")}
        if ((skew > 0.43) & (skew <= 1.3)) {return("Negativa")}
        if (skew > 1.3) {return("Muito negativa")}
        }
      }
    }
  if (method=="folk"|method=="mcA"|method=="mcB"|method=="trask"|method=="otto"){
    if (output == "phi"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if ((skew > -1) & (skew <= -0.3)) {return("Very negative")}
        if ((skew > -0.3) & (skew <= -0.1)) {return("Negative")}
        if ((skew > -0.1) & (skew <= 0.1)) {return("Approximately symmetrical")}
        if ((skew > 0.1) & (skew <= 0.3)) {return("Positive")}
        if ((skew > 0.3) & (skew <= 1)) {return("Very positive")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if ((skew > -1) & (skew <= -0.3)) {return("Muito negativa")}
        if ((skew > -0.3) & (skew <= -0.1)) {return("Negativa")}
        if ((skew > -0.1) & (skew <= 0.1)) {return("Aproximadamente sim\u00E9trica")}
        if ((skew > 0.1) & (skew <= 0.3)) {return("Positiva")}
        if ((skew > 0.3) & (skew <= 1)) {return("Muito positiva")}
        }
      }
    if (output == "metric"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if ((skew > -1) & (skew <= -0.3)) {return("Very positive")}
        if ((skew > -0.3) & (skew <= -0.1)) {return("Positive")}
        if ((skew > -0.1) & (skew <= 0.1)) {return("Approximately symmetrical")}
        if ((skew > 0.1) & (skew <= 0.3)) {return("Negative")}
        if ((skew > 0.3) & (skew <= 1)) {return("Very negative")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if ((skew > -1) & (skew <= -0.3)) {return("Muito positiva")}
        if ((skew > -0.3) & (skew <= -0.1)) {return("Positiva")}
        if ((skew > -0.1) & (skew <= 0.1)) {return("Aproximadamente sim\u00E9trica")}
        if ((skew > 0.1) & (skew <= 0.3)) {return("Negativa")}
        if ((skew > 0.3) & (skew <= 1)) {return("Muito negativa")}
        }
      }
    }
  }

