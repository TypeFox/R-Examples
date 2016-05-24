class.kurt<-function(kurt, method, output, lang){
  if (method == "moment"){
    if (output == "phi"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (kurt <= 1.70) {return("Very platykurtic")}
        if ((kurt > 1.70) & (kurt <= 2.55)) {return("Platykurtic")}
        if ((kurt > 2.55) & (kurt <= 3.70)) {return("Mesokurtic")}
        if ((kurt > 3.70) & (kurt <= 7.40)) {return("Leptokurtic")}
        if ((kurt > 7.40) & (kurt <= 15)) {return("Very leptokurtic")}
        if (kurt > 15) {return("Extremely leptokurtic")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (kurt <= 1.70) {return("Muito platic\u00FArtica")}
        if ((kurt > 1.70) & (kurt <= 2.55)) {return("Platic\u00FArtica")}
        if ((kurt > 2.55) & (kurt <= 3.70)) {return("Mesoc\u00FArtica")}
        if ((kurt > 3.70) & (kurt <= 7.40)) {return("Leptoc\u00FArtica")}
        if ((kurt > 7.40) & (kurt <= 15)) {return("Muito leptoc\u00FArtica")}
        if (kurt > 15) {return("Extremamente leptoc\u00FArtica")}
        }
      }
    if (output == "metric"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (kurt <= 1.70) {return("Very platykurtic")}
        if ((kurt > 1.70) & (kurt <= 2.55)) {return("Platykurtic")}
        if ((kurt > 2.55) & (kurt <= 3.70)) {return("Mesokurtic")}
        if ((kurt > 3.70) & (kurt <= 7.40)) {return("Leptokurtic")}
        if ((kurt > 7.40) & (kurt <= 15)) {return("Very leptokurtic")}
        if (kurt > 15) {return("Extremely leptokurtic")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (kurt <= 1.70) {return("Muito platic\u00FArtica")}
        if ((kurt > 1.70) & (kurt <= 2.55)) {return("Platic\u00FArtica")}
        if ((kurt > 2.55) & (kurt <= 3.70)) {return("Mesoc\u00FArtica")}
        if ((kurt > 3.70) & (kurt <= 7.40)) {return("Leptoc\u00FArtica")}
        if ((kurt > 7.40) & (kurt <= 15)) {return("Muito leptoc\u00FArtica")}
        if (kurt > 15) {return("Extremamente leptoc\u00FArtica")}
        }
      }
    }
  if (method=="folk"|method=="mcA"|method=="mcB"|method=="trask"|method=="otto"){
    if (output == "phi"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (kurt <= 0.67) {return("Very platykurtic")}
        if ((kurt > 0.67) & (kurt <= 0.9)) {return("Platykurtic")}
        if ((kurt > 0.9) & (kurt <= 1.11)) {return("Mesokurtic")}
        if ((kurt > 1.11) & (kurt <= 1.5)) {return("Leptokurtic")}
        if ((kurt > 1.5) & (kurt <= 3)) {return("Very leptokurtic")}
        if (kurt > 3) {return("Extremely leptokurtic")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (kurt <= 0.67) {return("Muito platic\u00FArtica")}
        if ((kurt > 0.67) & (kurt <= 0.9)) {return("Platic\u00FArtica")}
        if ((kurt > 0.9) & (kurt <= 1.11)) {return("Mesoc\u00FArtica")}
        if ((kurt > 1.11) & (kurt <= 1.5)) {return("Leptoc\u00FArtica")}
        if ((kurt > 1.5) & (kurt <= 3)) {return("Muito leptoc\u00FArtica")}
        if (kurt > 3) {return("Extremamente leptoc\u00FArtica")}
        }
      }
    if (output == "metric"){
      if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
        if (kurt <= 0.67) {return("Very platykurtic")}
        if ((kurt > 0.67) & (kurt <= 0.9)) {return("Platykurtic")}
        if ((kurt > 0.9) & (kurt <= 1.11)) {return("Mesokurtic")}
        if ((kurt > 1.11) & (kurt <= 1.5)) {return("Leptokurtic")}
        if ((kurt > 1.5) & (kurt <= 3)) {return("Very leptokurtic")}
        if (kurt > 3) {return("Extremely leptokurtic")}
        }
      if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
        if (kurt <= 0.67) {return("Muito platic\u00FArtica")}
        if ((kurt > 0.67) & (kurt <= 0.9)) {return("Platic\u00FArtica")}
        if ((kurt > 0.9) & (kurt <= 1.11)) {return("Mesoc\u00FArtica")}
        if ((kurt > 1.11) & (kurt <= 1.5)) {return("Leptoc\u00FArtica")}
        if ((kurt > 1.5) & (kurt <= 3)) {return("Muito leptoc\u00FArtica")}
        if (kurt > 3) {return("Extremamente leptoc\u00FArtica")}
        }
      }
    }
  }
