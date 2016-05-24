class.sort<-function(sort, method, output, lang){
  if (output == "phi"){
    if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
      if (sort <= 0.35) {return("Very well sorted")}
      if ((sort > 0.35) & (sort <= 0.5)) {return("Well sorted")}
      if ((sort > 0.5) & (sort <= 0.7)) {return("Moderately well sorted")}
      if ((sort > 0.7) & (sort <= 1)) {return("Moderately sorted")}
      if ((sort > 1) & (sort <= 2)) {return("Poorly sorted")}
      if ((sort > 2) & (sort <= 4)) {return("Very poorly sorted")}
      if (sort > 4) {return("Extremely poorly sorted")}
      }
    if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
      if (sort <= 0.35) {return("Muito bem selecionado")}
      if ((sort > 0.35) & (sort <= 0.5)) {return("Bem selecionado")}
      if ((sort > 0.5) & (sort <= 0.7)) {return("Moderadamente bem selecionado")}
      if ((sort > 0.7) & (sort<=1)) {return("Moderadamente selecionado")}
      if ((sort > 1) & (sort<=2)) {return("Pobremente selecionado")}
      if ((sort > 2) & (sort<=4)) {return("Muito pobremente selecionado")}
      if (sort > 4) {return("Extremamente mal selecionado")}
      }
    }
  if (output == "metric"){
    if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e"){
      if (sort <= 1.27) {return("Very well sorted")}
      if ((sort > 1.27) & (sort <= 1.41)) {return("Well sorted")}
      if ((sort > 1.41) & (sort <= 1.62)) {return("Moderately well sorted")}
      if ((sort > 1.62) & (sort <= 2)) {return("Moderately sorted")}
      if ((sort > 2) & (sort <= 4)) {return("Poorly sorted")}
      if ((sort > 4) & (sort <= 16)) {return("Very poorly sorted")}
      if (sort > 16) {return("Extremely poorly sorted")}
      }
    if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p"){
      if (sort <= 1.27) {return("Muito bem selecionado")}
      if ((sort > 1.27) & (sort <= 1.41)) {return("Bem selecionado")}
      if ((sort > 1.41) & (sort <= 1.62)) {return("Moderadamente bem selecionado")}
      if ((sort > 1.62) & (sort <= 2)) {return("Moderadamente selecionado")}
      if ((sort > 2) & (sort <= 4)) {return("Pobremente selecionado")}
      if ((sort > 4) & (sort <= 16)) {return("Muito pobremente selecionado")}
      if (sort > 16) {return("Extremamente pobremente selecionado")}
      }
    }
  }

