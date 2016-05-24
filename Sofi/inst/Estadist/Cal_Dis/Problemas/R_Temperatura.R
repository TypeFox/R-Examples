renderUI({
  if (!input$Ayuda_visible) return()
  
  if(Bus_Area %in% c("lower","upper")) {
    withMathJax(
      helpText('Es simple, lo primero es solo sustituir los valores de la siguiente forma: 
               $$Z=\\frac { X-\\mu  }{ \\sigma  }=\\frac { ',Cues_Val,'-\\',Temper,'  }{ \\',DS_Temper,' }=
               ',signif((Cues_Val-Temper)/DS_Temper,4),' $$'),
      helpText('Ahora la temperatura sigue una distribución normal con \\(\\mu\\) =0 y \\(\\sigma\\) =1')
      )
  }
  else {
    withMathJax(
      helpText('Es simple, lo primero es solo sustituir los valores de la siguiente forma \\(\\ { Z }_{ 1 }\\): 
               $${ Z }_{ 1 }=\\frac { X-\\mu  }{ \\sigma  }=\\frac { ',a_val,'-\\',Temper,'  }{ \\',DS_Temper,' }=
               ',signif((a_val-Temper)/DS_Temper,4),' $$ Y para \\(\\ { Z }_{ 2 }\\): $${ Z }_{ 2 }=\\frac { X-\\mu  }
               { \\sigma  }=\\frac { ',b_val,'-\\',Temper,'  }{ \\',DS_Temper,' }=
               ',signif((b_val-Temper)/DS_Temper,4),' $$')
      )
  }
  
  
  
  })