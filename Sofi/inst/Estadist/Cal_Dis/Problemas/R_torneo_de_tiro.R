renderUI({
  if (!input$Ayuda_visible) return()
  
  if(Bus_Area %in% c("lower","upper","equal")) {
    withMathJax(
      helpText('A modo de ejemplo, tenemos a \\(\\ k=\\) ',Aciertos,' \\(\\ n=\\) ',tiros,' y \\(\\ p=\\) ',probabi,'lo primero es sustituir los valores de la siguiente forma: 
               $$p(X=k)=\\left( \\begin{matrix} n \\\\ k \\end{matrix} \\right) { p }^{ k }\\bullet { q }^{ n-k }
                =p(X=',Aciertos,')=\\left( \\begin{matrix} ',tiros,' \\\\ ',Aciertos,' \\end{matrix} \\right) { ',probabi,' }^{ ',Aciertos,' }\\bullet { ',1-probabi,' }^{ ',tiros,'-',Aciertos,' }
               =',signif(factorial(tiros)/(factorial(Aciertos)*factorial(tiros-Aciertos))*probabi^Aciertos*(1-probabi)^(tiros-Aciertos),4),' $$')
      )
  }
  else {
    withMathJax(
      helpText('A modo de ejemplo, tenemos a \\(\\ k=\\) ',a_val,' \\(\\ n=\\) ',tiros,' y \\(\\ p=\\) ',probabi,'lo primero es sustituir los valores de la siguiente forma para \\(\\ { p }_{ 1 }\\): 
               $${ p }_{ 1 }(X=k)=\\left( \\begin{matrix} n \\\\ k \\end{matrix} \\right) { p }^{ k }\\bullet { q }^{ n-k }
               ={ p }_{ 1 }(X=',a_val,')=\\left( \\begin{matrix} ',tiros,' \\\\ ',a_val,' \\end{matrix} \\right) { ',probabi,' }^{ ',a_val,' }\\bullet { ',1-probabi,' }^{ ',tiros,'-',a_val,' }
               =',signif(factorial(tiros)/(factorial(a_val)*factorial(tiros-a_val))*probabi^a_val*(1-probabi)^(tiros-a_val),4),' $$ '),
      helpText( 'Y para \\(\\ k=\\) ',b_val,' tenemos que \\(\\ { p }_{ 2 }\\):
                  $${ p }_{ 2 }(X=k)=\\left( \\begin{matrix} n \\\\ k \\end{matrix} \\right) { p }^{ k }\\bullet { q }^{ n-k }
                ={ p }_{ 2 }(X=',b_val,')=\\left( \\begin{matrix} ',tiros,' \\\\ ',b_val,' \\end{matrix} \\right) { ',probabi,' }^{ ',b_val,' }\\bullet { ',1-probabi,' }^{ ',tiros,'-',b_val,' }
                =',signif(factorial(tiros)/(factorial(b_val)*factorial(tiros-b_val))*probabi^b_val*(1-probabi)^(tiros-b_val),4),' $$')
      )
  }
  })

