renderText({
  
  L = NULL ; U = NULL ; error = FALSE
  
  if (input$dist == "runif"){
    L = input$min ; U = input$max
    if (L > U){
      error = TRUE
    }
  }
  
  if (error)
    paste0()
  
  else{
    pop = parent()
    m_pop =  round(mean(pop),2)
    s_pop = round(sd(pop),2)
    
    n = input$n
    se=round(s_pop/sqrt(n),2)
    paste("De acuerdo con el teorema del límite central, la distribución de las medias de la muestra 
          (la distribución de muestreo) debe ser casi normal. La media de la distribución de muestreo 
          debe ser aproximadamente igual a la media de la población (", m_pop, ") 
          y el error estándar (la desviación estándar de la media de muestra) debe ser aproximadamente 
          igual a la SD de la población dividida por la raíz cuadrada del tamaño de la muestra (", s_pop,
          "/sqrt(",n, ") =", se,").")
  }
})