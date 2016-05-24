renderText({
  Temper<<- round(runif(1, 16.5, 18.5),digits=1)
  DS_Temper<<-round(runif(1, 4, 6),digits=1)
  #Cues_Val<-round(rnorm(1, mean=Temper, sd=DS_Temper),digits=1)
  Bus_Area<<-sample(c("lower", "upper", "both", "middle"), 1)
  
  if (Bus_Area=="lower") {
    Cues_Val<<-round(rnorm(1, mean=Temper, sd=DS_Temper),digits=0)
    Cuestion<-paste("media sea menor a ",Cues_Val," grados?")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Temper, sd_CalDis=DS_Temper, a_CalDis=Cues_Val))
  }
  else if(Bus_Area=="upper") {
    Cues_Val<<-round(rnorm(1, mean=Temper, sd=DS_Temper),digits=0)
    Cuestion<-paste("media sea mayor a ",Cues_Val," grados?")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Temper, sd_CalDis=DS_Temper, a_CalDis=Cues_Val))
  }
  else if(Bus_Area=="both") {
    Cues_Val<<-round(rnorm(20, mean=Temper, sd=DS_Temper),digits=0)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    Cuestion<-paste("media sea menos de ",a_val,"grados o más de ",b_val,"grados?")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Temper, sd_CalDis=DS_Temper, a_CalDis=a_val,b_CalDis=b_val))
  }
  else {
    Cues_Val<<-round(rnorm(3, mean=Temper, sd=DS_Temper),digits=0)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    Cuestion<-paste("media se encuentre entre ",a_val,"grados y ",b_val,"grados como máximo?")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Temper, sd_CalDis=DS_Temper, a_CalDis=a_val,b_CalDis=b_val))
  }
  
  #Valor_a<-input$a_CalDis
  paste("Se estima que en el estado de Aguascalientes las temperaturas siguen una distribución normal con una media  de ",Temper," grados centígrados y una 
        desviación típica de ",DS_Temper," grados centígrados, ¿Cuál es la probabilidad de que en un día cualquiera la 
        temperatura ",Cuestion)
})