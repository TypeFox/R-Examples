renderText({
  Peso<<- round(runif(1, 60, 80),digits=1)
  DS_Peso<<-round(runif(1, 2, 8),digits=1)
  #Cues_Val<-round(rnorm(1, mean=Peso, sd=DS_Peso),digits=1)
  Bus_Area<<-sample(c("lower", "upper", "both", "middle"), 1)
  
  if (Bus_Area=="lower") {
    Cues_Val<<-round(rnorm(1, mean=Peso, sd=DS_Peso),digits=0)
    Cuestion<-paste("un estudiante pese menos de ",Cues_Val,"kg.")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Peso, sd_CalDis=DS_Peso, a_CalDis=Cues_Val))
  }
  else if(Bus_Area=="upper") {
    Cues_Val<<-round(rnorm(1, mean=Peso, sd=DS_Peso),digits=0)
    Cuestion<-paste("el peso de un estudiante sea mayor a ",Cues_Val,"kg.")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Peso, sd_CalDis=DS_Peso, a_CalDis=Cues_Val))
  }
  else if(Bus_Area=="both") {
    Cues_Val<-round(rnorm(25, mean=Peso, sd=DS_Peso),digits=0)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    Cuestion<-paste("un estudiante pese menos de ",a_val,"kg o más de ",b_val,"kg.")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Peso, sd_CalDis=DS_Peso, a_CalDis=a_val,b_CalDis=b_val))
  }
  else {
    Cues_Val<-round(rnorm(3, mean=Peso, sd=DS_Peso),digits=0)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    Cuestion<-paste("el peso de un estudiante se encuentre entre ",a_val,"kg y ",b_val,"kg como máximo.")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis="rnorm", tail_CalDis=Bus_Area, mu_CalDis=Peso, sd_CalDis=DS_Peso, a_CalDis=a_val,b_CalDis=b_val))
  }
  
  #Valor_a<-input$a_CalDis
  paste("En una escuela se tienen los registros de los pesos de los estudiantes, de los que se tiene que la media de los pesos 
es ", Peso,"kg y la desviación típica ",DS_Peso,"kg suponiendo que los pesos se distribuyen normalmente, encontrar la probabilidad de que ",Cuestion)
  })