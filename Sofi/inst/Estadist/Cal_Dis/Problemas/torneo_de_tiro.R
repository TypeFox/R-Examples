renderText({
  Distribu<-"rbinom"
  tiros<<- round(runif(1, 5, 45),digits=0)
  probabi<<-round(runif(1, .4, .99),digits=2)
  lower_bound<-sample(c("open", "closed"), 1)
  upper_bound<-sample(c("open", "closed"), 1)
  
  #Cues_Val<-round(rnorm(1, mean=Peso, sd=DS_Peso),digits=1)
  Bus_Area<<-sample(c("lower", "upper", "both", "middle", "equal"), 1)

  if (Bus_Area=="lower") {
    Aciertos<<-rbinom(1, tiros, probabi)
    if (lower_bound=="closed"){
      Cuestion<-paste("exactamente en ",Aciertos," ocasiones o menos? ")
    }
    else{
      Cuestion<-paste("menos de ",Aciertos," ocasiones? ")
    }
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis=Distribu, tail_CalDis=Bus_Area, a_CalDis=Aciertos, n_CalDis=tiros, 
                                      p_CalDis=probabi, lower_bound_CalDis=lower_bound, upper_bound_CalDis=upper_bound))
  }
  else if(Bus_Area=="upper") {
    Aciertos<<-rbinom(1, tiros, probabi)
    if (lower_bound=="closed"){
      Cuestion<-paste("exactamente en ",Aciertos," ocasiones o mas? ")
    }
    else{
      Cuestion<-paste("mas de ",Aciertos," veces? ")
    }
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis=Distribu, tail_CalDis=Bus_Area, a_CalDis=Aciertos, n_CalDis=tiros, 
                                      p_CalDis=probabi, lower_bound_CalDis=lower_bound, upper_bound_CalDis=upper_bound))
  }
  else if(Bus_Area=="both") {
    
    Cues_Val<<-rbinom(3, tiros, probabi)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    if (lower_bound=="closed" & upper_bound=="closed"){
      Cuestion<-paste("exactamente en ",a_val," o menos ocasiones a la vez que pudiera acertar ",b_val," ocasiones exactamente o mas?")
    }
    else if (lower_bound=="closed" & upper_bound=="open"){
      Cuestion<-paste("exactamente en ",a_val," o menos ocasiones a la vez que pudiera acertar mas de ",b_val," ocasiones?")
    }
    else if (lower_bound=="open" & upper_bound=="closed"){
      Cuestion<-paste("menos de ",a_val," ocasiones a la vez que pudiera acertar ",b_val," ocasiones exactamente o mas?")
    }
    else
      {Cuestion<-paste("menos de ",a_val,"veces o más de ",b_val,"veces.")}
    
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis=Distribu, tail_CalDis=Bus_Area, a_CalDis=a_val, b_CalDis=b_val, n_CalDis=tiros, 
                                      p_CalDis=probabi, lower_bound_CalDis=lower_bound, upper_bound_CalDis=upper_bound))
  }
  else if(Bus_Area=="middle") {
    Cues_Val<<-rbinom(3, tiros, probabi)
    a_val<<-min(Cues_Val)
    b_val<<-max(Cues_Val)
    if (lower_bound=="closed" & upper_bound=="closed"){
      Cuestion<-paste("entre exactamente ",a_val,"  y ",b_val," ?")
    }
    else if (lower_bound=="closed" & upper_bound=="open"){
      Cuestion<-paste("mas o exactamente ",a_val," veces y menos de ",b_val," veces?")
    }
    else if (lower_bound=="open" & upper_bound=="closed"){
      Cuestion<-paste("mas de ",a_val," ocasiones y exactamente ",b_val," o menos ocasiones?")
    }
    else
      {Cuestion<-paste("mas de ",a_val," veces y menos de ",b_val,"veces?")}
    
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis=Distribu, tail_CalDis=Bus_Area, a_CalDis=a_val, b_CalDis=b_val, n_CalDis=tiros, 
                                      p_CalDis=probabi, lower_bound_CalDis=lower_bound, upper_bound_CalDis=upper_bound))
    
  }
  else {
    #Cues_Val<-round(rnorm(3, mean=Peso, sd=DS_Peso),digits=0)
    #a_val<-min(Cues_Val)
    #b_val<-max(Cues_Val)
    Aciertos<<-rbinom(1, tiros, probabi)
    Cuestion<-paste("exactamente en ",Aciertos," ocasiones?")
    Valor_Cues<<-as.numeric(Valor_Fin(dist_CalDis=Distribu, tail_CalDis=Bus_Area, a_CalDis=Aciertos, n_CalDis=tiros, 
                                      p_CalDis=probabi, lower_bound_CalDis=lower_bound, upper_bound_CalDis=upper_bound))
  }
  
  #Valor_a<-input$a_CalDis
  paste("Un  hombre que práctica para un torneo de tiro con arco dispara ",tiros,"veces, sus estadísticas dicen que tiene ",probabi,
        " de probabilidades de acierto.¿Cuál es la probabilidad de que acierte ",Cuestion)
})