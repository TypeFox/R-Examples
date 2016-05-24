library(shiny)

#source('./helper/chiTail.R')
#source('./helper/FTail.R')
#source('./helper/normTail.R')
source('./helper/Tails.R') #Calculadora de Distribución
# set mirror
options(repos=structure(c(CRAN="http://cran.rstudio.com")))
#options(encoding = 'UTF-8')

if (!("shiny" %in% names(installed.packages()[,"Package"]))) {install.packages("shiny")}
suppressMessages(library(shiny, quietly = TRUE))

if (!("openintro" %in% names(installed.packages()[,"Package"]))) {install.packages("openintro")}
suppressMessages(library(openintro, quietly = TRUE))

#if (!("knitr" %in% names(installed.packages()[,"Package"]))) {install.packages("knitr")}
#suppressMessages(library(knitr, quietly = TRUE))

defaults = list("tail_CalDis" = "lower",
                "lower_bound_CalDis" = "open",
                "upper_bound_CalDis" = "open")


generateResponse = function(response){
  if (response==1){
    print(sample(list("¡Correcto!","¡En el clavo!","¡Acertaste!","¡Felicidades!","Bien"),1)[[1]])
  }
  else if (response ==2){
    print(sample(list("Casi.","Cerca.","Sólo un poco fuera.."),1)[[1]])
  }
  else if (response == 3){
    print(sample(list("Estas frío...","Algo lejos...","Inténtalo de nuevo.", "Estas  perdido"),1)[[1]])
  }
}

shinyServer(function(input, output, session,updateNum)
{ 
  
  Valor_Cues <- 1 #current correlation
  score <- 0 #user's score
  answered <- FALSE # an indicator for whether question has been answered
  
##CalDis####
  output$tail_CalDis = renderUI(
  {
    #print("tail_CalDis")
    if (input$dist_CalDis == "rbinom")
    {
      selectInput(inputId = "tail_CalDis",
                  label = "Encuentra Área:",
                  choices = c("Cola inferior"="lower", 
                              "Cola superior"="upper", 
                              "Ambas colas"="both",
                              "Intermedio"="middle",
                              "Igualdad"="equal"),
                  selected = "lower")
    }
    else if (input$dist_CalDis == "rf" | input$dist_CalDis == "rchisq"){
      selectInput(inputId = "tail_CalDis",
                  label = "Encuentra Área:",
                  choices = c("Cola superior"="upper"),
                  selected = "upper")
    }
    else
    {
      selectInput(inputId = "tail_CalDis",
                  label = "Encuentra Área:",
                  choices = c("Cola inferior"="lower", 
                              "Cola superior"="upper", 
                              "Ambas colas"="both",
                              "Intermedio"="middle"),
                  selected = "lower")
    }
  })

  output$lower_bound_CalDis = renderUI(
  {
    #print("lower bound")

    if (input$dist_CalDis == "rbinom")
    {
      if (is.null(input$tail_CalDis))
      {
        shiny:::flushReact()
        return()
      }

      if (input$tail_CalDis %in% c("both","middle"))
      {
        selectInput(inputId = "lower_bound_CalDis",
                    label = "Límite inferior:",
                    choices = c("<" = "open", 
                                "\u2264" = "closed"),
                    selected = "open")
      }
      else if (input$tail_CalDis == "lower")
      {
        selectInput(inputId = "lower_bound_CalDis",
                    label = "Límite:",
                    choices = c("<" = "open", 
                                "\u2264" = "closed"),
                    selected = "open")
      }
      else if (input$tail_CalDis == "upper")
      {
        selectInput(inputId = "lower_bound_CalDis",
                    label = "Límite:",
                    choices = c(">" = "open", 
                                "\u2265" = "closed"),
                    selected = "open")
      }
    }
  })

  output$upper_bound_CalDis = renderUI(
  {
    #print("upper bound")

    if (input$dist_CalDis == "rbinom")
    {
      if (is.null(input$tail_CalDis))
      {
        shiny:::flushReact()
        return()
      }

      if (input$tail_CalDis == "middle")
      {
        selectInput(inputId = "upper_bound_CalDis",
                    label = "Límite superior:",
                    choices = c("<" = "open", 
                                "\u2264" = "closed"),
                    selected = "open")
      }
      else if (input$tail_CalDis == "both")
      {
        selectInput(inputId = "upper_bound_CalDis",
                    label = "Límite superior:",
                    choices = c(">" = "open", 
                                "\u2265" = "closed"),
                    selected = "open")
      }
    }
  })

  get_model_text = reactive(
  {
    if (is.null(input$tail_CalDis)){
      shiny:::flushReact()
      return()
    }

    low_less = "<"
    low_greater = ">"

    up_less = "<"
    up_greater = ">"

    if (input$dist_CalDis == "rbinom" & input$tail_CalDis != "equal")
    {
      if (is.null(input$lower_bound_CalDis))
      {
        shiny:::flushReact()
        return()
      }

      if (input$lower_bound_CalDis == "closed")
      {
        low_less = "\u2264"
        low_greater = "\u2265"
      }

      if (input$tail_CalDis %in% c("middle","both"))
      { 
        if (is.null(input$upper_bound_CalDis)){
          shiny:::flushReact()
          return()
        }

        if (input$upper_bound_CalDis == "closed")
        {
          up_less = "\u2264"
          up_greater = "\u2265"
        }
      }
    }

    text = ""
    if (length(input$tail_CalDis) != 0)
    {
      if (input$tail_CalDis == "lower")
      {
        # P(X < a)
        text = paste0("P(X ", low_less, " a)")
      }
      else if (input$tail_CalDis == "upper")
      {
        # P(X > a)
        text = paste0("P(X ", low_greater, " a)")
      }
      else if (input$tail_CalDis == "middle")
      {
        # P(a < X < b)
        text = paste0("P(a ", low_less, " X ", up_less, " b)")
      }
      else if (input$tail_CalDis == "both")
      {
        # P(X < a o X > b)
        text = paste0("P(X ", low_less, " a o X ", up_greater, " b)")
      }
      else if (input$tail_CalDis == "equal")
      {
        # P(X = a)
        text = paste0("P(X = a)")
      }
    }

    return(text)
  })

  output$model_CalDis = renderText(
  {
    #print("model_CalDis")

    get_model_text()
  })

  ####
  # Normal distribution #
  ####

  output$mean_CalDis = renderUI(
  {
    #print("mean")
    if (input$dist_CalDis == "rnorm")
    {
      numericInput("mu",
                  "Media",
                  value = 0,
                  min = -1000,
                  max = 1000,
                  step=0.1)
    }
  })
    
  output$sd_CalDis = renderUI(
  {
    #print("sd_CalDis")
    if (input$dist_CalDis == "rnorm")
    {
      numericInput("sd_CalDis",
                  "Desviacion estandar",
                  value = 1,
                  min = 0.1,
                  max = 20,
                  step=0.01)
    }
  })
  
  ###
  # t, F, X^2 distribution #
  ###

  output$df1_CalDis = renderUI(
  {
    #print("df1_CalDis")
    if (input$dist_CalDis %in% c("rt","rchisq","rf"))
    {
      numericInput(ifelse(input$dist_CalDis %in% c("rt","rchisq"), "df","df1_CalDis"),
                  "Grados de libertad",
                  value = 10,
                  min = 1,
                  max = 500)
    }
  })
  
  output$df2_CalDis = renderUI(
  {
    #print("df2_CalDis")
    if (input$dist_CalDis == "rf")
    {
      numericInput("df2_CalDis",
                  "Grados de libertad (2)",
                  value = 10,
                  min = 1,
                  max = 500)
    }
  })


  ####
  # Binomial distribution #
  ###

  output$n_CalDis = renderUI(
  {
    #print("n_CalDis")
    if (input$dist_CalDis == "rbinom")
    {
      numericInput("n_CalDis",
                  "n",
                  value = 10,
                  min = 1,
                  max = 1000,
                  step = 1)
    }
  })

  output$p_CalDis = renderUI(
  {
    #print("p_CalDis")
    if (input$dist_CalDis == "rbinom")
    {
      sliderInput("p_CalDis",
                  "p",
                  value = 0.5,
                  min = 0,
                  max = 1,
                  step = .01)
    }
  })
  



  output$a_CalDis = renderUI(
  {
    #print("a_CalDis")

    value = 1
    min = 0
    max = 1
    step = 1

    if (input$dist_CalDis == "rnorm")
    {
      find_normal_step = function(sd_CalDis)
      {
        10^round(log(7*sd_CalDis/100,10))
      }

      if (is.null(input$mu) | is.null(input$sd_CalDis)){
        shiny:::flushReact()
        return()
      }

      mu = input$mu
      sd_CalDis = input$sd_CalDis
      if (is.null(mu)) mu = 0
      if (is.null(sd_CalDis)) sd_CalDis = 1

      value = mu - 1.96 * sd_CalDis
      min   = mu - 4 * sd_CalDis
      max   = mu + 4 * sd_CalDis
      step  = 0.01#find_normal_step(sd_CalDis)
      #if (mu == 0 & sd_CalDis == 1) {step = .01}
    }
    else if (input$dist_CalDis == "rt")
    {
      value = -1.96 
      min   = -6
      max   = 6
      step  = 0.01
    }
    else if (input$dist_CalDis == "rf")
    {
      value = round(qf(.95,as.numeric(input$df1_CalDis),as.numeric(input$df2_CalDis)),digits=2)
      min   = 0
      max   = round(qf(.995,as.numeric(input$df1_CalDis),as.numeric(input$df2_CalDis))*1.05,digits=2)
      step  = 0.01
    }
    else if (input$dist_CalDis == "rchisq")
    {
      value = round(qchisq(.95,as.numeric(input$df)),digits=2)
      min   = 0
      max   = round(qchisq(.995,as.numeric(input$df)),digits=2)
      step  = 0.01
    }
    else if (input$dist_CalDis == "rbinom")
    {
      if (is.null(input$n_CalDis)){
        shiny:::flushReact()
        return()
      }

      value = round(input$n_CalDis/4)
      min = 0
      max = input$n_CalDis
      step = 1
    }

    numericInput("a_CalDis", "a",
                value = value,
                min   = min,
                max   = max,
                step  = step)
  })

  output$b_CalDis = renderUI(
  {
    #print("b_CalDis")
     
    if (is.null(input$tail_CalDis))
    {
      shiny:::flushReact()
      return()
    }
    
    if (input$tail_CalDis %in% c("middle","both"))
    {
      value = 1
      min = 0
      max = 1
      step = 1

      if (input$dist_CalDis == "rnorm")
      {
        find_normal_step = function(sd_CalDis)
        {
          10^round(log(7*sd_CalDis/100,10))
        }

        if (is.null(input$mu) | is.null(input$sd_CalDis)){
          shiny:::flushReact()
          return()
        }

        mu = input$mu
        sd_CalDis = input$sd_CalDis
        if (is.null(mu)) mu = 0
        if (is.null(sd_CalDis)) sd_CalDis = 1

        value = mu + 1.96 * sd_CalDis
        min   = mu - 4 * sd_CalDis
        max   = mu + 4 * sd_CalDis
        step  = 0.01#find_normal_step(sd_CalDis)
      }
      else if (input$dist_CalDis == "rt")
      {
        value = 1.96 
        min   = -6
        max   = 6
        step  = 0.01
      }
      else if (input$dist_CalDis == "rbinom")
      {
        if (is.null(input$n_CalDis)){
          shiny:::flushReact()
          return()
        }

        value = round(input$n_CalDis*3/4)
        min = 0
        max = input$n_CalDis
        step = 1
      }

      numericInput("b_CalDis", "b",
                  value = value,
                  min   = min,
                  max   = max,
                  step  = step)
    }
  })  


  ##
  # Plotting #
  ##
  
  output$plot_CalDis = renderPlot(
  { 
    #print("plot_CalDis")

    if (is.null(input$tail_CalDis) | is.null(input$a_CalDis))
    {
      shiny:::flushReact()
      return()
    }

    L = NULL
    U = NULL

    error = FALSE

    if (input$tail_CalDis == "lower" | input$tail_CalDis == "equal")
    {
      L = input$a_CalDis 
    }
    else if (input$tail_CalDis == "upper")
    {
      U = input$a_CalDis 
    }
    else if (input$tail_CalDis %in% c("both","middle"))
    {
      if (is.null(input$b_CalDis)){
        shiny:::flushReact()
        return()
      }
      
      L = input$a_CalDis
      U = input$b_CalDis

      if (L > U)
        error = TRUE
    }

    if (error)
    {
      plot(0,0,type='n',axes=FALSE,xlab="",ylab="",mar=c(1,1,1,1))
      text(0,0,"Error: Límite inferior mayor que el límite superior.",col="red",cex=2)
    }
    else
    {
      if (input$dist_CalDis == "rnorm" | input$dist_CalDis == "rt") 
      {
        M = NULL
        if (input$tail_CalDis == "middle")
        {
          M = c(L,U)
          L = NULL
          U = NULL
        }

        if(input$dist_CalDis == "rnorm")
        {
          if(is.null(input$mu) | is.null(input$sd_CalDis))
          {
            shiny:::flushReact()
            return()
          }

          normTail(m=input$mu, s=input$sd_CalDis, L=L, U=U, M=M, axes=3, cex.axis=1.5)
          title(main="Distribución Normal")
        }
        else if (input$dist_CalDis == "rt")
        {
          if(is.null(input$df))
          {
            shiny:::flushReact()
            return()
          }

          normTail(m=0, s=1, df=input$df, L=L, U=U, M=M, axes=3, cex.axis=1.5)
          title(main="Distribución t")
        }
      }
        else if (input$dist_CalDis == "rchisq")
        {
          if(is.null(input$df))
          {
            shiny:::flushReact()
            return()
          }
          M = NULL
          if (input$tail_CalDis == "middle")
          {
            M = c(L,U)
            L = NULL
            U = NULL
          }
          
          chiTail(U=U, df=input$df, xlim = c(0,round(qchisq(.995,input$df),digits=2)+1))
          title(main="Distribución Chi^2")
        }
        else if (input$dist_CalDis == "rf")
        {        
          if(is.null(input$df1_CalDis) | is.null(input$df2_CalDis))
          {
            shiny:::flushReact()
            return()
          }
        
          M = NULL
          if (input$tail_CalDis == "middle")
          {
            M = c(L,U)
            L = NULL
            U = NULL
          }
                   
          FTail(U=U,df_n=input$df1_CalDis, df_d=input$df2_CalDis)
          title(main="Distribución F")
        }
      else if (input$dist_CalDis == "rbinom")
      {
        if(  is.null(input$n_CalDis)
           | is.null(input$p_CalDis)
           | is.null(input$lower_bound_CalDis))
        {
          shiny:::flushReact()
          return()
        }

        if(input$tail_CalDis %in% c("both","middle") & is.null(input$upper_bound_CalDis))
        {
          shiny:::flushReact()
          return()
        }

        d = dbinom(0:input$n_CalDis,input$n_CalDis,input$p_CalDis)

        plot(0,0,type='n',xlim=c(-0.5,input$n_CalDis+0.5),ylim=c(0,max(d)),
             xlab="",ylab="", axes=FALSE)
        axis(1, cex.axis=1.5)
        axis(2, cex.axis=1.5)
        title(main=paste("Distribución Binomial"))



        for (k in 1:length(d)) 
        {
            col = NA

            if (input$tail_CalDis == "lower")
            {
              if (input$lower_bound_CalDis == "open"   & k-1 <  L) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & k-1 <= L) col = "#569BBD"
            }
            else if (input$tail_CalDis == "upper")
            {
              if (input$lower_bound_CalDis == "open"   & k-1 >  U) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & k-1 >= U) col = "#569BBD"
            }
            else if (input$tail_CalDis == "equal")
            {
              if (k-1 == L) col = "#569BBD"
            }
            else if (input$tail_CalDis == "both")
            {
              if (input$lower_bound_CalDis == "open"   & input$upper_bound_CalDis == "open"   & (k-1 <  L | k-1 >  U)) col = "#569BBD"
              if (input$lower_bound_CalDis == "open"   & input$upper_bound_CalDis == "closed" & (k-1 <  L | k-1 >= U)) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & input$upper_bound_CalDis == "open"   & (k-1 <= L | k-1 >  U)) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & input$upper_bound_CalDis == "closed" & (k-1 <= L | k-1 >= U)) col = "#569BBD"
            }
            else if (input$tail_CalDis == "middle")
            {
              if (input$lower_bound_CalDis == "open"   & input$upper_bound_CalDis == "open"   & k-1 >  L & k-1 <  U) col = "#569BBD"
              if (input$lower_bound_CalDis == "open"   & input$upper_bound_CalDis == "closed" & k-1 >  L & k-1 <= U) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & input$upper_bound_CalDis == "open"   & k-1 >= L & k-1 <  U) col = "#569BBD"
              if (input$lower_bound_CalDis == "closed" & input$upper_bound_CalDis == "closed" & k-1 >= L & k-1 <= U) col = "#569BBD"
            }

            p = matrix(c(-1.5+k,0, -0.5+k,0, -0.5+k,d[k], -1.5+k,d[k], -1.5+k,0),ncol=2,byrow=TRUE)
          
            polygon(p, col=col)
        }
      }
    }
  })

  ###
  # Calculations #
  ###

  Valor=reactive({
    #if (input$dist_CalDis == "rnorm")
    #{
    #  Val<-as.numeric(Valor_Fin(input$dist_CalDis, input$tail_CalDis, input$mu, input$sd_CalDis, input$a_CalDis, input$b_CalDis))
    #}
    #else if (input$dist_CalDis == "rbinom"){
    #  Val<-as.numeric(Valor_Fin(input$dist_CalDis, input$tail_CalDis, n_CalDis=input$n_CalDis, p_CalDis=input$p_CalDis))
    #}
    #else{
      Val<-as.numeric(Valor_Fin(input$dist_CalDis, input$tail_CalDis, input$mu, input$sd_CalDis, input$a_CalDis, input$b_CalDis, 
                                      input$df, input$df1_CalDis, input$df2_CalDis, input$n_CalDis, input$p_CalDis, input$lower_bound_CalDis, 
                                      input$upper_bound_CalDis))
    #  }
    
    return(Val)
  })
  
#  Error_Peso= reactive({
#    EP<-abs(Valor()-Valor_Cues)
#    return(EP)
#  })
  
  output$area_CalDis = renderText({
    #if (is.null(Valor()) | is.null(get_model_text())) return(NULL)
    #else{
      text = paste(get_model_text(),"=",signif(Valor(),3))
      #as.numeric(Valor_Final())
      if (is.null(input$a_CalDis) | is.null(text)) return(paste(signif(Valor(),3)))
      else {
      text = sub("a",input$a_CalDis,text)
        #if (is.null(input$b_CalDis) | is.null(text)) return(NULL)
        if (input$tail_CalDis %in% c("both","middle")){
          if (is.null(input$b_CalDis) | is.null(text)) return(NULL)
          text = sub("b",input$b_CalDis,text)
        }
        return(text)
        }
      
    #}
    
  })

  #output$Error_CalDis = renderText({
  #    if (is.null(Valor()) | is.null(Valor_Cues)) return(NULL)
  #    else{
  #      text1 = paste("Tu error es = ",signif(abs(Valor()-Valor_Cues),3))
  #      #as.numeric(Valor_Final())
  #      errr<-abs(Valor()-Valor_Cues)
  #      if (length(errr)==0) {return(NULL)}
  #      else if (errr<0.005) {text1 = "Felicidades"}
  #      
  #      
  #      #text = sub("a",input$a_CalDis,text1)
  #      #if (input$tail_CalDis %in% c("both","middle")) 
  #      #  text = sub("b",input$b_CalDis,text1)
  #      
  #      return(text1)
  #    }
  #    
  # })
  
  
  
 
  observe({
    #this observer monitors when input$newplot is invalidated
    #or when input$difficulty is invalidated
    #and generates a new plot
    if (input$sal == 1) stopApp()
    input$newdat
    if (input$Ejem_Dis == "Peso_Est"){
      output$Doc_Peso <- source('./Problemas/Peso_Estudiantes.R',local=T,encoding="UTF-8")$value
      output$Ayuda <- source('./Problemas/R_Peso_Estudiantes.R',local=T,encoding="UTF-8")$value
    }
    else if (input$Ejem_Dis == "Tiro_Arc"){
      output$Doc_Tiro = source('./Problemas/torneo_de_tiro.R',local=T,encoding="UTF-8")$value
      output$Ayuda = source('./Problemas/R_torneo_de_tiro.R',local=T,encoding="UTF-8")$value
    }
    else if (input$Ejem_Dis == "Temp_Est"){
      output$Doc_Temp = source('./Problemas/Temperatura.R',local=T,encoding="UTF-8")$value
      output$Ayuda <- source('./Problemas/R_Temperatura.R',local=T,encoding="UTF-8")$value
    }
    #display text
    output$status1 <- renderText({"Escribe tu respuesta y haz clic en 'Enviar'"})
    output$status2 <- renderText({""})
    output$status3 <- renderText({""})
    
    #reset answered status
    answered<<-FALSE

    
  })
  
  
  observe({
    
    input$newEje
    Ejemplos<-c("Peso_Est","Tiro_Arc","Temp_Est")
    
    isolate({
      EjeMenos<-Ejemplos[Ejemplos!=input$Ejem_Dis]
      Tomar_ejem<-sample(EjeMenos, 1)
      updateRadioButtons(session, "Ejem_Dis", selected = Tomar_ejem)
      #Valor()<-NULL
      })
    #Valor_Cues<-NULL
    #errr<-NULL
    #display text
    output$status1 <- renderText({"Escribe tu respuesta y haz clic en 'Enviar'"})
    output$status2 <- renderText({""})
    output$status3 <- renderText({""})
    
    #reset answered status
    answered<<-FALSE
  })
  
  observe({
    if (is.null(Valor())) return(NULL)
    else{
  updateNumericInput(session, "Res_Cuest", value = round(Valor(),4))}
    
  })
  
  observe({
    
    input$Resp_Ejem
    isolate({
      errr<-abs(Valor_Cues-input$Res_Cuest)
    })
    
    if(errr<0.01){
      output$status1 <- renderText({""})
      output$status2 <- renderText({paste(generateResponse(1))})
      output$status3 <- renderText({""})
      if(!answered){
        score <<- score+2
        output$score <- renderText({sprintf("Tu puntuación : %d",score)})
        answered <<- TRUE
      }
    }
    else if(errr<0.05){
      output$status1 <- renderText({""})
      output$status2 <- renderText({""})
      output$status3 <- renderText({generateResponse(2)})
      if(!answered){
        score <<- score+1
        output$score <- renderText({sprintf("Tu puntuación : %d",score)})
        answered <<- TRUE
      }
    }
    else{
      output$status1 <- renderText({""})
      output$status2 <- renderText({""})
      output$status3 <- renderText({paste(generateResponse(3), sprintf("(La respuesta correcta es: %.4f)",Valor_Cues))})
      #score <<- score-1
      answered <<- TRUE
    }
    
    
    #if (is.null(Valor()) | is.null(input$Res_Cuest)) return(NULL)
    #if (input$Ejem_Dis == "Peso_Est"){
    #  output$R_Doc_Peso = source('./Problemas/R_Peso_Estudiantes.R',local=T,encoding="UTF-8")$value
    #}
    #else if (input$Ejem_Dis == "Tiro_Arc"){
    #  output$R_Doc_Tiro = source('./Problemas/R_torneo_de_tiro.R',local=T,encoding="UTF-8")$value
    #}
    #else if (input$Ejem_Dis == "Temp_Est"){
    #  output$R_Doc_Temp = source('./Problemas/R_Temperatura.R',local=T,encoding="UTF-8")$value
    #}
    #Valor_Cues<-NULL
    
  })
  
  #Resp <- 0.5

  
  
##Reporte
#regFormula <- reactive({
#  as.formula(paste('mpg ~', input$x))
#})

#output$report <- renderUI({
#  src <- normalizePath('report.Rmd')
#  
  # temporarily switch to the temp dir, in case you do not have write
  # permission to the current working directory
#  owd <- setwd(tempdir())
#  on.exit(setwd(owd))
  #library(knitr)
#  opts_knit$set(root.dir = owd)#knitr::opts_knit$set(root.dir = owd)
  #options(encoding = 'UTF-8')
#  tagList(
#    HTML(knit2html(text = readLines(src), fragment.only = TRUE)),
    #HTML(knitr::knit2html(text = readLines(src), fragment.only = TRUE)),
    # typeset LaTeX math
#    tags$script(HTML('MathJax.Hub.Queue(["Typeset", MathJax.Hub]);')),
    # syntax highlighting
#    tags$script(HTML("if (hljs) $('#report pre code').each(function(i, e) {
#                     hljs.highlightBlock(e)
#});"))
#    )
#})


})
