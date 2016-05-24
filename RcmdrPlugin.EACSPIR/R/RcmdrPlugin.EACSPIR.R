# Codigo para el desarrollo del plugin de R-Commander para el manual EACSPIR #

.onAttach <- function(libname, pkgname){
  if (!interactive()) return()
  Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins
  if (!pkgname %in% plugins) {
    Rcmdr$plugins <- c(plugins, pkgname)
    options(Rcmdr=Rcmdr)
    if("package:Rcmdr" %in% search()) {
      if(!getRcmdr("autoRestart")) {
        closeCommander(ask=FALSE, ask.save=TRUE)
        Commander()
      }
    }
    else {
      Commander()
    }
  }
}

# Pendiente de introducir funcion findGlobals() desarrollada por J. Fox #
if (getRversion() >= '2.15.1') globalVariables(c('dataTab', 'optionsTab', 'tkselect', 'notebook',
                                                 'echocodigoVariable', 'creahtmlVariable', 'alternativaVariable', 'pruebaVariable', 'UMW', 'statisticsTab',
                                                 'alternativaFrame', 'pruebaFrame', 'opsFrame', 'descripVariable', 'graficoVariable', 'SCuadradosVariable',
                                                 '.prueba.avarMR', '.baseDatosActiva', 'SCuadradosFrame', 'statsFrame', 'tablaVariable', 'subsetVariable',
                                                 'subsetFrame', 'tablaFrame', 'statistics2Tab', 'porcentajesVariable', 'frecEspVariable',
                                                 'jicuadradoVariable', 'jiComponentesVariable', 'phiPearsonVariable', 'contingPearsonVariable', 'sakodaVariable',
                                                 'chuprovVariable', 'VCramerVariable', 'yuleVariable', 'lambdaVariable', 'tauVariable', 'theilVariable',
                                                 '.indices.bc', 'porcentajesFrame', 'esperadasFrame', 'errorpredFrame', 'CovVariable', 'pearsonVariable',
                                                 'spearmanVariable', 'kendallVariable', 'determVariable', '.indices.bn', 'cond1', 'cond2', 'cond3', 'cond4',
                                                 'gammaVariable', 'sommersVariable', 'wilsonVariable', '.indices.bo', 'dtm', 'identificaVariable', 'cajaFrame',
                                                 'cond', 'variable', 'escalaVariable', 'normalsupVariable', 'escalaFrame', 'histogramaFrame', 'oddsVariable',
                                                 'modaVariable', 'RVVariable', 'blauVariable', 'IVQVariable', 'teachmanVariable', '.indices.cat', 'numcatFrame',
                                                 'pruebaZ.apuntamiento', 'pruebaZ.forma', 'TW', 'adVariable', 'ksVariable', 'shapiroVariable',
                                                 'graficosVariable', '.groups', 'tab', '.norm.test', 'groupsFrame', '.TablaResumen', 'selectodasVariable',
                                                 'mediaVariable', 'medianaVariable', 'mediageomVariable', 'trimediaVariable', 'promcuarVariable', 'midRVariable',
                                                 'medrecVariable', 'varianciaVariable', 'dtVariable', 'CVVariable', 'dtgeomVariable', 'desvmedVariable',
                                                 'rangoVariable', 'IQRVariable', 'desvcuarVariable', 'madVariable', 'CVRVariable', 'ACentVariable', 'minVariable',
                                                 'maxVariable', 'Q1Variable', 'Q2Variable', 'Q3Variable', 'percentVariable', 'H1Variable', 'H3Variable',
                                                 'K2Variable', 'K3Variable', 'beta1Variable', 'gamma1Variable', 'beta2Variable', 'gamma2Variable', '.indices.num',
                                                 'tcFrame', 'tcRFrame', 'dispFrame', 'dispRFrame', 'statistics3Tab', 'posicFrame', 'posicRFrame',
                                                 'statistics4Tab', 'formaFrame', 'formaRFrame', '.indices.ord', 'ywttest'))

# Resumen descriptivo univariante para variables categoricas: indicadores y graficos #

catResumen <- function(x){
  ni <- table(x,useNA='always')
  fi <- ni/sum(ni)
  Ni <- cumsum(ni)
  Fi <- cumsum(fi)
  res <- as.data.frame(matrix(round(cbind(ni,Ni,fi,Fi),2),nrow=nrow(ni), byrow=FALSE,
                              dimnames=list(c(levels(x),'NAs'),c('ni','Ni','fi','Fi'))))
  res
}

resumen.categoricas <- function(){
  defecto <- list(x.inicial=NULL,echo.inicial="0",creahtml.inicial="0",
                  tab.inicial=0)
  dialog.valores <- getDialog("resumen.categoricas",defecto)
  initializeDialog(title=gettextRcmdr("Resumen distribucion de frecuencias"),use.tabs=TRUE)
  listaVar <- variableListBox(dataTab, Factors(), selectmode="multiple",
                              title=gettextRcmdr("Variables (escoja una o mas)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"factor"))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function()
  {
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    .BaseDatosActiva <- ActiveDataSet()
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("resumen.categoricas",list(x.inicial=x,echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))
    if (length(x) == 0){
      errorCondition(recall=resumen.categoricas, 
                     message=gettextRcmdr("Debe escoger una variable."))
      return()
    }    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html") 
    }
    for (variable in x)
    {
      instruccion1 <- paste(".TablaResumen <- catResumen(", .BaseDatosActiva, "$", variable, ")", sep="")
      if (echocodigo == 1)
      {
        logger(instruccion1)
      }
      justDoIt(instruccion1)
      doItAndPrint(paste(".TablaResumen # Resumen distribucion de frecuencias para", variable))
      if (creahtml == 1)
      {
        titulo <- paste("Distribucion frecuencias para variable ", 
                        variable,sep="")
        HTML(as.title(titulo),file=.archivo)
        HTML(.TablaResumen, file=.archivo)
        HTMLhr(file = .archivo)                   
      }
    }
    closeDialog()
    if (echocodigo == 1) logger("remove('.TablaResumen')") 
    remove(.TablaResumen, envir=.GlobalEnv)  
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="table",reset="resumen.categoricas",apply="resumen.categoricas")
  tkgrid(getFrame(listaVar), sticky="nw")
  tkgrid(opsFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tab.names=c("Datos","Opciones"))
}

catindices <- function(vars,statistics,kVariable=NULL){
  vars <- as.data.frame(vars)
  variables <- names(vars)
  .ni <- lapply(vars,table)
  odds <- function(x){
    res <- round(x/(sum(x,na.rm=TRUE)-x), 2)
    res
  }  
  moda <- function(vars){
    .ni <- table(vars)
    if (is.numeric(vars)) res <- as.numeric(names(.ni[which(.ni==max(.ni,na.rm=TRUE))]))
    else res <- names(.ni)[which(.ni==max(.ni))]
    res
  }
  RV <- function(x){
    res <- round(1 - max(x,na.rm=TRUE)/sum(x,na.rm=TRUE),2)
    res
  }
  blau <- function(x){
    res <- round(1-sum(prop.table(x)^2),2)
    res
  }
  IVQ <- function(x,kval){
    .ni <- lapply(x,table)
    res <- array()
    for (i in 1:length(x)){
      resi <- round((1-sum(prop.table(.ni[[i]])^2))/((kval[i]-1)/kval[i]),2)
      res[i] <- resi
    }
    res  
  }
  teachman <- function(x){
    res <- round(-sum(prop.table(x)[prop.table(x)!=0]*log(prop.table(x)[prop.table(x)!=0])),2)
    res
  }
  res <- list()
  res[[1]] <- variables
  res[[2]] <- statistics
  i <- 2 
  if ("Odds" %in% statistics) {
    .odds <- lapply(.ni,odds)
    i <- i + 1
    res[[i]] <- .odds
  }
  if ("Moda" %in% statistics) {
    .moda <- lapply(vars,moda)
    i <- i + 1
    res[[i]] <- .moda 
  }
  numfilas <- length(vars)
  numcolumnas <- sum(c("RV","Blau","IVQ","Teachman") %in% statistics)
  .tablaRes <- as.data.frame(matrix(nrow=numfilas,ncol=numcolumnas))   
  j <- 0
  if ("RV" %in% statistics) {
    .RV <- unlist(lapply(.ni,RV),use.names=T)
    j <- j + 1
    .tablaRes[,j]<-.RV
  }
  if ("Blau" %in% statistics) {
    .blau <- unlist(lapply(.ni,blau),use.names=T)
    j <- j + 1
    .tablaRes[,j]<-.blau
  }
  if ("IVQ" %in% statistics) {
    if (length(kVariable)==1) {
      if (kVariable == "<auto>") kVariable <- unlist(lapply(.ni,length),use.names=F)}
    if (length(kVariable) > 1) { 
      for (k in 1:length(vars)){
        if (kVariable[k] < length(.ni[[k]])) kVariable[k] <- length(.ni[[k]])
        if (is.na(kVariable[k]) | kVariable[k] < 1) kVariable[k] <- length(.ni[[k]])
        if (kVariable[k] == "<auto>") kVariable[k] <- length(.ni[[i]])
      }
    }
    .IVQ <- IVQ(vars,kVariable)
    j <- j + 1
    .tablaRes[,j]<-.IVQ
  }
  if ("Teachman" %in% statistics) {
    .teachman <- unlist(lapply(.ni,teachman),use.names=T) 
    j <- j + 1
    .tablaRes[j]<-.teachman
  }
  if (numcolumnas > 0){
    i <- i + 1
    res[[i]] <- .tablaRes
  }
  class(res) <- "catindices"
  res
}

print.catindices <- function(x,...){
  j<-3
  if ("Odds" %in% x[[2]]) {
    for (i in 1:length(x[[1]])) {cat(paste("# Odds para ", x[[1]][i],": ",sep= ''),"\n")
                                 names(x[[j]][i]) <- NULL
                                 print(x[[j]][i])
                                 cat("\n\n")
    }
    j <- j + 1
  }
  if ("Moda" %in% x[[2]]) {
    for (i in 1:length(x[[1]])) {cat(paste("# Moda para ", x[[1]][i],": ",sep= ''),"\n\n")
                                 names(x[[j]][[i]]) <- NULL
                                 print(x[[j]][[i]])
                                 cat("\n\n")
    }
    j <- j + 1
  }  
  if (sum(c("RV","Blau","IVQ","Teachman") %in% x[[2]]) > 0){
    names(x[[j]]) <- NULL
    .tablaResultados <- as.data.frame(x[[j]])
    rownames(.tablaResultados) <- x[[1]]
    colnames(.tablaResultados) <- c("RV","Blau","IVQ","Teachman")[c("RV","Blau","IVQ","Teachman") %in% x[[2]]]
    cat("# Otros indices descriptivos:","\n\n")
    print(.tablaResultados)  
    cat("\n\n")
  }
  invisible(x)  
}
    
indices.categoricas <- function(){
    defecto <- list(x.inicial=NULL,odds.inicial="0",moda.inicial="0",RV.inicial="0",
                    blau.inicial="0",IVQ.inicial="0",teachman.inicial="0",k.inicial="<auto>",
                    echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
    dialog.valores <- getDialog("indices.categoricas",defecto) 
    initializeDialog(title=gettextRcmdr("Otros indices para variables categoricas"),use.tabs=TRUE,
                                        tabs=c('dataTab','statisticsTab','optionsTab'))
    listaVar <- variableListBox(dataTab, Factors(), selectmode="multiple",
                                title=gettextRcmdr("Variables (escoja una o mas)"),
                                initialSelection=varPosn(dialog.valores$x.inicial,"factor"))
    onOK <- function(){
        tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1      
        x <- getSelection(listaVar)
        if (length(x) == 0){
            errorCondition(recall=indices.categoricas, message=gettextRcmdr("Debes escoger una variable."))
            return()
            }
        .BaseDatosActiva <- ActiveDataSet()
        echocodigo <- tclvalue(echocodigoVariable)
        oddsval <- tclvalue(oddsVariable)
        modaval <- tclvalue(modaVariable)
        RVval <- tclvalue(RVVariable)
        blauval <- tclvalue(blauVariable)
        IVQval <- tclvalue(IVQVariable)
        teachmanval <- tclvalue(teachmanVariable)
        creahtml <- tclvalue(creahtmlVariable)
        selec <- as.numeric(modaval) + as.numeric(oddsval) + as.numeric(RVval) +
          as.numeric(blauval) + as.numeric(IVQval) + 
          as.numeric(teachmanval)
        if (selec == 0){
          errorCondition(recall=indices.categoricas, 
                         message=gettextRcmdr("Debe escoger algun indicador."))
          return()
        }
        putDialog("indices.categoricas",list(x.inicial=x,odds.inicial=oddsval,moda.inicial=modaval,
                                             RV.inicial=RVval,blau.inicial=blauval,IVQ.inicial=IVQval,
                                             k.inicial=tclvalue(kVariable),teachman.inicial=teachmanval,
                                             echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))        
        if (IVQval == 1)
        {
          kval <- vector()
          opts <- options(warn=-1)
          options(opts)
          if (tclvalue(kVariable) != gettextRcmdr("<auto>")) {
            k <- c(gsub(" ", ",",gsub(", ", ",",tclvalue(kVariable))))
            k <- unlist(strsplit(k,","))
            if (length(k) < length(x))
            {
              tclvalue(kVariable) <- gettextRcmdr("<auto>")
              kval <- as.character(tclvalue(kVariable))
              Message(message=gettextRcmdr("Vector de categorias, invalido se utilizara '<auto>'."),
                      type="warning")
            }
            else {
              j <- 0
              for (variable in x){
                j <- j + 1
                kval[j] <- as.numeric(k[j]) 
                if (is.na(kval[j]) || kval[j] < 1)
                {
                  errorCondition(recall=indices.categoricas,
                                 message=gettextRcmdr
                                 ("El numero de categorias k debe ser un numero positivo"))
                  return()
                }  
              }
              kval <- paste("c(",paste(kval,collapse=",",sep=''),")",sep='')
            }
          } else kval <- tclvalue(kVariable)
        }

        vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
        else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
        if (length(x) == 1) variables <- paste(.BaseDatosActiva, "[", vars, "]", sep="")
        else variables <- paste(.BaseDatosActiva, "[,", vars, "]", sep="")
        stats <- paste("c(",
                       paste(c('"Odds"', '"Moda"', '"RV"', '"Blau"', '"IVQ"', '"Teachman"')
                             [c(oddsval, modaval, RVval, blauval, IVQval , teachmanval) == 1], 
                             collapse=", "), ")", sep="")
        if (IVQval == 1) {
          if (kval!="<auto>")
            instruccion1 <- paste(".indices.cat <- catindices(vars=",variables, ", statistics=", stats,", kVariable=", kval, ")", sep="")
          else
            instruccion1 <- paste(".indices.cat <- catindices(vars=",variables, ", statistics=", stats,", kVariable='", kval, "')", sep="")
        }
          else instruccion1 <- paste(".indices.cat <- catindices(vars=",variables, ", statistics=", stats,")", sep="")
        justDoIt(instruccion1)
        if (echocodigo == 1)
        {
          logger(instruccion1)
        }
          doItAndPrint(".indices.cat  # Indices para variables categoricas ")          
        if (creahtml == 1)
        {
          if (!file.exists("Informe de Resultados.html"))
            .archivo <- HTMLInitFile(file.path(getwd()),
                                     "Informe de Resultados", BackGroundColor="#FFFFCC")
          else
            .archivo <- file.path(getwd(), "Informe de Resultados.html")
            titulo <- "Descriptivos para variables categoricas"
            HTML(as.title(titulo),file=.archivo)
          i <- 1
          while (i <= length(x)){
            j <- 3
            if (oddsval ==1){
              HTML(paste("Odds para:",.indices.cat[[1]][i]), file=.archivo)
              .odds <- as.data.frame(.indices.cat[[j]][i])[2]
              colnames(.odds) <- 'Odds'
              rownames(.odds) <- unlist(as.data.frame(.indices.cat[[j]][i])[1],use.names=F)             
              HTML(.odds, file=.archivo)
              j <- j + 1
            }
            if (modaval ==1){
              HTML(paste("Moda para:",.indices.cat[[1]][i]), file=.archivo)
              .moda <- as.data.frame(.indices.cat[[j]][i])[1]
              colnames(.moda) <- 'Moda'   
              rownames(.moda) <- NULL 
              HTML(.moda, file=.archivo)
              j <- j + 1
            }
            i <- i +1
          }
          selec <- sum(as.numeric(RVval) + as.numeric(blauval) + as.numeric(IVQval) + as.numeric(teachmanval))
          if ( selec>0 ){
            .TablaRes <- as.data.frame(.indices.cat[[5]])
            rownames(.TablaRes) <- .indices.cat[[1]]
            colnames(.TablaRes) <- c("RV","Blau","IVQ","Teachman")[c("RV","Blau","IVQ","Teachman") %in% .indices.cat[[2]]]
            HTML("Otros indicadores: ", file=.archivo)
            HTML(.TablaRes, file=.archivo)
            HTMLhr(file = .archivo) 
          }
        }          

        remove(.indices.cat, envir=.GlobalEnv)
        closeDialog()
        if (echocodigo == 1) logger("remove(.indices.cat)")
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="RcmdrPlugin.EACSPIR",reset="indices.categoricas",apply="indices.categoricas")
    checkBoxes(statisticsTab,frame="numcatFrame",boxes=c("odds","moda","RV","blau","IVQ","teachman"),
               initialValues=c(dialog.valores$odds.inicial,dialog.valores$moda.inicial,
                               dialog.valores$RV.inicial,dialog.valores$blau.inicial,
                               dialog.valores$IVQ.inicial,dialog.valores$teachman.inicial),
               labels=gettextRcmdr(c("Odds ","Moda ","Razon de Variacion ","Indice de Diversidad de Blau ",
                                     "Indice de Variacion Cualitativa","Indice de Diversidad de Teachman ")), 
               title = gettextRcmdr("Indices"))    
    checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
               initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
               labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
               title = gettextRcmdr("Opciones"))
    rightFrame <- tkframe(statisticsTab)
    kVariable <- tclVar(dialog.valores$k.inicial)
    kField <- ttkentry(rightFrame, width="8", textvariable=kVariable)
    tkgrid(getFrame(listaVar), sticky="nw")
    tkgrid(labelRcmdr(rightFrame,text="  "))
    tkgrid(labelRcmdr(rightFrame,text="  "))
    tkgrid(labelRcmdr(rightFrame,text="  "))
    tkgrid(labelRcmdr(rightFrame,text="  ")) 
    tkgrid(labelRcmdr(rightFrame,text="  "))  
    tkgrid(labelRcmdr(rightFrame,text=gettextRcmdr("Categorias: k ="),fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
           kField,sticky="w")
    tkgrid(numcatFrame, labelRcmdr(statisticsTab,text="  "), rightFrame, sticky="nw")
    tkgrid(opsFrame, sticky="w")
    dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
                 tab.names=c("Datos","Estadisticos","Opciones"))
    }

pareto <- function(x,tablaVar){
  x <- as.data.frame(x)
  variable <- names(x)
  if (tablaVar == "ni") .tabla <- table(x)
  else .tabla <- table(x)/sum(table(x))
  .tabla <- .tabla[order(-.tabla)]
  titulo <- paste('Diagrama de Pareto para ', variable, sep="")
  tituloy <- if (tablaVar == "ni") "Frecuencias Absolutas"
  else "Frecuencias Relativas"
  par(mar=c(5,4,4,4))
  .x <- barplot(.tabla, main=titulo,ylab=tituloy,
                ylim=c(0,sum(.tabla)*1.05),col=heat.colors(length(.tabla)))
  lines(.x[1:length(.tabla)],cumsum(.tabla),type='b')
  box()
  axis(4,at=seq(0,max(cumsum(.tabla)),length=5),
       labels=paste(seq(0,1,length=5)*100,'%',sep=''))
  mtext('Porcentaje Acumulado', 4, line=2.5, las=3)
}

grafico.Pareto <- function(){
  defecto <- list(x.inicial=NULL,tabla.inicial="ni",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("grafico.Pareto",defecto) 
  initializeDialog(title=gettextRcmdr("Diagrama de Pareto"),use.tabs=TRUE)
  listaVar <- variableListBox(dataTab, Factors(),
                              title=gettextRcmdr("Variables (escoja una)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"factor"))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
    onOK <- function()
    {
      tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
      x <- getSelection(listaVar)
      .BaseDatosActiva <- ActiveDataSet()
      echocodigo <- tclvalue(echocodigoVariable)
      creahtml <- tclvalue(creahtmlVariable)
      tabVariable <- as.character(tclvalue(tablaVariable))
      putDialog("grafico.Pareto",list(x.inicial=x,tabla.inicial=tabVariable,
                                      echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))
      if (length(x) == 0){
        errorCondition(recall=grafico.Pareto, 
                       message=gettextRcmdr("Debe escoger una variable."))
        return()
      }
      vars <- paste('"', x, '"', sep="") 
      variables <- paste(.BaseDatosActiva, "[,", vars, "]", sep="")      
      instruccion1 <- paste("pareto(x=", .BaseDatosActiva, "[", vars,"], tablaVar='",tabVariable, "')", sep="")
      if (echocodigo == 1)
      {
        logger(instruccion1)
      }
      justDoIt(instruccion1)
      if (creahtml == 1)
      {
        if (!file.exists("Informe de Resultados.html"))
          .archivo <- HTMLInitFile(file.path(getwd()),
                                   "Informe de Resultados", BackGroundColor="#FFFFCC")
        else
          .archivo <- file.path(getwd(), "Informe de Resultados.html")
        titulo <- paste("Diagrama de Pareto para variable ",variable,sep="")
        HTML(as.title(titulo),file=.archivo)
        nombre.archivo <- paste("ParetoR",gsub(":","",substr(Sys.time(),12,19)),
                                ".jpg",sep="")
        dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                  width=500, height=500)
        HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
        HTMLhr(file = .archivo)
      }
      closeDialog()        
      tkfocus(CommanderWindow())
    }
  OKCancelHelp(helpSubject="table",reset="grafico.Pareto",apply="grafico.Pareto")
  tkgrid(getFrame(listaVar), sticky="nw") 
  radioButtons(dataTab,name = "tabla", buttons = c("niButton","fiButton"), 
               values = c("ni","fi"),
               labels = gettextRcmdr(c("Frecuencias Absolutas", "Frecuencias Relativas")), 
               initialValue = dialog.valores$tabla.inicial, 
               title = gettextRcmdr("Tablas basadas en:"))
  tkgrid(tablaFrame, sticky="w") 
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tab.names=c("Datos","Opciones"))
}

# Resumen descriptivo univariante para variables ordinales: indicadores y graficos #

ordindices <- function(vars,statistics,rec=NULL,propdat=NULL,percentil=NULL){
  
  vars <- as.data.frame(vars)
  variables <- names(vars)
  mediana <- function(x){
    res <- median(as.numeric(x),na.rm=TRUE)
    res
  }  
  moda <- function(x){
    .ni <- table(x)
    if (is.numeric(x)) res <- as.numeric(names(.ni[which(.ni==max(.ni,na.rm=TRUE))]))
    else res <- names(.ni)[which(.ni==max(.ni))]  
    res
  }
  trimedia <- function(x){
    .Finf <- fivenum(as.numeric(x))[2]
    .Md <- fivenum(as.numeric(x))[3]
    .Fsup <- fivenum(as.numeric(x))[4]
    res <- round((.Finf+2*.Md+.Fsup)/4,2)  
    res
  }
  
  promcuar <- function(x){
    .Q1 <- quantile(as.numeric(x),na.rm=TRUE)[2]
    .Q3 <- quantile(as.numeric(x),na.rm=TRUE)[4]
    res <- round((.Q1+.Q3)/2,2)
    names(res) <- NULL
    res
  }
  midR <- function(x){
    .min <- quantile(as.numeric(x),na.rm=TRUE)[1]
    .max <- quantile(as.numeric(x),na.rm=TRUE)[5]
    res <- round((.min+.max)/2,2)
    names(res) <- NULL    
    res
  }
  medrec <- function(x,rec){
    res <- round(mean(as.numeric(x),trim=rec,na.rm=TRUE),2)
    res
  }
  rango <- function(x){
    .min <- quantile(as.numeric(x),na.rm=TRUE)[1]
    .max <- quantile(as.numeric(x),na.rm=TRUE)[5]
    res <- round(.max-.min,2)
    names(res) <- NULL        
    res
  }
  iqr <- function(x){
    .Q1 <- quantile(as.numeric(x),na.rm=TRUE)[2]
    .Q3 <- quantile(as.numeric(x),na.rm=TRUE)[4]
    res <- round(.Q3-.Q1,2)
    names(res) <- NULL        
    res
  }
  mad <- function(x){
    .mediana <- median(as.numeric(x),na.rm=TRUE)
    res <- round(median(abs(as.numeric(x)-.mediana),na.rm=TRUE),2)
    res
  }
  cvr <- function(x){
    .Finf <- fivenum(as.numeric(x))[2]
    .Fsup <- fivenum(as.numeric(x))[4]
    res <- round((.Fsup-.Finf)/(.Finf+.Fsup), 2)
    res
  }
  desvcuar <- function(x){
    .Q1 <- quantile(as.numeric(x),na.rm=TRUE)[2]
    .Q3 <- quantile(as.numeric(x),na.rm=TRUE)[4]
    res <- round((.Q3-.Q1)/2,2)
    names(res) <- NULL    
    res
  }
  acent <- function(x,propdat){
    .ACinf <- quantile(as.numeric(x),probs=(1-propdat)/2,na.rm=TRUE)
    .ACsup <- quantile(as.numeric(x),probs=1-(1-propdat)/2,na.rm=TRUE)
    res <- round(.ACsup-.ACinf,2)
    names(res) <- NULL       
    res
  }
  q1 <- function(x){
    res <- round(quantile(as.numeric(x),na.rm=TRUE)[2],2) 
    names(res) <- NULL        
    res
  }
  q2 <- function(x){
    res <- round(quantile(as.numeric(x),na.rm=TRUE)[3],2) 
    names(res) <- NULL        
    res
  }
  q3 <- function(x){
    res <- round(quantile(as.numeric(x),na.rm=TRUE)[4],2) 
    names(res) <- NULL        
    res
  }
  pct <- function(x,percentil){
    res <- round(quantile(as.numeric(x),probs=percentil,na.rm=TRUE),2)
    names(res) <- NULL       
    res
  }
  h1 <- function(x) {
    .Finf <- fivenum(as.numeric(x))[2]
    .Md <- fivenum(as.numeric(x))[3]
    .Fsup <- fivenum(as.numeric(x))[4]
    res <- round((.Finf+.Fsup-2*.Md)/(2*.Md),2)
    res
  }
  h3 <- function(x) {
    .AC90 <- quantile(as.numeric(x),probs=0.9,na.rm=TRUE)
    .AC10 <- quantile(as.numeric(x),probs=0.1,na.rm=TRUE)
    .Md <- fivenum(as.numeric(x))[3]
    res <- round((.AC90+.AC10-2*.Md)/(2*.Md),2)
    names(res) <- NULL        
    res
  }
  k2 <- function(x){
    .AC90 <- quantile(as.numeric(x),probs=0.9,na.rm=TRUE)
    .AC10 <- quantile(as.numeric(x),probs=0.1,na.rm=TRUE)
    .Q1 <- quantile(as.numeric(x),na.rm=TRUE)[2]
    .Q3 <- quantile(as.numeric(x),na.rm=TRUE)[4]
    res <- round((.AC90-.AC10)/(1.9*(.Q3-.Q1)),2)
    names(res) <- NULL        
    res
  }
  k3 <- function(x){
    .Einf <- quantile(as.numeric(x),probs=0.125,na.rm=TRUE)
    .Esup <- quantile(as.numeric(x),probs=0.875,na.rm=TRUE)
    .Finf <- fivenum(as.numeric(x))[2]
    .Fsup <- fivenum(as.numeric(x))[4]
    res <- round((.Esup-.Einf)/(1.7*(.Fsup-.Finf)),2)
    names(res) <- NULL        
    res
  }
  res <- list()
  res[[1]] <- variables
  res[[2]] <- statistics
  res[[3]] <- rec
  res[[4]] <- propdat
  res[[5]] <- percentil
  i <- 5
  
  numfilas <- length(vars)
  TCcol <- sum(c("Mediana","Moda","Trimedia","PromCuar","midR","MediaRec") %in% statistics)
  dispcol <- sum(c("Rango", "IQR", "DesvCuar", "MAD", "CVR", "ACent") %in% statistics)
  posiccol <- sum(c("Min", "Max", "Q1", "Q2", "Q3") %in% statistics) + ("Pct" %in% statistics)*length(percentil)
  formacol <- sum(c("H1", "H3", "K2", "K3") %in% statistics)
  .tablaTC <- as.data.frame(matrix(nrow=numfilas,ncol=TCcol))
  .tabladisp <- as.data.frame(matrix(nrow=numfilas,ncol=dispcol))
  .tablaposic <- as.data.frame(matrix(nrow=numfilas,ncol=posiccol))
  .tablaforma <- as.data.frame(matrix(nrow=numfilas,ncol=formacol))
  
  if (TCcol > 0){
    j <- 0
    if ("Mediana" %in% statistics) {
      .mediana <- unlist(lapply(vars,mediana),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.mediana
    }
    if ("Moda" %in% statistics) {
      .moda <- lapply(vars,moda)
      
      for (variable in length(.moda))
      {
        if (length(.moda[[variable]]) > 1) {
          #.moda[[variable]] <- .moda[[variable]][1]
          Message(message=gettextRcmdr(paste("Variable ",variables[variable]," tiene mas de una moda: ",
                                             "En documento HTML se muestra un solo valor.", sep="")),
                  type="warning")        
        }
      }
      if (length(.moda[which(lapply(.moda,length)>1,useNames=FALSE)])>0)
        .moda[which(lapply(.moda,length)>1,useNames=FALSE)]<-.moda[[which(lapply(.moda,length)>1,useNames=FALSE)]][1]
      .moda <- unlist(.moda,use.names=F)
      j <- j + 1
      .tablaTC[,j] <- .moda
    }
    if ("Trimedia" %in% statistics) {
      .trimedia <- unlist(lapply(vars,trimedia),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.trimedia
    }
    if ("PromCuar" %in% statistics) {
      .promcuar <- unlist(lapply(vars,promcuar),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.promcuar
    }
    if ("midR" %in% statistics) {
      .midr <- unlist(lapply(vars,midR),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.midr
    }
    if ("MediaRec" %in% statistics) {
      .mediarec <- unlist(lapply(vars,medrec,rec),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.mediarec
    }    
    i <- i + 1
    res[[i]] <- .tablaTC
  }
  
  if (dispcol > 0){
    j <- 0
    if ("Rango" %in% statistics) {
      .rango <- unlist(lapply(vars,rango),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.rango
    }
    if ("IQR" %in% statistics) {
      .iqr <- unlist(lapply(vars,iqr),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.iqr
    }
    if ("DesvCuar" %in% statistics) {
      .dq <- unlist(lapply(vars,desvcuar),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.dq
    }
    if ("MAD" %in% statistics) {
      .mad <- unlist(lapply(vars,mad),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.mad
    }
    if ("CVR" %in% statistics) {
      .cvr <- unlist(lapply(vars,cvr),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.cvr
    }
    if ("ACent" %in% statistics) {
      .acent <- unlist(lapply(vars,acent,propdat),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.acent
    }    
    i <- i + 1
    res[[i]] <- .tabladisp
  }
  
  if (posiccol > 0){
    j <- 0
    if ("Min" %in% statistics) {
      .min <- unlist(lapply(lapply(vars,as.numeric),min,na.rm=TRUE),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.min
    }
    if ("Max" %in% statistics) {
      .max <- unlist(lapply(lapply(vars,as.numeric),max,na.rm=TRUE),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.max
    }
    if ("Q1" %in% statistics) {
      .q1 <- unlist(lapply(vars,q1),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q1
    }
    if ("Q2" %in% statistics) {
      .q2 <- unlist(lapply(vars,q2),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q2
    }
    if ("Q3" %in% statistics) {
      .q3 <- unlist(lapply(vars,q3),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q3
    }
    if ("Pct" %in% statistics) {
      .pct <- matrix(unlist(lapply(vars,pct,percentil),use.names=F),
                     nrow=numfilas,ncol=length(percentil),byrow=T)
      .tablaposic[,(j+1):(j+length(percentil))]<-.pct
    }    
    i <- i + 1
    res[[i]] <- .tablaposic
  }
  
  if (formacol > 0){
    j <- 0
    if ("H1" %in% statistics) {
      .h1 <- unlist(lapply(vars,h1),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.h1
    }
    if ("H3" %in% statistics) {
      .h3 <- unlist(lapply(vars,h3),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.h3
    }
    if ("K2" %in% statistics) {
      .k2 <- unlist(lapply(vars,k2),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.k2
    }
    if ("K3" %in% statistics) {
      .k3 <- unlist(lapply(vars,k3),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.k3
    }  
    i <- i + 1
    res[[i]] <- .tablaforma
  }  
  class(res) <- "ordindices"
  res
}

print.ordindices <- function(x,...){
  j<-6
  if (sum(c("Mediana","Moda","Trimedia","PromCuar","midR","MediaRec") %in% x[[2]]) > 0){
    .tablaTC <- as.data.frame(x[[j]])
    rownames(.tablaTC) <- x[[1]]
    colnames(.tablaTC) <- c("Mediana","Moda","Trimedia","PromCuar","midR",
                            paste("MediaRec(",x[[3]]*100,"%)",sep=''))[c("Mediana",
                                                                         "Moda","Trimedia","PromCuar","midR","MediaRec") %in% x[[2]]]
    cat("# Indices de tendencia central: \n\n")
    print(.tablaTC)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("Rango", "IQR", "DesvCuar", "MAD", "CVR", "ACent") %in% x[[2]]) > 0){
    .tabladisp <- as.data.frame(x[[j]])
    rownames(.tabladisp) <- x[[1]]
    colnames(.tabladisp) <- c("Rango", "IQR", "DesvCuar", "MAD", "CVR",
                              paste("ACent(",x[[4]]*100,"%)",sep=''))[c("Rango",
                                                                        "IQR","DesvCuar", "MAD", "CVR", "ACent") %in% x[[2]]]
    cat("# Indices de dispersion:","\n\n")
    print(.tabladisp)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("Min", "Max", "Q1", "Q2", "Q3") %in% x[[2]]) + ("Pct" %in% x[[2]])*length(x[[5]]) > 0){
    .tablaposic <- as.data.frame(x[[j]])
    rownames(.tablaposic) <- x[[1]]
    colnames(.tablaposic) <- c("Min", "Max", "Q1", "Q2", "Q3",
                               paste(if (!is.null(x[[5]])) x[[5]]*100 else "Pct",'%',sep=''))[c("Min",
                                                                                                "Max", "Q1", "Q2", "Q3",rep("Pct",if (!is.null(x[[5]])) length(x[[5]])
                                                                                                                            else 1)) %in% x[[2]]]
    cat("# Indices de posicion: \n\n")
    print(.tablaposic)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("H1", "H3", "K2", "K3") %in% x[[2]]) > 0){
    .tablaforma <- as.data.frame(x[[j]])
    rownames(.tablaforma) <- x[[1]]
    colnames(.tablaforma) <- c("H1", "H3", "K2", "K3")[c("H1", "H3", "K2", "K3") %in% x[[2]]]
    cat("# Indices de forma: \n\n")
    print(.tablaforma)  
    cat("\n\n")
    j <- j+1
  }
  
  invisible(x)  
}

resumen.ordinales <- function(){
  defecto <- list(x.inicial=NULL,mediana.inicial="0",moda.inicial="0",trimedia.inicial="0",
                  promcuar.inicial="0",midR.inicial="0",medrec.inicial="0",trim.inicial="0.05",
                  rango.inicial="0",IQR.inicial="0",desvcuar.inicial="0",mad.inicial="0",
                  CVR.inicial="0",ACent.inicial="0",AC.inicial="0.9",min.inicial="0",
                  max.inicial="0",Q1.inicial="0",Q2.inicial="0",Q3.inicial="0",
                  percent.inicial="0",percentil.inicial="0, .25, .5, .75, 1",
                  H1.inicial="0",H3.inicial="0",K2.inicial="0",K3.inicial="0",
                  echo.inicial="0",creahtml.inicial="0",selectodas.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("resumen.ordinales",defecto) 
  initializeDialog(title=gettextRcmdr("Indices para variables ordinales"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','statistics2Tab','statistics3Tab','statistics4Tab','optionsTab'))
  listaVar <- variableListBox(dataTab, Factors(), selectmode="multiple",
                              title=gettextRcmdr("Variables (escoja una o mas)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"factor"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    if (length(x) == 0){
      errorCondition(recall=resumen.ordinales, message=gettextRcmdr("Debes escoger una variable."))
      return()
    }
    for (variable in x)
    {          
      justDoIt(paste("cond <- !is.ordered(",paste(ActiveDataSet(),"$",variable,sep=""),")",sep=""))
      if (cond){
        errorCondition(recall=resumen.ordinales, message=gettextRcmdr(paste("Variable ",variable, " no es ordinal.",sep='')))
        return()
      }
      remove("cond", envir=.GlobalEnv)
    }        
    .BaseDatosActiva <- ActiveDataSet()
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    selectodas <- tclvalue(selectodasVariable)
    if (selectodas == 1)
    {
      medianaval = modaval = trimediaval = promcuarval = midRval = 
        medrecval = rangoval = IQRval = desvcuarval =
        madval = CVRval = ACentval = minval = maxval = Q1val = Q2val =
        Q3val = percentval = H1val = H3val = K2val = K3val = TRUE
    }
    else
    {
      medianaval <- tclvalue(medianaVariable)
      modaval <- tclvalue(modaVariable)
      trimediaval <- tclvalue(trimediaVariable)
      promcuarval <- tclvalue(promcuarVariable)
      midRval <- tclvalue(midRVariable)
      medrecval <- tclvalue(medrecVariable)
      rangoval <- tclvalue(rangoVariable)
      IQRval <- tclvalue(IQRVariable)
      desvcuarval <- tclvalue(desvcuarVariable)
      madval <- tclvalue(madVariable)
      CVRval <- tclvalue(CVRVariable)
      ACentval <- tclvalue(ACentVariable)  
      minval <- tclvalue(minVariable)
      maxval <- tclvalue(maxVariable)
      Q1val <- tclvalue(Q1Variable)
      Q2val <- tclvalue(Q2Variable)
      Q3val <- tclvalue(Q3Variable)
      percentval <- tclvalue(percentVariable)
      H1val <- tclvalue(H1Variable)
      H3val <- tclvalue(H3Variable)
      K2val <- tclvalue(K2Variable)
      K3val <- tclvalue(K3Variable)
    }        
    selec <- as.numeric(medianaval) + as.numeric(modaval) + 
      as.numeric(trimediaval) + as.numeric(promcuarval) +
      as.numeric(midRval) + as.numeric(medrecval)
    selec2 <- as.numeric(rangoval) +
      as.numeric(IQRval) + as.numeric(desvcuarval) +        
      as.numeric(madval) + as.numeric(CVRval) +
      as.numeric(ACentval)
    selec3 <- as.numeric(minval) + as.numeric(maxval) + 
      as.numeric(Q1val) + as.numeric(Q2val) + 
      as.numeric(Q3val) + as.numeric(percentval)
    selec4 <- as.numeric(H1val) + as.numeric(H3val) + 
      as.numeric(K2val) + as.numeric(K3val)
    seleccion <- selec + selec2 + selec3 + selec4
    if (seleccion == 0){
      errorCondition(recall=resumen.ordinales, 
                     message=gettextRcmdr("Debes escoger algun indicador."))
      return()
    }
    putDialog("resumen.ordinales",list(x.inicial=x,mediana.inicial=medianaval,moda.inicial=modaval,
                                       trimedia.inicial=trimediaval,promcuar.inicial=promcuarval,midR.inicial=midRval,
                                       medrec.inicial=medrecval,trim.inicial=tclvalue(trimVariable),
                                       rango.inicial=rangoval,IQR.inicial=IQRval,desvcuar.inicial=desvcuarval,
                                       mad.inicial=madval,CVR.inicial=CVRval,ACent.inicial=ACentval,
                                       AC.inicial=tclvalue(ACVariable),min.inicial=minval,max.inicial=maxval,
                                       Q1.inicial=Q1val,Q2.inicial=Q2val,Q3.inicial=Q3val,percent.inicial=percentval,
                                       percentil.inicial=tclvalue(percentilVariable),H1.inicial=H1val,
                                       H3.inicial=H3val,K2.inicial=K2val,K3.inicial=K3val,
                                       selectodas.inicial=selectodas,echo.inicial=echocodigo,
                                       creahtml.inicial=creahtml,tab.inicial=tab))   
    if (percentval == 1)
    {
      pct <- c(gsub(" ", ",",gsub(", ", ",",tclvalue(percentilVariable))))
      pct1 <- as.numeric(unlist(strsplit(pct,",")))
      if ( is.na(pct1) || (sum(pct1<0.0)>0) || (sum(pct1>1.0)>0) || (sum(!is.numeric(pct1))>0) )
      {
        pct <- paste(seq(0.,1.,.25),collapse=",")
        Message(message=gettextRcmdr("Vector de percentiles invalido. Se utilizara vector por defecto."),
                type="warning")
      }
    } else pct <- NULL
    if (medrecval == 1)
    {
      rec <- as.numeric(tclvalue(trimVariable))
      if ( rec < .0 || rec > .5 || !is.numeric(rec) )
      {
        rec <- 0.05
        Message(message=gettextRcmdr("Proporcion de recorte invalida se utilizara valor por defecto."),
                type="warning")              
      } 
    } else rec <- NULL
    if (ACentval == 1)
    {
      propdat <- as.numeric(tclvalue(ACVariable))
      if ( propdat < .0 || propdat > 1. || !is.numeric(propdat) )
      {
        prop.dat <- 0.9
        Message(message=gettextRcmdr("Proporcion de datos invalida se utilizara valor por defecto."),
                type="warning")              
      }
    } else propdat <- NULL
    
    vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
    else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
    if (length(x) == 1) variables <- paste(.BaseDatosActiva, "[", vars, "]", sep="")
    else variables <- paste(.BaseDatosActiva, "[,", vars, "]", sep="")
    stats <- paste("c(",
                   paste(c('"Mediana"', '"Moda"', '"Trimedia"', '"PromCuar"', '"midR"', '"MediaRec"',
                           '"Rango"', '"IQR"', '"DesvCuar"', '"MAD"', '"CVR"', '"ACent"',
                           '"Min"', '"Max"', '"Q1"', '"Q2"', '"Q3"', '"Pct"',
                           '"H1"', '"H3"', '"K2"', '"K3"')
                         [c(medianaval, modaval, trimediaval, promcuarval, midRval, medrecval,
                            rangoval, IQRval, desvcuarval, madval, CVRval, ACentval,
                            minval, maxval, Q1val, Q2val, Q3val, percentval,
                            H1val, H3val, K2val, K3val) == 1], 
                         collapse=", "), ")", sep="")
    if (percentval ==1)
      instruccion1 <- paste(".indices.ord <- ordindices(vars=",variables, ", statistics=", stats, if (!is.null(rec)){paste(",rec=",rec)},
                            if(!is.null(propdat)){paste(",propdat=",propdat)},",percentil=c(",pct,")",")", sep="")
    else
      instruccion1 <- paste(".indices.ord <- ordindices(vars=",variables, ", statistics=", stats, if (!is.null(rec)){paste(",rec=",rec)},
                            if(!is.null(propdat)){paste(",propdat=",propdat)},")", sep="")  
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(".indices.ord  # Indices para variables ordinales ")  
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- "Indicadores descriptivos para variables ordinales"
      HTML(as.title(titulo),file=.archivo)
      j <- 6
      if ( selec>0 ){
        .TablaRes <- as.data.frame(.indices.ord[[j]])
        rownames(.TablaRes) <- .indices.ord[[1]]
        colnames(.TablaRes) <- c("Mediana","Moda","Trimedia","PromCuar","midR",
                                 paste("MediaRec(",.indices.ord[[3]]*100,"%)",sep=''))[c("Mediana", "Moda", "Trimedia",
                                                                                         "PromCuar", "midR", "MediaRec") %in% .indices.ord[[2]]]
        HTML("Indicadores de tendencia central: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }
      if ( selec2>0 ){
        .TablaRes <- as.data.frame(.indices.ord[[j]])
        rownames(.TablaRes) <- .indices.ord[[1]]
        colnames(.TablaRes) <- c("Rango", "IQR", "DesvCuar", "MAD", "CVR",
                                 paste("ACent(",.indices.ord[[4]]*100,"%)",sep=''))[c("Rango", "IQR",
                                                                                      "DesvCuar", "MAD", "CVR", "ACent") %in% .indices.ord[[2]]]
        HTML("Indicadores de dispersion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      } 
      if ( selec3>0 ){
        .TablaRes <- as.data.frame(.indices.ord[[j]])
        rownames(.TablaRes) <- .indices.ord[[1]]
        colnames(.TablaRes) <- c("Min", "Max", "Q1", "Q2", "Q3",
                                 paste(.indices.ord[[5]]*100,'%',sep=''))[c("Min", "Max", "Q1", 
                                                                            "Q2", "Q3",rep("Pct",length(.indices.ord[[5]]))) %in% .indices.ord[[2]]]
        HTML("Indicadores de posicion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }   
      if ( selec4>0 ){
        .TablaRes <- as.data.frame(.indices.ord[[j]])
        rownames(.TablaRes) <- .indices.ord[[1]]
        colnames(.TablaRes) <- c("H1", "H3", "K2", "K3")[c("H1", "H3", "K2", "K3")
                                                         %in% .indices.ord[[2]]]
        HTML("Indicadores de forma: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }       
    }
    
    remove(.indices.ord, envir=.GlobalEnv)
    closeDialog()
    if (echocodigo == 1) logger("remove(.indices.ord)")
    tkfocus(CommanderWindow())    
  }
  OKCancelHelp(helpSubject="RcmdrPlugin.EACSPIR",reset="resumen.ordinales",apply="resumen.ordinales")
  tkgrid(getFrame(listaVar), sticky="nw")
  checkBoxes(statisticsTab,frame="tcFrame",boxes=c("mediana","moda","trimedia"),
             initialValues=c(dialog.valores$mediana.inicial,dialog.valores$moda.inicial,
                             dialog.valores$trimedia.inicial),
             labels=gettextRcmdr(c("Mediana ","Moda ","Trimedia ")), 
             title = gettextRcmdr("Indices de tendencia central"))
  checkBoxes(statisticsTab,frame="tcRFrame",boxes=c("promcuar","midR","medrec"),
             initialValues=c(dialog.valores$promcuar.inicial,
                             dialog.valores$midR.inicial,dialog.valores$medrec.inicial),
             labels=gettextRcmdr(c("Promedio de cuartiles ","Rango medio ","Media Recortada ")), 
             title = gettextRcmdr(" "))  
  trimFrame <- tkframe(statisticsTab)
  trimVariable <- tclVar(dialog.valores$trim.inicial)
  trimField <- ttkentry(trimFrame, width="8", textvariable=trimVariable)
  tkgrid(labelRcmdr(trimFrame,text="  "))
  tkgrid(labelRcmdr(trimFrame,text="  "))
  tkgrid(labelRcmdr(trimFrame,text="  ")) 
  tkgrid(labelRcmdr(trimFrame,text=gettextRcmdr("Proporcion datos recortados = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         trimField,sticky="w")
  tkgrid(tcFrame, labelRcmdr(statisticsTab,text="  "),tcRFrame, labelRcmdr(statisticsTab,text="  "), 
         trimFrame, sticky="nw")
  
  checkBoxes(statistics2Tab,frame="dispFrame",boxes=c("rango","IQR","desvcuar"),
             initialValues=c(dialog.valores$rango.inicial,dialog.valores$IQR.inicial,
                             dialog.valores$desvcuar.inicial),
             labels=gettextRcmdr(c("Amplitud ","Amplitud intercuartil (IQR) ","Desviacion Cuartil")), 
             title = gettextRcmdr("Indices de dispersion")) 
  checkBoxes(statistics2Tab,frame="dispRFrame",boxes=c("mad","CVR","ACent"),
             initialValues=c(dialog.valores$mad.inicial,
                             dialog.valores$CVR.inicial,dialog.valores$ACent.inicial),
             labels=gettextRcmdr(c("Mediana Desviaciones Absolutas (MAD) ","Coeficiente Variacion Robusto ",
                                   "Desviacion centilica ")), 
             title = gettextRcmdr("  "))  
  ACFrame <- tkframe(statistics2Tab)
  ACVariable <- tclVar(dialog.valores$AC.inicial)
  ACField <- ttkentry(ACFrame, width="4", textvariable=ACVariable)
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text=gettextRcmdr("Proporcion datos utilizados = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         ACField,sticky="w")
  tkgrid(dispFrame, labelRcmdr(statistics2Tab,text="  "),
         dispRFrame, labelRcmdr(statistics2Tab,text="  "), ACFrame, sticky="nw")
  
  checkBoxes(statistics3Tab,frame="posicFrame",boxes=c("Q1","Q2","Q3"),
             initialValues=c(dialog.valores$Q1.inicial,dialog.valores$Q2.inicial,
                             dialog.valores$Q3.inicial),
             labels=gettextRcmdr(c("Primer Cuartil",
                                   "Segundo Cuartil ","Tercer Cuartil ")), 
             title = gettextRcmdr("Indices de posicion"))
  checkBoxes(statistics3Tab,frame="posicRFrame",boxes=c("min","max","percent"),
             initialValues=c(dialog.valores$min.inicial,dialog.valores$max.inicial,
                             dialog.valores$percent.inicial),
             labels=gettextRcmdr(c("Minimo","Maximo","Cuantilas")), 
             title = gettextRcmdr("  "))  
  percentFrame <- tkframe(statistics3Tab)
  percentilVariable <- tclVar(dialog.valores$percentil.inicial)
  percentField <- ttkentry(percentFrame, width="15", textvariable=percentilVariable)
  tkgrid(labelRcmdr(percentFrame,text="  "))
  tkgrid(labelRcmdr(percentFrame,text="  "))
  tkgrid(labelRcmdr(percentFrame,text="  ")) 
  tkgrid(labelRcmdr(percentFrame,text=gettextRcmdr("Seleccione cuantilas = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         percentField,sticky="w")
  tkgrid(posicFrame, labelRcmdr(statistics3Tab,text="  "),
         posicRFrame, labelRcmdr(statistics3Tab,text="  "), percentFrame, sticky="nw")  
  
  checkBoxes(statistics4Tab,frame="formaFrame",boxes=c("H1","H3"),
             initialValues=c(dialog.valores$H1.inicial,dialog.valores$H3.inicial),
             labels=gettextRcmdr(c("Coef. Asimetria H1 ","Coef. Asimetria H3 ")), 
             title = gettextRcmdr("Indices de forma")) 
  checkBoxes(statistics4Tab,frame="formaRFrame",boxes=c("K2","K3"),
             initialValues=c(dialog.valores$K2.inicial,dialog.valores$K3.inicial),
             labels=gettextRcmdr(c("Coef. Apuntamiento K2 ","Coef. Apuntamiento K3 ")), 
             title = gettextRcmdr("   ")) 
  tkgrid(formaFrame, labelRcmdr(statistics4Tab,text="  "), formaRFrame, sticky="w")
  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml","selectodas"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial,
                             dialog.valores$selectodas.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ",
                                   "Calcular todos los indices ")), 
             title = gettextRcmdr("Opciones"))
  tkgrid(opsFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','statistics2Tab',
                                                      'statistics3Tab','statistics4Tab','optionsTab'),
               tab.names=c("Datos","Tend. Central","Dispersion","Posicion","Forma","Opciones"))
  
}

# Resumen descriptivo univariante para variables cuantitativas: indicadores y graficos #

numindices <- function(vars,statistics,rec=NULL,propdat=NULL,percentil=NULL){
  
  vars <- as.data.frame(vars)
  variables <- names(vars)
  
  media <- function(x){
    res <- round(mean(x,na.rm=TRUE),2)
    res
  }    
  mediana <- function(x){
    res <- round(median(x,na.rm=TRUE),2)
    res
  }  
  moda <- function(x){
    .ni <- table(x)
    res <- as.numeric(names(.ni[which(.ni==max(.ni,na.rm=TRUE))]))
    if (length(res) > 1) res <- res[[1]]
    res
  }
  moda2 <- function(x){
    .ni <- table(x)
    res <- as.numeric(names(.ni[which(.ni==max(.ni,na.rm=TRUE))]))
    res
  }  
  mediageom <- function(x){
    res <- round(mean(log(x)[x!=0],na.rm=TRUE),2)
    res
  }
  trimedia <- function(x){
    .Finf <- fivenum(x)[2]
    .Md <- fivenum(x)[3]
    .Fsup <- fivenum(x)[4]
    res <- round((.Finf+2*.Md+.Fsup)/4,2)  
    res
  }
  
  promcuar <- function(x){
    .Q1 <- quantile(x,na.rm=TRUE)[2]
    .Q3 <- quantile(x,na.rm=TRUE)[4]
    res <- round((.Q1+.Q3)/2,2)
    names(res) <- NULL
    res
  }
  midR <- function(x){
    .min <- quantile(x,na.rm=TRUE)[1]
    .max <- quantile(x,na.rm=TRUE)[5]
    res <- round((.min+.max)/2,2)
    names(res) <- NULL    
    res
  }
  medrec <- function(x,rec){
    res <- round(mean(x,trim=rec,na.rm=TRUE),2)
    res
  }
  variancia <- function(x){
    res <- round(var(x,na.rm=TRUE),2)
    res
  }  
  dt <- function(x){
    res <- round(sd(x,na.rm=TRUE),2)
    res
  }  
  dtgeom <- function(x){
    res <- round(sd(log(x)[x!=0],na.rm=TRUE),2)
    res
  }    
  desvmed <- function(x){
    .media <- mean(x,na.rm=TRUE)
    res <- round(sum(abs(x - .media),na.rm=TRUE)/length(na.omit(x)),2)
    res
  }  
  CV <- function(x){
    .media <- mean(x,na.rm=TRUE)
    .dt <- sd(x,na.rm=TRUE)
    res <- round(.dt/.media,2)
    res
  }
  rango <- function(x){
    .min <- quantile(x,na.rm=TRUE)[1]
    .max <- quantile(x,na.rm=TRUE)[5]
    res <- round(.max-.min,2)
    names(res) <- NULL        
    res
  }
  iqr <- function(x){
    .Q1 <- quantile(x,na.rm=TRUE)[2]
    .Q3 <- quantile(x,na.rm=TRUE)[4]
    res <- round(.Q3-.Q1,2)
    names(res) <- NULL        
    res
  }
  mad <- function(x){
    .mediana <- median(x,na.rm=TRUE)
    res <- round(median(abs(x-.mediana),na.rm=TRUE),2)
    res
  }
  cvr <- function(x){
    .Finf <- fivenum(x)[2]
    .Fsup <- fivenum(x)[4]
    res <- round((.Fsup-.Finf)/(.Finf+.Fsup), 2)
    res
  }
  desvcuar <- function(x){
    .Q1 <- quantile(x,na.rm=TRUE)[2]
    .Q3 <- quantile(x,na.rm=TRUE)[4]
    res <- round((.Q3-.Q1)/2,2)
    names(res) <- NULL    
    res
  }
  acent <- function(x,propdat){
    .ACinf <- quantile(x,probs=(1-propdat)/2,na.rm=TRUE)
    .ACsup <- quantile(x,probs=1-(1-propdat)/2,na.rm=TRUE)
    res <- round(.ACsup-.ACinf,2)
    names(res) <- NULL       
    res
  }
  q1 <- function(x){
    res <- round(quantile(x,na.rm=TRUE)[2],2) 
    names(res) <- NULL        
    res
  }
  q2 <- function(x){
    res <- round(quantile(x,na.rm=TRUE)[3],2) 
    names(res) <- NULL        
    res
  }
  q3 <- function(x){
    res <- round(quantile(x,na.rm=TRUE)[4],2) 
    names(res) <- NULL        
    res
  }
  pct <- function(x,percentil){
    res <- round(quantile(x,probs=percentil,na.rm=TRUE),2)
    names(res) <- NULL       
    res
  }
  h1 <- function(x) {
    .Finf <- fivenum(x)[2]
    .Md <- fivenum(x)[3]
    .Fsup <- fivenum(x)[4]
    res <- round((.Finf+.Fsup-2*.Md)/(2*.Md),2)
    res
  }
  h3 <- function(x) {
    .AC90 <- quantile(x,probs=0.9,na.rm=TRUE)
    .AC10 <- quantile(x,probs=0.1,na.rm=TRUE)
    .Md <- fivenum(x)[3]
    res <- round((.AC90+.AC10-2*.Md)/(2*.Md),2)
    names(res) <- NULL        
    res
  }
  k2 <- function(x){
    .AC90 <- quantile(x,probs=0.9,na.rm=TRUE)
    .AC10 <- quantile(x,probs=0.1,na.rm=TRUE)
    .Q1 <- quantile(x,na.rm=TRUE)[2]
    .Q3 <- quantile(x,na.rm=TRUE)[4]
    res <- round((.AC90-.AC10)/(1.9*(.Q3-.Q1)),2)
    names(res) <- NULL        
    res
  }
  k3 <- function(x){
    .Einf <- quantile(x,probs=0.125,na.rm=TRUE)
    .Esup <- quantile(x,probs=0.875,na.rm=TRUE)
    .Finf <- fivenum(x)[2]
    .Fsup <- fivenum(x)[4]
    res <- round((.Esup-.Einf)/(1.7*(.Fsup-.Finf)),2)
    names(res) <- NULL        
    res
  }
  beta1 <- function(x){
    .media <- mean(x,na.rm=TRUE)
    .n <- length(na.omit(x))
    res <- round((sum((x-.media)^3,na.rm=TRUE)/.n)^2/(sum((x-.media)^2,na.rm=TRUE)/.n)^3,2)
    res
  }
  gamma1 <- function(x){
    .media <- mean(x,na.rm=TRUE)
    .dt <- sd(x,na.rm=TRUE)
    .n <- length(na.omit(x))
    res <- round(.n*sum((x-.media)^3,na.rm=TRUE)/((.n-1)*(.n-2))/(.dt^3),2)
    res
  }  
  beta2 <- function(x){
    .media <- mean(x,na.rm=TRUE)
    .n <- length(na.omit(x))
    res <- round(sum((x-.media)^4,na.rm=TRUE)/.n/(sum((x-.media)^2,na.rm=TRUE)/.n)^2,2)
    res
  }
  gamma2 <- function(x){
    .media <- mean(x,na.rm=TRUE)
    .dt <- sd(x,na.rm=TRUE)
    .n <- length(na.omit(x))
    res <- round((.n*(.n+1)*sum((x-.media)^4,na.rm=TRUE)/((.n-1)*(.n-2)*(.n-3))-3*
                    sum((x-.media)^2,na.rm=TRUE)^2/((.n-2)*(.n-3)))/(.dt^4),2)
    res
  }   
  
  res <- list()
  res[[1]] <- variables
  res[[2]] <- statistics
  res[[3]] <- rec
  res[[4]] <- propdat
  res[[5]] <- percentil
  i <- 5
  
  numfilas <- length(vars)
  TCcol <- sum(c("Media", "Mediana", "Moda", "MediaGeom", "Trimedia", "PromCuar",
                 "midR", "MediaRec") %in% statistics)
  dispcol <- sum(c("Variancia", "DT", "DTGeom", "DesvMed", "CV", "Rango",
                   "IQR", "DesvCuar", "MAD", "CVR", "ACent") %in% statistics)
  posiccol <- sum(c("Min", "Max", "Q1", "Q2", "Q3") %in% statistics) + ("Pct" %in% statistics)*length(percentil)
  formacol <- sum(c("H1", "H3", "K2", "K3", "beta1", "gamma1", "beta2", "gamma2") %in% statistics)
  .tablaTC <- as.data.frame(matrix(nrow=numfilas,ncol=TCcol))
  .tabladisp <- as.data.frame(matrix(nrow=numfilas,ncol=dispcol))
  .tablaposic <- as.data.frame(matrix(nrow=numfilas,ncol=posiccol))
  .tablaforma <- as.data.frame(matrix(nrow=numfilas,ncol=formacol))
  
  if (TCcol > 0){
    j <- 0
    if ("Media" %in% statistics) {
      .media <- unlist(lapply(vars,media),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.media
    }    
    if ("Mediana" %in% statistics) {
      .mediana <- unlist(lapply(vars,mediana),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.mediana
    }
    if ("Moda" %in% statistics) {
      .moda <- unlist(lapply(vars,moda),use.names=F)
      .moda2 <- lapply(vars,moda2)
      for (k in length(.moda2))
      {
        if (length(.moda2[[k]]) > 1) {
          Message(message=gettextRcmdr(paste("Variable ",variables[k]," tiene mas de una moda: ",
                                             "En documento HTML se muestra un solo valor.", sep="")),
                  type="warning") 
        }
      }     
      j <- j + 1
      .tablaTC[,j]<-.moda
    }
    if ("MediaGeom" %in% statistics) {
      .mediageom <- unlist(lapply(vars,mediageom),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.mediageom
    }    
    if ("Trimedia" %in% statistics) {
      .trimedia <- unlist(lapply(vars,trimedia),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.trimedia
    }
    if ("PromCuar" %in% statistics) {
      .promcuar <- unlist(lapply(vars,promcuar),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.promcuar
    }
    if ("midR" %in% statistics) {
      .midr <- unlist(lapply(vars,midR),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.midr
    }
    if ("MediaRec" %in% statistics) {
      .mediarec <- unlist(lapply(vars,medrec,rec),use.names=F)
      j <- j + 1
      .tablaTC[,j]<-.mediarec
    }    
    i <- i + 1
    res[[i]] <- .tablaTC
  }
  
  if (dispcol > 0){
    j <- 0
    if ("Variancia" %in% statistics) {
      .var <- unlist(lapply(vars,variancia),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.var
    }
    if ("DT" %in% statistics) {
      .dt <- unlist(lapply(vars,dt),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.dt
    }    
    if ("DTGeom" %in% statistics) {
      .dtgeom <- unlist(lapply(vars,dtgeom),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.dtgeom
    }    
    if ("DesvMed" %in% statistics) {
      .desvmed <- unlist(lapply(vars,desvmed),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.desvmed
    }
    if ("CV" %in% statistics) {
      .cv <- unlist(lapply(vars,CV),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.cv
    }    
    if ("Rango" %in% statistics) {
      .rango <- unlist(lapply(vars,rango),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.rango
    }
    if ("IQR" %in% statistics) {
      .iqr <- unlist(lapply(vars,iqr),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.iqr
    }
    if ("DesvCuar" %in% statistics) {
      .dq <- unlist(lapply(vars,desvcuar),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.dq
    }
    if ("MAD" %in% statistics) {
      .mad <- unlist(lapply(vars,mad),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.mad
    }
    if ("CVR" %in% statistics) {
      .cvr <- unlist(lapply(vars,cvr),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.cvr
    }
    if ("ACent" %in% statistics) {
      .acent <- unlist(lapply(vars,acent,propdat),use.names=F)
      j <- j + 1
      .tabladisp[,j]<-.acent
    }    
    i <- i + 1
    res[[i]] <- .tabladisp
  }
  
  
  if (posiccol > 0){
    j <- 0
    if ("Min" %in% statistics) {
      .min <- unlist(lapply(lapply(vars,as.numeric),min,na.rm=TRUE),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.min
    }
    if ("Max" %in% statistics) {
      .max <- unlist(lapply(lapply(vars,as.numeric),max,na.rm=TRUE),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.max
    }
    if ("Q1" %in% statistics) {
      .q1 <- unlist(lapply(vars,q1),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q1
    }
    if ("Q2" %in% statistics) {
      .q2 <- unlist(lapply(vars,q2),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q2
    }
    if ("Q3" %in% statistics) {
      .q3 <- unlist(lapply(vars,q3),use.names=F)
      j <- j + 1
      .tablaposic[,j]<-.q3
    }
    if ("Pct" %in% statistics) {
      .pct <- matrix(unlist(lapply(vars,pct,percentil),use.names=F),
                     nrow=numfilas,ncol=length(percentil),byrow=T)
      .tablaposic[,(j+1):(j+length(percentil))]<-.pct
    }    
    i <- i + 1
    res[[i]] <- .tablaposic
  }
  
  if (formacol > 0){
    j <- 0
    if ("H1" %in% statistics) {
      .h1 <- unlist(lapply(vars,h1),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.h1
    }
    if ("H3" %in% statistics) {
      .h3 <- unlist(lapply(vars,h3),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.h3
    }
    if ("K2" %in% statistics) {
      .k2 <- unlist(lapply(vars,k2),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.k2
    }
    if ("K3" %in% statistics) {
      .k3 <- unlist(lapply(vars,k3),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.k3
    }  
    
    if ("beta1" %in% statistics) {
      .beta1 <- unlist(lapply(vars,beta1),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.beta1
    }
    if ("gamma1" %in% statistics) {
      .gamma1 <- unlist(lapply(vars,gamma1),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.gamma1
    }  
    if ("beta2" %in% statistics) {
      .beta2 <- unlist(lapply(vars,beta2),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.beta2
    }  
    if ("gamma2" %in% statistics) {
      .gamma2 <- unlist(lapply(vars,gamma2),use.names=F)
      j <- j + 1
      .tablaforma[,j]<-.gamma2
    }    
    i <- i + 1
    res[[i]] <- .tablaforma   
  }
  
  class(res) <- "numindices"
  res
}

print.numindices <- function(x,...){
  j<-6
  if (sum(c("Media", "Mediana", "Moda", "MediaGeom", "Trimedia", "PromCuar",
            "midR", "MediaRec") %in% x[[2]]) > 0){
    .tablaTC <- as.data.frame(x[[j]])
    rownames(.tablaTC) <- x[[1]]
    colnames(.tablaTC) <- c("Media", "Mediana", "Moda", "MediaGeom", "Trimedia", "PromCuar",
                            "midR", paste("MediaRec(",x[[3]]*100,"%)",sep=''))[c("Media", 
                                                                                 "Mediana", "Moda", "MediaGeom", "Trimedia", "PromCuar",
                                                                                 "midR", "MediaRec") %in% x[[2]]]
    cat("# Indices de tendencia central: \n\n")
    print(.tablaTC)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("Variancia", "DT", "DTGeom", "DesvMed", "CV", "Rango",
            "IQR", "DesvCuar", "MAD", "CVR", "ACent") %in% x[[2]]) > 0){
    .tabladisp <- as.data.frame(x[[j]])
    rownames(.tabladisp) <- x[[1]]
    colnames(.tabladisp) <- c("Variancia", "DT", "DTGeom", "DesvMed", "CV", "Rango",
                              "IQR", "DesvCuar", "MAD", "CVR",
                              paste("ACent(",x[[4]]*100,"%)",sep=''))[c("Variancia", 
                                                                        "DT", "DTGeom", "DesvMed", "CV", "Rango", "IQR", "DesvCuar",
                                                                        "MAD", "CVR", "ACent") %in% x[[2]]]
    cat("# Indices de dispersion:","\n\n")
    print(.tabladisp)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("Min", "Max", "Q1", "Q2", "Q3") %in% x[[2]]) + ("Pct" %in% x[[2]])*length(x[[5]]) > 0){
    .tablaposic <- as.data.frame(x[[j]])
    rownames(.tablaposic) <- x[[1]]
    colnames(.tablaposic) <- c("Min", "Max", "Q1", "Q2", "Q3",
                               paste(if (!is.null(x[[5]])) x[[5]]*100 else "Pct",'%',sep=''))[c("Min",
                                                                                                "Max", "Q1", "Q2", "Q3",rep("Pct",if (!is.null(x[[5]])) length(x[[5]])
                                                                                                                            else 1)) %in% x[[2]]]
    cat("# Indices de posicion: \n\n")
    print(.tablaposic)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("H1", "H3", "K2", "K3", "beta1", "gamma1", "beta2", "gamma2") %in% x[[2]]) > 0){
    .tablaforma <- as.data.frame(x[[j]])
    rownames(.tablaforma) <- x[[1]]
    colnames(.tablaforma) <- c("H1", "H3", "K2", "K3", "beta1", "gamma1", "beta2", "gamma2")[c("H1", 
                                                                                               "H3", "K2", "K3", "beta1", "gamma1", "beta2", "gamma2") %in% x[[2]]]
    cat("# Indices de forma: \n\n")
    print(.tablaforma)  
    cat("\n\n")
  }
  
  invisible(x)  
}

resumen.numericas <- function(){
  defecto <- list(x.inicial=NULL,media.inicial="0",mediana.inicial="0",moda.inicial="0",
                  mediageom.inicial="0",trimedia.inicial="0",promcuar.inicial="0",
                  midR.inicial="0",medrec.inicial="0",trim.inicial="0.05",variancia.inicial="0",
                  dt.inicial="0",dtgeom.inicial="0",desvmed.inicial="0",CV.inicial="0",
                  rango.inicial="0",IQR.inicial="0",desvcuar.inicial="0",mad.inicial="0",
                  CVR.inicial="0",ACent.inicial="0",AC.inicial="0.9",min.inicial="0",
                  max.inicial="0",Q1.inicial="0",Q2.inicial="0",Q3.inicial="0",
                  percent.inicial="0",percentil.inicial="0, .25, .5, .75, 1",
                  H1.inicial="0",H3.inicial="0",K2.inicial="0",K3.inicial="0",
                  beta1.inicial="0",gamma1.inicial="0",beta2.inicial="0",gamma2.inicial="0",
                  echo.inicial="0",creahtml.inicial="0",selectodas.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("resumen.numericas",defecto) 
  initializeDialog(title=gettextRcmdr("Indices para variables cuantitativas"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','statistics2Tab','statistics3Tab','statistics4Tab','optionsTab'))
  listaVar <- variableListBox(dataTab, Numeric(), selectmode="multiple",
                              title=gettextRcmdr("Variables (escoja una o mas)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"numeric"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    if (length(x) == 0){
      errorCondition(recall=resumen.ordinales, message=gettextRcmdr("Debes escoger una variable."))
      return()
    }        
    .BaseDatosActiva <- ActiveDataSet()
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    selectodas <- tclvalue(selectodasVariable)
    if (selectodas == 1)
    {
      mediaval = medianaval = modaval = mediageomval = trimediaval =
        promcuarval = midRval = medrecval = varianciaval = dtval = CVval =
        dtgeomval = desvmedval = rangoval = IQRval = desvcuarval =
        madval = CVRval = ACentval = minval = maxval = Q1val = Q2val =
        Q3val = percentval = H1val = H3val = K2val = K3val = beta1val =
        gamma1val = beta2val = gamma2val = TRUE
    }
    else
    {
      mediaval <- tclvalue(mediaVariable)
      medianaval <- tclvalue(medianaVariable)
      modaval <- tclvalue(modaVariable)
      mediageomval <- tclvalue(mediageomVariable)
      trimediaval <- tclvalue(trimediaVariable)
      promcuarval <- tclvalue(promcuarVariable)
      midRval <- tclvalue(midRVariable)
      medrecval <- tclvalue(medrecVariable)
      varianciaval <- tclvalue(varianciaVariable)
      dtval <- tclvalue(dtVariable)
      CVval <- tclvalue(CVVariable)
      dtgeomval <- tclvalue(dtgeomVariable)
      desvmedval <- tclvalue(desvmedVariable)
      rangoval <- tclvalue(rangoVariable)
      IQRval <- tclvalue(IQRVariable)
      desvcuarval <- tclvalue(desvcuarVariable)
      madval <- tclvalue(madVariable)
      CVRval <- tclvalue(CVRVariable)
      ACentval <- tclvalue(ACentVariable)  
      minval <- tclvalue(minVariable)
      maxval <- tclvalue(maxVariable)
      Q1val <- tclvalue(Q1Variable)
      Q2val <- tclvalue(Q2Variable)
      Q3val <- tclvalue(Q3Variable)
      percentval <- tclvalue(percentVariable)
      H1val <- tclvalue(H1Variable)
      H3val <- tclvalue(H3Variable)
      K2val <- tclvalue(K2Variable)
      K3val <- tclvalue(K3Variable)
      beta1val <- tclvalue(beta1Variable)        
      gamma1val <- tclvalue(gamma1Variable)        
      beta2val <- tclvalue(beta2Variable)        
      gamma2val <- tclvalue(gamma2Variable)
    }        
    selec <- as.numeric(mediaval) + as.numeric(medianaval) + 
      as.numeric(modaval) + as.numeric(mediageomval) + 
      as.numeric(trimediaval) + as.numeric(promcuarval) +
      as.numeric(midRval) + as.numeric(medrecval)
    selec2 <- as.numeric(varianciaval) + as.numeric(dtval) + 
      as.numeric(CVval) + as.numeric(dtgeomval) + 
      as.numeric(desvmedval) + as.numeric(rangoval) +
      as.numeric(IQRval) + as.numeric(desvcuarval) +        
      as.numeric(madval) + as.numeric(CVRval) +
      as.numeric(ACentval)
    selec3 <- as.numeric(minval) + as.numeric(maxval) + 
      as.numeric(Q1val) + as.numeric(Q2val) + 
      as.numeric(Q3val) + as.numeric(percentval)
    selec4 <- as.numeric(H1val) + as.numeric(H3val) + 
      as.numeric(K2val) + as.numeric(K3val) + 
      as.numeric(beta1val) + as.numeric(gamma1val) +
      as.numeric(beta2val) + as.numeric(gamma2val)
    seleccion <- selec + selec2 + selec3 + selec4
    if (seleccion == 0){
      errorCondition(recall=resumen.numericas, 
                     message=gettextRcmdr("Debes escoger algun indicador."))
      return()
    }
    putDialog("resumen.numericas",list(x.inicial=x,media.inicial=mediaval,mediana.inicial=medianaval,
                                       moda.inicial=modaval,mediageom.inicial=mediageomval,
                                       trimedia.inicial=trimediaval,promcuar.inicial=promcuarval,
                                       midR.inicial=midRval,medrec.inicial=medrecval,
                                       trim.inicial=tclvalue(trimVariable),variancia.inicial=varianciaval,
                                       dt.inicial=dtval,dtgeom.inicial=dtgeomval,desvmed.inicial=desvmedval,
                                       CV.inicial=CVval,rango.inicial=rangoval,IQR.inicial=IQRval,
                                       desvcuar.inicial=desvcuarval,mad.inicial=madval,CVR.inicial=CVRval,
                                       ACent.inicial=ACentval,AC.inicial=tclvalue(ACVariable),min.inicial=minval,
                                       max.inicial=maxval,Q1.inicial=Q1val,Q2.inicial=Q2val,Q3.inicial=Q3val,
                                       percent.inicial=percentval,percentil.inicial=tclvalue(percentilVariable),
                                       H1.inicial=H1val,H3.inicial=H3val,K2.inicial=K2val,K3.inicial=K3val,
                                       beta1.inicial=beta1val,gamma1.inicial=gamma1val,beta2.inicial=beta2val,
                                       gamma2.inicial=gamma2val,selectodas.inicial=selectodas,
                                       echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))   
    if (percentval == 1)
    {
      pct <- c(gsub(" ", ",",gsub(", ", ",",tclvalue(percentilVariable))))
      pct1 <- as.numeric(unlist(strsplit(pct,",")))
      if ( is.na(pct1) || (sum(pct1<0.0)>0) || (sum(pct1>1.0)>0) || (sum(!is.numeric(pct1))>0) )
      {
        pct <- paste(seq(0.,1.,.25),collapse=",")
        Message(message=gettextRcmdr("Vector de percentiles invalido. Se utilizara vector por defecto."),
                type="warning")
      }
    } else pct <- NULL
    if (medrecval == 1)
    {
      rec <- as.numeric(tclvalue(trimVariable))
      if ( rec < .0 || rec > .5 || !is.numeric(rec) )
      {
        rec <- 0.05
        Message(message=gettextRcmdr("Proporcion de recorte invalida se utilizara valor por defecto."),
                type="warning")              
      } 
    } else rec <- NULL
    if (ACentval == 1)
    {
      propdat <- as.numeric(tclvalue(ACVariable))
      if ( propdat < .0 || propdat > 1. || !is.numeric(propdat) )
      {
        prop.dat <- 0.9
        Message(message=gettextRcmdr("Proporcion de datos invalida se utilizara valor por defecto."),
                type="warning")              
      }
    } else propdat <- NULL
    
    vars <- if (length(x) == 1) paste('"', x, '"', sep="") 
    else paste("c(", paste('"', x, '"', collapse=", ", sep=""), ")", sep="")
    if (length(x) == 1) variables <- paste(.BaseDatosActiva, "[", vars, "]", sep="")
    else variables <- paste(.BaseDatosActiva, "[,", vars, "]", sep="")
    stats <- paste("c(",
                   paste(c('"Media"', '"Mediana"', '"Moda"', '"MediaGeom"', '"Trimedia"', '"PromCuar"',
                           '"midR"', '"MediaRec"', '"Variancia"', '"DT"', '"DTGeom"', '"DesvMed"', '"CV"',
                           '"Rango"', '"IQR"', '"DesvCuar"', '"MAD"', '"CVR"', '"ACent"',
                           '"Min"', '"Max"', '"Q1"', '"Q2"', '"Q3"', '"Pct"',
                           '"H1"', '"H3"', '"K2"', '"K3"', '"beta1"', '"gamma1"', '"beta2"', '"gamma2"')
                         [c(mediaval, medianaval, modaval, mediageomval, trimediaval, promcuarval, midRval, 
                            medrecval, varianciaval, dtval, dtgeomval, desvmedval, CVval, rangoval, IQRval, 
                            desvcuarval, madval, CVRval, ACentval, minval, maxval, Q1val, Q2val, Q3val,
                            percentval, H1val, H3val, K2val, K3val, beta1val, gamma1val, beta2val,
                            gamma2val) == 1], 
                         collapse=", "), ")", sep="")
    if (percentval ==1)
      instruccion1 <- paste(".indices.num <- numindices(vars=",variables, ", statistics=", stats, if (!is.null(rec)){paste(",rec=",rec)},
                            if(!is.null(propdat)){paste(",propdat=",propdat)},",percentil=c(",pct,")",")", sep="")
    else
      instruccion1 <- paste(".indices.num <- numindices(vars=",variables, ", statistics=", stats, if (!is.null(rec)){paste(",rec=",rec)},
                            if(!is.null(propdat)){paste(",propdat=",propdat)},")", sep="")  
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(".indices.num  # Indices para variables cuantitativas ")  
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- "Indicadores descriptivos para variables cuantitativas"
      HTML(as.title(titulo),file=.archivo)
      j <- 6
      if ( selec>0 ){
        .TablaRes <- as.data.frame(.indices.num[[j]])
        rownames(.TablaRes) <- .indices.num[[1]]
        colnames(.TablaRes) <- c("Media","Mediana","Moda", "MediaGeom", "Trimedia","PromCuar","midR",
                                 paste("MediaRec(",.indices.num[[3]]*100,"%)",sep=''))[c("Media","Mediana","Moda", "MediaGeom",
                                                                                         "Trimedia","PromCuar","midR", "MediaRec") %in% .indices.num[[2]]]
        HTML("Indicadores de tendencia central: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }
      if ( selec2>0 ){
        .TablaRes <- as.data.frame(.indices.num[[j]])
        rownames(.TablaRes) <- .indices.num[[1]]
        colnames(.TablaRes) <- c("Variancia", "DT", "DTGeom", "DesvMed", "CV","Rango", "IQR",
                                 "DesvCuar", "MAD", "CVR",paste("ACent(",.indices.num[[4]]*100,"%)",sep=''))[c("Variancia",
                                                                                                               "DT", "DTGeom", "DesvMed", "CV","Rango", "IQR", "DesvCuar", "MAD", "CVR", "ACent") %in% .indices.num[[2]]]
        HTML("Indicadores de dispersion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      } 
      if ( selec3>0 ){
        .TablaRes <- as.data.frame(.indices.num[[j]])
        rownames(.TablaRes) <- .indices.num[[1]]
        colnames(.TablaRes) <- c("Min", "Max", "Q1", "Q2", "Q3",
                                 paste(.indices.num[[5]]*100,'%',sep=''))[c("Min", "Max", "Q1", 
                                                                            "Q2", "Q3",rep("Pct",length(.indices.num[[5]]))) %in% .indices.num[[2]]]
        HTML("Indicadores de posicion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }   
      if ( selec4>0 ){
        .TablaRes <- as.data.frame(.indices.num[[j]])
        rownames(.TablaRes) <- .indices.num[[1]]
        colnames(.TablaRes) <- c("H1", "H3", "K2", "K3","beta1", "gamma1", "beta2", "gamma2")[c("H1",
                                                                                                "H3", "K2", "K3","beta1", "gamma1", "beta2", "gamma2") %in% .indices.num[[2]]]
        HTML("Indicadores de forma: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }       
    }
    remove(.indices.num, envir=.GlobalEnv)
    closeDialog()
    if (echocodigo == 1) logger("remove(.indices.num)")
    tkfocus(CommanderWindow())    
  }
  OKCancelHelp(helpSubject="RcmdrPlugin.EACSPIR",reset="resumen.numericas",apply="resumen.numericas")
  tkgrid(getFrame(listaVar), sticky="nw")
  checkBoxes(statisticsTab,frame="tcFrame",boxes=c("media","mediana","mediageom","moda"),
             initialValues=c(dialog.valores$media.inicial,dialog.valores$mediana.inicial,
                             dialog.valores$mediageom.inicial,dialog.valores$moda.inicial),
             labels=gettextRcmdr(c("Media ","Mediana ","Media Geometrica ","Moda ")), 
             title = gettextRcmdr("Indices de tendencia central"))
  checkBoxes(statisticsTab,frame="tcRFrame",boxes=c("trimedia","promcuar","midR","medrec"),
             initialValues=c(dialog.valores$trimedia.inicial,dialog.valores$promcuar.inicial,
                             dialog.valores$midR.inicial,dialog.valores$medrec.inicial),
             labels=gettextRcmdr(c("Trimedia ","Promedio de cuartiles ","Rango medio ","Media Recortada ")), 
             title = gettextRcmdr(" "))  
  trimFrame <- tkframe(statisticsTab)
  trimVariable <- tclVar(dialog.valores$trim.inicial)
  trimField <- ttkentry(trimFrame, width="8", textvariable=trimVariable)
  tkgrid(labelRcmdr(trimFrame,text="  "))
  tkgrid(labelRcmdr(trimFrame,text="  "))
  tkgrid(labelRcmdr(trimFrame,text="  ")) 
  tkgrid(labelRcmdr(trimFrame,text="  ")) 
  tkgrid(labelRcmdr(trimFrame,text=gettextRcmdr("Proporcion datos recortados = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         trimField,sticky="w")
  tkgrid(tcFrame, labelRcmdr(statisticsTab,text="  "),tcRFrame, labelRcmdr(statisticsTab,text="  "), 
         trimFrame, sticky="nw")
  
  checkBoxes(statistics2Tab,frame="dispFrame",boxes=c("variancia","dt","dtgeom","desvmed","rango","IQR"),
             initialValues=c(dialog.valores$variancia.inicial,dialog.valores$dt.inicial,
                             dialog.valores$dtgeom.inicial,dialog.valores$desvmed.inicial,
                             dialog.valores$rango.inicial,dialog.valores$IQR.inicial),
             labels=gettextRcmdr(c("Variancia ", "Desv. Tipica ", "Desviacion Geometrica ", "Desv. Media ", "Amplitud ",
                                   "Amplitud intercuartil (IQR) ")), 
             title = gettextRcmdr("Indices de dispersion")) 
  checkBoxes(statistics2Tab,frame="dispRFrame",boxes=c("CV","desvcuar","mad","CVR","ACent"),
             initialValues=c(dialog.valores$CV.inicial,dialog.valores$desvcuar.inicial,
                             dialog.valores$mad.inicial,dialog.valores$CVR.inicial,
                             dialog.valores$ACent.inicial),
             labels=gettextRcmdr(c("Coeficiente Variacion ","Desviacion Cuartil ", 
                                   "Mediana Desviaciones Absolutas (MAD) ","Coeficiente Variacion Robusto ",
                                   "Desviacion centilica ")), 
             title = gettextRcmdr("  "))  
  ACFrame <- tkframe(statistics2Tab)
  ACVariable <- tclVar(dialog.valores$AC.inicial)
  ACField <- ttkentry(ACFrame, width="4", textvariable=ACVariable)
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text="  "))
  tkgrid(labelRcmdr(ACFrame,text="  ")) 
  tkgrid(labelRcmdr(ACFrame,text="  "))   
  tkgrid(labelRcmdr(ACFrame,text=gettextRcmdr("Proporcion datos utilizados = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         ACField,sticky="w")
  tkgrid(dispFrame, labelRcmdr(statistics2Tab,text="  "),
         dispRFrame, labelRcmdr(statistics2Tab,text="  "), ACFrame, sticky="nw")
  
  checkBoxes(statistics3Tab,frame="posicFrame",boxes=c("Q1","Q2","Q3"),
             initialValues=c(dialog.valores$Q1.inicial,dialog.valores$Q2.inicial,
                             dialog.valores$Q3.inicial),
             labels=gettextRcmdr(c("Primer Cuartil",
                                   "Segundo Cuartil ","Tercer Cuartil ")), 
             title = gettextRcmdr("Indices de posicion"))
  checkBoxes(statistics3Tab,frame="posicRFrame",boxes=c("min","max","percent"),
             initialValues=c(dialog.valores$min.inicial,dialog.valores$max.inicial,
                             dialog.valores$percent.inicial),
             labels=gettextRcmdr(c("Minimo ","Maximo ","Cuantilas ")), 
             title = gettextRcmdr("  "))  
  percentFrame <- tkframe(statistics3Tab)
  percentilVariable <- tclVar(dialog.valores$percentil.inicial)
  percentField <- ttkentry(percentFrame, width="15", textvariable=percentilVariable)
  tkgrid(labelRcmdr(percentFrame,text="  "))
  tkgrid(labelRcmdr(percentFrame,text="  "))
  tkgrid(labelRcmdr(percentFrame,text="  ")) 
  tkgrid(labelRcmdr(percentFrame,text=gettextRcmdr("Seleccione cuantilas = "),
                    fg=getRcmdr("title.color"),font="RcmdrTitleFont"),
         percentField,sticky="w")
  tkgrid(posicFrame, labelRcmdr(statistics3Tab,text="  "),
         posicRFrame, labelRcmdr(statistics3Tab,text="  "), percentFrame, sticky="nw")  
  
  checkBoxes(statistics4Tab,frame="formaFrame",boxes=c("H1","H3","beta1","gamma1"),
             initialValues=c(dialog.valores$H1.inicial,dialog.valores$H3.inicial,
                             dialog.valores$beta1.inicial,dialog.valores$gamma1.inicial),
             labels=gettextRcmdr(c("Coef. Asimetria H1 ","Coef. Asimetria H3 ",
                                   "Coef. Asimetria Pearson ","Coef. Asimetria Fisher ")), 
             title = gettextRcmdr("Indices de forma")) 
  checkBoxes(statistics4Tab,frame="formaRFrame",boxes=c("K2","K3","beta2","gamma2"),
             initialValues=c(dialog.valores$K2.inicial,dialog.valores$K3.inicial,
                             dialog.valores$beta2.inicial,dialog.valores$gamma2.inicial),
             labels=gettextRcmdr(c("Coef. Apuntamiento K2 ","Coef. Apuntamiento K3 ",
                                   "Coef. Apuntamiento Pearson ","Coef. Apuntamiento Fisher ")), 
             title = gettextRcmdr("   ")) 
  tkgrid(formaFrame, labelRcmdr(statistics4Tab,text="  "), formaRFrame, sticky="w")
  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml","selectodas"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial,
                             dialog.valores$selectodas.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ",
                                   "Calcular todos los indices ")), 
             title = gettextRcmdr("Opciones"))
  tkgrid(opsFrame, sticky="w")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','statistics2Tab',
                                                      'statistics3Tab','statistics4Tab','optionsTab'),
               tab.names=c("Datos","Tend. Central","Dispersion","Posicion","Forma","Opciones"))
  
}

histfun <- function(data,tscale,nsup,interval){
  variable <- unlist(strsplit(deparse(substitute(data)), "[$]"))[2]
  titulop <- paste("Histograma para ",variable, sep="")
  if (tscale == "frequency")
  {
    tituloy <- "Frecuencias absolutas"
    frecuencia  <- TRUE
  }
  if (tscale == "density")
  {
    tituloy <- "Densidades"
    frecuencia <- FALSE
  }
  titulox <- "Intervalos" 
  if (nsup == TRUE)
  {
    .xfit <- seq(min(data,na.rm=TRUE),max(data,na.rm=TRUE),length=1000)
    .yfit <- dnorm(x=.xfit,mean=mean(data,na.rm=TRUE),sd=sd(data,na.rm=TRUE))
    .h <- hist(data,breaks=interval,plot=FALSE)
    if (frecuencia==TRUE) .yfit <- .yfit*diff(.h$mids[1:2])*length(na.omit(data))
    if (frecuencia==TRUE) maxlim <- max(max(.h$counts),max(.yfit))
    else maxlim <- max(max(.h$density),max(.yfit))
    .h <- hist(data,freq=frecuencia,breaks=interval,main=titulop,ylab=tituloy,xlab=titulox,col='red',
               ylim=c(0,maxlim))
    lines(.xfit,.yfit,col='blue',lwd=2)
  }
  else
    .h <- hist(data,freq=frecuencia,breaks=interval,main=titulop,ylab=tituloy,xlab=titulox,col='red')
  box()
}

# Funcion creada a partir de la funcion Hist de John Fox incluida en R-Commander #
histograma <- function(){
  defecto <- list(x.inicial=NULL,escala.inicial="frequency",intervalos.inicial="<auto>",
                  normalsup.inicial="0",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("histograma",defecto)
  initializeDialog(title=gettextRcmdr("Histograma"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','optionsTab'))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    .BaseDatosActiva <- ActiveDataSet()
    interv <- tclvalue(intervalosVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    escala <- tclvalue(escalaVariable)
    normalsup <- tclvalue(normalsupVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    if (length(x) == 0){
      errorCondition(recall=histograma, message=gettextRcmdr("Debe escoger una variable."))
      return()
    }
    putDialog("histograma",list(x.inicial=x,escala.inicial=escala,intervalos.inicial=interv,
                                normalsup.inicial=normalsup,echo.inicial=echocodigo,
                                creahtml.inicial=creahtml,tab.inicial=tab)) 
    opts <- options(warn=-1)
    interv <- if (interv == gettextRcmdr("<auto>")) '"Sturges"'
    else as.numeric(interv)
    vars <- x 
    bd <- paste(.BaseDatosActiva, "$", vars, sep="")
    normalsup <- as.logical(as.numeric(normalsup))
    options(opts)
    instruccion1 <- paste("histfun(data=",bd, ", tscale='", escala, "', nsup=", normalsup,", interval=", interv, ")", sep="")
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    justDoIt(instruccion1)    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Histograma para variable ",x,sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("HistogramaR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()        
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="hist",reset="histograma",apply="histograma")
  listaVar <- variableListBox(dataTab, Numeric(),
                              title=gettextRcmdr("Variables (escoja una)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"numeric"))
  radioButtons(statisticsTab,name = "escala", buttons = c("frequency","density"), values = c("frequency","density"), 
               labels = gettextRcmdr(c("Frecuencias Absolutas","Densidades")),
               title = gettextRcmdr("Escala eje ordenadas"),
               initialValue = dialog.valores$escala.inicial)
  checkBoxes(statisticsTab,frame="histogramaFrame",boxes=c("normalsup"),
             initialValues=c(dialog.valores$normalsup.inicial),
             labels=gettextRcmdr(c("Superponer curva normal ")), 
             title = gettextRcmdr("Configuracion del histograma")) 
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  intervalosFrame <- tkframe(statisticsTab)
  intervalosVariable <- tclVar(dialog.valores$intervalos.inicial)
  intervalosField <- ttkentry(intervalosFrame, width="8", textvariable=intervalosVariable)  
  tkgrid(getFrame(listaVar), sticky="nw")
  tkgrid(escalaFrame, sticky="w")
  tkgrid(histogramaFrame, sticky="w") 
  tkgrid(labelRcmdr(intervalosFrame, text = gettextRcmdr("Numero de intervalos: ")), 
         intervalosField, sticky = "w")
  tkgrid(intervalosFrame, sticky = "w")
  tkgrid.configure(intervalosField, sticky = "e")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))       
}
           
bpfun <- function(data, intervinf, intervsup, ident){
  variable <- unlist(strsplit(deparse(substitute(data)), "[$]"))[2] 
  titulop <- paste("Diagrama de caja para ",variable, sep="")
  .bxp1 <- boxplot.stats(as.numeric(data),coef=intervinf)  
  .bxp2 <- boxplot.stats(as.numeric(data),coef=intervsup)
  boxplot(as.numeric(data),main=titulop,col='red',outpch=NA)
  .selec <- .bxp1$out %in% .bxp2$out
  .anom <- .bxp1$out
  .anom[.selec] <- NA
  points(rep(1, length(.anom)), .anom, pch = 1, col = 'blue')
  .extrem <- .bxp2$out
  points(rep(1, length(.extrem)), .extrem, pch = 8, col = 'red')
  if (ident == TRUE)
  {
    identify(rep(1,length(data)),as.numeric(data),rownames(data.frame(data)))
  } 
}

diagrama.caja.ord <- function (){
  defecto <- list(x.inicial=NULL,intervinf.inicial="1.5",intervsup.inicial="3.0",
                  identifica.inicial="0",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("diagrama.caja.ord",defecto)
  initializeDialog(title=gettextRcmdr("Diagrama de Caja"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','optionsTab'))
  listaVar <- variableListBox(dataTab, Factors(),
                              title=gettextRcmdr("Variables (escoja una)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"factor"))
  checkBoxes(statisticsTab,frame="cajaFrame",boxes=c("identifica"),
             initialValues=c(dialog.valores$identifica.inicial),
             labels=gettextRcmdr(c("Identificar sujetos (por numero de fila ocupado) ")), 
             title = gettextRcmdr("Configuracion del grafico")) 
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  intervalosFrame <- tkframe(statisticsTab)
  intervinfVariable <- tclVar(dialog.valores$intervinf.inicial)
  intervinfField <- ttkentry(intervalosFrame, width="6", textvariable=intervinfVariable)
  intervsupVariable <- tclVar(dialog.valores$intervsup.inicial)
  intervsupField <- ttkentry(intervalosFrame, width="6", textvariable=intervsupVariable)
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    .BaseDatosActiva <- ActiveDataSet()
    intervinf <- as.numeric(tclvalue(intervinfVariable))
    intervsup <- as.numeric(tclvalue(intervsupVariable))
    if ( is.na(intervinf) || (intervinf<0) || (!is.numeric(intervinf)) ||
           is.na(intervsup) || (intervsup<0) || (!is.numeric(intervsup)) ||
           (intervsup < intervinf) )
    {
      intervinf <- 1.5
      intervsup <- 3.0
      Message(message=gettextRcmdr("Coeficientes inferior y superior no validos. Se utilizara valores por defecto."),
              type="warning")
    }   
    opts <- options(warn=-1)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    identif <- as.logical(as.numeric(tclvalue(identificaVariable)))
    putDialog("diagrama.caja.ord",list(x.inicial=x,intervinf.inicial=intervinf,intervsup.inicial=intervsup,
                                       identifica.inicial=identif,echo.inicial=echocodigo,
                                       creahtml.inicial=creahtml,tab.inicial=tab)) 
    if (length(x) == 0){
      errorCondition(recall=diagrama.caja.ord, message=gettextRcmdr("Debe escoger una variable."))
      return()
    }
    bd <- paste(.BaseDatosActiva, "$", x, sep="")
    options(opts)       
    justDoIt(paste("cond <- !is.ordered(",paste(.BaseDatosActiva,"$",x,sep=""),")",sep=""))
    if (cond){
      errorCondition(recall=diagrama.caja.ord, message=gettextRcmdr(paste("Variable ",x, " no es ordinal.",sep='')))
      return()
    } 
    remove("cond", envir=.GlobalEnv)
    instruccion1 <- paste("bpfun(data=",bd, ", intervinf=", intervinf, ", intervsup=", intervsup,
                          ", ident=", identif, ")", sep="")
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    justDoIt(instruccion1)    
    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Diagrama de caja para variable ",x,sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("DiagramaCajaR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()        
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="boxplot",reset="diagrama.caja.ord",apply="diagrama.caja.ord")
  tkgrid(getFrame(listaVar), sticky="nw")
  tkgrid(cajaFrame, sticky="w") 
  tkgrid(labelRcmdr(intervalosFrame, text = gettextRcmdr("Coeficiente Limite inferior: ")), 
         intervinfField, sticky = "w")
  tkgrid(labelRcmdr(intervalosFrame, text = gettextRcmdr("Coeficiente Limite superior: ")), 
         intervsupField, sticky = "w")
  tkgrid(intervalosFrame, sticky = "w")
  tkgrid.configure(intervinfField, sticky = "e")
  tkgrid.configure(intervsupField, sticky = "e")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))    
}       

diagrama.caja <- function (){
  defecto <- list(x.inicial=NULL,intervinf.inicial="1.5",intervsup.inicial="3.0",
                  identifica.inicial="0",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("diagrama.caja",defecto)
  initializeDialog(title=gettextRcmdr("Diagrama de Caja"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','optionsTab'))
  listaVar <- variableListBox(dataTab, Numeric(),
                              title=gettextRcmdr("Variables (escoja una)"),
                              initialSelection=varPosn(dialog.valores$x.inicial,"numeric"))
  checkBoxes(statisticsTab,frame="cajaFrame",boxes=c("identifica"),
             initialValues=c(dialog.valores$identifica.inicial),
             labels=gettextRcmdr(c("Identificar sujetos (por numero de fila ocupado) ")), 
             title = gettextRcmdr("Configuracion del grafico")) 
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  intervalosFrame <- tkframe(statisticsTab)
  intervinfVariable <- tclVar(dialog.valores$intervinf.inicial)
  intervinfField <- ttkentry(intervalosFrame, width="6", textvariable=intervinfVariable)
  intervsupVariable <- tclVar(dialog.valores$intervsup.inicial)
  intervsupField <- ttkentry(intervalosFrame, width="6", textvariable=intervsupVariable)
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    x <- getSelection(listaVar)
    .BaseDatosActiva <- ActiveDataSet()
    intervinf <- as.numeric(tclvalue(intervinfVariable))
    intervsup <- as.numeric(tclvalue(intervsupVariable))
    if ( is.na(intervinf) || (intervinf<0) || (!is.numeric(intervinf)) ||
           is.na(intervsup) || (intervsup<0) || (!is.numeric(intervsup)) ||
           (intervsup < intervinf) )
    {
      intervinf <- 1.5
      intervsup <- 3.0
      Message(message=gettextRcmdr("Coeficientes inferior y superior no validos. Se utilizara valores por defecto."),
              type="warning")
    }   
    opts <- options(warn=-1)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    identif <- as.logical(as.numeric(tclvalue(identificaVariable)))
    putDialog("diagrama.caja",list(x.inicial=x,intervinf.inicial=intervinf,intervsup.inicial=intervsup,
                                   identifica.inicial=identif,echo.inicial=echocodigo,
                                   creahtml.inicial=creahtml,tab.inicial=tab)) 
    if (length(x) == 0){
      errorCondition(recall=diagrama.caja, message=gettextRcmdr("Debe escoger una variable."))
      return()
    }
    bd <- paste(.BaseDatosActiva, "$", x, sep="")
    options(opts)
    instruccion1 <- paste("bpfun(data=",bd, ", intervinf=", intervinf, ", intervsup=", intervsup,
                          ", ident=", identif, ")", sep="")
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    justDoIt(instruccion1)    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Diagrama de caja para variable ",x,sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("DiagramaCajaR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()        
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="boxplot",reset="diagrama.caja",apply="diagrama.caja")
  tkgrid(getFrame(listaVar), sticky="nw")
  tkgrid(cajaFrame, sticky="w") 
  tkgrid(labelRcmdr(intervalosFrame, text = gettextRcmdr("Coeficiente Limite inferior: ")), 
         intervinfField, sticky = "w")
  tkgrid(labelRcmdr(intervalosFrame, text = gettextRcmdr("Coeficiente Limite superior: ")), 
         intervsupField, sticky = "w")
  tkgrid(intervalosFrame, sticky = "w")
  tkgrid.configure(intervinfField, sticky = "e")
  tkgrid.configure(intervsupField, sticky = "e")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

# Resumen descriptivo bivariante para variables categoricas: indicadores y graficos #

bivcat <- function(rowvar,colvar,statistics=NULL,tables=NULL,subset=NULL,pcttable=NULL,data){
  res <- list()
  res[[1]] <- rowvar
  res[[2]] <- colvar
  res[[3]] <- statistics
  res[[4]] <- tables
  res[[5]] <- pcttable
  i <- 6
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))
  res[[i]] <- .Tabla
  MArow <- sum(c("JiCuadrado", "Phi", "CoefCont", "Sakoda", "Chuprov", "VCramer") %in% statistics) +
    + ("Yule" %in% statistics)*3
  MEProw <- sum(c("Lambda","Theil") %in% statistics)*3 + ("Tau" %in% statistics)*2
  .tablaMA <- array(dim=c(MArow,1))
  .tablaMEP <- array(dim=c(MEProw,1))
  if ("Porcentajes" %in% tables){
    if (pcttable == "fila")
    {
      .porcentajes <- rowPercents(.Tabla)
    }
    if (pcttable == "columna")
    {
      .porcentajes <- colPercents(.Tabla)
    }
    if (pcttable == "total")
    {
      .porcentajes <- totPercents(.Tabla)
    }
    i <- i + 1
    res[[i]] <- round(.porcentajes,2)
  }
  if ("FrecEsperadas" %in% tables)
  {
    .esperadas <- chisq.test(.Tabla, correct=FALSE)$expected
    mensAviso <- NULL
    if (0 < (emq1 <- sum(.esperadas < 1))) mensAviso <- paste(emq1,
                                                              gettextRcmdr("frecuencias esperadas menores que 1"))
    if (0 < (emq5 <- sum(.esperadas < 5))) mensAviso <- paste(mensAviso, "\n", emq5,
                                                              gettextRcmdr(" frecuencias esperadas menores que 5"), sep="")
    if (!is.null(mensAviso)) Message(message=mensAviso,
                                     type="warning")
    i <- i + 1
    res[[i]] <- round(.esperadas,2)
  }
  if ("JiComponentes" %in% tables)
  {
    .Componentes <- round(chisq.test(.Tabla, correct=FALSE)$residuals^2,2)
    i <- i + 1
    res[[i]] <- round(.Componentes,2)
  }  
  if (MArow > 0){
    j <- 0
    if ( (sum(c("JiCuadrado", "Phi", "CoefCont", "Sakoda", "Chuprov", "VCramer") %in% statistics) +
            "JiComponentes" %in% tables) > 0 )
    {
      .Jicuadrado <- chisq.test(.Tabla, correct=FALSE)$statistic
      names(.Jicuadrado)<-NULL
    }
    if ("JiCuadrado" %in% statistics)
    {
      j <- j + 1
      .tablaMA[j] <- .Jicuadrado
    }
    if ("Phi" %in% statistics)
    {
      .phi <- sqrt(.Jicuadrado/sum(.Tabla))
      j <- j + 1
      .tablaMA[j] <- .phi
    }
    if ("CoefCont" %in% statistics)
    {
      .Coef.Contingencia <- sqrt(.Jicuadrado/(sum(.Tabla)+.Jicuadrado))
      j <- j + 1
      .tablaMA[j] <- .Coef.Contingencia
    }
    if ("Sakoda" %in% statistics)
    {
      .sakoda<- sqrt(min(dim(.Tabla))*.Jicuadrado/((min(dim(.Tabla))-1)*(sum(.Tabla)+.Jicuadrado)))
      j <- j + 1
      .tablaMA[j]<-.sakoda
    }  
    if ("Chuprov" %in% statistics)
    {
      .chuprov <- sqrt(.Jicuadrado/(sum(.Tabla)*(dim(.Tabla)[1]-1)*(dim(.Tabla)[2]-1)))
      j <- j + 1
      .tablaMA[j]<-.chuprov
    }
    if ("VCramer" %in% statistics)
    {
      .VCramer <- sqrt(.Jicuadrado/((min(dim(.Tabla))-1)*sum(.Tabla)))
      j <- j + 1
      .tablaMA[j]<-.VCramer
    }
    if ("Yule" %in% statistics)
    {
      if ( dim(.Tabla)[1] != 2 || dim(.Tabla)[2] != 2){
        Message(message=gettextRcmdr("La tabla de contingencia no es de tamano 2x2: \n Indices de Yule no se calcularan"),
                type="warning")
        return()
      }
      .a <- .Tabla[1,1]
      .b <- .Tabla[1,2]
      .c <- .Tabla[2,1]
      .d <- .Tabla[2,2]
      .Q <- (.a*.d-.b*.c)/(.a*.d+.b*.c)
      .Y <- (sqrt(.a*.d)-sqrt(.b*.c))/(sqrt(.a*.d)+sqrt(.b*.c))
      .V <- (.a*.d-.b*.c)/((.a+.b)*(.a+.c)*(.b+.d)*(.c+.d))
      j <- j + 1
      .tablaMA[j:(j+2)] <- c(.Q,.Y,.V)
    }
    i <- i + 1
    res[[i]] <- round(.tablaMA,2)
  }
  if (MEProw > 0){
    j <- 0
    if ("Lambda" %in% statistics)
    {
      lambda.a.b <- (sum(apply(.Tabla,2,max)/sum(.Tabla)) - max(rowSums(.Tabla))/sum(.Tabla))/(1 - 
                                                                                                 max(rowSums(.Tabla))/sum(.Tabla))
      lambda.b.a <- (sum(apply(.Tabla,1,max)/sum(.Tabla)) - max(colSums(.Tabla))/sum(.Tabla))/(1 - 
                                                                                                 max(colSums(.Tabla))/sum(.Tabla))
      .lambda <- (lambda.a.b + lambda.b.a)/2
      j <- j + 1
      .tablaMEP[j:(j+2)] <- c(lambda.a.b,lambda.b.a,.lambda)
      j <- j+2
    }    
    if ("Tau" %in% statistics)
    {
      tau.a.b <- (sum((.Tabla/sum(.Tabla))^2/matrix(colSums(.Tabla)[col(.Tabla)]/sum(.Tabla),
                                                    nrow=nrow(.Tabla)))-sum((rowSums(.Tabla)/sum(.Tabla))^2))/(1-sum((rowSums(.Tabla)/sum(.Tabla))^2))
      tau.b.a <- (sum((.Tabla/sum(.Tabla))^2/matrix(rowSums(.Tabla)[row(.Tabla)]/sum(.Tabla),
                                                    nrow=nrow(.Tabla)))-sum((colSums(.Tabla)/sum(.Tabla))^2))/(1-sum((colSums(.Tabla)/sum(.Tabla))^2))
      j <- j+1
      .tablaMEP[j:(j+1)] <- c(tau.a.b,tau.b.a)
      j <- j+1
      
    }
    if ("Theil" %in% statistics)
    {
      H.a.b <- -sum(.Tabla/sum(.Tabla)*log(.Tabla/sum(.Tabla)),na.rm=TRUE)
      H.a <- -sum(rowSums(.Tabla)/sum(.Tabla)*log(rowSums(.Tabla)/sum(.Tabla)),na.rm=TRUE)
      H.b <- -sum(colSums(.Tabla)/sum(.Tabla)*log(colSums(.Tabla)/sum(.Tabla)),na.rm=TRUE)
      theil.a.b <- (H.a + H.b - H.a.b)/H.a
      theil.b.a <- (H.a + H.b - H.a.b)/H.b
      .theil <- 2*(H.a + H.b - H.a.b)/(H.a+H.b)
      j <- j+1
      .tablaMEP[j:(j+2)] <- c(theil.a.b,theil.b.a,.theil)
    }       
    i <- i + 1
    res[[i]] <- round(.tablaMEP,2)    
  }
  class(res) <- "bivcat"
  res  
}

print.bivcat <- function(x,...){
  cat("Tabla de contingencia para ",x[[1]], " y ",x[[2]],": \n\n",sep='')
  print(x[[6]])
  cat("\n\n")
  j<-7
  if ("Porcentajes" %in% x[[4]]){
    msg <- switch(x[[5]],
                  'total'='Tabla porcentajes respecto al total: ',
                  'fila'='Tabla porcentajes respecto a marginales fila: ',
                  'columna'='Tabla porcentajes respecto a marginales columna: ')
    cat(paste(msg,"\n\n"))
    print(x[[j]])
    cat("\n\n")
    j <- j + 1
  }
  if ("FrecEsperadas" %in% x[[4]]){
    cat(paste("Frecuencias Esperadas: \n\n"))
    print(x[[j]])
    cat("\n\n")
    j <- j + 1
  }
  if ("JiComponentes" %in% x[[4]]){
    cat(paste("Descomposicion del estadistico Ji cuadrado: \n\n"))
    print(x[[j]])
    cat("\n\n")
    j <- j + 1
  }  
  if (sum(c("JiCuadrado", "Phi", "CoefCont", "Sakoda", "Chuprov", "VCramer", "Yule") %in% x[[3]]) > 0){
    .tablaMA <- as.data.frame(x[[j]])
    rownames(.tablaMA) <- c("Ji Cuadrado", "Phi de Pearson", "Coef. Contingencia de Pearson", "Transf. Sakoda", 
                            "Coef. Chuprov", "V de Cramer",c("Q de Yule","Y de Yule","V de Yule"))[
                              c("JiCuadrado", "Phi", "CoefCont", "Sakoda", "Chuprov", "VCramer",rep("Yule",3)) %in% x[[3]]]
    colnames(.tablaMA) <- "Valores"
    cat("Coeficientes de asociacion: \n\n")
    print(.tablaMA)  
    cat("\n\n")
    j <- j+1
  }
  if (sum(c("Lambda", "Tau", "Theil") %in% x[[3]]) > 0){
    .tablaMEP <- as.data.frame(x[[j]])
    rownames(.tablaMEP) <-  c(c("Lambda A/B", "Lambda B/A","Lambda (simetrica)"), c("Tau A/B", 
                                                                                    "Tau B/A"), c("Theil A/B","Theil B/A","Theil (simetrica)"))[
                                                                                      c(rep("Lambda",3), rep("Tau",2), rep("Theil",3)) %in% x[[3]]]
    colnames(.tablaMEP) <- "Valores"
    cat("Medidas del error de prediccion: \n\n")
    print(.tablaMEP)  
    cat("\n\n")
    j <- j+1
  }
  invisible(x)  
}

bivariante.categoricas <- function(){
  defecto <- list(fila.inicial=NULL,columna.inicial=NULL,porcentajes.inicial="ninguno",
                  frecEsp.inicial="0",jicuadrado.inicial="0",
                  phiPearson.inicial="0",sakoda.inicial="0",VCramer.inicial="0",
                  jicomponentes.inicial="0",contingPearson.inicial="0",chuprov.inicial="0",
                  yule.inicial="0",lambda.inicial="0",tau.inicial="0",theil.inicial="0",
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("bivariante.categoricas",defecto) 
  initializeDialog(title=gettextRcmdr("Descripcion tablas de contingencia"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','statistics2Tab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  filaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                             initialSelection=varPosn(dialog.valores$fila.inicial,"factor")) 
  columnaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                                initialSelection=varPosn(dialog.valores$columna.inicial,"factor"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  radioButtons(statisticsTab,name = "porcentajes", buttons = c("fila","columna", "total", "ninguno"),
               values = c("fila", "columna","total", "ninguno"), 
               labels = gettextRcmdr(c("Porcentajes respecto a marginales fila",
                                       "Porcentajes respecto a marginales columna","Porcentajes respecto al total",
                                       "Ningun porcentaje")),
               title = gettextRcmdr("Calcular Porcentajes"),
               initialValue = dialog.valores$porcentajes.inicial)
  checkBoxes(statistics2Tab,frame="esperadasFrame",boxes=c("frecEsp"),
             initialValues=c(dialog.valores$frecEsp.inicial),
             labels=gettextRcmdr(c("Calcular frecuencias esperadas ")), 
             title = gettextRcmdr("Frecuencias Esperadas")) 
  checkBoxes(statistics2Tab,frame="statsFrame",boxes=c("jicuadrado","phiPearson","sakoda","VCramer"),
             initialValues=c(dialog.valores$jicuadrado.inicial,dialog.valores$phiPearson.inicial,
                             dialog.valores$sakoda.inicial,dialog.valores$VCramer.inicial),
             labels=gettextRcmdr(c("Ji Cuadrado ","Phi de Pearson ","Transformacion de Sakoda ",
                                   "Coef. Contingencia Cramer ")), 
             title = gettextRcmdr("Coeficientes de Asociacion"))
  rightFrame <- tkframe(statistics2Tab)
  checkBoxes(statistics2Tab,frame="rightFrame",boxes=c("jiComponentes","contingPearson","chuprov","yule"),
             initialValues=c(dialog.valores$jicomponentes.inicial,dialog.valores$contingPearson.inicial,
                             dialog.valores$chuprov.inicial,dialog.valores$yule.inicial),
             labels=gettextRcmdr(c("Descomposicion Ji cuadrado de Pearson ","Coef. Contingencia Pearson ",
                                   "Coef. Contingencia Chuprov ","Coeficientes Yule (Tablas 2x2) ")), 
             title = gettextRcmdr(" "))
  checkBoxes(statistics2Tab,frame="errorpredFrame",boxes=c("lambda","tau","theil"),
             initialValues=c(dialog.valores$lambda.inicial,dialog.valores$tau.inicial,
                             dialog.valores$theil.inicial),
             labels=gettextRcmdr(c("Lambda Goodman-Kruskal ","Tau Goodman-Kruskal ",
                                   "Coef. Incertidumbre Theil ")), 
             title = gettextRcmdr("Medidas Error de Prediccion"))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))    
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    fila <- getSelection(filaVar)
    columna <- getSelection(columnaVar)
    if (length(fila) == 0 || length(columna) == 0){
      errorCondition(recall=bivariante.categoricas, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (fila == columna) {
      errorCondition(recall=bivariante.categoricas, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    porcentajes <- as.character(tclvalue(porcentajesVariable))
    esperadas <- tclvalue(frecEspVariable)
    jicuadrado <- tclvalue(jicuadradoVariable)
    jicomponentes <- tclvalue(jiComponentesVariable)
    phival <- tclvalue(phiPearsonVariable)
    contingval <- tclvalue(contingPearsonVariable)
    sakodaval <- tclvalue(sakodaVariable)
    chuprovval <- tclvalue(chuprovVariable)
    VCramerval <- tclvalue(VCramerVariable)
    yuleval <- tclvalue(yuleVariable)
    lambdaval <- tclvalue(lambdaVariable)
    tauval <- tclvalue(tauVariable)
    theilval <- tclvalue(theilVariable)
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("bivariante.categoricas",list(fila.inicial=fila,columna.inicial=columna,
                                            porcentajes.inicial=porcentajes,frecEsp.inicial=esperadas,
                                            jicuadrado.inicial=jicuadrado,jicomponentes.inicial=jicomponentes,
                                            phiPearson.inicial=phival,contingPearson.inicial=contingval,
                                            sakoda.inicial=sakodaval,chuprov.inicial=chuprovval,
                                            VCramer.inicial=VCramerval,yule.inicial=yuleval,lambda.inicial=lambdaval,
                                            tau.inicial=tauval, theil.inicial=theilval,subconjunto.inicial=subconjunto,
                                            echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))     
    selec <- as.numeric(jicuadrado) + as.numeric(phival) + as.numeric(contingval) + 
      as.numeric(sakodaval) + as.numeric(chuprovval) + as.numeric(VCramerval) +
      as.numeric(yuleval)*3
    selec2 <- as.numeric(lambdaval)*3 + as.numeric(tauval)*2 + as.numeric(theilval)*3 
    if (selec+selec2 >0)
    stats <- paste("c(",
                   paste(c('"JiCuadrado"', '"Phi"', '"CoefCont"', '"Sakoda"', '"Chuprov"',
                           '"VCramer"', '"Yule"', '"Lambda"', '"Tau"', '"Theil"')
                         [c(jicuadrado,phival,contingval,sakodaval,chuprovval,VCramerval,yuleval,
                            lambdaval,tauval,theilval) == 1], 
                         collapse=", "), ")", sep="") 
    else stats <- 'NULL'
    if ((porcentajes!='ninguno') || (as.numeric(esperadas)+as.numeric(jicomponentes))>0)
    tabs <- paste("c(",
                  paste(c('"Porcentajes"', '"FrecEsperadas"', '"JiComponentes"')
                        [c((porcentajes!="ninguno"),c(esperadas,jicomponentes) == 1)], 
                        collapse=", "), ")", sep="")
    else tabs <- 'NULL'
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste(".indices.bc <- bivcat(rowvar='",fila, "', colvar='",columna, "', statistics=", stats,
                            ", tables=",tabs, if(porcentajes!="ninguno"){paste(", pcttable='",porcentajes,"'",sep='')},
                            ", data='",ActiveDataSet(),"')", sep="")
    else {
      instruccion1 <- paste(".indices.bc <- bivcat(rowvar='",fila, "', colvar='",columna, "', statistics=", stats,
                            ", tables=",tabs, if(porcentajes!="ninguno"){paste(", pcttable='",porcentajes,"'",sep='')},
                            ", subset='",subconjunto, "', data='",ActiveDataSet(),"')", sep="")
    }
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(paste(".indices.bc # Descripcion bivariante de: ",fila, 
                       " y ", columna, sep=""))    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")          
      titulo <- paste("Descripcion bivariante de datos categoricos: ",fila, 
                      " y ", columna, sep="")     
      HTML(as.title(titulo),file=.archivo)
      HTML('Tabla de contigencia: ', file=.archivo)
      HTML(.indices.bc[[6]], file=.archivo)
      j <- 7
      if ("Porcentajes" %in% .indices.bc[[4]]){
        msg <- switch(.indices.bc[[5]],
                      'total'='Tabla porcentajes respecto al total: ',
                      'fila'='Tabla porcentajes respecto a marginales fila: ',
                      'columna'='Tabla porcentajes respecto a marginales columna: ')
        HTML(msg, file=.archivo)        
        HTML(.indices.bc[[j]], file=.archivo)
        j <- j + 1
      }
      if ("FrecEsperadas" %in% .indices.bc[[4]]){
        msg <- 'Frecuencias esperadas:'
        HTML(msg, file=.archivo)        
        HTML(.indices.bc[[j]], file=.archivo)
        j <- j + 1
      }
      if ("JiComponentes" %in% .indices.bc[[4]]){
        msg <- 'Descomposicion del estadistico Ji cuadrado::'
        HTML(msg, file=.archivo)        
        HTML(.indices.bc[[j]], file=.archivo)
        j <- j + 1
      }        
      if ( selec>0 ){
        .TablaRes <- as.data.frame(.indices.bc[[j]])
        rownames(.TablaRes) <- c("Ji Cuadrado", "Phi de Pearson", "Coef. Contingencia de Pearson", "Transf. Sakoda", 
                                 "Coef. Chuprov", "V de Cramer",c("Q de Yule","Y de Yule","V de Yule"))[
                                   c("JiCuadrado", "Phi", "CoefCont", "Sakoda", "Chuprov", "VCramer",rep("Yule",3)) %in% .indices.bc[[3]]]
        colnames(.TablaRes) <- "Valores"
        HTML("Coeficientes de asociacion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      } 
      if ( selec2>0 ){
        .TablaRes <- as.data.frame(.indices.bc[[j]])
        rownames(.TablaRes) <-  c(c("Lambda A/B", "Lambda B/A","Lambda (simetrica)"), c("Tau A/B", 
                                                                                        "Tau B/A"), c("Theil A/B","Theil B/A","Theil (simetrica)"))[
                                                                                          c(rep("Lambda",3), rep("Tau",2), rep("Theil",3)) %in% .indices.bc[[3]]]
        colnames(.TablaRes) <- "Valores"
        HTML("Medidas del error de prediccion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }
    }
    remove(.indices.bc, envir=.GlobalEnv) 
    closeDialog()
    tkfocus(CommanderWindow())    
  }
  OKCancelHelp(helpSubject="xtabs",reset="bivariante.categoricas",apply="bivariante.categoricas")
  tkgrid(getFrame(filaVar), labelRcmdr(variablesFrame, text="    "), getFrame(columnaVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(porcentajesFrame, sticky="w")  
  tkgrid(esperadasFrame, sticky="w")
  tkgrid(statsFrame, rightFrame, sticky="w") 
  tkgrid(errorpredFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab',
                                                      'statistics2Tab','optionsTab'),
               tab.names=c("Datos","Tablas","Estadisticos","Opciones")) 
}

barragfun <- function(rowvar,colvar,tabVariable='ni',subset=NULL,data){
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))  
  if (tabVariable == "fi")
  {  
    .Tabla <- .Tabla/sum(.Tabla)
  }
  titulo <- paste("Barras agrupadas para ",rowvar," y ",colvar, sep="")
  tituloy <- if (tabVariable == "ni") "Frecuencias Absolutas"
  else "Frecuencias Relativas"
  eval(parse(text=paste('barplot(.Tabla,beside=TRUE,main=titulo,ylab=tituloy,xlab="',colvar,
                        '",ylim=c(0,max(.Tabla)*1.05),col=heat.colors(length(levels(',data,'$',rowvar,
                        '))),legend.text=TRUE,args.legend=list(x="topright",title="',rowvar,'")',')',sep='')))
  box()   
}

barras.agrupadas <- function(){
  defecto <- list(fila.inicial=NULL,columna.inicial=NULL,tabla.inicial="ni",
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("barras.agrupadas",defecto) 
  initializeDialog(title=gettextRcmdr("Diagrama de barras agrupadas"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  filaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                             initialSelection=varPosn(dialog.valores$fila.inicial,"factor")) 
  columnaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                                initialSelection=varPosn(dialog.valores$columna.inicial,"factor"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  radioButtons(statisticsTab,name = "tabla", buttons = c("niBoton","fiBoton"), values = c("ni","fi"), 
               labels = gettextRcmdr(c("Frecuencias Absolutas", "Frecuencias Relativas")),
               title = gettextRcmdr("Tabla basada en:"),
               initialValue = dialog.valores$tabla.inicial)
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))     
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    fila <- getSelection(filaVar)
    columna <- getSelection(columnaVar)
    if (length(fila) == 0 || length(columna) == 0){
      errorCondition(recall=barras.agrupadas, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (fila == columna) {
      errorCondition(recall=barras.agrupadas, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    tabla <- as.character(tclvalue(tablaVariable))
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("barras.agrupadas",list(fila.inicial=fila,columna.inicial=columna,
                                      tabla.inicial=tabla,subconjunto.inicial=subconjunto,
                                      echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))         
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste("barragfun(rowvar='",fila, "', colvar='",columna, "', tabVariable='", tabla,
                            "', data='",ActiveDataSet(),"')", sep="")
    else {
      instruccion1 <- paste("barragfun(rowvar='",fila, "', colvar='",columna, "', tabVariable='", tabla,
                            "', subset='",subconjunto, "', data='",ActiveDataSet(),"')", sep="")
    }
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Diagrama de barras agrupadas para : ",fila, 
                      " y ", columna, sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("BarrasAgrupadasR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)      
    }
    closeDialog()        
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="barplot",reset="barras.agrupadas",apply="barras.agrupadas")
  tkgrid(getFrame(filaVar), labelRcmdr(variablesFrame, text="    "), getFrame(columnaVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(tablaFrame, sticky="w") 
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))        
}

mosaicofun <- function(rowvar,colvar,subset=NULL,data){
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))  
  titulo <- paste("Grafico de mosaico para ",rowvar," y ",colvar, sep="")
  mosaicplot(.Tabla,main=titulo,col=1:dim(.Tabla)[2])
  box()
}

grafico.mosaico <- function(){
  defecto <- list(fila.inicial=NULL,columna.inicial=NULL,
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("grafico.mosaico",defecto) 
  initializeDialog(title=gettextRcmdr("Grafico de Mosaico"),use.tabs=TRUE,
                   tabs=c('dataTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  filaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                             initialSelection=varPosn(dialog.valores$fila.inicial,"factor")) 
  columnaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                                initialSelection=varPosn(dialog.valores$columna.inicial,"factor"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))       
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    fila <- getSelection(filaVar)
    columna <- getSelection(columnaVar)
    if (length(fila) == 0 || length(columna) == 0){
      errorCondition(recall=grafico.mosaico, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (fila == columna) {
      errorCondition(recall=grafico.mosaico, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("grafico.mosaico",list(fila.inicial=fila,columna.inicial=columna,
                                     subconjunto.inicial=subconjunto,echo.inicial=echocodigo,
                                     creahtml.inicial=creahtml,tab.inicial=tab))         
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste("mosaicofun(rowvar='",fila, "', colvar='",columna, 
                            "', data='",ActiveDataSet(),"')", sep="")
    else {
      instruccion1 <- paste("mosaicofun(rowvar='",fila, "', colvar='",columna,
                            "', subset='",subconjunto, "', data='",ActiveDataSet(),"')", sep="")
    }
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Grafico de mosaico para : ",fila, 
                      " y ", columna, sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("GraficoMosaicoR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()        
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="mosaicplot",reset="grafico.mosaico",apply="grafico.mosaico")
  tkgrid(getFrame(filaVar), labelRcmdr(variablesFrame, text="    "), getFrame(columnaVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Datos","Opciones"))         
}

# Resumen descriptivo bivariante para variables ordinales: indicadores y graficos #

bivord <- function(rowvar,colvar,statistics=NULL,tables=FALSE,subset=NULL,pcttable=NULL,data){
  res <- list()
  res[[1]] <- rowvar
  res[[2]] <- colvar
  res[[3]] <- statistics
  res[[4]] <- pcttable
  i <- 5
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))
  res[[i]] <- .Tabla
  filas <- row(.Tabla)
  columnas <- col(.Tabla)
  n <- sum(.Tabla)
  q <- min(dim(.Tabla))
  C <- sum(.Tabla * mapply(function(f, c){sum(.Tabla[(filas > f) &
                                                       (columnas > c)])}, f = filas, c = columnas))
  D <- sum(.Tabla * mapply(function(f, c){sum(.Tabla[(filas > f) &
                                                       (columnas < c)])}, f = filas, c = columnas))
  E.X <- (sum(apply(.Tabla,1,sum)^2)-n)/2
  E.Y <- (sum(apply(.Tabla,2,sum)^2)-n)/2
  E.XY <- (sum(.Tabla^2)-n)/2
  tablares <- matrix(c(C,D,E.X,E.Y,E.XY),nrow=5,ncol=1)
  i <- i + 1
  res[[i]] <- tablares
  if (tables == TRUE){
    if (pcttable == "fila")
    {
      .porcentajes <- rowPercents(.Tabla)
    }
    if (pcttable == "columna")
    {
      .porcentajes <- colPercents(.Tabla)
    }
    if (pcttable == "total")
    {
      .porcentajes <- totPercents(.Tabla)
    }
    i <- i + 1
    res[[i]] <- round(.porcentajes,2)
  }
  MArow <- ("Gamma" %in% statistics) +  ("Tau" %in% statistics)*3 + ("Sommers" %in% statistics)*3 + 
    ("Wilson" %in% statistics)
  .tablaMA <- array(dim=c(MArow,1))  
  if (MArow > 0){
    j <- 0
    if ("Gamma" %in% statistics)
    {
      gammaGK <- (C-D)/(C+D)
      j <- j + 1
      .tablaMA[j] <- gammaGK
    }
    if ("Tau" %in% statistics)
    {
      tau.a <- (C-D)/choose(n,2)
      tau.b <- (C-D)/sqrt((C+D+E.X-E.XY)*(C+D+E.Y-E.XY))
      tau.c <- 2*q*(C-D)/(n^2*(q-1))
      j <- j + 1
      .tablaMA[j:(j+2)] <- c(tau.a,tau.b,tau.c)
      j <- j+2
    }
    if ("Sommers" %in% statistics)
    {
      sommers.x.y <- (C-D)/(C+D+E.X-E.XY)
      sommers.y.x <- (C-D)/(C+D+E.Y-E.XY)
      .sommers <- (C-D)/(C+D+(E.X+E.Y)/2-E.XY)
      j <- j + 1
      .tablaMA[j:(j+2)] <- c(sommers.x.y,sommers.y.x,.sommers)
      j <- j+2
    }  
    if ("Wilson" %in% statistics)
    {
      .wilson <- 2*(C-D)/(choose(n,2)-E.XY)
      j <- j + 1
      .tablaMA[j] <- .wilson
    }
    i <- i + 1
    res[[i]] <- round(.tablaMA,2)
  }
  class(res) <- "bivord"
  res  
}

print.bivord <- function(x,...){
  cat("Tabla de contingencia para ",x[[1]], " y ",x[[2]],": \n\n",sep='')
  print(x[[5]])
  cat("\n\n")
  cat("Descripcion tabla de contingencia: \n\n")
  tabladesc <- x[[6]]
  rownames(tabladesc) <- c("Pares concordantes","Pares discordantes",
                           paste("Empates",x[[1]]),paste("Empates",x[[2]]),
                           paste("Empates",x[[1]],"y",x[[2]]))
  colnames(tabladesc) <- "Valores"
  print(tabladesc)
  cat("\n\n")  
  j<-7
  if (!is.null(x[[4]])){
    msg <- switch(x[[4]],
                  'total'='Tabla porcentajes respecto al total: ',
                  'fila'='Tabla porcentajes respecto a marginales fila: ',
                  'columna'='Tabla porcentajes respecto a marginales columna: ')
    cat(paste(msg,"\n\n"))
    print(x[[j]])
    cat("\n\n")
    j <- j + 1
  }
  if (sum(c("Gamma", "Tau", "Sommers", "Wilson") %in% x[[3]]) > 0){
    .tablaMA <- as.data.frame(x[[j]])
    rownames(.tablaMA) <- c("Gamma de Goodman-Kruskal", c("Tau a de Kendall","Tau b de Kendall","Tau c de Kendall"),
                            c("d de Sommers X/Y", "d de Sommers Y/X","d de Sommers (simetrica)"),"e de Wilson")[
                              c("Gamma",rep("Tau",3),rep("Sommers",3),"Wilson") %in% x[[3]]]
    colnames(.tablaMA) <- "Valores"
    cat("Coeficientes de asociacion: \n\n")
    print(.tablaMA)  
    cat("\n\n")
    j <- j+1
  }
  invisible(x)  
}

bivariante.ordinales <- function(){
  defecto <- list(fila.inicial=NULL,columna.inicial=NULL,
                  porcentajes.inicial="ninguno",
                  gamma.inicial="0",tau.inicial="0",sommers.inicial="0",wilson.inicial="0",
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("bivariante.ordinales",defecto) 
  initializeDialog(title=gettextRcmdr("Descripcion bivariante datos ordinales"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','statistics2Tab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  
  filaVar <- variableListBox(variablesFrame, c(Factors(),Numeric()), title=gettextRcmdr("Variables (escoja una)"),
                             initialSelection=varPosn(dialog.valores$fila.inicial,"factor")) 
  columnaVar <- variableListBox(variablesFrame, c(Factors(),Numeric()), title=gettextRcmdr("Variables (escoja una)"),
                                initialSelection=varPosn(dialog.valores$columna.inicial,"factor"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  radioButtons(statisticsTab,name = "porcentajes", buttons = c("fila","columna", "total", "ninguno"),
               values = c("fila", "columna","total", "ninguno"), 
               labels = gettextRcmdr(c("Porcentajes respecto a marginales fila",
                                       "Porcentajes respecto a marginales columna","Porcentajes respecto al total",
                                       "Ningun porcentaje")),
               title = gettextRcmdr("Calcular Porcentajes"),
               initialValue = dialog.valores$porcentajes.inicial)
  checkBoxes(statistics2Tab,frame="statsFrame",boxes=c("gamma","tau","sommers","wilson"),
             initialValues=c(dialog.valores$gamma.inicial,dialog.valores$tau.inicial,
                             dialog.valores$sommers.inicial,dialog.valores$wilson.inicial),
             labels=gettextRcmdr(c("Coef. Gamma de Goodman-Kruskal ","Coefs. Tau de Kendall ","Indices d de Sommers ",
                                   "Indice e de Wilson ")), 
             title = gettextRcmdr("Coeficientes de Asociacion"))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))    
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    fila <- getSelection(filaVar)
    columna <- getSelection(columnaVar)
    if (length(fila) == 0 || length(columna) == 0){
      errorCondition(recall=bivariante.ordinales, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (fila == columna) {
      errorCondition(recall=bivariante.ordinales, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    justDoIt(paste("cond1 <- !is.ordered(",paste(ActiveDataSet(),"$",fila,sep=""),")",sep=""))
    justDoIt(paste("cond2 <- !is.numeric(",paste(ActiveDataSet(),"$",fila,sep=""),")",sep=""))
    justDoIt(paste("cond3 <- !is.ordered(",paste(ActiveDataSet(),"$",columna,sep=""),")",sep=""))
    justDoIt(paste("cond4 <- !is.numeric(",paste(ActiveDataSet(),"$",columna,sep=""),")",sep=""))        
    if (cond1 && cond2 || cond3 && cond4){
      errorCondition(recall=bivariante.ordinales, message=gettextRcmdr("Escoja variables ordinales"))
      return()
    }
    remove(list=c("cond1","cond2","cond3","cond4"), envir=.GlobalEnv)
    porcentajes <- as.character(tclvalue(porcentajesVariable))
    gammaval <- tclvalue(gammaVariable)
    tauval <- tclvalue(tauVariable)
    sommersval <- tclvalue(sommersVariable)
    wilsonval <- tclvalue(wilsonVariable)
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("bivariante.ordinales",list(fila.inicial=fila,columna.inicial=columna,
                                          porcentajes.inicial=porcentajes,
                                          gamma.inicial=gammaval,tau.inicial=tauval,sommers.inicial=sommersval,
                                          wilson.inicial=wilsonval,subconjunto.inicial=subconjunto,
                                          echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab)) 
    selec <- as.numeric(gammaval) + as.numeric(tauval)*3 + as.numeric(sommersval)*3 +
      as.numeric(wilsonval)
    if (selec>0)
      stats <- paste("c(",
                     paste(c('"Gamma"', '"Tau"', '"Sommers"', '"Wilson"')
                           [c(gammaval,tauval,sommersval,wilsonval) == 1], 
                           collapse=", "), ")", sep="") 
    else stats <- 'NULL'
    tabs <- if (porcentajes !='ninguno') TRUE else FALSE 
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste(".indices.bo <- bivord(rowvar='",fila, "', colvar='",columna,"', statistics=",stats,
                            ", tables=",tabs, if(porcentajes!="ninguno"){paste(", pcttable='",porcentajes,"'",sep='')},
                            ", data='",ActiveDataSet(),"')", sep="")
    else {
      instruccion1 <- paste(".indices.bo <- bivord(rowvar='",fila, "', colvar='",columna, "', statistics=", stats,
                            ", tables=",tabs, if(porcentajes!="ninguno"){paste(", pcttable='",porcentajes,sep='')},
                            "', subset='",subconjunto, "', data='",ActiveDataSet(),"')", sep="")
    }
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(paste(".indices.bo # Descripcion bivariante de: ",fila, 
                       " y ", columna, sep=""))
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Descripcion bivariante de datos ordinales: ",fila, 
                      " y ", columna, sep="")
      HTML(as.title(titulo),file=.archivo)      
      HTML('Tabla de contigencia', file=.archivo)
      HTML(.indices.bo[[5]], file=.archivo)
      HTML('Descripcion tabla de contigencia', file=.archivo)
      tabladesc <- .indices.bo[[6]]
      rownames(tabladesc) <- c("Pares concordantes","Pares discordantes",
                               paste("Empates",.indices.bo[[1]]),paste("Empates",.indices.bo[[2]]),
                               paste("Empates",.indices.bo[[1]],"y",.indices.bo[[2]]))
      colnames(tabladesc) <- "Valores"     
      HTML(tabladesc, file=.archivo)      
      j <- 7
      if (tabs == TRUE){
        msg <- switch(.indices.bo[[4]],
                      'total'='Tabla porcentajes respecto al total: ',
                      'fila'='Tabla porcentajes respecto a marginales fila: ',
                      'columna'='Tabla porcentajes respecto a marginales columna: ')
        HTML(msg, file=.archivo)        
        HTML(.indices.bo[[j]], file=.archivo)
        j <- j + 1
      }     
      if ( selec>0 ){
        .TablaRes <- as.data.frame(.indices.bo[[j]])
        rownames(.TablaRes) <- c("Gamma", c("Tau a","Tau b","Tau c"),c("d Sommers X/Y","d Sommers Y/X",
                                                                       "d Sommers (simetrica)"),"e Wilson")[
                                                                         c("Gamma", rep("Tau",3), rep("Sommers",3), "Wilson") %in% .indices.bo[[3]]]
        colnames(.TablaRes) <- "Valores"
        HTML("Coeficientes de asociacion: ", file=.archivo)
        HTML(.TablaRes, file=.archivo)
        HTMLhr(file = .archivo)
        j <- j+1
      }                 
    }
    remove(.indices.bo, envir=.GlobalEnv) 
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="xtabs",reset="bivariante.ordinales",apply="bivariante.ordinales")
  tkgrid(getFrame(filaVar), labelRcmdr(variablesFrame, text="    "), getFrame(columnaVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(porcentajesFrame, sticky="w")  
  tkgrid(statsFrame, sticky="w") 
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab',
                                                      'statistics2Tab','optionsTab'),
               tab.names=c("Datos","Tablas","Estadisticos","Opciones")) 
}

dispordfun <- function(rowvar,colvar,subset=NULL,data){
  .Tabla <- eval(parse(text=paste('xtabs(~',rowvar,'+',colvar, if (!is.null(subset)) {paste(',subset=',subset,'')},
                                  ',data=',data,')')))
  def.par <- par(no.readonly = TRUE)
  sup <- max(rowSums(.Tabla),colSums(.Tabla))
  rangox <- eval(parse(text=paste('c(min(as.numeric(',data,'$',rowvar,'),na.rm=TRUE),
                                  max(as.numeric(',data,'$',rowvar,'),na.rm=TRUE))',sep="")))
  rangoy <- eval(parse(text=paste('c(min(as.numeric(',data,'$',colvar,'),na.rm=TRUE),
                                  max(as.numeric(',data,'$',colvar,'),na.rm=TRUE))',sep="")))
  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
  par(mar=c(5,4,3,3))
  eval(parse(text=paste('plot(as.numeric(',data,'$',rowvar,'),
                        as.numeric(',data,'$',colvar,'),xlim=rangox,ylim=rangoy,
                        ,xlab="',rowvar,'",ylab="',colvar,'")',sep="")))
  z <- eval(parse(text=paste('with(',data,',merge(data.frame(',rowvar,',',
                             colvar,'),melt(table(',rowvar,',',colvar,')),sort =F)$value)',sep="")))
  z <- eval(parse(text=paste('1.5*z/length(',data,'$',rowvar,')',sep="")))
  eval(parse(text=paste("symbols(na.omit(data.frame(",data,"$",rowvar,",",data,"$",colvar,")),
                        circles=z,inches=F, bg='grey',fg=NA,add=T)",sep="")))
  ok <- eval(parse(text=paste('is.finite(',data,'$',rowvar,') & is.finite(',data,'$',colvar,')',sep="")))
  if (any(ok)) 
    eval(parse(text=paste('lines(stats::lowess(',data,'$',rowvar,'[ok],',data,'$',colvar,'[ok],
                          f = 2/3, iter = 3), col = "red")',sep="")))
  box()
  par(mar=c(0,3,1,1))
  eval(parse(text=paste('barplot(table(',data,'$',rowvar,'),axes=FALSE,ylim=c(0,', sup,'),space=0.1,col="grey")',sep=""))) 
  eval(parse(text=paste('title("Diagrama de puntos para ',rowvar,' y ',colvar,'")',sep="")))
  par(mar=c(3,0,1,1)) 
  eval(parse(text=paste('barplot(table(',data,'$',colvar,'),axes=FALSE,xlim=c(0,', sup,'),space=0.1,horiz=TRUE,col="grey")',sep="")))
  par(def.par)
}

dispersion.ordinales <- function(){
  defecto <- list(fila.inicial=NULL,columna.inicial=NULL,
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("dispersion.ordinales",defecto) 
  initializeDialog(title=gettextRcmdr("Diagrama de puntos para datos ordinales"),use.tabs=TRUE,
                   tabs=c('dataTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  filaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                             initialSelection=varPosn(dialog.valores$fila.inicial,"factor")) 
  columnaVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Variables (escoja una)"),
                                initialSelection=varPosn(dialog.valores$columna.inicial,"factor"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))       
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    fila <- getSelection(filaVar)
    columna <- getSelection(columnaVar)
    if (length(fila) == 0 || length(columna) == 0){
      errorCondition(recall=barras.agrupadas, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (fila == columna) {
      errorCondition(recall=barras.agrupadas, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    justDoIt(paste("cond1 <- !is.ordered(",paste(ActiveDataSet(),"$",fila,sep=""),")",sep=""))
    justDoIt(paste("cond2 <- !is.ordered(",paste(ActiveDataSet(),"$",columna,sep=""),")",sep=""))
    if (cond1 || cond2){
      errorCondition(recall=dispersion.ordinales, message=gettextRcmdr("Escoja variables ordinales"))
      return()
    }
    remove(list=c("cond1","cond2"), envir=.GlobalEnv)    
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("dispersion.ordinales",list(fila.inicial=fila,columna.inicial=columna,
                                          subconjunto.inicial=subconjunto,
                                          echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))         
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste("dispordfun(rowvar='",fila, "', colvar='",columna,"', data='",ActiveDataSet(),"')", 
                            sep="")
    else {
      instruccion1 <- paste("dispordfun(rowvar='",fila, "', colvar='",columna,"', subset='",subconjunto,
                            "', data='",ActiveDataSet(),"')", sep="")
    }
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Diagrama de puntos para datos ordinales: ",fila, 
                      " y ", columna, sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("DiagramaPuntosOrdinalesR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="plot",reset="dispersion.ordinales",apply="dispersion.ordinales")
  tkgrid(getFrame(filaVar), labelRcmdr(variablesFrame, text="    "), getFrame(columnaVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Datos","Opciones"))        
}

# Resumen descriptivo bivariante para variables cuantitativas I/R: indicadores y graficos #

bivnum <- function(datos,variables,stats,subconjunto=NULL){
  if (is.null(subconjunto)) datos <- eval(parse(text='datos[variables]'))
  else datos <- eval(parse(text=paste('subset(datos,subset=',subconjunto,')[variables]')))
  res <- list()
  res[[1]] <- variables
  res[[2]] <- stats
  res[[3]] <- subconjunto
  numfilas <- sum(c("Variancia","Pearson","","Spearman","Kendall") %in% stats) + ("CoefDeterm" %in% stats)*
    sum(c("Pearson","","Spearman","Kendall") %in% stats)
  .tablaAsoc <- as.data.frame(matrix(nrow=numfilas,ncol=1))
  j <- 0
  if ("Covariancia" %in% stats)
  {
    covariancia <- cov(datos,use='na.or.complete')[1,2]
    j <- j + 1
    .tablaAsoc[j,]<-covariancia
  }  
  if ("Pearson" %in% stats)
  {
    correlacion <- cor(datos,method='pearson',use='na.or.complete')[1,2]
    j <- j + 1
    .tablaAsoc[j,]<-correlacion
    if ("CoefDeterm" %in% stats)
    {
      R2pearson <- correlacion^2
      j <- j + 1
      .tablaAsoc[j,]<-R2pearson
    }
  }
  if ("Spearman" %in% stats)
  {
    correlacion <- cor(datos,method='spearman',use='na.or.complete')[1,2]
    j <- j + 1
    .tablaAsoc[j,]<-correlacion
    if ("CoefDeterm" %in% stats)
    {
      R2spearman <- correlacion^2
      j <- j + 1
      .tablaAsoc[j,]<-R2spearman
    }
  }
  if ("Kendall" %in% stats)
  {
    correlacion <- cor(datos,method='kendall',use='na.or.complete')[1,2]
    j <- j + 1
    .tablaAsoc[j,]<-correlacion
    if ("CoefDeterm" %in% stats)
    {
      R2kendall <- correlacion^2
      j <- j + 1
      .tablaAsoc[j,]<-R2kendall     
    }
  }
  res[[4]] <- .tablaAsoc
  class(res) <- "bivnum"
  res
}

print.bivnum <- function(x,...){
  cat("Coeficientes de asociacion para ",x[[1]][1], " y ",x[[1]][2],": \n\n",sep='')
  .tabla <- as.data.frame(round(x[[4]],3))
  rownames(.tabla) <- c(if ("Covariancia"%in% x[[2]]) "Covariancia",if ("Pearson"%in% x[[2]]) "r de Pearson",
                        if ("Pearson"%in% x[[2]] & "CoefDeterm" %in% x[[2]])"r Cuadrado (Pearson)",
                        if ("Spearman"%in% x[[2]]) "r de Spearman", 
                        if ("Spearman"%in% x[[2]] & "CoefDeterm" %in% x[[2]])"r Cuadrado (Spearman)",
                        if ("Kendall"%in% x[[2]]) "tau de Kendall", 
                        if ("Kendall"%in% x[[2]] & "CoefDeterm" %in% x[[2]]) "r Cuadrado (Kendall)")
  colnames(.tabla) <- "Valores"
  print(.tabla)  
  cat("\n\n")
  invisible(x)  
}

bivariante.numericas <- function(){
  defecto <- list(var1.inicial=NULL,var2.inicial=NULL,Cov.inicial="0",
                  pearson.inicial="0",spearman.inicial="0",kendall.inicial="0",determ.inicial="0",
                  subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("bivariante.numericas",defecto) 
  initializeDialog(title=gettextRcmdr("Descripcion bivariante datos cuantitativos"),use.tabs=TRUE,
                   tabs=c('dataTab','statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variable1Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var1.inicial,"numeric")) 
  variable2Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var2.inicial,"numeric"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  
  checkBoxes(statisticsTab,frame="statsFrame",boxes=c("Cov","pearson","spearman","kendall","determ"),
             initialValues=c(dialog.valores$Cov.inicial,dialog.valores$pearson.inicial,
                             dialog.valores$spearman.inicial,dialog.valores$kendall.inicial,
                             dialog.valores$determ.inicial),
             labels=gettextRcmdr(c("Coeficiente de covariancia ","Coeficiente de correlacion de Pearson ",
                                   "Coeficiente de correlacion de Spearman ","Coeficiente de correlacion de Kendall ",
                                   "Coeficiente de determinacion ")), 
             title = gettextRcmdr("Coeficientes de Asociacion"))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones")) 
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    var1 <- getSelection(variable1Var)
    var2 <- getSelection(variable2Var)
    if (length(var1) == 0 || length(var2) == 0){
      errorCondition(recall=bivariante.numericas, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (var1 == var2) {
      errorCondition(recall=bivariante.numericas, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    Covval <- tclvalue(CovVariable)
    pearsonval <- tclvalue(pearsonVariable)
    spearmanval <- tclvalue(spearmanVariable)
    kendallval <- tclvalue(kendallVariable)   
    determval <- tclvalue(determVariable)
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    selec <- as.numeric(Covval) + as.numeric(pearsonval) + as.numeric(spearmanval) +
      as.numeric(kendallval) + as.numeric(determval)*(as.numeric(pearsonval) + 
                                                        as.numeric(spearmanval) +as.numeric(kendallval))
    if (selec == 0){
      errorCondition(recall=bivariante.numericas, 
                     message=gettextRcmdr("Debe escoger algun indicador."))
      return()
    }
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("bivariante.numericas",list(var1.inicial=var1,var2.inicial=var2,
                                          Cov.inicial=Covval,
                                          pearson.inicial=pearsonval,spearman.inicial=spearmanval,
                                          kendall.inicial=kendallval,determ.inicial=determval,
                                          subconjunto.inicial=subconjunto,
                                          echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab)) 
    vars <- paste("c(", paste('"', c(var1,var2), '"', collapse=", ", sep=""), ")", sep="")
    stats <- paste("c(",
                   paste(c('"Covariancia"','"Pearson"','"Spearman"','"Kendall"','"CoefDeterm"')
                         [c(Covval,pearsonval,spearmanval,kendallval,determval) == 1], 
                         collapse=", "), ")", sep="")    
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste(".indices.bn <- bivnum(datos=",ActiveDataSet(),",variables=",vars,",stats=", stats,
                            ")", sep="")
    else {
      instruccion1 <- paste(".indices.bn <- bivnum(datos=",ActiveDataSet(),",variables=c(",vars,"),stats=", stats,
                            ", subconjunto='",subconjunto,"')", sep="")
    }
    logger(instruccion1)
    justDoIt(instruccion1)
    
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(paste(".indices.bn # Descripcion bivariante de: ",var1, 
                       " y ", var2, sep=""))        
    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Descripcion bivariante de datos cuantitativos: ",var1, 
                      " y ", var2, sep="")
      HTML(as.title(titulo),file=.archivo)
      .TablaRes <- as.data.frame(round(.indices.bn[[4]],3))
      rownames(.TablaRes) <- c(if ("Covariancia"%in% .indices.bn[[2]]) "Covariancia",if ("Pearson"%in% .indices.bn[[2]]) "r de Pearson",
                               if ("Pearson"%in% .indices.bn[[2]] & "CoefDeterm" %in% .indices.bn[[2]])"r Cuadrado (Pearson)",
                               if ("Spearman"%in% .indices.bn[[2]]) "r de Spearman", 
                               if ("Spearman"%in% .indices.bn[[2]] & "CoefDeterm" %in% .indices.bn[[2]])"r Cuadrado (Spearman)",
                               if ("Kendall"%in%.indices.bn[[2]]) "tau de Kendall", 
                               if ("Kendall"%in% .indices.bn[[2]] & "CoefDeterm" %in% .indices.bn[[2]]) "r Cuadrado (Kendall)")
      colnames(.TablaRes) <- "Valores"
      HTML("Coeficientes de asociacion: ", file=.archivo)
      HTML(.TablaRes, file=.archivo)
      HTMLhr(file = .archivo)      
      HTMLhr(file = .archivo)      
    }
    remove(.indices.bn, envir=.GlobalEnv) 
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="cov",reset="bivariante.numericas",apply="bivariante.numericas")
  tkgrid(getFrame(variable1Var), labelRcmdr(variablesFrame, text="    "), getFrame(variable2Var), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w") 
  tkgrid(statsFrame, sticky="w") 
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones")) 
}

dispnumfun <- function(datos,variables,subconjunto=NULL){
  var1 <- variables[1]
  var2 <- variables[2]
  def.par <- par(no.readonly = TRUE)
  if (is.null(subconjunto)) datos <- eval(parse(text='datos[variables]'))
  else datos <- eval(parse(text=paste('subset(datos,subset=',subconjunto,')[variables]')))
  xhist <- hist(datos[,1],plot=FALSE)
  yhist <- hist(datos[,2],plot=FALSE)
  sup <- max(c(xhist$counts,yhist$counts))
  rangox <- c(min(datos[,1],na.rm=TRUE),max(datos[,1],na.rm=TRUE))
  rangoy <- c(min(datos[,2],na.rm=TRUE),max(datos[,2],na.rm=TRUE))
  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
  par(mar=c(5,4,3,3))
  plot(datos[,1],datos[,2],xlim=rangox,ylim=rangoy,xlab=var1,ylab=var2)
  box()
  par(mar=c(0,3,1,1))
  barplot(xhist$counts,axes=FALSE,ylim=c(0, sup),space=0,col='grey')
  title(paste("Diagrama de dispersion para ",var1," y ",var2,sep=""))
  par(mar=c(3,0,1,1))
  barplot(yhist$counts,axes=FALSE,xlim=c(0, sup),space=0,col='grey',horiz=TRUE)
  par(def.par)
}

dispersion.numericas <- function(){
  defecto <- list(var1.inicial=NULL,var2.inicial=NULL,subconjunto.inicial=gettextRcmdr("<all valid cases>"),
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("dispersion.numericas",defecto) 
  initializeDialog(title=gettextRcmdr("Diagrama de dispersion para variables cuantitativas"),use.tabs=TRUE,
                   tabs=c('dataTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variable1Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var1.inicial,"numeric")) 
  variable2Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var2.inicial,"numeric"))  
  subsetBox(dataTab, subset.expression=dialog.valores$subconjunto.inicial)
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))  
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1    
    var1 <- getSelection(variable1Var)
    var2 <- getSelection(variable2Var)
    if (length(var1) == 0 || length(var2) == 0){
      errorCondition(recall=dispersion.numericas, message=gettextRcmdr("Debe seleccionar dos variables."))
      return()
    }
    if (var1 == var2) {
      errorCondition(recall=dispersion.numericas, message=gettextRcmdr("Debe seleccionar dos variables distintas."))
      return()
    }
    subconjunto <- tclvalue(subsetVariable)
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    putDialog("dispersion.numericas",list(var1.inicial=var1,var2.inicial=var2,
                                          subconjunto.inicial=subconjunto,
                                          echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab)) 
    vars <- paste("c(", paste('"', c(var1,var2), '"', collapse=", ", sep=""), ")", sep="")   
    if (trim.blanks(subconjunto) == gettextRcmdr("<all valid cases>")) 
      instruccion1 <- paste("dispnumfun(datos=",ActiveDataSet(),",variables=",vars,")", sep="")
    else {
      instruccion1 <- paste("dispnumfun(datos=",ActiveDataSet(),",variables=c(",vars,"), subconjunto='",
                            subconjunto,"')", sep="")
    }
    justDoIt(instruccion1)  
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }       
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Diagrama de puntos para datos cuantitativos: ",var1, 
                      " y ", var2, sep="")
      HTML(as.title(titulo),file=.archivo)
      nombre.archivo <- paste("DiagramaDispersionR",gsub(":","",substr(Sys.time(),12,19)),
                              ".jpg",sep="")
      dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                width=500, height=500)
      HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
      HTMLhr(file = .archivo)
    }
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="plot",reset="dispersion.numericas",apply="dispersion.numericas")
  tkgrid(getFrame(variable1Var), labelRcmdr(variablesFrame, text="    "), getFrame(variable2Var), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Datos","Opciones"))        
}

# Alguna funciones necesarias para llevar a cabo el test de Yuen-Welch, extraidas del paquete WRS #

winvar <- function(x,tr=.2,na.rm=FALSE){
  if(na.rm)x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  winvar<-var(y)
  winvar
}

yuen.test <- function(x,y,tr=.2,alpha=.05){
  if(tr==.5)stop("No se puede utilizar tr=0.5")
  if(tr>.25)print("Aviso: con tr>.25 el control sobre el error tipo I puede ser insuficiente")
  x<-x[!is.na(x)]
  y<-y[!is.na(y)]
  h1<-length(x)-2*floor(tr*length(x))
  h2<-length(y)-2*floor(tr*length(y))
  q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
  q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
  df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
  crit<-qt(1-alpha/2,df)
  dif<-mean(x,tr)-mean(y,tr)
  low<-dif-crit*sqrt(q1+q2)
  up<-dif+crit*sqrt(q1+q2)
  test<-abs(dif/sqrt(q1+q2))
  pval<-2*(1-pt(test,df))
  res <- list(ci=c(low,up),p.value=pval,dif=dif,se=sqrt(q1+q2),
              teststat=test,crit=crit,gl=df,ic=1-alpha,rec=tr)
  class(res)<-"yuen"
  res
}

print.yuen <- function (x,digits = max(4,getOption("digits") - 4),...)
{
  cat("\n")
  cat("\t",'Prueba t de Yuen-Welch',"\n")
  cat("\n")
  cat(paste("T de Yuen = ",round(x$teststat,4),", gl = ",round(x$gl,2),", valor p = ",
            round(x$p.value,4),sep=""),"\n")
  cat("Hipotesis Alternativa: Diferencia entre medias recortadas distinta de 0","\n")
  cat(paste("Intervalo de confianza ",x$ic*100,"%: ",sep=""),"\n")
  cat("\t",round(x$ci,4),"\n")
  cat(paste("Estimacion puntual de la diferencia de medias recortadas al ",x$rec*100,"%: ",sep=""),"\n")
  cat("\t",round(x$dif,4),"\n")  
  cat("\n\n")
  invisible(x)  
}

yuenWelch.tTest <- function(){
  defecto <- list(grupo.inicial=NULL,var.inicial=NULL,trim.inicial="0.2",confint.inicial=".95",
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("yuenWelch.tTest",defecto) 
  initializeDialog(title=gettextRcmdr("Prueba t de Yuen-Welch"),use.tabs=TRUE,tabs=c('dataTab',
                                                                                     'statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  grupoVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Grupos (escoja una)"),
                              initialSelection=varPosn(dialog.valores$grupo.inicial,"factor")) 
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                 initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))    
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1     
    grupo <- getSelection(grupoVar)
    if (length(grupo) == 0) {
      errorCondition(recall=yuenWelch.tTest, message=gettextRcmdr("Debe seleccionar una variable de agrupacion."))
      return()
    }
    varresp <- getSelection(variableVar)
    if (length(varresp) == 0) {
      errorCondition(recall=yuenWelch.tTest, message=gettextRcmdr("Debe seleccionar una variable de respuesta."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    rec <- as.numeric(tclvalue(trimVariable))
    if (rec > 1 | rec < 0)
    {
      rec <- 0.2
      Message(message=gettextRcmdr("Proporcion de recorte invalida, se utilizara valor por defecto."),
              type="warning")              
    }
    ic <- as.numeric(tclvalue(ConfIntVariable))
    if ( ic < .0 || ic > 1. || !is.numeric(ic) )
    {
      ic <- 0.95
      Message(message=gettextRcmdr("Nivel de confianza invalido, se utilizara valor por defecto."),
              type="warning")              
    }
    putDialog("yuenWelch.tTest",list(grupo.inicial=grupo,var.inicial=varresp,
                                     trim.inicial=rec,confint.inicial=ic,
                                     echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))     
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Prueba t de Yuen-Welch: ",varresp, 
                      " por ", grupo, sep="")
      HTML(as.title(titulo),file=.archivo)
    }
    closeDialog()
    instruccion1 <- paste("x <-", ActiveDataSet(),"$",varresp,"[",ActiveDataSet(),"$",grupo,
                          "==levels(",ActiveDataSet(),"$",grupo,")[1]]",sep="")
    justDoIt(instruccion1)
    instruccion2 <- paste("y <- ", ActiveDataSet(),"$",varresp,"[",ActiveDataSet(),"$",grupo,
                          "==levels(",ActiveDataSet(),"$",grupo,")[2]]",sep="")
    justDoIt(instruccion2)
    instruccion3 <- paste("ywttest <- yuen.test(x,y,tr=",rec,",alpha=1-",ic,")",sep="")
    justDoIt(instruccion3)   
    if (echocodigo == 1)
    {
      logger(instruccion1)
      logger(instruccion2)
      logger(instruccion3)
    }
    doItAndPrint(paste("ywttest # Prueba t de Yuen-Welch para ", varresp," segun ",grupo,sep=""))
    if (creahtml == 1)
    {
      
      HTML('Prueba t de Yuen-Welch', file=.archivo)
      HTML(paste("datos: ",ActiveDataSet(),"$",varresp," por ",ActiveDataSet(),"$",grupo,sep=""), file=.archivo)
      HTML(paste("T de Yuen = ",round(ywttest[[5]],4),", gl = ",round(ywttest[[7]],2),", valor p = ",
                 round(ywttest[[2]],4),sep=""), file=.archivo)
      HTML("Hipotesis Alternativa: Diferencia entre medias recortadas distinta de 0", file=.archivo)
      HTML(paste("Intervalo de confianza ",ic*100,"%: ",sep=""), file=.archivo)
      HTML(round(ywttest[[1]],4), file=.archivo)
      HTML(paste("Estimacion puntual de la diferencia de medias recortadas al ",rec*100,"%: ",sep=""), file=.archivo)
      HTML(round(ywttest[[3]],4), file=.archivo)
      HTMLhr(file = .archivo)
    }
    remove(list=c('x','y','ywttest'), envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t.test",reset="yuenWelch.tTest",apply="yuenWelch.tTest")
  tkgrid(getFrame(grupoVar), labelRcmdr(variablesFrame, text="    "), getFrame(variableVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  groupsLabel(groupsBox=grupoVar,dataTab)
  trimFrame <- tkframe(statisticsTab)
  trimVariable <- tclVar(dialog.valores$trim.inicial)
  trimField <- ttkentry(trimFrame, width="8", textvariable=trimVariable)
  tkgrid(labelRcmdr(trimFrame,text=gettextRcmdr("Proporcion datos recortados = "),font="RcmdrTitleFont"),
         trimField,sticky="w")
  ConfIntFrame <- tkframe(statisticsTab)
  ConfIntVariable <- tclVar(dialog.valores$confint.inicial)
  ConfIntField <- ttkentry(ConfIntFrame, width="8", textvariable=ConfIntVariable)
  tkgrid(labelRcmdr(ConfIntFrame,text=gettextRcmdr("Nivel de Confianza = "),font="RcmdrTitleFont"),
         ConfIntField,sticky="w")  
  tkgrid(trimFrame, sticky="nw")
  tkgrid(ConfIntFrame, sticky="nw")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

# Funcion adaptada de twoSampleWilcoxonTest de J. Fox #

UMWfun <- function(y,x,alternativa,conf.int,prueba,data){
  res <- list()
  res[[1]] <- x
  res[[2]] <- y
  res[[3]] <- data
  i <- 4
  meds <- eval(parse(text=paste("tapply(",data,"$", y,",",data,"$",x,
                                ", median, na.rm=TRUE)", sep="")))
  res[[4]] <- meds
  g1 <- eval(parse(text=paste("with(",data,",",y,"[as.numeric(factor(",x,"))==1])",sep="")))
  g2 <- eval(parse(text=paste("with(",data,",",y,"[as.numeric(factor(",x,"))==2])",sep="")))
  rangos <- rank(c(na.omit(g1),na.omit(g2)))
  R1 <- sum(rangos[1:length(g1)])
  R2 <- sum(rangos[(length(na.omit(g1))+1):length(rangos)])
  Rangos <- eval(parse(text=paste("matrix(round(c(R1/length(na.omit(g1)),R1,R2/length(na.omit(g2)),R2),2),
                                  nrow=2,byrow=T,dimnames=list(levels(",data,"$",x,"),c('Rango promedio','Suma Rangos')))",sep="")))
  res[[5]] <- Rangos
  if (prueba == "default"){
    UMW <- eval(parse(text=paste("wilcox.test(",y, " ~ ",x, ', alternative="', 
                                 alternativa, '",conf.int=TRUE,conf.level=',conf.int,",data=",data, ")", sep="")))
  }
  if (prueba == "normal"){
    UMW <- eval(parse(text=paste("wilcox.test(",y, " ~ ",x, ', alternative="', 
                                 alternativa, '",conf.int=TRUE,conf.level=',conf.int,",data=",data, "exact=FALSE,
                                 correct=FALSE)",sep="")))
  }
  if (prueba == "correct"){
    UMW <- eval(parse(text=paste("wilcox.test(",y, " ~ ",x, ', alternative="', 
                                 alternativa, '",conf.int=TRUE,conf.level=',conf.int,",data=",data, "exact=FALSE,
                                 correct=TRUE)",sep="")))
  }  
  else {
    UMW <- eval(parse(text=paste("wilcox.test(",y, " ~ ",x, ", alternative='", 
                                 alternativa, "', exact=", prueba=="exact", 
                                 ", correct=", prueba=="correct",", conf.int=TRUE, conf.level=",conf.int,", data=",
                                 data, ")", sep="")))
  }
  res[[6]] <- UMW
  class(res) <- "UMW"
  res
  }

print.UMW <- function (x,digits = max(4,getOption("digits") - 4),...)
{
  cat('Medianas para',x[[2]],"segun",x[[1]],"\n")
  cat("\n")
  print(x[[4]])  
  cat("\n\n")  
  cat("Resumen descriptivo por rangos \n")
  cat("\n")
  print(x[[5]])  
  cat("\n\n")  
  cat('Prueba U de Mann-Whitney para',x[[2]],"por",x[[1]],"\n")
  cat("\n")
  print(x[[6]])  
  cat("\n\n")
  invisible(x)  
}

U.Mann.Whitney <- function(){
  defecto <- list(grupo.inicial=NULL,var.inicial=NULL,alt.inicial="two.sided",confint.inicial=".95",
                  prueba.inicial="default",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("U.Mann.Whitney",defecto) 
  initializeDialog(title=gettextRcmdr("Prueba U de Mann-Whitney"),use.tabs=TRUE,tabs=c('dataTab',
                                                                                       'statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  grupoVar <- variableListBox(variablesFrame, Factors(), title=gettextRcmdr("Grupos (escoja una)"),
                              initialSelection=varPosn(dialog.valores$grupo.inicial,"factor")) 
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                 initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1     
    grupo <- getSelection(grupoVar)
    if (length(grupo) == 0) {
      errorCondition(recall=U.Mann.Whitney, message=gettextRcmdr("Debe seleccionar una variable de agrupacion."))
      return()
    }
    varresp <- getSelection(variableVar)
    if (length(varresp) == 0) {
      errorCondition(recall=U.Mann.Whitney, message=gettextRcmdr("Debe seleccionar una variable de respuesta."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    ic <- as.numeric(tclvalue(ConfIntVariable))
    if ( ic < .0 || ic > 1. || !is.numeric(ic) )
    {
      ic <- 0.95
      Message(message=gettextRcmdr("Nivel de confianza invalido, se utilizara valor por defecto."),
              type="warning")              
    }
    alternativa <- as.character(tclvalue(alternativaVariable))
    prueba <- as.character(tclvalue(pruebaVariable))
    putDialog("U.Mann.Whitney",list(grupo.inicial=grupo,var.inicial=varresp,
                                    alt.inicial=alternativa,confint.inicial=ic,prueba.inicial=prueba,
                                    echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))
    .baseDatosActiva <- ActiveDataSet()
    opts <- options(warn=-1)
    instruccion1 <- paste("UMW <- UMWfun(y='", varresp, "',x='", grupo, "', alternativa='", 
                          alternativa, "',conf.int=",ic,",prueba='",prueba,"',data='", .baseDatosActiva, "')",
                          sep="")
    justDoIt(instruccion1)        
    if (echocodigo==1) logger(instruccion1)
    doItAndPrint(paste("UMW # Prueba U de Mann-Whitney para ", varresp," segun ",grupo,sep="")) 
    closeDialog()
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Prueba U de Mann-Whitney: ",varresp, 
                      " por ", grupo, sep="")
      HTML(as.title(titulo),file=.archivo)
      HTML(paste('Medianas para ',varresp," segun ",grupo,sep=""), file=.archivo)
      .medstabla <- as.data.frame(matrix(UMW[[4]],2,1,dimnames=list(names(UMW[[4]]),"Mediana")))
      HTML(.medstabla,file=.archivo)
      HTML('Resumen Descriptivo',file=.archivo)
      HTML(as.data.frame(UMW[[5]]),file=.archivo)
      HTML(paste('Prueba U de Mann-Whitney',if(prueba=='exact')" (exacta)",
                 if(prueba=='correct')" (Normal con correccion por continuidad)",
                 if(prueba=='normal')" (Normal)",sep=""),file=.archivo)
      HTML(paste("datos: ",ActiveDataSet(),"$",varresp," segun ",ActiveDataSet(),"$",grupo,sep=""), file=.archivo)
      if (UMW[[6]][3] > 0.0001) pval<-paste(" = ",round(UMW[[6]][[3]],4),sep="") else pval<- " < 0.0001"
      HTML(paste("U = ",round(UMW[[6]][[1]],4),", valor p",pval,sep=""), file=.archivo)
      HTML("Hipotesis Alternativa: Diferencia entre distribuciones distinta de 0", file=.archivo)
      HTML(paste("Intervalo de confianza ",ic*100,"%: ",sep=""), file=.archivo)
      HTML(round(UMW[[6]][[8]],4), file=.archivo)
      HTML(paste("Estimacion puntual de la diferencia de medianas: ",sep=""), file=.archivo)
      meddif <- UMW[[6]][[9]]
      names(meddif) <- NULL
      HTML(round(meddif,4), file=.archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove('UMW')")
    remove('UMW', envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="wilcox.test",reset="U.Mann.Whitney",apply="U.Mann.Whitney")
  tkgrid(getFrame(grupoVar), labelRcmdr(variablesFrame, text="    "), getFrame(variableVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  groupsLabel(groupsBox=grupoVar,dataTab)
  statsFrame <- tkframe(statisticsTab)
  radioButtons(statsFrame,name="alternativa", buttons=c("twosided", "less", "greater"),
               values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Bilateral", "Diferencia < 0", "Diferencia > 0")),
               title=gettextRcmdr("Hipotesis alternativa"), initialValue= dialog.valores$alt.inicial)
  radioButtons(statsFrame,name="prueba", buttons=c("default", "exact", "normal", "correct"), 
               labels=gettextRcmdr(c("Por defecto", "Exacta", "Aproximacion normal",
                                     "Aprox. normal con\ncorrecc. continuidad")), 
               title=gettextRcmdr("Tipo de Prueba"), initialValue= dialog.valores$prueba.inicial)
  ConfIntFrame <- tkframe(statsFrame)
  ConfIntVariable <- tclVar(dialog.valores$confint.inicial)
  ConfIntField <- ttkentry(ConfIntFrame, width = "6", 
                           textvariable = ConfIntVariable)
  tkgrid(labelRcmdr(ConfIntFrame, text = gettextRcmdr("Nivel de confianza"), 
                    fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
  tkgrid(ConfIntField, sticky = "w")
  tkgrid(alternativaFrame, labelRcmdr(statsFrame, text = "    "), 
         ConfIntFrame, labelRcmdr(statsFrame, text = "    "), 
         pruebaFrame, sticky = "nw")
  tkgrid(statsFrame, sticky = "nw")  
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

TWfun <- function(var1,var2,alternativa,conf.int,prueba,data){
  res <- list()
  res[[1]] <- var1
  res[[2]] <- var2
  res[[3]] <- data
  meddif <- eval(parse(text=paste("median(", data, "$", var1, " - ", data, "$", var2, 
                                  ", na.rm=TRUE)", sep="")))
  res[[4]] <- meddif
  m1 <- eval(parse(text=paste(data,"$",var1,sep="")))
  m2 <- eval(parse(text=paste(data,"$",var2,sep="")))
  difs <- na.omit((m1-m2)[(m1-m2)!=0])
  Rpos <- sum(rank(abs(difs))[difs>0])
  difpos <- sum(difs>0)
  Rneg <- sum(rank(abs(difs))[difs<0])
  difneg <- sum(difs<0)
  empates <- sum((m1-m2)==0)
  Rangos <- matrix(round(c(difpos,if (difpos !=0) Rpos/difpos else 0,
                           Rpos,difneg,if (difneg !=0) Rneg/difneg else 0,Rneg,empates,0,0),2),nrow=3,byrow=T,
                   dimnames=list(c('Rangos positivos','Rangos negativos','Empates'),c('Frecuencia','Rango promedio','Suma Rangos')))
  res[[5]] <- Rangos
  if (prueba == "default"){
    TW <- eval(parse(text=paste("wilcox.test(",data,"$",var1, ",",data,"$",var2, ", alternative='", 
                                alternativa, "',conf.int=TRUE,conf.level=",conf.int,",paired=TRUE)", sep="")))
  }
  if (prueba == "normal"){
    TW <- eval(parse(text=paste("wilcox.test(",data,"$",var1, ",",data,"$",var2, ", alternative='", 
                                alternativa, "',conf.int=TRUE,conf.level=",conf.int,",exact=FALSE,
                                correct=FALSE,paired=TRUE)", sep="")))
  }
  if (prueba == "correct"){
    TW <- eval(parse(text=paste("wilcox.test(",data,"$",var1, ",",data,"$",var2, ", alternative='", 
                                alternativa, "',conf.int=TRUE,conf.level=",conf.int,",exact=FALSE,
                                correct=TRUE,paired=TRUE)", sep="")))
  }  
  else {
    TW <- eval(parse(text=paste("wilcox.test(",data,"$",var1, ",",data,"$",var2, ", alternative='", 
                                alternativa, "',conf.int=TRUE,conf.level=",conf.int,",exact=TRUE,
                                correct=FALSE,paired=TRUE)", sep="")))
  }
  res[[6]] <- TW
  class(res) <- "TW"
  res  
  }

print.TW <- function (x,digits = max(4,getOption("digits") - 4),...)
{
  cat('Diferencia de medianas para',x[[1]],"y",x[[2]],":",x[[4]],"\n")
  cat("\n\n")  
  cat("Resumen descriptivo por rangos \n")
  cat("\n")
  print(x[[5]])  
  cat("\n\n")  
  cat('Prueba T de Wilcoxon para',x[[1]],"y",x[[2]],"\n")
  cat("\n")
  print(x[[6]])  
  cat("\n\n")
  invisible(x)  
}

pruebaT.Wilcoxon <- function(){
  defecto <- list(var1.inicial=NULL,var2.inicial=NULL,alt.inicial="two.sided",confint.inicial=".95",
                  prueba.inicial="default",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("pruebaT.Wilcoxon",defecto) 
  initializeDialog(title=gettextRcmdr("Prueba T de Wilcoxon"),use.tabs=TRUE,tabs=c('dataTab',
                                                                                   'statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variable1Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variable 1 (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var1.inicial,"numeric")) 
  variable2Var <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variable 2 (escoja una)"),
                                  initialSelection=varPosn(dialog.valores$var2.inicial,"numeric"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1      
    var1 <- getSelection(variable1Var)
    if (length(var1) == 0) {
      errorCondition(recall=pruebaT.Wilcoxon, message=gettextRcmdr("Debe seleccionar una variable de agrupacion."))
      return()
    }
    var2 <- getSelection(variable2Var)
    if (length(var2) == 0) {
      errorCondition(recall=pruebaT.Wilcoxon, message=gettextRcmdr("Debe seleccionar una variable de respuesta."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    ic <- as.numeric(tclvalue(ConfIntVariable))
    if ( ic < .0 || ic > 1. || !is.numeric(ic) )
    {
      ic <- 0.95
      Message(message=gettextRcmdr("Nivel de confianza invalido, se utilizara valor por defecto."),
              type="warning")              
    }
    alternativa <- as.character(tclvalue(alternativaVariable))
    prueba <- as.character(tclvalue(pruebaVariable))
    putDialog("pruebaT.Wilcoxon",list(var1.inicial=var1,var2.inicial=var2,
                                      alt.inicial=alternativa,confint.inicial=ic,prueba.inicial=prueba,
                                      echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))
    opts <- options(warn=-1)
    .baseDatosActiva <- ActiveDataSet()
    instruccion1 <- paste("TW <- TWfun(var1='", var1, "',var2='", var2, "', alternativa='", 
                          alternativa, "',conf.int=",ic,",prueba='",prueba,"',data='", .baseDatosActiva, "')",
                          sep="")
    justDoIt(instruccion1)        
    if (echocodigo==1) logger(instruccion1)
    doItAndPrint(paste("TW # Prueba T de Wilcoxon para ", var1," y ",var2,sep=""))  
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Prueba de Wilcoxon (muestras relacionadas) para: ",var1, 
                      " y ", var2, sep="")
      HTML(as.title(titulo),file=.archivo)
      HTML(paste('Diferencia de medianas entre ',var1," y ",var2," : ",TW[[4]],sep=""), file=.archivo)
      HTML(paste('Resumen Descriptivo para ', var1, ' - ', var2,sep=""),file=.archivo)
      HTML(as.data.frame(TW[[5]]),file=.archivo)
      HTML(paste('Prueba T de Wilcoxon',if(prueba=='exact')" (exacta)",
                 if(prueba=='correct')" (Normal con correccion por continuidad)",
                 if(prueba=='normal')" (Normal)",sep=""),file=.archivo)
      HTML(paste("datos: ",ActiveDataSet(),"$",var1," y ",ActiveDataSet(),"$",var2,sep=""), file=.archivo)
      if (TW[[6]][[3]] > 0.0001) pval<-paste(" = ",round(TW[[6]][[3]],4),sep="") else pval<- " < 0.0001"
      HTML(paste("V = ",round(TW[[6]][[1]],4),", valor p",pval,sep=""), file=.archivo)
      HTML("Hipotesis Alternativa: Diferencia entre distribuciones distinta de 0", file=.archivo)
      HTML(paste("Intervalo de confianza ",ic*100,"%: ",sep=""), file=.archivo)
      HTML(round(TW[[6]][[8]],4), file=.archivo)
      HTML(paste("Estimacion puntual de la (pseudo)mediana de las diferencias: ",sep=""), file=.archivo)
      meddifest <- TW[[6]][[9]]
      names(meddifest) <- NULL
      HTML(round(meddifest,4), file=.archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove('TW')")
    remove('TW', envir=.GlobalEnv)
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="wilcox.test",reset="pruebaT.Wilcoxon",apply="pruebaT.Wilcoxon")
  tkgrid(getFrame(variable1Var), labelRcmdr(variablesFrame, text="    "), getFrame(variable2Var), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  statsFrame <- tkframe(statisticsTab)
  radioButtons(statsFrame,name="alternativa", buttons=c("twosided", "less", "greater"),
               values=c("two.sided", "less", "greater"),
               labels=gettextRcmdr(c("Bilateral", "Diferencia < 0", "Diferencia > 0")),
               title=gettextRcmdr("Hipotesis alternativa"), initialValue= dialog.valores$alt.inicial)
  ConfIntFrame <- tkframe(statsFrame)
  ConfIntVariable <- tclVar(dialog.valores$confint.inicial)
  ConfIntField <- ttkentry(ConfIntFrame, width = "6", 
                           textvariable = ConfIntVariable)  
  radioButtons(statsFrame,name="prueba", buttons=c("default", "exact", "normal", "correct"),
               values = c("default", "exact","normal", "correct"),
               labels=gettextRcmdr(c("Por defecto", "Exacta", "Aproximacion normal",
                                     "Aprox. normal con\ncorrecc. continuidad")), 
               title=gettextRcmdr("Tipo de Prueba"), initialValue= dialog.valores$prueba.inicial)
  tkgrid(labelRcmdr(ConfIntFrame, text = gettextRcmdr("Nivel de confianza"), 
                    fg = getRcmdr("title.color"), font="RcmdrTitleFont"), sticky = "w")
  tkgrid(ConfIntField, sticky = "w")
  tkgrid(alternativaFrame, labelRcmdr(statsFrame, text = "    "), 
         ConfIntFrame, labelRcmdr(statsFrame, text = "    "), 
         pruebaFrame, sticky = "nw")
  tkgrid(statsFrame, sticky = "nw")  
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

# Pruebas conformidad: forma y apuntamiento #

formafun <- function(data){
  variable <- unlist(strsplit(deparse(substitute(data)), "[$]"))[2]
  res <- list()
  res[[1]] <- variable
  n <- length(na.omit(data))
  sesgo <- (n*sum((data-mean(data,na.rm=TRUE))^3,na.rm=TRUE))/((n-1)*(n-2))/(sd(data,na.rm=TRUE)^3)
  res[[2]] <- round(sesgo,2)
  error.simetria <- sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)))
  res[[3]]<-round(error.simetria,2)
  sim.estandar <- sesgo/error.simetria
  res[[4]] <- round(sim.estandar,2)
  if (sign(sim.estandar) > 0) valor.p <- 2*pnorm(sim.estandar,lower.tail=FALSE)
  else valor.p <- 2*pnorm(sim.estandar)
  res[[5]] <- round(valor.p,3)
  class(res) <- "Zforma"
  res
}

print.Zforma <- function (x,digits = max(4,getOption("digits") - 4),...)
{
  cat('Coeficiente de asimetria para ',x[[1]],": ",x[[2]],"\n",sep='')
  cat("\n")
  cat('Error tipico del coeficiente de asimetria para ',x[[1]],": ",x[[3]],"\n",sep='')
  cat("\n")
  cat('Coeficiente de simetria estandarizado para ',x[[1]],": ",x[[4]],"\n",sep='')
  cat("\n")
  cat('Significacion estadistica (bilateral) asociada para ',x[[1]],": ",x[[5]],"\n",sep='')
  cat("\n")
  invisible(x)  
}

prueba.conformidad.forma <- function(){
  defecto <- list(var.inicial=NULL,echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("prueba.conformidad.forma",defecto) 
  initializeDialog(title=gettextRcmdr("Prueba conformidad parametro de forma"),
                   use.tabs=TRUE,tabs=c('dataTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                 initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))  
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1  
    var <- getSelection(variableVar)
    if (length(var) == 0) {
      errorCondition(recall=prueba.conformidad.forma, message=gettextRcmdr("Debe seleccionar una variable."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    opts <- options(warn=-1)
    vars <- var 
    bd <- paste(ActiveDataSet(), "$", vars, sep="")
    options(opts)
    putDialog("prueba.conformidad.forma",list(var.inicial=var,echo.inicial=echocodigo,creahtml.inicial=creahtml,
                                              tab.inicial=tab))  
    instruccion1 <- paste("pruebaZ.forma <- formafun(data=",bd,")", sep="")
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    justDoIt(instruccion1)
    doItAndPrint(paste("pruebaZ.forma # Prueba de conformidad respecto al parametro de forma para", var)) 
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Prueba de conformidad respecto al parametro de forma para: ",var,
                      sep="")
      HTML(as.title(titulo),file=.archivo)
      HTML('Prueba de Conformidad para el parametro de forma',file=.archivo)
      HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
      if (pruebaZ.forma[[5]] > 0.0001) pval<-paste(" = ",round(pruebaZ.forma[[5]],4),sep="") else pval<- " < 0.0001"
      HTML(paste("Coeficiente de simetria = ",round(pruebaZ.forma[[2]],4),", Error tipico = ",round(pruebaZ.forma[[3]],4),
                 sep=""),file=.archivo)
      HTML(paste("Coeficiente de simetria estandarizado = ",round(pruebaZ.forma[[4]],4),", valor p",pval,sep=""),
           file=.archivo)
      HTML("Hipotesis Alternativa: Parametro beta distinto de 0", file=.archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove('pruebaZ.forma')")
    remove('pruebaZ.forma', envir=.GlobalEnv)
    closeDialog()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RCmdrPlugin.EACSPIR",reset="prueba.conformidad.forma",apply="prueba.conformidad.forma")
  tkgrid(getFrame(variableVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Datos","Opciones"))
}

apuntfun <- function(data){
  variable <- unlist(strsplit(deparse(substitute(data)), "[$]"))[2]
  res <- list()
  res[[1]] <- variable
  n <- length(na.omit(data))
  curtosis <-c(n*(n+1)*sum((data-mean(data,na.rm=TRUE))^4,na.rm=TRUE)/((n-1)*(n-2)*(n-3))-3*sum((data-
                                                                                                   mean(data,na.rm=TRUE))^2,na.rm=TRUE)^2/((n-2)*(n-3)))/(sd(data,na.rm=TRUE)^4)
  res[[2]] <- round(curtosis,2)
  error.curtosis <- sqrt((24*n*(n-1)^2)/((n-3)*(n-2)*(n+3)*(n+5)))
  res[[3]]<-round(error.curtosis,2)
  curtosis.estandar <- curtosis/error.curtosis
  res[[4]] <- round(curtosis.estandar,2)
  if (sign(curtosis.estandar) > 0) valor.p <- 2*pnorm(curtosis.estandar,lower.tail=FALSE)
  else valor.p <- 2*pnorm(curtosis.estandar)
  res[[5]] <- round(valor.p,3)
  class(res) <- "Zapunt"
  res  
}

print.Zapunt <- function (x,digits = max(4,getOption("digits") - 4),...)
{
  cat('Coeficiente de apuntamiento para ',x[[1]],": ",x[[2]],"\n",sep='')
  cat("\n")
  cat('Error tipico del coeficiente de apuntamiento para ',x[[1]],": ",x[[3]],"\n",sep='')
  cat("\n")
  cat('Coeficiente de apuntamiento estandarizado para ',x[[1]],": ",x[[4]],"\n",sep='')
  cat("\n")
  cat('Significacion estadistica (bilateral) asociada para ',x[[1]],": ",x[[5]],"\n",sep='')
  cat("\n")
  invisible(x)  
}

prueba.conformidad.apuntamiento <- function(){
  defecto <- list(var.inicial=NULL,echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("prueba.conformidad.apuntamiento",defecto) 
  initializeDialog(title=gettextRcmdr("Prueba conformidad parametro de apuntamiento"),
                   use.tabs=TRUE,tabs=c('dataTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                 initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))   
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1  
    var <- getSelection(variableVar)
    if (length(var) == 0) {
      errorCondition(recall=prueba.conformidad.apuntamiento, message=gettextRcmdr("Debe seleccionar una variable."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    opts <- options(warn=-1)
    vars <- var 
    bd <- paste(ActiveDataSet(), "$", vars, sep="")
    options(opts)
    putDialog("prueba.conformidad.apuntamiento",list(var.inicial=var,echo.inicial=echocodigo,creahtml.inicial=creahtml,
                                                     tab.inicial=tab))  
    instruccion1 <- paste("pruebaZ.apuntamiento <- apuntfun(data=",bd,")", sep="")
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    justDoIt(instruccion1)
    doItAndPrint(paste("pruebaZ.apuntamiento # Prueba de conformidad respecto al parametro de apuntamiento para", var)) 
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Prueba de conformidad respecto al parametro de apuntamiento para: ",var,
                      sep="")
      HTML(as.title(titulo),file=.archivo)
    }
    closeDialog()
    
    if (creahtml == 1)
    {
      HTML('Prueba de Conformidad para el parametro de apuntamiento',file=.archivo)
      HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
      if (pruebaZ.apuntamiento[[5]] > 0.0001) pval<-paste(" = ",round(pruebaZ.apuntamiento[[5]],4),sep="") 
      else pval<- " < 0.0001"
      HTML(paste("Coeficiente de apuntamiento = ",round(pruebaZ.apuntamiento[[2]],4),", Error tipico = ",
                 round(pruebaZ.apuntamiento[[3]],4),
                 sep=""),file=.archivo)
      HTML(paste("Coeficiente de apuntamiento estandarizado = ",round(pruebaZ.apuntamiento[[4]],4),
                 ", valor p",pval,sep=""),file=.archivo)
      HTML("Hipotesis Alternativa: Parametro gamma distinto de 0", file=.archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove('pruebaZ.apuntamiento')")
    remove('pruebaZ.apuntamiento', envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RCmdrPlugin.EACSPIR",reset="prueba.conformidad.apuntamiento",
               apply="prueba.conformidad.apuntamiento")
  tkgrid(getFrame(variableVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Datos","Opciones"))
}

# Determinacion del tamano de la muestra: media y proporcion #

tam.muestra <-function(conf,e,est=c("proporcion","media"),inf=TRUE,N=NULL,pi=0.5,sigma=NULL){
    if ((conf >= 1) || (conf <= 0)) stop("Especifique el valor para el nivel de confianza");
    alfa <- 1-conf
    if (!is.numeric(e)) stop("Especifique el valor de precision");
    if ( (inf==FALSE) && !is.numeric(N) ) stop("Especifique el tamano de la poblacion")
    est <- match.arg(est)
    if (est == "proporcion"){
      if ((pi > 1) || (pi < 0)) stop("Especifique el valor para la proporcion poblacional");
      if (inf == TRUE) n <- qnorm(alfa/2)^2*pi*(1-pi)/e^2
      else n <- qnorm(alfa/2)^2*pi*(1-pi)*N/(e^2*(N-1)+qnorm(alfa/2)^2*pi*(1-pi))
      res <- list(parametro=est,prop=pi,precision=e,confianza=conf,muestra=n)
      class(res) <- "tmuestra"
    }
    if (est == "media"){
      if (!is.numeric(sigma)) stop("Especifique el valor para la desviacion tipica poblacional");
      if (inf == TRUE) n <- qnorm(alfa/2)^2*sigma^2/e^2
      else n <- qnorm(alfa/2)^2*sigma^2*N/(e^2*(N-1)+qnorm(alfa/2)^2*sigma^2)
      res <- list(parametro=est,desv=sigma,precision=e,confianza=conf,muestra=n)
      class(res) <- "tmuestra"
    }
    res
}

print.tmuestra <- function (x,digits = max(4,getOption("digits") - 4),...)
{
    if (x$parametro == 'proporcion')
    {
      cat("\n")
      cat("Parametro a estimar: ",x$parametro, " = ", x$prop, "\n")
      cat("Precision: ", x$precision, "\n")
      cat("Nivel de confianza: ", x$confianza, "\n")
      cat("Tamano de la muestra requerido: ", x$muestra, "observaciones", "\n")
    }  
    if (x$parametro == 'media')
    {
      cat("\n")      
      cat("Parametro a estimar: ",x$parametro, "\n")
      cat("Sigma = ", x$desv, "\n")
      cat("Precision: ", x$precision, "\n")
      cat("Nivel de confianza: ", x$confianza, "\n")
      cat("Tamano de la muestra requerido: ", x$muestra, "observaciones", "\n")
    }    
    cat("\n")
    invisible(x)  
}

determ.tam.proporcion <- function(){
  defecto <- list(prop.inicial="0.5",pobN.inicial="<auto>",precision.inicial="0.05",confint.inicial=".95",
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("determ.tam.proporcion",defecto) 
  initializeDialog(title=gettextRcmdr("Determinacion tamano muestra: proporciones"),
                   use.tabs=TRUE,tabs=c('statisticsTab','optionsTab'))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == statisticsTab$ID) 0 else 1         
    proporcion <- as.numeric(tclvalue(propVariable))
    pobN <- tclvalue(pobNVariable)
    precision <- as.numeric(tclvalue(precisionVariable))
    echocodigo <- tclvalue(echocodigoVariable)
    confianza <- as.numeric(tclvalue(ConfIntVariable))
    creahtml <- tclvalue(creahtmlVariable)
    if ((proporcion > 1) || (proporcion < 0)) {
      errorCondition(recall=determ.tam.proporcion,
                     message=gettextRcmdr("Especifique un valor del parametro proporcion."))
      return()
    }
    if (!is.numeric(precision)) {
      precision <- 0.05
      Message(message=gettextRcmdr("Valor de precision invalido, se utilizara valor por defecto."))
      return()
    }
    if (pobN == gettextRcmdr("<auto>")) N <- NULL
    else { 
      N <- as.numeric(pobN)
      if ( !is.numeric(N)) {
        N <- NULL
        Message(message=gettextRcmdr("Tamano de la poblacion invalido, se utilizara valor por defecto."))
        return()
      }
    }
    if ( confianza < .0 || confianza > 1. || !is.numeric(confianza) )
    {
      confianza <- 0.95
      Message(message=gettextRcmdr("Nivel de confianza invalido, se utilizara valor por defecto."),
              type="warning")              
    }
    putDialog("determ.tam.proporcion",list(prop.inicial=proporcion,pobN.inicial=if (is.null(N)) "<auto>" else N,
                                           precision.inicial=precision,confint.inicial=confianza,
                                           echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- "Determinacion del tamano de la muestra: proporciones"
      HTML(as.title(titulo),file=.archivo)
    }
    closeDialog()
    if (is.numeric(N))
      instruccion <- paste("tam.muestra(conf=",confianza,",e=",precision,",est='proporcion'",
                           ",inf=FALSE,N=",N,",pi=",proporcion,")",sep="")        
    else
      instruccion <- paste("dtm <- tam.muestra(conf=",confianza,",e=",precision,",est='proporcion'",
                           ",inf=TRUE,N=NULL,pi=",proporcion,")",sep="")                
    justDoIt(instruccion)
    if (echocodigo==1) logger(instruccion)
    doItAndPrint("dtm # Determinacion tamano de la muestra: proporcion ")        
    if (creahtml == 1)
    {
      HTML('Determinacion del tamano de la muestra: Caso proporcion',file=.archivo)
      HTML(dtm,file= .archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove(list=c('dtm'))")
    remove(list=c('dtm'), envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RcmdrPugin.EACSPIR",reset="determ.tam.proporcion",apply="determ.tam.proporcion")
  propFrame <- tkframe(statisticsTab)
  propVariable <- tclVar(dialog.valores$prop.inicial)
  propField <- ttkentry(propFrame, width="8", textvariable=propVariable)
  tkgrid(labelRcmdr(propFrame,text=gettextRcmdr("Parametro proporcion = "),font="RcmdrTitleFont"),
         propField,sticky="w")
  pobNFrame <- tkframe(statisticsTab)
  pobNVariable <- tclVar(dialog.valores$pobN.inicial)
  pobNField <- ttkentry(pobNFrame, width="8", textvariable=pobNVariable)
  tkgrid(labelRcmdr(pobNFrame,text=gettextRcmdr("Tamano poblacional = "),font="RcmdrTitleFont"),
         pobNField,sticky="w")
  precisionFrame <- tkframe(statisticsTab)
  precisionVariable <- tclVar(dialog.valores$precision.inicial)
  precisionField <- ttkentry(precisionFrame, width="8", textvariable=precisionVariable)
  tkgrid(labelRcmdr(precisionFrame,text=gettextRcmdr("Valor de precision = "),font="RcmdrTitleFont"),
         precisionField,sticky="w")  
  ConfIntFrame <- tkframe(statisticsTab)
  ConfIntVariable <- tclVar(dialog.valores$confint.inicial)
  ConfIntField <- ttkentry(ConfIntFrame, width="8", textvariable=ConfIntVariable)
  tkgrid(labelRcmdr(ConfIntFrame,text=gettextRcmdr("Nivel de Confianza = "),font="RcmdrTitleFont"),
         ConfIntField,sticky="w")  
  tkgrid(propFrame, sticky="nw")
  tkgrid(pobNFrame, sticky="nw")
  tkgrid(precisionFrame, sticky="nw")
  tkgrid(ConfIntFrame, sticky="nw")
  tkgrid(opsFrame, sticky="nw")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('statisticsTab','optionsTab'),
               tab.names=c("Estadisticos","Opciones"))  
}

determ.tam.media <- function(){
  defecto <- list(sigma.inicial="1.0",pobN.inicial="<auto>",precision.inicial="0.10",confint.inicial=".95",
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("determ.tam.media",defecto) 
  initializeDialog(title=gettextRcmdr("Determinacion tamano muestra: medias"),
                   use.tabs=TRUE,tabs=c('statisticsTab','optionsTab'))
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == statisticsTab$ID) 0 else 1    
    sigma <- as.numeric(tclvalue(sigmaVariable))
    pobN <- tclvalue(pobNVariable)
    precision <- as.numeric(tclvalue(precisionVariable))
    echocodigo <- tclvalue(echocodigoVariable)
    confianza <- as.numeric(tclvalue(ConfIntVariable))
    creahtml <- tclvalue(creahtmlVariable)
    if (!is.numeric(sigma)) {
      errorCondition(recall=determ.tam.media,
                     message=gettextRcmdr("Especifique un valor del parametro sigma."))
      return()
    }
    if (!is.numeric(precision)) {
      precision <- 0.05
      Message(message=gettextRcmdr("Valor de precision invalido, se utilizara valor por defecto."))
      return()
    }
    if (pobN == gettextRcmdr("<auto>")) N <- NULL
    else { 
      N <- as.numeric(pobN)
      if ( !is.numeric(N)) {
        N <- NULL
        Message(message=gettextRcmdr("Tamano de la poblacion invalido, se utilizara valor por defecto."))
        return()
      }
    }
    if ( confianza < .0 || confianza > 1. || !is.numeric(confianza) )
    {
      confianza <- 0.95
      Message(message=gettextRcmdr("Nivel de confianza invalido, se utilizara valor por defecto."),
              type="warning")              
    }
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- "Determinacion del tamano de la muestra: medias"
      HTML(as.title(titulo),file=.archivo)
    }
    putDialog("determ.tam.media",list(sigma.inicial=sigma,pobN.inicial=if (is.null(N)) "<auto>" else N,
                                      precision.inicial=precision,confint.inicial=confianza,
                                      echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))    
    closeDialog()
    if (is.numeric(N))
      instruccion <- paste("dtm <- tam.muestra(conf=",confianza,",e=",precision,",est='media'",
                           ",inf=FALSE,N=",N,",sigma=",sigma,")",sep="")        
    else
      instruccion <- paste("dtm <- tam.muestra(conf=",confianza,",e=",precision,",est='media'",
                           ",inf=TRUE,N=NULL,sigma=",sigma,")",sep="")                
    justDoIt(instruccion)
    if (echocodigo==1) logger(instruccion)
    doItAndPrint("dtm # Determinacion tamano de la muestra: media ")        
    if (creahtml == 1)
    {
      HTML('Determinacion del tamano de la muestra: Caso media',file=.archivo)
      HTML(dtm,file= .archivo)
      HTMLhr(file = .archivo)
    }
    if (echocodigo == 1) logger("remove(list=c('dtm'))")
    remove(list=c('dtm'), envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RcmdrPugin.EACSPIR",reset="determ.tam.media",apply="determ.tam.media")
  sigmaFrame <- tkframe(statisticsTab)
  sigmaVariable <- tclVar(dialog.valores$sigma.inicial)
  sigmaField <- ttkentry(sigmaFrame, width="8", textvariable=sigmaVariable)
  tkgrid(labelRcmdr(sigmaFrame,text=gettextRcmdr("Parametro sigma = "),font="RcmdrTitleFont"),
         sigmaField,sticky="w")
  pobNFrame <- tkframe(statisticsTab)
  pobNVariable <- tclVar(dialog.valores$pobN.inicial)
  pobNField <- ttkentry(pobNFrame, width="8", textvariable=pobNVariable)
  tkgrid(labelRcmdr(pobNFrame,text=gettextRcmdr("Tamano poblacional = "),font="RcmdrTitleFont"),
         pobNField,sticky="w")
  precisionFrame <- tkframe(statisticsTab)
  precisionVariable <- tclVar(dialog.valores$precision.inicial)
  precisionField <- ttkentry(precisionFrame, width="8", textvariable=precisionVariable)
  tkgrid(labelRcmdr(precisionFrame,text=gettextRcmdr("Valor de precision = "),font="RcmdrTitleFont"),
         precisionField,sticky="w")  
  ConfIntFrame <- tkframe(statisticsTab)
  ConfIntVariable <- tclVar(dialog.valores$confint.inicial)
  ConfIntField <- ttkentry(ConfIntFrame, width="8", textvariable=ConfIntVariable)
  tkgrid(labelRcmdr(ConfIntFrame,text=gettextRcmdr("Nivel de Confianza = "),font="RcmdrTitleFont"),
         ConfIntField,sticky="w")  
  tkgrid(sigmaFrame, sticky="nw")
  tkgrid(pobNFrame, sticky="nw")
  tkgrid(precisionFrame, sticky="nw")
  tkgrid(ConfIntFrame, sticky="nw")
  tkgrid(opsFrame, sticky="nw")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('statisticsTab','optionsTab'),
               tab.names=c("Estadisticos","Opciones"))
}

# Ajuste a la normal: pruebas de bondad de ajuste y graficos #

normfun <- function(variable,stats,group=NULL){
  res <- list()
  res[[1]] <- variable
  res[[2]] <- group
  res[[3]] <- stats
  i <- 4 
  opts <- options(warn=-1)  
  if ('Pearsontest' %in% stats)
  { 
    if (!is.null(group))
      .prueba.pearson <- eval(parse(text=paste("by(",variable,",",group,",pearson.test)",sep="")))
    else
      .prueba.pearson <- eval(parse(text=paste("pearson.test(",variable,")",sep="")))
    if (!is.null(group)){                   
      for (j in 1:dim(.prueba.pearson))
        .prueba.pearson[[j]]$data.name <- variable
    }
    res[[i]] <- .prueba.pearson
    i <- i + 1
  }
  if ('ADtest' %in% stats)
  { 
    if (!is.null(group)){
      .prueba.AD <- eval(parse(text=paste("by(",variable,",",group,",ad.test)",sep="")))
      for (j in 1:dim(.prueba.AD))
        .prueba.AD[[j]]$data.name <- variable      
    }
    else
      .prueba.AD <- eval(parse(text=paste("ad.test(",variable,")",sep="")))
    res[[i]] <- .prueba.AD
    i <- i + 1  
  }
  if ('KStest' %in% stats)
  { 
    if (!is.null(group)){
      .prueba.KS <- eval(parse(text=paste("by(",variable,","
                                          ,group,",function(x)ks.test(x,'pnorm',mean(x,na.rm=T),sd(x,na.rm=T)))"
                                          ,sep="")))
      for (j in 1:dim(.prueba.KS))
        .prueba.KS[[j]]$data.name <- variable        
    }
    else
      .prueba.KS <- eval(parse(text=paste("ks.test(",variable,",'pnorm',mean(",
                                          variable,",na.rm=T),sd(",variable,",na.rm=T))",sep="")))
    res[[i]] <- .prueba.KS
    i <- i + 1    
  }
  if ('Shapirotest' %in% stats)
  { 
    if (!is.null(group)){
      .prueba.shapiro <- eval(parse(text=paste(".prueba.shapiro <- by(",variable,",",group,",shapiro.test)",sep="")))
      for (j in 1:dim(.prueba.shapiro))
        .prueba.shapiro[[j]]$data.name <- variable        
    }  
    else
      .prueba.shapiro <- eval(parse(text=paste("shapiro.test(",variable,")",sep="")))
    res[[i]] <- .prueba.shapiro
    i <- i + 1
  }
  class(res) <- "prnorm"
  res
}

print.prnorm <- function(x,...){
  j<-4
  if ("Pearsontest" %in% x[[3]]){
    cat("# Prueba de ajuste a la normal: Ji-Cuadrado de Pearson : \n\n")
    print(x[[j]])  
    cat("\n\n")
    j <- j+1
  }
  if ("ADtest" %in% x[[3]]){
    cat("# Prueba de ajuste a la normal: Anderson-Darling : \n\n")
    print(x[[j]])  
    cat("\n\n")
    j <- j+1
  }
  if ("KStest" %in% x[[3]]){
    cat("# Prueba de ajuste a la normal: Kolmogorov-Smirnov : \n\n")
    print(x[[j]])  
    cat("\n\n")
    j <- j+1
  }
  if ("Shapirotest" %in% x[[3]]){
    cat("# Prueba de ajuste a la normal: Shapiro-Wilk : \n\n")
    print(x[[j]])  
    cat("\n\n")
    j <- j+1
  }
  
  invisible(x)  
}

normgrfun <- function(datos){
  par(mfrow=c(2,2))
  variable <- unlist(strsplit(deparse(substitute(datos)), "[$]"))[2]
  datos <- as.vector(na.omit(datos))
  plot(density(datos),main=paste('Densidad suavizada para ',variable,sep=""),ylab='Densidades')
  box()
  x <- datos
  h <- hist(x,freq=TRUE,plot=FALSE)
  plot(h,col='red', xlab='Intervalos', ylab='Frecuencias',main=paste('Histograma para ',variable,sep=""))
  xfit<-seq(min(x,na.rm=T),max(x,na.rm=T),length=1000)
  yfit<-dnorm(xfit,mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))
  yfit<-yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col='blue', lwd=2)
  box()
  z <- as.numeric(scale(x))
  qqnorm(z,xlab='Cuantilas teoricas',ylab='Cuantilas empiricas',main=paste('Grafico QQ para ',variable,sep=""))
  abline(0,1)
  box()
  plot(sort(z),pnorm(sort(z)),type='l',col='red',
       main=paste("Grafico de cuantilas para ",variable,sep=""),xlab='Puntuaciones Z',
       ylab='Probabilidad Acumulada')
  plot(ecdf(z),add=TRUE)
  box()
}

pruebas.normalidad <- function(){
  defecto <- list(var.inicial=NULL,pearson.inicial="0",ks.inicial="0",ad.inicial="0",shapiro.inicial="0",
                  grupo.inicial=NULL,graficos.inicial="0",echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("pruebas.normalidad",defecto)
  initializeDialog(title=gettextRcmdr("Pruebas de Ajuste Distribucion Normal"),
                   use.tabs=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Variables (escoja una)"),
                                 initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))
  groupsBox(recall=pruebas.normalidad, label=gettextRcmdr("Pruebas segun:"), 
            initialLabel=if (is.null(dialog.valores$grupo.inicial)) gettextRcmdr("Pruebas segun grupos") 
            else paste(gettextRcmdr("Pruebas segun:"), dialog.valores$grupo.inicial), 
            initialGroup=dialog.valores$grupo.inicial, window = dataTab)  
  checkBoxes(statisticsTab,frame="statsFrame",boxes=c("pearson","ks","ad","shapiro","graficos"),
             initialValues=c(dialog.valores$pearson.inicial,dialog.valores$ks.inicial,
                             dialog.valores$ad.inicial,dialog.valores$shapiro.inicial,
                             dialog.valores$graficos.inicial),
             labels=gettextRcmdr(c("Prueba Ji-Cuadrado de Pearson ",
                                   "Prueba de Kolmogorov-Smirnov ",
                                   "Prueba Anderson-Darling ",
                                   "Prueba Shapiro-Wilk ",
                                   "Representaciones graficas ")), 
             title = gettextRcmdr("Pruebas de normalidad"))  
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))     
  onOK <- function(){
    var <- getSelection(variableVar)
    if (length(var) == 0) {
      errorCondition(recall=pruebas.normalidad, message=gettextRcmdr("Debe seleccionar una variable."))
      return()
    }
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    pearsonval <- tclvalue(pearsonVariable)
    adval <- tclvalue(adVariable)
    ksval <-tclvalue(ksVariable)
    shapiroval <- tclvalue(shapiroVariable)
    grafval <- tclvalue(graficosVariable)
    seleccion <- as.numeric(pearsonval) + as.numeric(adval) + 
      as.numeric(ksval) + as.numeric(shapiroval) + as.numeric(grafval)
    if (seleccion == 0){
      errorCondition(recall=pruebas.normalidad, 
                     message=gettextRcmdr("Debe escoger algun indicador."))
      return()
    }
    putDialog("pruebas.normalidad",list(grupo.inicial=if (.groups != FALSE) .groups else NULL,var.inicial=var,
                                        pearson.inicial=pearsonval,ad.inicial=adval,ks.inicial=ksval,
                                        shapiro.inicial=shapiroval,graficos.inicial=grafval,
                                        echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))
    .BaseDatosActiva <- ActiveDataSet()
    if (.groups != FALSE) {
      grupos <- paste(.BaseDatosActiva, "$", .groups, sep="")
    } 
    variable <- paste(.BaseDatosActiva, "$", var, sep="")    
    stats <- paste("c(",
                   paste(c('"Pearsontest"', '"ADtest"', '"KStest"', '"Shapirotest"')
                         [c(pearsonval, adval, ksval, shapiroval) == 1],collapse=", "), ")", sep="")
    if (.groups != FALSE)
      instruccion1 <- paste(".norm.test <- normfun(var='",variable, "', stats=", stats, ", group='",grupos, "')", sep="")
    else
      instruccion1 <- paste(".norm.test <- normfun(var='",variable, "', stats=", stats, ")", sep="")  
    justDoIt(instruccion1)
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(".norm.test  # Pruebas de normalidad ")    
    
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- paste("Pruebas de ajuste a la distribucion normal para: ",var,
                      sep="")
      HTML(as.title(titulo),file=.archivo)
      j <- 4
      if (pearsonval == 1)
      { 
        if (.groups == FALSE){
          HTML('Prueba de ajuste a la normal: Ji-Cuadrado de Pearson ',file=.archivo)
          HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
          if (.norm.test[[j]][[2]] > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[2]],4),sep="") else pval<- " < 0.0001"
          HTML(paste("Ji Cuadrado de Pearson = ",round(.norm.test[[j]][[1]],4),", valor p",pval,sep=""),
               file=.archivo)
          j <- j+1
          HTMLhr(file = .archivo) 
        }
        else {
          nombres <- eval(parse(text=paste("levels(",.norm.test[[2]][[1]],")",sep='')))
          HTML('Prueba de ajuste a la normal: Ji-Cuadrado de Pearson ',file=.archivo)
          for (i in 1:dim(.norm.test[[j]])){
            HTML(paste("datos: ",ActiveDataSet(),"$",var," y ",grupos," = ",nombres[i],sep=""), file=.archivo)
            if (.norm.test[[j]][[i]]$p.value > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[i]]$p.value,4),sep="") 
            else pval<- " < 0.0001"
            HTML(paste("Ji Cuadrado de Pearson = ",round(.norm.test[[j]][[i]]$statistic,4),", valor p",pval,sep=""),
                 file=.archivo)
            HTML("----------------------------------------------------------------------------------",file=.archivo) 
          }
          j <- j+1
          HTMLhr(file = .archivo)        
        }
      }
      if (adval == 1){
        if (.groups == FALSE){
          HTML('Prueba de ajuste a la normal: Anderson-Darling ',file=.archivo)
          HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
          if (.norm.test[[j]][[2]] > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[2]],4),sep="") else pval<- " < 0.0001"
          HTML(paste("A = ",round(.norm.test[[j]][[1]],4),", valor p",pval,sep=""),
               file=.archivo)
          j <- j+1
          HTMLhr(file = .archivo) 
        }
        else {
          nombres <- eval(parse(text=paste("levels(",.norm.test[[2]][[1]],")",sep='')))
          HTML('Prueba de ajuste a la normal: Anderson-Darling ',file=.archivo)
          for (i in 1:dim(.norm.test[[j]])){
            HTML(paste("datos: ",ActiveDataSet(),"$",var," y ",grupos," = ",nombres[i],sep=""), file=.archivo)
            if (.norm.test[[j]][[i]]$p.value > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[i]]$p.value,4),sep="") 
            else pval<- " < 0.0001"
            HTML(paste("A = ",round(.norm.test[[j]][[i]]$statistic,4),", valor p",pval,sep=""),
                 file=.archivo)
            HTML("----------------------------------------------------------------------------------",file=.archivo) 
          }
          j <- j+1
          HTMLhr(file = .archivo)        
        }
      }
      if (ksval == 1){
        if (.groups == FALSE){
          HTML('Prueba de ajuste a la normal: Kolmogorov-Smirnov ',file=.archivo)
          HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
          if (.norm.test[[j]][[2]] > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[2]],4),sep="") else pval<- " < 0.0001"
          HTML(paste("D = ",round(.norm.test[[j]][[1]],4),", valor p",pval,sep=""),
               file=.archivo)
          j <- j+1
          HTMLhr(file = .archivo) 
        }
        else {
          nombres <- eval(parse(text=paste("levels(",.norm.test[[2]][[1]],")",sep='')))
          HTML('Prueba de ajuste a la normal: Kolmogorov-Smirnov ',file=.archivo)
          for (i in 1:dim(.norm.test[[j]])){
            HTML(paste("datos: ",ActiveDataSet(),"$",var," y ",grupos," = ",nombres[i],sep=""), file=.archivo)
            if (.norm.test[[j]][[i]]$p.value > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[i]]$p.value,4),sep="") 
            else pval<- " < 0.0001"
            HTML(paste("D = ",round(.norm.test[[j]][[i]]$statistic,4),", valor p",pval,sep=""),
                 file=.archivo)
            HTML("----------------------------------------------------------------------------------",file=.archivo) 
          }
          j <- j+1
          HTMLhr(file = .archivo)        
        }
      }
      if (shapiroval == 1){
        if (.groups == FALSE){
          HTML('Prueba de ajuste a la normal: Shapiro-Wilk ',file=.archivo)
          HTML(paste("datos: ",ActiveDataSet(),"$",var,sep=""), file=.archivo)
          if (.norm.test[[j]][[2]] > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[2]],4),sep="") else pval<- " < 0.0001"
          HTML(paste("W = ",round(.norm.test[[j]][[1]],4),", valor p",pval,sep=""),
               file=.archivo)
          j <- j+1
          HTMLhr(file = .archivo) 
        }
        else {
          nombres <- eval(parse(text=paste("levels(",.norm.test[[2]][[1]],")",sep='')))
          HTML('Prueba de ajuste a la normal: Shapiro-Wilk ',file=.archivo)
          for (i in 1:dim(.norm.test[[j]])){
            HTML(paste("datos: ",ActiveDataSet(),"$",var," y ",grupos," = ",nombres[i],sep=""), file=.archivo)
            if (.norm.test[[j]][[i]]$p.value > 0.0001) pval<-paste(" = ",round(.norm.test[[j]][[i]]$p.value,4),sep="") 
            else pval<- " < 0.0001"
            HTML(paste("W = ",round(.norm.test[[j]][[i]]$statistic,4),", valor p",pval,sep=""),
                 file=.archivo)
            HTML("----------------------------------------------------------------------------------",file=.archivo) 
          }
          j <- j+1
          HTMLhr(file = .archivo)        
        }
      }      
    }
    closeDialog()
    
    if (.groups != FALSE & grafval ==1){
      grafval <- 0.95
      Message(message=gettextRcmdr("Analisis segun variable de agrupacion: No se realizaran graficos agrupados."),
              type="warning")            
    }
    if (grafval == 1){
      instruccion1 <- paste("normgrfun(datos=",variable,")", sep="")
      justDoIt(instruccion1)
      if (echocodigo == 1)
      {
        logger(instruccion1)
      }    
      if (creahtml == 1)
      {
        titulo <- paste("Graficas para variable ",var,sep="")
        HTML(as.title(titulo),file=.archivo)
        nombre.archivo <- paste("GraficasNormR",gsub(":","",substr(Sys.time(),12,19)),
                                ".jpg",sep="")
        dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                  width=500, height=500)
        HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
        HTMLhr(file = .archivo)
      }      
    }
    
    remove(.norm.test, envir=.GlobalEnv)    
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RcmdrPlugin.EACSPIR",reset="pruebas.normalidad",apply="pruebas.normalidad")
  tkgrid(getFrame(variableVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(statsFrame, sticky="nw")
  tkgrid(opsFrame, sticky="w")
  tkgrid(groupsFrame, sticky = "w", padx=6)  
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

avarMRfun <- function(.baseDatosActiva, idvar, intra, mr, sctype,descrip){
  res <- list()
  res[[1]] <- .baseDatosActiva
  res[[2]] <- descrip
  datos.ANOVAMR <- eval(parse(text=paste("ezANOVA(.baseDatosActiva,dv=.(",noquote(mr),"),wid=.(",noquote(idvar),"),
                                         within=.(",noquote(intra),"),type=",sctype,",detailed=TRUE)")))
  res[[3]] <- datos.ANOVAMR[[1]]
  j <- 4
  if (length(levels(.baseDatosActiva[[2]])) >= 3){  
    res[[4]] <- datos.ANOVAMR[[2]]
    res[[5]] <- datos.ANOVAMR[[3]]
    j <- 6
  }
  if (descrip == 1)
  {
    descriptivo.ANOVAMR <- eval(parse(text=paste("ezStats(.baseDatosActiva,dv=.(",noquote(mr),"),
                                                 wid=.(",noquote(idvar),"),within=.(",noquote(intra),"))")))
    res[[j]] <- descriptivo.ANOVAMR     
  }
  class(res) <- "avarMR"
  res
}

print.avarMR <- function(x,...){
  cat('Resumen ANOVA: \n\n')
  print(x[[3]])
  cat('\n\n')
  j <- 4
  if (length(levels(x[[1]][[2]])) >= 3){
    cat('Prueba Esfericidad de Mauchly: \n\n')
    print(x[[4]])
    cat('\n\n')
    cat('Correcciones Esfericidad: \n\n')
    print(x[[5]])  
    cat("\n\n")
    j <- 6
  }
  
  if (x[[2]] == TRUE){
    cat('Descriptivos ANOVA: \n\n')
    print(x[[j]])
    cat("\n\n")    
  }
  invisible(x) 
}

avarMRgraf <- function(.baseDatosActiva,idvar,intra,mr){
  grafico.Medias <- eval(parse(text=paste("ezPlot(.baseDatosActiva,dv=.(",noquote(mr),"),
                                          wid=.(",noquote(idvar),"),within=.(",noquote(intra),"),
                                          x = .(",noquote(intra),"), do_lines = FALSE,x_lab = '",intra,
                                          "',y_lab = '",mr,"')")))
  print(grafico.Medias)
}

avar.MR <- function(){
  defecto <- list(var.inicial=NULL,ID.inicial=NULL,intra.inicial="f.intrasujeto",medrep.inicial="medida.rep",
                  descrip.inicial="0",grafico.inicial="0",SC.inicial="3",
                  echo.inicial="0",creahtml.inicial="0",tab.inicial=0)
  dialog.valores <- getDialog("avar.MR",defecto)
  initializeDialog(title=gettextRcmdr("ANOVA Medidas Repetidas"),
                   use.tabs=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'))
  variablesFrame <- tkframe(dataTab)
  variableVar <- variableListBox(variablesFrame, Numeric(), title=gettextRcmdr("Medidas repetidas (escoja dos o mas)"),
                                 selectmode="multiple", initialSelection=varPosn(dialog.valores$var.inicial,"numeric"))
  IDVar <- variableListBox(variablesFrame, 
                           title=gettextRcmdr("Variable de identificacion (Opcional)"),
                           initialSelection=varPosn(dialog.valores$ID.inicial))   
  checkBoxes(statisticsTab,frame="statsFrame",boxes=c("descrip","grafico"),
             initialValues=c(dialog.valores$descrip.inicial,dialog.valores$grafico.inicial),
             labels=gettextRcmdr(c("Mostrar descriptivos del ANOVA ",
                                   "Mostrar grafico de medias ")), 
             title = gettextRcmdr("Pruebas Descriptivas"))
  radioButtons(statisticsTab, name="SCuadrados", buttons=c("SC1Boton", "SC2Boton","SC3Boton"),
               values=c(1,2,3),labels=gettextRcmdr(c("Tipo I", "Tipo II", "Tipo III")),
               title=gettextRcmdr("Suma de Cuadrados:"),initialValue=dialog.valores$SC.inicial)                                                    
  checkBoxes(optionsTab,frame="opsFrame",boxes=c("echocodigo","creahtml"),
             initialValues=c(dialog.valores$echo.inicial,dialog.valores$creahtml.inicial),
             labels=gettextRcmdr(c("Mostrar en pantalla el codigo de R ejecutado ","Generar informe de resultados ")), 
             title = gettextRcmdr("Opciones"))
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1       
    var <- getSelection(variableVar)
    if (length(var) < 2) {
      errorCondition(recall=avar.MR, message=gettextRcmdr("Debe seleccionar 2 variables como minimo."))
      return()
    }
    idvar <-getSelection(IDVar)
    if (length(idvar)==0)
      var <- paste('"', var, '"', sep="")
    else
    {
      var <- paste('"', c(idvar,var), '"', sep="")          
    }       
    echocodigo <- tclvalue(echocodigoVariable)
    creahtml <- tclvalue(creahtmlVariable)
    descripval <- tclvalue(descripVariable)
    grafval <- tclvalue(graficoVariable)
    intra <- tclvalue(intraVariable)
    mr <- tclvalue(medrepVariable)
    sctype <- as.numeric(tclvalue(SCuadradosVariable))
    putDialog("avar.MR",list(var.inicial=getSelection(variableVar),
                             ID.inicial=if (length(idvar) != 0) idvar else NULL,
                             intra.inicial=intra,medrep.inicial=mr,
                             descrip.inicial=descripval,
                             grafico.inicial=grafval,SC.inicial=sctype,
                             echo.inicial=echocodigo,creahtml.inicial=creahtml,tab.inicial=tab))    
    if (length(idvar) == 0){
      instruccion <- paste(".baseDatosActiva <- na.omit(",ActiveDataSet(),"[,c(", paste(var, collapse=","),")])",sep="")
      justDoIt(instruccion)
      instruccion2 <- "ObsNumero <- as.factor(1:nrow(.baseDatosActiva))"
      justDoIt(instruccion2)
      instruccion3 <- ".baseDatosActiva <- cbind(.baseDatosActiva,ObsNumero)"
      justDoIt(instruccion3)
      instruccion4 <- paste(".baseDatosActiva <- melt.data.frame(.baseDatosActiva, id.vars=c('ObsNumero'), variable_name='",
                            intra,"')",sep='')
      justDoIt(instruccion4)
      idvar <- "ObsNumero"
      justDoIt(paste("colnames(.baseDatosActiva)[3] <- '",mr,"'",sep=''))        
    }
    else{
      instruccion <- paste(".baseDatosActiva <- na.omit(",ActiveDataSet(),"[,c(", paste(var, collapse=","),")])",sep="")
      justDoIt(instruccion)
      instruccion2 <- paste(".baseDatosActiva <- melt.data.frame(.baseDatosActiva, id.vars='",idvar,"',variable_name='",
                            intra,"')",sep='')
      justDoIt(instruccion2)
      justDoIt(paste("colnames(.baseDatosActiva)[3] <- '",mr,"'",sep=''))          
    } 
    instruccion1 <- paste(" .prueba.avarMR <- avarMRfun(.baseDatosActiva, idvar='", idvar, "',intra='",intra,
                          "',mr='",mr,"',sctype=",sctype,",descrip=",as.logical(as.numeric(descripval)),")", sep="")  
    justDoIt(instruccion1)
    if (echocodigo == 1)
    {
      logger(instruccion1)
    }
    doItAndPrint(".prueba.avarMR  # Analisis de la varianza de medidas repetidas ")
    if (grafval == 1){
      instruccion1 <- paste(" avarMRgraf(.baseDatosActiva, idvar='", idvar, "',intra='",intra,
                            "',mr='",mr,"')", sep="")  
      justDoIt(instruccion1)
      if (echocodigo == 1)
      {
        logger(instruccion1)
      }    
      
    }
    closeDialog()
    if (creahtml == 1)
    {
      if (!file.exists("Informe de Resultados.html"))
        .archivo <- HTMLInitFile(file.path(getwd()),
                                 "Informe de Resultados", BackGroundColor="#FFFFCC")
      else
        .archivo <- file.path(getwd(), "Informe de Resultados.html")
      titulo <- "ANOVA de Medidas Repetidas"
      HTML(as.title(titulo),file=.archivo)
      HTML("Resumen ANOVA",file=.archivo)
      HTML(as.data.frame(.prueba.avarMR[[3]]),file=.archivo)
      j <- 4
      if (length(levels(.baseDatosActiva[[2]])) >= 3){
        HTML("Prueba Esfericidad de Mauchly",file=.archivo)
        HTML(as.data.frame(.prueba.avarMR[[4]]),file=.archivo)
        HTML("Correcciones Esfericidad",file=.archivo)
        HTML(as.data.frame(.prueba.avarMR[[5]]),file=.archivo)
        j <- 6
      }
      
      if (descripval == 1)
      {
        HTML("Descripcion ANOVA",file=.archivo)
        HTML(as.data.frame(.prueba.avarMR[[j]]),file=.archivo)
      }
      if (grafval == 1)
      {
        HTML("Grafica de medias ",file=.archivo)
        nombre.archivo <- paste("GraficaMediasANOVAMR",gsub(":","",substr(Sys.time(),12,19)),
                                ".jpg",sep="")
        dev.print(jpeg, filename=paste(getwd(),"/",nombre.archivo,sep=""),
                  width=500, height=500)
        HTMLInsertGraph(nombre.archivo,file=.archivo,append=TRUE)
        HTMLhr(file = .archivo)
      }
      if (grafval != 1) HTMLhr(file = .archivo)          
    }
    remove(list=c('.baseDatosActiva','.prueba.avarMR'), envir=.GlobalEnv)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="RcmdrPlugin.EACSPIR",reset="avar.MR",apply="avar.MR")
  tkgrid(getFrame(variableVar), labelRcmdr(variablesFrame, text="    "), getFrame(IDVar), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  intraFrame <- tkframe(dataTab)
  intraVariable <- tclVar(dialog.valores$intra.inicial)
  intraField <- ttkentry(intraFrame, width="12", textvariable=intraVariable)
  tkgrid(labelRcmdr(intraFrame,text=gettextRcmdr("Nombre del factor Intrasujetos:      "),font="RcmdrTitleFont"),
         intraField,sticky="w")
  medrepFrame <- tkframe(dataTab)
  medrepVariable <- tclVar(dialog.valores$medrep.inicial)
  medrepField <- ttkentry(medrepFrame, width="12", textvariable=medrepVariable)
  tkgrid(labelRcmdr(medrepFrame,text=gettextRcmdr("Nombre de la variable cuantitativa:        "),font="RcmdrTitleFont"),
         medrepField,sticky="w")  
  tkgrid(intraFrame, sticky="nw")
  tkgrid(medrepFrame, sticky="nw")
  tkgrid(SCuadradosFrame, sticky="nw")
  tkgrid(statsFrame, sticky="nw")  
  tkgrid(opsFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','statisticsTab','optionsTab'),
               tab.names=c("Datos","Estadisticos","Opciones"))
}

ayuda.RcmdrPlugin.EACSPIR <- function(){
    doItAndPrint("help(\"RcmdrPlugin.EACSPIR\")")
    invisible(NULL)
}
