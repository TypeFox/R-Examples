#' @title avalia Ajuste
#' @description this function evaluates the quality of the adjustment of the statistical model, rom observed data and those estimated by the model, observed
#' @param dataFrame dataFrane with information observed, estimated
#' @param variavelObservados vector of values observed.
#' @param variavelEstimados vector of values estimated.
#' @param linear boolean is linear model
#' @param nParametros number of parameters used in the adjusted model
#' @param intercepto if you model is no-intercepto use FALSE
#' @param plot Vector graphic information
#' @param modelo the name of the adjusted model
#' @param resumo if you want summary information, use TRUE
#' @param emf to save the graphic in the format emf use TRUE
#' @importFrom "stats" "scatter.smooth"
#' @importFrom "grDevices" "dev.off" "png"
#' @importFrom "graphics" "abline" "hist" "mtext" "par"
#' @importFrom "stats" "cor"
#' @importFrom "devEMF" "emf"
#' @importFrom "png" "readPNG"
#' @export
avaliaAjuste <- function(dataFrame, variavelObservados, variavelEstimados, linear=TRUE, nParametros=NA, intercepto = TRUE, plot=NA, modelo=NA, resumo=FALSE, emf=TRUE) {
  if (is.na(nParametros) || nParametros < 0) {
    stop("Enter the number of parameters used in the adjusted model.")
  }
  if (is.na(modelo) || !is.character(modelo)) {
    stop("Enter the name of the adjusted model.")
  }
  if (!is.data.frame(dataFrame)) {
    stop("Enter a valid data.frame!")
  }
  if (!is.character(variavelObservados) || !is.character(variavelEstimados)) {
    stop("t must enter the name of the observed variables (variavelObservados) and estimated (variavelEstimados) ")
  }
  if (variavelObservados == variavelEstimados) {
    stop("he variables observed and estimated can not the same.")
  }
  cat(iconv("\nCalculando estatisticas", from="UTF-8", to="LATIN1"))

  if (!is.na(plot[1])) {
    dir.create(plot[1], showWarnings = F, recursive=T)
  }


  ####################################################
  strObservados = paste0("dataFrame$", variavelObservados)
  strEstimados = paste0("dataFrame$", variavelEstimados)
  strResiduos = paste0("dataFrame$residuo_", variavelObservados)
  strResiduosPerc = paste0("dataFrame$residuoPerc_", variavelObservados)
  strCalculaResiduos = paste0(strResiduos, " = ", strObservados, " - ", strEstimados)
  eval(parse(text=strCalculaResiduos))

  strCalculaResiduosPerc = paste0(strResiduosPerc," = residuoPerc(observados = ", strObservados, ", estimados = ",strEstimados,")")

  ####################################################
  eval(parse(text=strCalculaResiduosPerc))

  resultado = list()

  if (linear) {
    if (intercepto) {
      resultado$r2 = R21a(observados = eval(parse(text=strObservados)), estimados = eval(parse(text=strEstimados)), k=nParametros)
    } else {
      resultado$r2 = R29a(observados = eval(parse(text=strObservados)), estimados = eval(parse(text=strEstimados)), k=nParametros)
    }
  }

  resultado$correlacao = cor(eval(parse(text=strEstimados)), eval(parse(text=strObservados)))
  resultado$cv = syx(eval(parse(text=strObservados)), eval(parse(text=strEstimados)), n=nrow(dataFrame), p=nParametros)

  resultado$cvPerc =  syxPerc(resultado$cv, eval(parse(text=strObservados)))

  resultado$rmse =  rmse( eval(parse(text=strObservados)), eval(parse(text=strEstimados)))
  resultado$rmsePerc = calculaPerc(resultado$rmse,eval(parse(text=strObservados)))

  resultado$bias = bias(observados = eval(parse(text=strObservados)), estimados = eval(parse(text=strEstimados)))
  resultado$biasPerc = calculaPerc(resultado$bias,eval(parse(text=strObservados)))

  resultado$mae  =  mae( eval(parse(text=strObservados)), eval(parse(text=strEstimados)))
  resultado$rrmse = rrmse( eval(parse(text=strObservados)), eval(parse(text=strEstimados)))
  resultado$ce = ce(eval(parse(text=strObservados)), eval(parse(text=strEstimados)))

  if (resumo) {
    resultado$cv = NULL
    resultado$cvPerc = NULL
    resultado$mae = NULL
    resultado$rrmse = NULL
    resultado$ce = NULL
  }


  if (any(!is.na(plot))) {
    cat(iconv("\nGerando grafico...", from="UTF-8", to="LATIN1"))

    strVariavelXResiduo = paste0("dataFrame$", plot[5])
    strVariavel_Residuo = paste0("dataFrame$", plot[6])

    strVariavel_ResiduoPerc = paste0("dataFrame$", plot[7])

    Sys.setlocale('Portuguese_Brazil.1252', category="LC_ALL")

    file = NULL
    str_file=paste0(plot[1], modelo)
    if (emf) {
      file=paste0(str_file, ".emf")
      emf(file, bg = "transparent", family = "Times", pointsize = 10)
    } else {
      file = paste0(str_file, ".png")
      png(file, width = 800, height = 800)
    }

    par(mfrow=c(2,2))

    #scatter.smooth
    scatter.smooth(eval(parse(text=strObservados)),eval(parse(text=strEstimados)), pch=18, col=2,
                   ylab=eval(parse(text=plot[2])),
                   xlab=eval(parse(text=plot[3])))
    abline(0,1)
    mtext(modelo, side = 1, line = -66, outer =TRUE)

    scatter.smooth(eval(parse(text=strVariavelXResiduo)), eval(parse(text=strVariavel_Residuo)), col=2, pch=18,
                   ylab=iconv("Residuo absoluto", from="UTF-8", to="LATIN1"),
                   xlab=eval(parse(text=plot[4])))
    abline(0,0)

    scatter.smooth(eval(parse(text=strVariavelXResiduo)), eval(parse(text=strVariavel_ResiduoPerc)), col=2, pch=18,
                   ylab=iconv("Residuo (%)", from="UTF-8", to="LATIN1"),
                   xlab=eval(parse(text=plot[4])))
    abline(0,0)

    hist(eval(parse(text=strVariavel_Residuo)), xlab=iconv("Residuos", from="UTF-8", to="LATIN1"), breaks=100, main=iconv("Histograma de residuos", from="UTF-8", to="LATIN1"))
    dev.off()
    par(mfrow=c(1,1))


    dev.off()
    if (!emf) {
      img <- readPNG(paste0(str_file, ".png"))
      grid::grid.raster(img)
    }
    cat(paste0(iconv("\nGrafico gerado: '", from="UTF-8", to="LATIN1"), str_file, ".png'"))
  }

  cat("\nDone!")
  return(list(resultado=data.frame(modelo=modelo, resultado), dataFrame = dataFrame))
}
