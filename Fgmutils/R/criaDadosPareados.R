##' @title Create Date Paired
##' @description paired a dataframe
##' @param dataFrame dataframe that you want to pair dataFrame must contain columns cod_id, ANO_MEDICAO1, ANO_MEDICAO2, DAP1, DAP2, HT1, HT2, ID_PROJETO
##' @param campoChave character the column that will be paired
##' @param campoComparacao character the field used to compare the period of change
##' @param camposPareados vector the fields that will be paired exemple CampoesPareados=c(dap,ht)
##' @param camposNaoPareados the fields he wants to be present without the paired
##' @return will be returned a dataframe containing columns
##' cod_id, ANO_MEDICAO1, ANO_MEDICAO2, DAP1, DAP2, HT1, HT2, ID_PROJETO
##' @import data.table
##' @import sqldf
##' @importFrom "utils" "setTxtProgressBar" "txtProgressBar"
##' @export
criaDadosPareados <- function (dataFrame, campoChave, campoComparacao, camposPareados, camposNaoPareados)
{
  ini = Sys.time()

  if (!campoChave %in% names(dataFrame)) {
    stop("O campo ",campoChave," nao existe no dataframe informado")
  }

  if (!campoComparacao %in%  names(dataFrame)) {
    stop("O campo ",campoComparacao," nao existe no dataframe informado")
  }

  naoExistemPareados = camposPareados[which((camposPareados %in% names(dataFrame)) %in% FALSE)]
  if (length(naoExistemPareados>0)) {
    stop("Os seguintes campos nao existem no dataframe informado: ", naoExistemPareados)
  }

  naoExistemNaoPareados = camposNaoPareados[which((camposNaoPareados %in% names(dataFrame)) %in% FALSE)]
  if (length(naoExistemNaoPareados>0)) {
    stop("Os seguintes campos nao existem no dataframe informado:\n", naoExistemNaoPareados)
  }

  ini = Sys.time()
  cat("\f")
  cat(paste("\nMontando dados..."))

  sql = "SELECT * from pdataFrame order by campoChave, campoComparacao"
  sql = gsub("pdataFrame", toString(substitute(dataFrame)), sql)
  sql = gsub("campoChave", campoChave, sql)
  sql = gsub("campoComparacao", campoComparacao, sql)
  dataFrame = sqldf(sql)
  dataFrame = data.table(dataFrame)

  ## Translate characters in character vectors, from upper to lower case
  setnames(dataFrame,names(dataFrame),tolower(names(dataFrame)))
  camposPareados = tolower(camposPareados)
  campoComparacao = tolower(campoComparacao)
  campoChave = tolower(campoChave)
  camposNaoPareados = tolower(camposNaoPareados)

  cont = 1;

  ##################Generation of string to assemble the dfr.
  ret = NULL
  dfrString =  "dfr <- data.frame("
  eval(parse(text=paste0("ret=verificaTipoColuna(dataFrame$", campoChave,")")))
  aux = paste0(campoChave, "=", ret)
  dfrString = paste0(dfrString,aux)
  eval(parse(text=paste0("ret=verificaTipoColuna(dataFrame$", campoComparacao,")")))
  aux1 = paste0(campoComparacao, "1=", ret)
  aux2 = paste0(campoComparacao, "2=", ret)
  aux = paste0(',', aux1,", " ,aux2)
  dfrString = paste0(dfrString,aux)


  #Creating the matched data

  for (i in 1:length(camposPareados))
  {
    eval(parse(text=paste0("ret=verificaTipoColuna(dataFrame$", camposPareados[i],")")))
    aux1 = paste0(camposPareados[i], "1=", ret)
    aux2 = paste0(camposPareados[i], "2=", ret)
    aux = paste0(',', aux1,", " ,aux2)
    dfrString = paste0(dfrString, aux)
  }
  #Creating unmatched data
  for (i in 1:length(camposNaoPareados))
  {
    eval(parse(text=paste0("ret=verificaTipoColuna(dataFrame$", camposNaoPareados[i],")")))
    aux = paste0(camposNaoPareados[i], "=", ret)
    dfrString = paste0(dfrString,", " ,aux)
  }

  #Finishing the dfr
  dfrString = paste0(dfrString, ",stringsAsFactors=FALSE) ")
  eval(parse(text=paste0(dfrString)))
  dfr = data.table(dfr)
  n = nrow(dataFrame)
  pb = txtProgressBar(min=1, max=n, style=3)

  for (i in 2:n) {
    setTxtProgressBar(pb, i)
    #catching the field campoComparacao

    #campoComparacao
    eval(parse(text=paste0(campoComparacao, "1 = dataFrame$", campoComparacao,"[i-1]")))
    eval(parse(text=paste0(campoComparacao, "2 = dataFrame$", campoComparacao,"[i]")))
    #campoChave
    eval(parse(text=paste0(campoChave, "1 = dataFrame$", campoChave,"[i-1]")))
    eval(parse(text=paste0(campoChave, "2 = dataFrame$", campoChave,"[i]")))

    stringIf = paste0("if((", campoComparacao,"2 > ", campoComparacao,"1) && (", campoChave, "1 == ", campoChave,"2 )){")

    #if is true
    #catching the value of the matched fields
    for (j in 1:length(camposPareados)) {
      aux1 = paste0(camposPareados[j], "1=retornaValor(dataFrame$",camposPareados[j],"[i-1])");
      aux2 = paste0(camposPareados[j], "2=retornaValor(dataFrame$",camposPareados[j],"[i])");
      aux = paste(aux1,"; ",   aux2, ";")
      stringIf = paste0(stringIf, aux);
    }
    #catching the values of unpaired fields
    for (j in 1:length(camposNaoPareados)) {
      aux = paste0(camposNaoPareados[j], "=retornaValor(dataFrame$",camposNaoPareados[j],"[i])");
      stringIf = paste0(stringIf, aux,";");
    }




    #create a linha
    linha = NULL
    stringLinha = "linha = list("
    aux = paste0(campoChave, "=", campoChave,"1,")
    stringLinha = paste0(stringLinha, aux);
    aux1 = paste0(campoComparacao, "1 = ",campoComparacao, "1");
    aux2 = paste0(campoComparacao, "2 = ",campoComparacao, "2");
    aux = paste(aux1,", ",   aux2, ", ")
    stringLinha = paste0(stringLinha, aux);


    #paired fields
    if(length(camposPareados)==1){
      aux1 = paste0(camposPareados[1], "1=",camposPareados[1], "1");
      aux2 = paste0(camposPareados[1], "2=",camposPareados[1], "2");
      aux = paste(aux1,", ",   aux2, ", ")
      stringLinha = paste0(stringLinha, aux);

    }
    else{
      for (j in 1:(length(camposPareados)-1)) {
        aux1 = paste0(camposPareados[j], "1=",camposPareados[j], "1");
        aux2 = paste0(camposPareados[j], "2=",camposPareados[j], "2");
        aux = paste(aux1,", ",   aux2, ", ")
        stringLinha = paste0(stringLinha, aux);
      }
      aux1 = paste0(camposPareados[j+1], "1=",camposPareados[j+1], "1");
      aux2 = paste0(camposPareados[j+1], "2=",camposPareados[j+1], "2");
      aux = paste(aux1,", ",   aux2, ", ")
      stringLinha = paste0(stringLinha, aux);
    }
    #Fields unpaired
    if(length(camposNaoPareados)==1){
      aux = paste0(camposNaoPareados[1], "=",camposNaoPareados[1]);
      stringLinha = paste0(stringLinha, aux, ');');

    }
    else
    {
      for (j in 1:(length(camposNaoPareados)-1))
      {
        aux = paste0(camposNaoPareados[j], "=",camposNaoPareados[j]);
        stringLinha = paste0(stringLinha, aux, ',');
      }
      aux = paste0(camposNaoPareados[j+1], "=",camposNaoPareados[j+1]);
      stringLinha = paste0(stringLinha, aux, ');');

    }

    stringIf = paste0(stringIf, stringLinha)
    stringIf = paste0(stringIf, "; dfr = rbindlist(list(dfr, linha)) ; cont = cont + 1} else{next}")
    eval(parse(text=paste0(stringIf)))


  }
  fim = Sys.time()
  tempo = fim-ini
  cat("\n\n")
  print(tempo)
  remove(ini, fim, cont, linha, pb, tempo)

  return (dfr)
}
