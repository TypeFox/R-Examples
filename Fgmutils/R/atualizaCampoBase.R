##' @title updated base field
##' @description this function update certain fields in a dataframe, based on the provided key
##' @param camposAtualizar is the vector you want to update
##' @param baseAgrupada It is the database that contains the data you want to update on dataframe
##' @param baseAtualizar It is dataframe that you want to change fields
##' @param keys are the keys of the table that will be used in the compare
##' @return baseAtualizar with the updated fields according to baseAgrupada
##' @import data.table
##' @export
atualizaCampoBase <- function (camposAtualizar, baseAgrupada, baseAtualizar, keys){

  ini = Sys.time()

  #Picking up the columns of the database to be updated.
  baseAtualizar = data.table(baseAtualizar)  # Base to be updated.
  baseAgrupada = data.table(baseAgrupada) # Base with grouped data
  
  for (i in 1:length(camposAtualizar)) {
   comando = paste0("baseAtualizar$",camposAtualizar[i],"=-999")
   eval(parse(text=comando))
   remove(comando)
  }	  

  tipo = NULL
  tipo1 = NULL
  tipo2 = NULL
  #check that the column already exists in the base, failing that, will create
  for(i in 1:(length(camposAtualizar))){

    if(!(camposAtualizar[i] %in% names(baseAtualizar))){


      (novaColuna <- as.numeric(rep(-1, nrow(baseAtualizar))))

      baseAtualizar <- cbind(baseAtualizar, novaColuna)
      id_col = which(names(baseAtualizar) == "novaColuna") # Take the column identifier in the database.
      setnames(baseAtualizar, "novaColuna", camposAtualizar[i])
    }
    eval(parse(text=paste0("tipo = typeof(baseAgrupada$",camposAtualizar[i],")")))
    #to hit the type



    if(tipo == "character"){
      eval(parse(text=paste0("baseAtualizar$",camposAtualizar[i]," = as.character('-1')")))
    }
    else
      eval(parse(text=paste0("baseAtualizar$",camposAtualizar[i]," = as.numeric(-1)")))
  }

  #Putting the keys in the same type
  for(i in 1:(length(keys))){
    eval(parse(text=paste0("tipo1 = typeof(baseAgrupada$",keys[i],")")))
    eval(parse(text=paste0("tipo2 = typeof(baseAtualizar$",keys[i],")")))
    if (tipo1 != tipo2){
      eval(parse(text=paste0("baseAtualizar$",keys[i]," = as.character(baseAtualizar$",keys[i],")")))
      eval(parse(text=paste0("baseAgrupada$",keys[i]," = as.character(baseAgrupada$",keys[i],")")))

    }

  }


  #Generating the keys to setkey
  chavesSet = ""


  if((length(keys))>1){



    for(i in 1:(length(keys)-1)){
      chavesSet = paste0(chavesSet, keys[i], ", ")
    }
    chavesSet = paste0(chavesSet, keys[i+1])

    eval(parse(text=paste0("setkey(baseAtualizar, ", chavesSet,")")))
    eval(parse(text=paste0("setkey(baseAgrupada, ", chavesSet,")")))





    cat("\n")
    for (i in 1:nrow(baseAgrupada)) {
      chaves=""
      for(j in 1:(length(keys)-1)){

        eval(parse(text=paste0(keys[j],"=retornaValor(baseAgrupada$", keys[j], "[i])")))
        eval(parse(text=paste0("chaves=paste0(chaves, ", keys[j],", ', ')" )))
      }
      eval(parse(text=paste0(keys[j+1],"=retornaValor(baseAgrupada$", keys[j+1], "[i])")))
      eval(parse(text=paste0("chaves=paste0(chaves, ", keys[j+1],")" )))
      eval(parse(text=paste0("setkey(baseAtualizar, ", chavesSet,")")))
      cat(".")
      varNovoCampo = NULL
      for(j in 1:(length(camposAtualizar))){
        eval(parse(text=paste0("varNovoCampo = retornaValor(baseAgrupada$", camposAtualizar[j], "[i])")))


        eval(parse(text=paste0("baseAtualizar[list(", chaves,"), ",camposAtualizar[j], " := ", varNovoCampo,"]")))
      }



    }
  }
  else {

    chavesSet = paste0(chavesSet, keys[1])

    eval(parse(text=paste0("setkey(baseAtualizar, ", chavesSet,")")))
    eval(parse(text=paste0("setkey(baseAgrupada, ", chavesSet,")")))


    for (i in 1:nrow(baseAgrupada)) {
      chaves=""
      aux = ""
      eval(parse(text=paste0(keys[1],"=retornaValor(baseAgrupada$", keys[1], "[i])")))
      eval(parse(text=paste0("chaves=paste0(chaves, ", keys[1],")" )))
      eval(parse(text=paste0("setkey(baseAtualizar, ", chavesSet,")")))
      cat(".")

      for(j in 1:(length(camposAtualizar))){

        eval(parse(text=paste0("varNovoCampo = retornaValor(baseAgrupada$", camposAtualizar[j], "[i])")))
        eval(parse(text=paste0("baseAtualizar[list(", chaves,"), ",camposAtualizar[j], " := ", varNovoCampo,"]")))
      }



    }




  }
  fim = Sys.time()
  tempo = fim-ini
  remove(fim, ini)
  cat("\n")
  print(tempo)


  return (baseAtualizar)

}
