setMethod("initialize", "yuima.characteristic",
           function(.Object, equation.number, time.scale){
             if(equation.number==length(time.scale)){
               .Object@equation.number <- equation.number
               .Object@time.scale <- time.scale
             }else if(length(time.scale)==1){
               time.scale <- rep(time.scale, equation.number)
               .Object@equation.number <- equation.number
               .Object@time.scale <- time.scale
             }else{
               yuima.warn("Dimension missmatch")
               return(NULL)
             }
             return(.Object)
           })


setCharacteristic <-
  function(equation.number=1, time.scale=1){
    return(new("yuima.characteristic", equation.number, time.scale=time.scale))
  }
