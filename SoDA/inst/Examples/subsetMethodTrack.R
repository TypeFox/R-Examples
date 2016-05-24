setMethod("[",
      signature(x = "trackNumeric", i = "ANY", j = "missing"),
      function(x, i) {
            x@.Data[i]
      })
