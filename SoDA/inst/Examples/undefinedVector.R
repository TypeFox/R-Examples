setMethod("undefined<-",
          signature(x="vector", value = "vector"),
          function(x, codes = NULL, value) {
              if(length(codes) > 0)
                  x[x %in% codes] <- NA
              x[is.na(x)] <- value
              x
          })
