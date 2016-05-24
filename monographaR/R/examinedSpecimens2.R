examinedSpecimens2 <-
function (data) 
{
  colnames(data) <- c("spp", "col", "cn", "h", "hn", "cntr", "state", "city")
  data$specimens <- paste(data$col, data$cn)
  miss.col <- which(is.na(data$col))
  miss.cn <- which(is.na(data$cn))
  missing <- c(miss.col, miss.cn)
  missing <- unique(missing)
  data$specimens[missing] <- paste(data$h[missing], data$hn[missing])
  filename = ""
    #cat("Examine", file = filename, fill = T, append = T)
    sp.data <- data
    cntrs <- levels(as.factor(as.character(sp.data$cntr)))
    for (c in 1:length(cntrs)) {
      cat(cntrs[c], ". ", sep = "", file = filename, append = T)
      cntr.data <- sp.data[sp.data$cntr == cntrs[c], ]
      states <- levels(as.factor(as.character(cntr.data$state)))
      for (st in 1:length(states)) {
        cat(states[st], ": ", sep = "", file = filename, 
            append = T)
        state.data <- cntr.data[cntr.data$state == states[st], 
                                ]
        city <- levels(as.factor(as.character(state.data$city)))
        if (length(city) == 1) {
          cat(city, ", ", sep = "", file = filename, 
              append = T)
          city.data <- state.data[state.data$city == 
                                    city, ]
          spec.levs <- levels(as.factor(as.character(city.data$specimens)))
          sep.cities <- rep(",", length(spec.levs))
          sep.cities[length(spec.levs)] <- "."
          for (co in 1:length(spec.levs)) {
            cat(spec.levs[co], " (", sep = "", file = filename, 
                append = T)
            spec.data <- city.data[city.data$specimens == 
                                     spec.levs[co], ]
            herbs <- levels(as.factor(as.character(spec.data$h)))
            if (length(herbs) == 1) {
              cat(herbs, ")", sep = "", file = filename, 
                  append = T)
            }
            else {
              for (h in 1:(length(herbs) - 1)) {
                cat(herbs[h], ", ", sep = "", file = filename, 
                    append = T)
              }
              cat(herbs[length(herbs)], ")", sep = "", 
                  file = filename, append = T)
            }
            cat(sep.cities[co], " ", sep = "", file = filename, 
                append = T)
          }
        }
        else {
          for (ci in 1:(length(city) - 1)) {
            cat(city[ci], ", ", sep = "", file = filename, 
                append = T)
            city.data <- state.data[state.data$city == 
                                      city[ci], ]
            spec.levs <- levels(as.factor(as.character(city.data$specimens)))
            sep.cities <- rep(",", length(spec.levs))
            sep.cities[length(spec.levs)] <- ";"
            for (co in 1:length(spec.levs)) {
              cat(spec.levs[co], " (", sep = "", file = filename, 
                  append = T)
              spec.data <- city.data[city.data$specimens == 
                                       spec.levs[co], ]
              herbs <- levels(as.factor(as.character(spec.data$h)))
              if (length(herbs) == 1) {
                cat(herbs, ")", sep = "", file = filename, 
                    append = T)
              }
              else {
                for (h in 1:(length(herbs) - 1)) {
                  cat(herbs[h], ", ", sep = "", file = filename, 
                      append = T)
                }
                cat(herbs[length(herbs)], ")", sep = "", 
                    file = filename, append = T)
              }
              cat(sep.cities[co], " ", sep = "", file = filename, 
                  append = T)
            }
          }
          cat(city[length(city)], ", ", sep = "", file = filename, 
              append = T)
          city.data <- state.data[state.data$city == 
                                    city[length(city)], ]
          spec.levs <- levels(as.factor(as.character(city.data$specimens)))
          sep.cities <- rep(",", length(spec.levs))
          sep.cities[length(spec.levs)] <- "."
          for (co in 1:length(spec.levs)) {
            cat(spec.levs[co], " (", sep = "", file = filename, 
                append = T)
            spec.data <- city.data[city.data$specimens == 
                                     spec.levs[co], ]
            herbs <- levels(as.factor(as.character(spec.data$h)))
            if (length(herbs) == 1) {
              cat(herbs, ")", sep = "", file = filename, 
                  append = T)
            }
            else {
              for (h in 1:(length(herbs) - 1)) {
                cat(herbs[h], ", ", sep = "", file = filename, 
                    append = T)
              }
              cat(herbs[length(herbs)], ")", sep = "", 
                  file = filename, append = T)
            }
            cat(sep.cities[co], " ", sep = "", file = filename, 
                append = T)
          }
        }
      }
    }
    cat("\n", file = filename, append=T)
}
