examinedSpecimens <-
function(data, filename="examined.txt") {
  if (class(data) != "data.frame") {
    stop("data must be a data.frame")
  }
  if (ncol(data) != 8) {
    stop("data must have 8 columns, see help(\"examinedSpecimens\")")
  }
  colnames(data) <- c("spp", "col", "cn", "h", "hn", "cntr", "state", "city")
  message("Assuming the columns are ordered as: species, collector name, collector number, herbarium acronym, herbarium number, country, state, municipality")
  paste(data$col, data$cn) -> data$specimens
  which(is.na(data$col)) -> miss.col
  which(is.na(data$cn)) -> miss.cn
  c(miss.col, miss.cn) -> missing
  unique(missing) -> missing
  paste(data$h[missing], data$hn[missing]) -> data$specimens[missing]
  levels(as.factor(data$spp)) -> spp 
  cat("Specimens Examined", file=filename)
  for (s in 1:length(spp)) {
    cat("\n", file=filename, fill=T, append=T)
    cat(spp[s], file=filename, fill=T, append=T)
    cat("", file=filename, fill=T, append=T)
    data[data$spp == spp[s],] -> sp.data
    levels(as.factor(as.character(sp.data$cntr))) -> cntrs
    for (c in 1:length(cntrs)) {
      cat(cntrs[c],". ", sep="", file=filename, append=T)
      sp.data[sp.data$cntr == cntrs[c],] -> cntr.data
      levels(as.factor(as.character(cntr.data$state))) -> states
      for (st in 1:length(states)) {
        cat(states[st],": ", sep="", file=filename, append=T)
        cntr.data[cntr.data$state == states[st],] -> state.data
        levels(as.factor(as.character(state.data$city))) -> city
        if (length(city) == 1) {
          cat(city,", ", sep="", file=filename, append=T)
          state.data[state.data$city == city,] -> city.data
          levels(as.factor(as.character(city.data$specimens))) -> spec.levs
          sep.cities <- rep(",", length(spec.levs))
          sep.cities[length(spec.levs)] <- "."
          for (co in 1:length(spec.levs)) {
            cat(spec.levs[co]," (", sep="", file=filename, append=T)
            city.data[city.data$specimens == spec.levs[co],] -> spec.data
            levels(as.factor(as.character(spec.data$h))) -> herbs
            if (length(herbs) == 1) {
              cat(herbs,")", sep="", file=filename, append=T)
            } else {
              for (h in 1:(length(herbs)-1)) {
                cat(herbs[h],", ", sep="", file=filename, append=T)
              }
              cat(herbs[length(herbs)],")", sep="", file=filename, append=T)
            }
            cat(sep.cities[co], " ", sep="", file=filename, append=T)
          }
        } else {
          for (ci in 1:(length(city)-1)) {
            cat(city[ci],", ", sep="", file=filename, append=T)
            state.data[state.data$city == city[ci],] -> city.data
            levels(as.factor(as.character(city.data$specimens))) -> spec.levs
            sep.cities <- rep(",", length(spec.levs))
            sep.cities[length(spec.levs)] <- ";"          
            for (co in 1:length(spec.levs)) {
              cat(spec.levs[co]," (", sep="", file=filename, append=T)
              city.data[city.data$specimens == spec.levs[co],] -> spec.data
              levels(as.factor(as.character(spec.data$h))) -> herbs
              if (length(herbs) == 1) {
                cat(herbs,")", sep="", file=filename, append=T)
              } else {
                for (h in 1:(length(herbs)-1)) {
                  cat(herbs[h],", ", sep="", file=filename, append=T)
                }
                cat(herbs[length(herbs)],")", sep="", file=filename, append=T)
              }
              cat(sep.cities[co], " ", sep="", file=filename, append=T)
            }
          }
          cat(city[length(city)],", ", sep="", file=filename, append=T)
          state.data[state.data$city == city[length(city)],] -> city.data
          levels(as.factor(as.character(city.data$specimens))) -> spec.levs
          sep.cities <- rep(",", length(spec.levs))
          sep.cities[length(spec.levs)] <- "."  
          for (co in 1:length(spec.levs)) {
            cat(spec.levs[co]," (", sep="", file=filename, append=T)
            city.data[city.data$specimens == spec.levs[co],] -> spec.data
            levels(as.factor(as.character(spec.data$h))) -> herbs
            if (length(herbs) == 1) {
              cat(herbs,")", sep="", file=filename, append=T)
            } else {
              for (h in 1:(length(herbs)-1)) {
                cat(herbs[h],", ", sep="", file=filename, append=T)
              }
              cat(herbs[length(herbs)],")", sep="", file=filename, append=T)
            }
            cat(sep.cities[co], " ", sep="", file=filename, append=T)
          }
        }
      }
    }
  }
  if (nchar(filename) > 0) {
    cat("The examined specimens list was saved in:")
    cat("\n", getwd())
  }
}
