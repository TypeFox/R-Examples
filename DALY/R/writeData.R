## Helper function for 'saveDALYdata'
## Write 'DALY' database data for one outcome

writeData <-
function(txt, dist, strat, data, file){
  write(c(txt, tclvalue(dist), tclvalue(strat)),
        file = file, ncolumns = 6, sep = "\t", append = TRUE)
  for (i in seq(5))
    write(getData(data, "data")[i,],
          file = file, ncolumns = 6, sep = "\t", append = TRUE)
}