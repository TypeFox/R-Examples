## Sets LE table to 'Standard Life Expectancy'
## Update 'LE' and '.LE'

setStdLE <-
function(){
  for (i in seq(21)){
    DALYassign("LE", DALYget("stdM")[i], i, 1)
    DALYassign("LE", DALYget("stdF")[i], i, 2)
  }
  DALYupdate(".LE")
}