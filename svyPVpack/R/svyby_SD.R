svyby_SD <-
function(formula, design, deff=FALSE, na.rm=TRUE)
{
  # formel in einen variablennamen umwandeln --> hier werden also die PV rausgefischt
  NAM1 <- gsub("\\s*~\\s*","",formula)[2]
  # übergibt hier die richtige funktion
  withReplicates(design, weigh.SD, vari=NAM1, na.rm=na.rm)  
}
