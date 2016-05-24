###########################################################################################
## Función que divide un objeto 'data.frame' en una lista dividida en variables          ##
##   cuantitativas y cualitativas 	                                                 ##
##											 ##
##											 ##
## Elaborado por: Pedro Cesar del Campo Neira	
## Revisado por: Campo Elías Pardo					 ##
## Universidad Nacional de Colombia							 ##
## 											 ##
## requiere:ade4      library(ade4)							 ##
##											 ##
## Fac.Num  ( tabla     := objeto 'data frame' )	        			 ##
## 											 ##
###########################################################################################

Fac.Num <- function (tabla) 
{
    if (!is.data.frame(tabla)){return(cat("Only 'data.frame' \n"))}
    
    c.var <- dim(tabla)[2]
    clase <- NULL
    
    for (i in 1:c.var) {
      if ( class(tabla[, i])=="character" ) tabla[,i] <- as.factor(tabla[,i])
      clase <- cbind(clase, class(tabla[, i]))
    }
        
    SALIDA <- NULL
    SALIDA <- as.list(SALIDA)
    SALIDA$numeric <- tabla[(1:c.var)[clase == "numeric"]]
    SALIDA$factor  <- tabla[(1:c.var)[clase == "factor" ]]
    SALIDA$integer <- tabla[(1:c.var)[clase == "integer"]]
    
    if (dim(SALIDA$numeric)[2] == 0) {SALIDA$numeric <- NULL }
    if (dim(SALIDA$factor)[2]  == 0) {SALIDA$factor  <- NULL }
    if (dim(SALIDA$integer)[2] == 0) {SALIDA$integer <- NULL }
    
    return(SALIDA)
}

 