
##################################################
# converts a dataframe in longformat (produced
# in jomo) into a list of datasets

jomo2datlist <- function( jomo.dataframe , variable="Imputation" ){
    dat <- jomo.dataframe
    M1 <- max( unique( dat$Imputation ) )
    datlist <- as.list( 1:M1 )
    ind <- which( colnames(dat) == variable )
    for ( mm in 1:M1){
        datlist[[mm]] <- dat[ dat$Imputation == mm , - ind ]
                    }
    return(datlist)
            }

