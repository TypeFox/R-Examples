AT.FLUKA.read.USRBIN.regs <- function (exp.name, number.of.runs, unit, data.source = "local", 
    vol.cm3 = NULL, density.g.cm3 = NULL)                        
{                                                                
    for (cur.run in 1:number.of.runs) {                          
        file.name <- paste(exp.name, AT.add.leading.zeros(cur.run, 
            3), "_fort.", unit, sep = "")                          
        if (data.source == "condor") {                             
            file.name <- paste(exp.name, "_node", AT.add.leading.zeros(cur.run - 
                1, 4), "001_fort.", unit, sep = "")                              
        }                                                                        
        if (data.source == "condor_cleaned") {                                   
            file.name <- paste(exp.name, AT.add.leading.zeros(cur.run -          
                1, 5), "_fort.", unit, sep = "")                                 
        }                                                                        
        input <- scan(file = file.name, what = "character", strip.white = TRUE,  
            sep = "")                                                            
        no.of.particles <- as.numeric(gsub(",", "", input[grep("followed",       
            input) + 1]))                                                        
        outputs <- grep("Region", input)                                         
        if (cur.run == 1) {                                                      
            names <- gsub(" ", "", input[outputs + 4])                           
            bins <- as.numeric(input[outputs + 10])                              
            index <- c(1, cumsum(bins)[-length(bins)] + 1)                       
            df <- data.frame(idx = 1:sum(bins), name = character(sum(bins)),     
                index = numeric(sum(bins)), bin = numeric(sum(bins)),            
                E.dep.GeV = numeric(sum(bins)), E2.dep.GeV2 = numeric(sum(bins)))
            class(df$name) <- "character"                                        
            for (i in 1:length(names)) {                                         
                ii <- df$idx >= index[i] & df$idx < (index[i] +                  
                  bins[i])                                                       
                df$name[ii] <- as.character(rep(names[i], sum(ii)))              
                df$bin[ii] <- 1:sum(ii)                                          
                tmp <- as.numeric(input[outputs[i] + 66 + (1:bins[i]) -          
                  1])                                                            
                df$E.dep.GeV <- tmp                                              
                df$E2.dep.GeV2 <- tmp^2                                          
            }                                                                    
        }                                                                        
        else {                                                                   
            for (i in 1:length(names)) {                                         
                ii <- df$idx >= index[i] & df$idx < (index[i] +                  
                  bins[i])                                                       
                tmp <- as.numeric(input[outputs[i] + 66 + (1:bins[i]) -          
                  1])                                                            
                df$E.dep.GeV <- df$E.dep.GeV + tmp
                df$E2.dep.GeV2 <- df$E2.dep.GeV2 + tmp^2
            }
        }
    }

    if (!(is.null(vol.cm3) | is.null(density.g.cm3))) {
        if (length(vol.cm3) == 1) {
            df$vol.cm3 <- rep(vol.cm3, nrow(df))
        }
        else {
            df$vol.cm3 <- vol.cm3
        }
        if (length(density.g.cm3) == 1) {
            df$density.g.cm3 <- rep(density.g.cm3, nrow(df))
        }
        else {
            df$density.g.cm3 <- density.g.cm3
        }
        E.dep.GeV.cm3 <- df$E.dep.GeV/(number.of.runs * df$vol.cm3)
        stdev.E.dep.GeV.cm3 <- sqrt(df$E2.dep.GeV2/(number.of.runs *
            df$vol.cm3) - E.dep.GeV.cm3^2)
        sterr.E.dep.GeV.cm3 <- stdev.E.dep.GeV.cm3/sqrt(number.of.runs)
        df$D.Gy <- E.dep.GeV.cm3 * 0.0000001602176462/df$density.g.cm3
        df$sterr.D.Gy <- sterr.E.dep.GeV.cm3 * 0.0000001602176462/df$density.g.cm3
        df$E.dep.GeV <- NULL
        df$E2.dep.GeV2 <- NULL
        df$vol.cm3 <- NULL
        df$density.g.cm3 <- NULL
    }else{
        df$E.dep.GeV       <- df$E.dep.GeV/(number.of.runs)
        df$stdev.E.dep.GeV <- sqrt(df$E2.dep.GeV2/number.of.runs - df$E.dep.GeV^2)
        df$sterr.E.dep.GeV <- df$stdev.E.dep.GeV/sqrt(number.of.runs)
	df$E2.dep.GeV2 <- NULL
    }
    
    df$idx <- NULL
    df$index <- NULL
    return(df)
}