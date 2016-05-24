AT.FLUKA.read.USRBIN.mesh <- function(exp.name, number.of.runs, unit, data.source = 'local', density.g.cm3 = 1.0)
{

    for (cur.run in 1:number.of.runs){
        # DEBUG: cur.run <- 1
        file.name            <-    paste(    exp.name,                
                                    AT.add.leading.zeros(cur.run, 3),
                                    "_fort.",
                                    unit,
                                    sep = "")
        if(data.source == "condor"){                                                         # condor naming in case cleaning hasn't been run
            file.name            <-    paste(    exp.name,
                                        "_node",
                                        AT.add.leading.zeros(cur.run - 1, 4),
                                        "001_fort.",
                                        unit,
                                        sep = "")
        }
        if(data.source == "condor_cleaned"){
            file.name            <-    paste(    exp.name,                                    # condor
                                        AT.add.leading.zeros(cur.run - 1, 5),
                                        "_fort.",
                                        unit,
                                        sep = "")
        }
        input               <-    scan(file = file.name, what = "character", strip.white = TRUE, sep = "")
        
        # find no. of particles
        no.of.particles    <-    as.numeric(gsub(",", "", input[grep("followed", input) + 1]))

        # find usrbin outputs
        outputs            <-    grep("Cartesian", input) 

        if (cur.run == 1){
            names               <- gsub(" ", "", input[outputs + 4])
            bins                <- as.numeric(input[outputs + 43])
            index               <- c(1, cumsum(bins)[-length(bins)] + 1)

            # build data.frame
            df    <-    data.frame( idx           = 1:sum(bins),
                                    name          = character(sum(bins)),
                                    bin           = numeric(sum(bins)),
                                    E.GeV.cm3     = numeric(sum(bins)),
                                    E2.GeV2.cm6   = numeric(sum(bins)))

            class(df$name)        <-    "character"

            for (i in 1:length(names)){
                #i <- 1
                ii                   <- df$idx >= index[i] & df$idx < (index[i] + bins[i])
                df$name[ii]          <- as.character(rep(names[i], sum(ii)))
                df$bin[ii]           <- 1:sum(ii)
                tmp                  <- as.numeric(input[outputs[i] + 63 + (1:bins[i]) - 1])
                df$E.GeV.cm3         <- tmp
                df$E2.GeV2.cm6       <- tmp^2
            }
        }else{
            for (i in 1:length(names)){
                #i <- 1
                ii                   <- df$idx >= index[i] & df$idx < (index[i] + bins[i])
                tmp                  <- as.numeric(input[outputs[i] + 63 + (1:bins[i]) - 1])
                df$E.GeV.cm3         <- df$E.GeV.cm3 + tmp
                df$E2.GeV2.cm6       <- df$E2.GeV2.cm6 + tmp^2
            }
        }
    }

    df$E.GeV.cm3         <- df$E.GeV.cm3 / number.of.runs 
    stdev.E.GeV.cm3      <- sqrt(df$E2.GeV2.cm6 / number.of.runs  - df$E.GeV.cm3^2)				# valid only for large (>10) numbers of number.of.runs, when 1/N ~ 1/(N-1) for estimating the stdev
    sterr.E.GeV.cm3	     <- stdev.E.GeV.cm3 / sqrt(number.of.runs)
    df$D.Gy              <- df$E.GeV.cm3 * 1.602176462e-7 / density.g.cm3
    df$sterr.D.Gy        <- sterr.E.GeV.cm3 * 1.602176462e-7 / density.g.cm3

    df$E.GeV.cm3         <- NULL
    df$E2.GeV2.cm6       <- NULL
    df$idx               <- NULL

    return(df)
}
