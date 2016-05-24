AT.FLUKA.read.USRTRACK <- function(exp.name, number.of.runs, unit, data.source = 'local', compress = TRUE)
{
    for (cur.run in 1:number.of.runs){
        #cur.run <- 1
        file.name  <- paste( exp.name,                                                      # default: local
                         AT.add.leading.zeros(cur.run, 3),
                                        "_fort.",
                                        unit,
                                        sep = "")
        if(data.source == "condor"){                                                        # condor naming in case cleaning hasn't been run
            file.name            <-    paste(    exp.name,
                                            "_node",
                                            AT.add.leading.zeros(cur.run - 1, 4),
                                            "001_fort.",
                                            unit,
                                            sep = "")
        }
        if(data.source == "condor_cleaned"){
            file.name            <-    paste(    exp.name,                                   # condor
                                            AT.add.leading.zeros(cur.run - 1, 5),
                                            "_fort.",
                                            unit,
                                            sep = "")
        }
        
        input            <- scan(file = file.name, what = "character", strip.white = TRUE, sep = "")

        # find no. of particles
        no.of.particles  <- as.numeric(gsub(",", "", input[grep("followed", input) + 1]))

        # find USRTRACK outputs
        outputs          <- grep("Track", input) 

        # first run only (assuming all number.of.runs are absolutely parallel, except for no.of.particles)
        if (cur.run == 1){
            names           <- gsub(" ", "", input[outputs + 3])
            regions         <- as.numeric(input[outputs + 11])
            bins            <- as.numeric(input[outputs + 34])
            E.low.GeV       <- as.numeric(input[outputs + 30])
            E.high.GeV      <- as.numeric(input[outputs + 32])
            binning         <- input[outputs + 26]
            ii              <- binning == "linear"
            E.step.GeV      <- 0
            E.step.GeV[ii]  <- as.numeric(input[outputs + 37][ii])  		        # linear binning
            E.step.GeV[!ii] <- as.numeric(gsub(")", "", input[outputs + 38][!ii]))  # log binning

            # reformat particle.name to macht libamtrack format        
            particle.names  <- AT.FLUKA.particle.name.to.libamtrack.particle.name( FLUKA.particle.names = as.character(names))
            particle.no     <- AT.particle.no.from.particle.name(particle.names)

            # build data.frame
            df           <- data.frame( particle.no      = rep(particle.no, bins),
                                        region           = rep(regions, bins),
                                        stringsAsFactors = FALSE)
            
            
            # Read in fluence values (as matrix for sake of speed)
            values       <- matrix(ncol = 4, nrow = sum(bins))
            values[,1]   <- 1:sum(bins)
            index        <- c(1, cumsum(bins)[-length(bins)] + 1)
            for (i in 1:length(names)){
                #i <- 1
                idx            <- index[i]:(index[i] + bins[i] - 1)
                if(binning[i] == "linear"){
                    values[idx,2]   <- E.step.GeV[i] * ((1:bins[i]) - 0.5)
                    values[idx,3]   <- rep(E.step.GeV[i], bins[i])
                    values[idx,4]   <- as.numeric(input[outputs[i] + 48 + (1:bins[i]) - 1])
                }else{
                    low             <- E.low.GeV[i] * E.step.GeV[i]^((1:bins[i]) - 1)
                    high            <- E.low.GeV[i] * E.step.GeV[i]^(1:bins[i])
                    values[idx,2]   <- sqrt(low * high)
                    values[idx,3]   <- high - low
                    values[idx,4]   <- as.numeric(input[outputs[i] + 47 + (1:bins[i]) - 1])
                }
            }
        }else{                # for all later number.of.runs, just add fluences and fluences^2 (for stdev)
            for (i in 1:length(names)){
                #i <- 1
                idx            <- index[i]:(index[i] + bins[i] - 1)
                values[idx,4]  <- values[idx,4] + as.numeric(input[outputs[i] + 47 + (1:bins[i]) - 1])
             }
            
        }
    }

    df$E.MeV.u         <- values[,2] * 1000   # convert energy scale GeV/u -> MeV/u
    df$DE.MeV.u        <- values[,3] * 1000
    df$fluence.cm2     <- values[,4] * values[,3] / number.of.runs

    # reduce size
    if(compress == TRUE){
      ii <- df$fluence.cm2 == 0
      df <- df[!ii,]
    }

    return(df)
}
