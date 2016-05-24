make.shortnames.coefficients <- function(thenames,
                                    formula,
                                    method="MLE",
                                    model="additive",
                                    Spline="b-spline",
                                    knots.Bh,
                                    degree.Bh=3,
                                    intercept.Bh=TRUE,
                                    log.Bh=FALSE,
                                    name.runningtime=".t"){

  # the baseline
  # ncoef_baseline = length(Knots.bh) + degree.h + intercept.bh + log.bh  
  thenames[1:(length(knots.Bh) + degree.Bh + intercept.Bh)] <- paste("Baseline hazard:", 1:(length(knots.Bh) + degree.Bh + intercept.Bh), sep="")

  if(log.Bh){
  thenames[1+ length(knots.Bh) + degree.Bh + intercept.Bh] <- "Bh:log(t)"
  }

  if(length(thenames) > 0){
    for (i in 1:length(thenames)) {
      if (substr(thenames[i],1,7)=="NPHNLL(") {
        if (method == "MLE" | (method == "glm" & model == "additive")) {
          
        # get the call to NPHNLL
          
          thecall <- NULL

          if(grepl(":XxXxXXxXxX", thenames[i])){
                  # alpha(x)
            indx <- as.integer(strsplit(thenames[i],":XxXxXXxXxX ")[[1]][2])
            thecall <-strsplit(thenames[i],":XxXxXXxXxX")[[1]][1]
            thecall <- match.call(NPHNLL, parse(text=thecall))
            xtvar <- as.character(thecall[[2]])
          }
          else if (grepl(":TtTtTTtTtT", thenames[i])){
                  # beta(t)
            indx <- strsplit(thenames[i],":TtTtTTtTtT ")[[1]][2]
            thecall <-strsplit(thenames[i],":TtTtTTtTtT")[[1]][1]
            thecall <- match.call(NPHNLL, parse(text=thecall))
            xtvar <- (formula[[2]])[[2]]
          }
          thenames[i] <- paste("NPHNLL(", as.character(thecall[[2]]), ", ", (formula[[2]])[[2]], ")", xtvar, ":", indx, sep="")
        }
        else if (method == "glm" & model == "multiplicative"){
          # replace name.runningtime by the timevar in the initial formula : (formula[[2]])[[2]]
          thenames[i] <- gsub(name.runningtime, (formula[[2]])[[2]], thenames[i])
        }
      }

      # numbre of ")" in the name
      nbr<-length(strsplit(thenames[i],")")[[1]])
      
# NLL effect 
      if (substr(thenames[i],1,4)=="NLL(") {
        thecall <- NULL
        indx <- as.integer(strsplit(thenames[i],")")[[1]][nbr])
        thecall <-paste(paste(strsplit(thenames[i],")")[[1]][1:(nbr-1)], collapse=")"), ")", sep="")
        thecall <- match.call(NLL, parse(text=thecall))
        thenames[i] <- paste("NLL(", as.character(thecall[[2]]), "):", indx, sep="")
      }
        
# NPH effect 
      if (substr(thenames[i],1,4)=="NPH(") {
        thecall <- NULL
        indx <- as.integer(strsplit(thenames[i],")")[[1]][nbr])
        thecall <-paste(paste(strsplit(thenames[i],")")[[1]][1:(nbr-1)], collapse=")"), ")", sep="")
        thecall <- match.call(NPH, parse(text=thecall))
        thenames[i] <- paste("NPH(", as.character(thecall[[2]]), ", ", (formula[[2]])[[2]], "):", indx, sep="")
      }
    }
  }
  return(thenames)

}

