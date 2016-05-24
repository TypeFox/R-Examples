#'Visualise the best climate window
#'
#'Create a scatterplot showing the fit of the best climate window model through
#'the biological data.
#'@param dataset A dataframe containing information on all fitted climate 
#'  windows. Output from \code{\link{climatewin}}.
#'@param bestmodel A model object. The strongest climate window model. Output 
#'  from \code{\link{singlewin}} or \code{\link{climatewin}}.
#'@param bestmodeldata A dataframe with the data used to 
#'  fit the strongest climate window model. Output from \code{\link{singlewin}} 
#'  or \code{\link{climatewin}}.
#'@return Returns a scatterplot with a fitted line to show the fit of the best 
#'  model through the data.
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'# Visualise the best climate window from the datasets Mass and MassClimate
#'
#'data(MassOutput)
#'data(Mass)
#'data(MassClimate)
#'
#'single <- singlewin(xvar = list(Temp = MassClimate$Temp), 
#'                    cdate = MassClimate$Date, bdate = Mass$Date, 
#'                    baseline = lm(Mass ~ 1, data = Mass),
#'                    furthest = 72, closest = 15, 
#'                    stat = "mean", func = "lin", 
#'                    type = "fixed", cutoff.day = 20, cutoff.month = 5, 
#'                    cmissing = FALSE, cinterval = "day")
#'            
#'plotbest(dataset = MassOutput, bestmodel = single$BestModel,
#'         bestmodeldata = single$BestModelData)
#'              
#'@import ggplot2
#'@export

plotbest <- function(dataset, bestmodel, bestmodeldata){
  names(bestmodeldata)[1] <- "Yvar"
  
  if(is.null(bestmodeldata$WGdev) == FALSE){
    with(bestmodeldata, {
      ggplot(bestmodeldata, aes(y = Yvar, x = climate))+
        geom_point(size = 1, alpha = 0.5)+
        geom_abline(intercept = coef(bestmodel)[1], slope = coef(bestmodel)[2], colour = "red")+
        geom_abline(intercept = coef(bestmodel)[1], slope = coef(bestmodel)[3], colour = "blue")+
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16)) +
        ggtitle("Output of best model") +
        ylab("Biological variable")
    })
  } else {
    if (dataset$Function[1] == "log"){
      names(bestmodeldata)[(ncol(bestmodeldata)-1)] <- "climate"
    }
    if (dataset$Function[1] == "inv"){
      names(bestmodeldata)[ncol(bestmodeldata) - 1]   <- "climate"
      class(bestmodeldata[, ncol(bestmodeldata) - 1]) <- class(bestmodeldata[, ncol(bestmodeldata) - 1])[-match("AsIs", class(bestmodeldata[, ncol(bestmodeldata) - 1]))]
      #WHEN WE USE INVERSE FUNCTION 'climate' becomes class AsIs which the graphs can't deal with
      #With this class change, we turn the 'climate' value in to a basic numeric.
    }
    
    #TEST IF THERE ARE MODEL WEIGHTS
    if(is.null(weights(bestmodel)) == TRUE || sum(weights(bestmodel)) == nrow(bestmodeldata)){
      #TEST IF THERE ARE ADDITIONAL COVARIATES IN THE MODEL
      if(ncol(bestmodeldata) == 2 || dataset$Function[1] == "quad" & ncol(bestmodeldata) == 3 || dataset$Function[1] == "cub" & ncol(bestmodeldata) == 4){
        with(bestmodeldata, {
          ggplot(bestmodeldata, aes(x = climate, y = Yvar), environment = environment()) +
            geom_point(size = 1, alpha = 1) +
            geom_line(data = cbind(bestmodeldata, pred = predict(bestmodel, type = "response", allow.new.levels = TRUE)), aes(y = pred)) +
            theme_classic() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(size = 0.25, colour = "black"),
                  plot.title = element_text(size = 16)) +
            ggtitle("Output of best model") +
            ylab("Biological variable") +    
            if (dataset$Function[1] == "log"){
              xlab("Log of climate variable")
            } else if (dataset$Function[1] == "inv"){
              xlab("Inverse of climate variable")
            } else {
              xlab("Climate variable")
            }
        }
        )       
      } else {
        col <- 1
        if(dataset$Function[1] == "quad"){
          col <- 2
        } 
        if(dataset$Function[1] == "cub"){
          col <- 3
        }
        xval <- seq(from = min(bestmodeldata$climate), to = max(bestmodeldata$climate),
                    by = (max(bestmodeldata$climate) - min(bestmodeldata$climate)) / (nrow(bestmodeldata)))
        #When we have additional covariates, we need to integrate them in to the dataset for the predictions
        #However, we need to determine the mean (or reference category) for each of these variables
        newdat <- matrix(ncol = ncol(bestmodeldata) - col, nrow = nrow(bestmodeldata) + 1)
        #Create a matrix which has columns for all variables (bar Yvar because this will be calculated with predict)
        #nrow is 1 larger than predicted due to the length of xval
        newdat <- as.data.frame(newdat)
        newdat[, 1] <- xval
        #The first column of the matrix will always be the same as the xval
        for(cols in 2:(ncol(bestmodeldata) - col)){ #This will go through every column except for Yvar
          if(is.character(bestmodeldata[, cols]) == FALSE){ #If the variable is not categorical then take the mean
            newdat[, cols] <- mean(bestmodeldata[, cols])
          } else { #If it is categorical, simply take the first category
            newdat[, cols] = bestmodeldata[1, cols]          
          }
        }
        names(newdat) <- c("climate", names(bestmodeldata)[2:(ncol(bestmodeldata) - col)])  #Then change the names so they match what would be in the model
        pred          <- predict(bestmodel, newdata = newdat, type = "response", allow.new.levels = TRUE)
        with(bestmodeldata, {
          ggplot(bestmodeldata, aes(x = climate, y = Yvar), environment = environment()) +
            geom_point(size = 1, alpha = 1) +
            geom_line(data = newdat, aes(y = pred)) +
            theme_classic() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(size = 0.25, colour = "black"),
                  plot.title = element_text(size = 16)) +
            ggtitle("Output of best model") +
            ylab("Biological variable") +    
            if (dataset$Function[1] == "log"){
              xlab("Log of climate variable")
            } else if (dataset$Function[1] == "inv"){
              xlab("Inverse of climate variable")
            } else {
              xlab("Climate variable")
            }
        }
        )  
      }
    } else {
      if(ncol(bestmodeldata) == 3 || dataset$Function[1] == "quad" & ncol(bestmodeldata) == 4 || dataset$Function[1] == "cub" & ncol(bestmodeldata) == 5){ 
        if (dataset$Function[1] == "log" || dataset$Function[1] == "inv"){
          names(bestmodeldata)[ncol(bestmodeldata) - 1] <- "climate"  
        }
        with(bestmodeldata, {
          ggplot(bestmodeldata, aes(x = climate, y = Yvar), environment = environment()) +
            geom_point(size = 1, alpha = 1) +
            geom_line(data = cbind(bestmodeldata, pred = predict(bestmodel, type = "response", allow.new.levels = TRUE)), aes(y = pred)) +
            theme_classic() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(size = 0.25, colour = "black"),
                  plot.title = element_text(size = 16)) +
            ggtitle("Output of best model") +
            ylab("Biological variable") +    
            if (dataset$Function[1] == "log"){
              xlab("Log of climate variable")
            } else if (dataset$Function[1] == "inv"){
              xlab("Inverse of climate variable")
            } else {
              xlab("Climate variable")
            }
        }
        )
      } else {
        col <- 1
        if(dataset$Function[1] == "quad"){
          col <- 2
        } 
        if(dataset$Function[1] == "cub"){
          col <- 3
        }
        xval <- seq(from = min(bestmodeldata$climate), to = max(bestmodeldata$climate),
                    by = (max(bestmodeldata$climate) - min(bestmodeldata$climate)) / (nrow(bestmodeldata)))
        #When we have additional covariates, we need to integrate them in to the dataset for the predictions
        #However, we need to determine the mean (or reference category) for each of these variables
        newdat <- matrix(ncol = ncol(bestmodeldata) - col, nrow = nrow(bestmodeldata) + 1)
        #Create a matrix which has columns for all variables (bar Yvar because this will be calculated with predict)
        #nrow is 1 larger than predicted due to the length of xval
        newdat      <- as.data.frame(newdat)
        newdat[, 1] <- xval
        #The first column of the matrix will always be the same as the xval
        for(cols in 2:(ncol(bestmodeldata) - col)){ #This will go through every column except for Yvar
          if(is.character(bestmodeldata[, cols]) == FALSE){ #If the variable is not categorical then take the mean
            newdat[, cols] <- mean(bestmodeldata[, cols])
          } else { #If it is categorical, simply take the first category
            newdat[, cols] = bestmodeldata[1, cols]          
          }
        }
        names(newdat) <- c("climate", names(bestmodeldata)[2:(ncol(bestmodeldata) - col)])  #Then change the names so they match what would be in the model
        pred <- predict(bestmodel, newdata = newdat,  type = "response", allow.new.levels = TRUE)
        with(bestmodeldata, {
          ggplot(bestmodeldata, aes(x = climate, y = Yvar), environment = environment()) +
            geom_point(size = 1, alpha = 1) +
            geom_line(data = newdat, aes(y = pred)) +
            theme_classic() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(size = 0.25, colour = "black"),
                  plot.title = element_text(size = 16)) +
            ggtitle("Output of best model") +
            ylab("Biological variable") +    
            if (dataset$Function[1] == "log"){
              xlab("Log of climate variable")
            } else if (dataset$Function[1] == "inv"){
              xlab("Inverse of climate variable")
            } else {
              xlab("Climate variable")
            }
        })
      }
    } 
  }
}