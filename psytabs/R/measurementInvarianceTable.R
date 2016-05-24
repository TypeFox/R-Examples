measurementInvarianceTable <-
function (measurement.invariance) {
  # which `tests' do we have?
  scaled <- FALSE
  TESTS <- unlist(lapply(measurement.invariance$fit.configural@Fit@test, "[", "test"))
  if(any(c("satorra.bentler", "yuan.bentler") %in% TESTS)) {
    scaled <- TRUE
  }
  
  if(!scaled) {
    out.configural <- lavaan::fitMeasures(measurement.invariance$fit.configural, c("chisq", "df", "cfi", "rmsea", "bic"))
    out.configural <- data.frame(t(as.matrix(out.configural)))
    out.loadings <- lavaan::fitMeasures(measurement.invariance$fit.loadings, c("chisq", "df", "cfi", "rmsea", "bic"))
    out.loadings <- data.frame(t(as.matrix(out.loadings)))
    out.intercepts <- lavaan::fitMeasures(measurement.invariance$fit.intercepts, c("chisq", "df", "cfi", "rmsea", "bic"))
    out.intercepts <- data.frame(t(as.matrix(out.intercepts)))
    
    if(!is.null(measurement.invariance$fit.residuals)) { 
      out.residuals <- lavaan::fitMeasures(measurement.invariance$fit.residuals, c("chisq", "df", "cfi", "rmsea", "bic"))
      out.residuals <- data.frame(t(as.matrix(out.residuals)))
      out.means <- lavaan::fitMeasures(measurement.invariance$fit.means, c("chisq", "df", "cfi", "rmsea", "bic"))
      out.means <- data.frame(t(as.matrix(out.means)))
      
      out1 <- rbind(out.configural, out.loadings, out.intercepts, out.residuals, out.means)
      
      out.cm1 <- compareModels(out.loadings, out.configural)
      out.cm2 <- compareModels(out.intercepts, out.loadings)
      out.cm3 <- compareModels(out.residuals, out.intercepts)
      out.cm4 <- compareModels(out.means, out.residuals)
      out2 <- rbind(NA, out.cm1, out.cm2, out.cm3, out.cm4)
      
    } else {
      out.means <- lavaan::fitMeasures(measurement.invariance$fit.means, c("chisq", "df", "cfi", "rmsea", "bic"))
      out.means <- data.frame(t(as.matrix(out.means)))
      
      out1 <- rbind(out.configural, out.loadings, out.intercepts, out.means)
      
      out.cm1 <- compareModels(out.loadings, out.configural)
      out.cm2 <- compareModels(out.intercepts, out.loadings)
      out.cm3 <- compareModels(out.means, out.intercepts)
      out2 <- rbind(NA, out.cm1, out.cm2, out.cm3)
      
    }
    
    
    out <- cbind(out1, out2)
    table <- out[c("chisq", "df", "delta.chisq", "delta.df", "delta.p.value", "cfi", "delta.cfi", "rmsea", "delta.rmsea", "bic", "delta.bic")]
    
    table <- plyr::colwise(makeNumeric)(table)
    table[,c(1, 3, 10, 11)] <- round(table[,c(1, 3, 10, 11)], 1)
    table[,c(5, 6, 7, 8, 9)] <- round(table[,c(5, 6, 7, 8, 9)], 3)
    names(table) <- c("chi2", "df", "Dchi2", "Ddf", "Dpval", "cfi", "Dcfi", "rmsea", "Drmsea", "bic", "Dbic")
      
  } else { #Scaled!
    out.configural <- lavaan::fitMeasures(measurement.invariance$fit.configural, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
    out.configural <- data.frame(t(as.matrix(out.configural)))
    out.loadings <- lavaan::fitMeasures(measurement.invariance$fit.loadings, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
    out.loadings <- data.frame(t(as.matrix(out.loadings)))
    out.intercepts <- lavaan::fitMeasures(measurement.invariance$fit.intercepts, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
    out.intercepts <- data.frame(t(as.matrix(out.intercepts)))
    
    if(!is.null(measurement.invariance$fit.residuals)) { 
      out.residuals <- lavaan::fitMeasures(measurement.invariance$fit.residuals, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
      out.residuals <- data.frame(t(as.matrix(out.residuals)))
      out.means <- lavaan::fitMeasures(measurement.invariance$fit.means, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
      out.means <- data.frame(t(as.matrix(out.means)))
      
      out1 <- rbind(out.configural, out.loadings, out.intercepts, out.residuals, out.means)
      
      out.loadings.scaling = measurement.invariance$fit.loadings@Fit@test[[2]]$scaling.factor
      out.configural.scaling = measurement.invariance$fit.configural@Fit@test[[2]]$scaling.factor
      out.intercepts.scaling = measurement.invariance$fit.intercepts@Fit@test[[2]]$scaling.factor
      out.residuals.scaling = measurement.invariance$fit.residuals@Fit@test[[2]]$scaling.factor
      out.means.scaling = measurement.invariance$fit.means@Fit@test[[2]]$scaling.factor
      
      out.cm1 <- compareModels(out.loadings, out.configural, scaled = TRUE, out.loadings.scaling, out.configural.scaling)
      out.cm2 <- compareModels(out.intercepts, out.loadings, scaled = TRUE, out.intercepts.scaling, out.loadings.scaling)
      out.cm3 <- compareModels(out.residuals, out.intercepts, scaled = TRUE, out.residuals.scaling, out.intercepts.scaling)
      out.cm4 <- compareModels(out.means, out.residuals, scaled = TRUE, out.means.scaling, out.residuals.scaling)
      out2 <- rbind(NA, out.cm1, out.cm2, out.cm3, out.cm4)
      
    } else {
      out.means <- lavaan::fitMeasures(measurement.invariance$fit.means, c("chisq", "df.scaled", "cfi.scaled", "rmsea.scaled", "bic", "chisq.scaled"))
      out.means <- data.frame(t(as.matrix(out.means)))
      
      out1 <- rbind(out.configural, out.loadings, out.intercepts, out.means)
      
      out.loadings.scaling = measurement.invariance$fit.loadings@Fit@test[[2]]$scaling.factor
      out.configural.scaling = measurement.invariance$fit.configural@Fit@test[[2]]$scaling.factor
      out.intercepts.scaling = measurement.invariance$fit.intercepts@Fit@test[[2]]$scaling.factor
      out.means.scaling = measurement.invariance$fit.means@Fit@test[[2]]$scaling.factor
      
      out.cm1 <- compareModels(out.loadings, out.configural, scaled = TRUE, out.loadings.scaling, out.configural.scaling)
      out.cm2 <- compareModels(out.intercepts, out.loadings, scaled = TRUE, out.intercepts.scaling, out.loadings.scaling)
      out.cm3 <- compareModels(out.means, out.intercepts, scaled = TRUE, out.means.scaling, out.intercepts.scaling)
      out2 <- rbind(NA, out.cm1, out.cm2, out.cm3)
      
    }
    
    out <- cbind(out1, out2)
    table <- out[c("chisq.scaled", "df.scaled", "delta.chisq.scaled", "delta.df.scaled", "p.value.scaled", "cfi.scaled", "delta.cfi.scaled", "rmsea.scaled", "delta.rmsea.scaled", "bic", "delta.bic")]
    asNumeric <- function(x) x <- as.numeric(x)
    table <- plyr::colwise(asNumeric)(table)
    table[,c(1, 3, 10, 11)] <- round(table[,c(1, 3, 10, 11)], 1)
    table[,c(5, 6, 7, 8, 9)] <- round(table[,c(5, 6, 7, 8, 9)], 3)
    names(table) <- c("chi2S", "df", "Dchi2S", "Ddf", "Dpval", "cfi", "Dcfi", "rmsea", "Drmsea", "bic", "Dbic")
    
  }
  
  if(dim(table)[1] > 4) {
    rownames(table) <- c("Configural", "Metric", "Scalar", "Residual", "Mean") 
  } else {
    rownames(table) <- c("Configural", "Metric", "Scalar", "Mean")
  }
  
  table
}