#####
#####
###     Function "opscale" 
###     Version 1.1.
###     July 30, 2014
###

###     William G. Jacoby
###     E-mail: jacoby@msu.edu
###

#####
#####
###     The following five functions, "init.os.1",
###     "init.os.2", "create.u.1", "create.u.2", 
###     and "test.mono", are used in the "opscale"
###     function.
#####
#####

###
###   Function creating initial OS values
###   for discrete data   
###

   init.os.1 <- function (x1, x2) {
      os <- t(as.matrix(tapply(x2, x1, mean, na.rm = T)))
      os
   }
   
###
###   Function creating initial OS values
###   for continuous data   
###

   init.os.2 <- function (x1, x2) {
      for (i in sort(unique(x1))) {
         if (i == min(x1, na.rm = TRUE)) {os <- c()
                          orig <- c() 
         }
     {os <- c(os, sort(unique(x2[x1 == i])))
      orig <- c(orig, rep(i, length(sort(unique(x2[x1 == i])))))
     }
     }
   os <- rbind(os, orig)   
   os
   }

###
###   Function creating indicator matrix,
###   u, for discrete data
###

   create.u.1 <- function (x1) {
      u <- x1 %*% matrix(1,1, length(sort(unique(x1)))) ==
        (matrix(1, length(x1), 1) %*% t(as.matrix(sort(unique(x1)))))
      u
   }

###
###   Function creating indicator matrix,
###   u, for continuous data
###

   create.u.2 <- function (x1, x2, os) {
      upart1 <- x1 %*% matrix(1,1, length(os[1,])) ==
             matrix(1, length(x1), 1) %*% t(as.matrix(os[2,]))
      upart2 <- x2 %*% matrix(1,1, length(os[1,])) ==
             matrix(1, length(x2), 1) %*% t(as.matrix(os[1,]))
      u <- upart1 * upart2
      u
   }

###
###   Function for comparing adjacent values
###   of OS data. Average downward 
###   where necessary, to make OS
###   values weakly monotonic to
###   input qualitative values
###

   test.mono <- function (u, os) {
     
     categ.n <- apply(u,2,sum, na.rm = T)

     for (i in 2:length(os[1,])) {
       if (os[1,i] < os[1,i-1])    {
         for (i2 in 1:(i-1))     {
         val <- sum(os[1,(i-i2):i] * categ.n[(i-i2):i]) /
sum(categ.n[(i-i2):i])
         os[1,(i-i2):i] <- rep(val, i2 + 1)
         if (((i-i2) == 1) || (os[1,(i-i2)] >= os[1,i-i2-1])) break
         }
       }
     }
    os
   }


###
###     The "opscale" function rescales a vector of qualitative
###     data values to a least-squares fit with a vector of
###     quantitative data values, subject to constraints imposed
###     by the measurement level and process of the qualitative data.
###
###     The arguments are as follows:
###
###          x.qual   A vector of data values, assumed nominal or
###                   ordinal level measurement
###         x.quant   A vector of quantitative values, probably 
###                   obtained as model estimates
###           level   Measurement level of qualitative vector, with
###                   1 = nominal and 2 = ordinal                  
###         process   Measurement process of qualitative vector, with
###                   1 = discrete and 2 = continuous  
###       na.impute   FALSE (default) = leave missing values in qualitative
###                   vector as missing in optimally scaled vector, 
###                   TRUE = assign quantitative vector values to 
###                   optimally scaled vector for missing entries in 
###                   qualitative vector
###       na.assign   TRUE (default) = if quantitative value is missing, assign
###                   mean of optimally scaled values for the corresponding
###                   qualitative value to the optimally scaled vector,
###                   FALSE = if quantitative value is missing, leave optimally
###                   scaled value missing, even if qualitative value is present
###         rescale   TRUE (default) = rescale optimally scaled values to the same 
###                   mean and standard deviation as the values in the original 
###                   qualitative vector.  FALSE = do not rescale optimally scaled
###                   vector
#####
#####


opscale <- function (x.qual, 
                     x.quant=seq(1:length(x.qual)), 
                     level=1, 
                     process=1, 
                     na.impute=FALSE, 
                     na.assign=TRUE,
                     rescale=TRUE
                     ) {


###
###   Check input data
###

   if (length(x.qual) != length(x.quant))
     stop ("Input vectors must be same length")
   
   if (!is.numeric(x.quant))
     stop ("Quantitative vector must be numeric")
   
   if (!(level==1 | level==2))
      stop ("Level argument must be 1 or 2")
      
   if (!(process==1 | process==2))
      stop ("Process argument must be 1 or 2")

###
###   Assign name of qualitative vector to object "varname"
###
      
   varname <- deparse(substitute(x.qual))
   
###
###   Test qualitative vector for character entries.
###   If necessary, convert to numeric entries
###

   x.qual.char <- c()
   
   if (!is.numeric(x.qual)) {
      x.qual.char <- c(x.qual.char, x.qual)
      x.qual <- as.numeric(as.factor(x.qual))
      }

###
###   Obtain mean and standard deviation of
###   qualitative values, for rescaling
###   the OS vector
###

   if (rescale == TRUE) {
      qual.mean <- mean(x.qual, na.rm = TRUE)
      qual.sd <- sd(x.qual, na.rm = TRUE)
      }

###
###   Test for qualitative values 
###   that have only missing quantitative
###   values in corresponding entries. When
###   this occurs, assign missing value to 
###   entry in qualitative vector
###

   x.qual[x.qual == sort(unique(x.qual))[
           tapply(x.quant, x.qual, length) ==
           tapply(x.quant, x.qual, function (x) {sum(is.na(x))})]
           ] <- NA

###
###   Create initial set of OS values
###   and indicator matrix
###

     if (level == 1 | process == 1) {
        os <- init.os.1 (x.qual, x.quant)
         u <- create.u.1 (x.qual)
        }
     else {
        os <- init.os.2 (x.qual, x.quant)
         u <- create.u.2 (x.qual, x.quant, os)
        }

###
###   Compare adjacent OS values for
###   ordinal data, and average, as
###   necessary, to maintain monotonicity
###

     if (level == 2) {
        os <- test.mono (u, os)
        }

###
###   Assign final OS values to observations
###   in a new vector, os.data
###

   os.data <- as.vector(u %*% as.matrix(os[1,]))

###
###   If qualitative data are continuous-
###   nominal, repeat analysis, treating as
###   continuous ordinal, relative to OS
###   values already calculated
###

   if (level == 1 & process == 2) {
      os <- init.os.2(os.data, x.quant)
       u <- create.u.2 (os.data, x.quant, os)
      os <- test.mono (u, os)
      os.data <- as.vector(u %*% as.matrix(os[1,]))
      }

###
###   If specified, assign quantitative values to
###   entries that were missing in qualitative vector
###

   if (na.impute == TRUE) {
      os.data[which(is.na(x.qual))] <- x.quant[which(is.na(x.qual))]
      }

###
###   If specified, assign mean optimally scaled value for
###   value of qualitative vector to missing quantitative
###   vector values in optimally scaled vector
###

   if (na.assign == TRUE & sum(is.na(x.quant)) > 0) {
     os.data[which(is.na(x.quant))] <-
     (x.qual[which(is.na(x.quant))] %*%
matrix(1,1,length(sort(unique(x.qual)))) ==
     (matrix(1, length(which(is.na(x.quant))), 1) %*% 
     t(as.matrix(sort(unique(x.qual)))))) %*% 
     as.matrix(tapply(os.data, x.qual, mean, na.rm = T))
     }
   
   if (na.assign == FALSE) {os.data[which(is.na(x.quant))] <- NA}

###
###   If specified, rescale optimally scaled values
###   to same mean and standard deviation as vector
###   of qualitative data values. Collect mean and
###   standard deviation of "raw" optimally scaled
###   values for inclusion in OS object.
###

   os.raw.mean <- mean(os.data, na.rm = TRUE)
   os.raw.sd <- sd(os.data, na.rm = TRUE)
   
   if (rescale == TRUE) {
      os.data <- (as.vector(scale(os.data)) * qual.sd) +
      qual.mean
      }


###
###   Return list, opscaled
###

   if (length(x.qual.char) > 0) {x.qual <- x.qual.char}
   
   if (level == 1) {measlev <- "Nominal"}
      else {measlev <- "Ordinal"}
      
   if (process == 1) {measproc <- "Discrete"}
      else {measproc <- "Continuous"}

   opscaled <- list(qual = x.qual, 
                    quant = x.quant,
                    os = os.data,
                    varname = varname,
                    measlev = measlev,
                    measproc = measproc,
                    rescale = rescale,
                    os.raw.mean = os.raw.mean,
                    os.raw.sd = os.raw.sd
                    )

   class(opscaled) <- "opscale"
   
   opscaled

}

###
###   Function to
###   plot the optimally scaled values
###   against the original data values
###

   os.plot <- function (x.qual, os.data,
      main.title = "Plot of optimal transformation") {

   if (length(x.qual) != length(os.data))
     stop ("Input vectors must be same length")
   
   if (!is.numeric(os.data))
     stop ("Optimally scaled vector must be numeric")
   
      x.qual.values <- c()
         if (!is.numeric(x.qual)) {
            x.qual.values <- levels(as.factor(x.qual))
            x.qual <- as.numeric(as.factor(x.qual))
         }

      x <- init.os.2(x.qual, os.data)
      low <- min(x) - ((max(x) - min(x)) * .05)
      high <- max(x) + ((max(x) - min(x)) * .05)

      if (length(x.qual.values) == 0) {
         osplot <- xyplot(x[1,] ~ x[2,],
            aspect = 1,
            xlim = c(low, high),
            ylim = c(low, high),
            type = "b",
            col = "black",
            xlab = "Original data values",
            ylab = "Optimally scaled values",
            main = main.title
            )
         }
      else {
         osplot <- xyplot(x[1,] ~ x[2,],
            aspect = 1,
            xlim = c(low, high),
            ylim = c(low, high),
            type = "b",
            col = "black",
            xlab = "Original data values",
            ylab = "Optimally scaled values",
            main = main.title,
            scales = list(x = list(
               at = as.numeric(as.factor(x.qual.values)),
               labels = x.qual.values,
               rot = 45))
            )
       }      
      osplot
      }


###
###   Function to 
###   produce a Shepard diagram
###   for the optimally scaled vector
###

   shep.plot <- function (x.quant, os.data,
      main.title = "Shepard Diagram") {

   if (length(x.quant) != length(os.data))
     stop ("Input vectors must be same length")
   
   if (!is.numeric(x.quant))
     stop ("Quantitative vector must be numeric")
   
   if (!is.numeric(os.data))
     stop ("Optimally scaled vector must be numeric")

      shepplot <- xyplot(x.quant ~ os.data,
         aspect = 1,
         col = "black",
         xlab = "Optimally scaled values",
         ylab = "Quantitative data values",
         main = main.title
         )
      shepplot
      }  
   

###
###   Function to 
###   calculate Stress statistics
###

   calc.stress <- function (quant, os, 
      rescale = FALSE, 
      os.raw.mean = mean(os, na.rm = TRUE), 
      os.raw.sd = sd(os, na.rm = TRUE))  {

   if (length(quant) != length(os))
     stop ("Input vectors must be same length")
   
   if (!is.numeric(quant))
     stop ("Quantitative vector must be numeric")
   
   if (!is.numeric(os))
     stop ("Optimally scaled vector must be numeric")

      if (rescale == TRUE) {
      os <- (as.vector(scale(os)) * os.raw.sd) + os.raw.mean
      }
      resids <- quant - os
      raw.stress <- sum(resids^2, na.rm=T)
      stress1 <- (raw.stress / sum(quant^2, na.rm=T))^.5
      stress2 <- (raw.stress / sum((quant - mean(quant, na.rm=T))^2,
na.rm=T))^.5
      scoeffs <- c(stress.1 = stress1, stress.2 = stress2, raw.stress =
raw.stress)
      scoeffs
      }

###
###     Create generic print, summary, and
###     plot methods for objects
###     of class "opscale". Create summary object,
###     and print.summary method for object
###

   print.opscale <- function (x, ...) {
      cat("Optimal Scaling Object for Variable", x$varname,":","\n\n")
      print(data.frame(x[1:3]))
      invisible(x)
      }   

   summary.opscale <- function (object, ...) {
      x.qual.values <- c()
         if (!is.numeric(object$qual)) {
            x.qual.values <- levels(as.factor(object$qual))
            object$qual <- as.numeric(as.factor(object$qual))
         }

      values <- as.data.frame(t(init.os.2(object$qual, object$os)) %*%
matrix(c(0,1,1,0),2,2))
      if (length(x.qual.values > 0)) {
         if (object$measproc == "Discrete") {
         values[,1] <- x.qual.values
         }
         else {
         values[,1] <- rep(x.qual.values, apply(create.u.1(values[,1]), 2,
sum))
         }
      }
      colnames(values) <- c("initial.values", "OS.values")
      rownames(values) <- NULL
      results <- list(values = values, varname = object$varname,
         measlev = object$measlev, measproc = object$measproc)
      
      class(results) <- "summary.opscale"
      results
      }
      
   print.summary.opscale <- function(x, ...) {
      cat("Optimal Transformation for Variable", x$varname,":","\n\n",
         "     Measurement Level: ", x$measlev, "\n",
         "   Measurement Process: ", x$measproc, "\n\n")
      print(x$values)
      invisible(x)
      }
   
   plot.opscale <- function (x, ...) {
      os.plot(x$qual, x$os, 
      main.title = paste("Optimal Transformation for Variable:", x$varname))
      }

###
###   Create methods for calculating
###   stress coefficients and for producing
###   the Shepard diagram from an opscale object
###
    
    stress <- function (x, ...) {
      if (class(x) != "opscale")
         stop("Object must be of class 'opscale'")
       cat("Stress Values for Variable:", x$varname, "\n\n")
       calc.stress(quant = x$quant, 
          os = x$os, 
          rescale = x$rescale, 
          os.raw.mean = x$os.raw.mean, 
          os.raw.sd = x$os.raw.sd)
       }
      

   shepard <- function (x, ...) {
      if (class(x) != "opscale")
         stop("Object must be of class 'opscale'")
      shep.plot(x$quant, x$os,
      main.title = paste("Shepard Diagram for Variable:", x$varname))
      }

