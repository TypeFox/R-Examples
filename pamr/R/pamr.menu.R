pamr.menu <- function(data) {
  done <- FALSE
  junk.train <- NULL
  junk.results <- NULL
  while(!done) {
    cat("", fill = TRUE)
    switch(menu(c("pamr.train", "pamr.cv", "pamr.plotcv", 
                  "pamr.plotcen", "pamr.confusion", 
                  "pamr.plotcvprob", "pamr.geneplot", 
                  "pamr.listgenes", 
                  "pamr.train with heterogeneity analysis", 
                  "Exit")),
           junk.train <- pamr.train(data),
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               junk.results <- pamr.cv(junk.train, data)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               pamr.plotcv(junk.results)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.plotcen(junk.train, data, threshold = 
                            threshold)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.confusion(junk.results, threshold = 
                              threshold)
             }
           }
           ,
           {
             if(is.null(junk.results)) {
               cat("Error: need to run pamr.cv first", fill
                   = TRUE)
             }
             if(!is.null(junk.results)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.plotcvprob(junk.results, data, threshold
                               = threshold)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.geneplot(junk.train, data, threshold = 
                             threshold)
             }
           }
           ,
           {
             if(is.null(junk.train)) {
               cat("Error: need to run pamr.train first", 
                   fill = TRUE)
             }
             if(!is.null(junk.train)) {
               cat("threshold?")
               threshold <- scan("", nlines = 1)
               pamr.listgenes(junk.train, data, threshold = 
                              threshold)
             }
           }
           ,
           {
             junkk.train <- NULL
             cat("Normal class?", fill = TRUE)
             normal <- scan("", nlines = 1, what = "")
             junk.train <- pamr.train(data, hetero = normal)
           }
           ,
           done <- TRUE)
  }
  cat("Done\n")
}

pamr.pairscore <-function(x, pair.ind=NULL) {
}

