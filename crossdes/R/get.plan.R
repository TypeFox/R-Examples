get.plan <- function(trt, k = trt, maxsub = 1000) {
  
  if (trt < 2) {
    stop("Number of treatments must be at least 2")
  }
  if ((k < 2) || (k > trt)) {
    stop("Number of periods must be larger than one and no larger than the number of treatments.")
  }
  
  primep100 <- c(2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 
    37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97)
  
  choi <- choices(trt, k, maxsub)  # Get possible methods of construction for the specified parameters
  if (!sum(choi[[1]])) {
    # Get new parameters, until there is at least one feasible method
    while (!sum(choi[[1]])) {
      dummy <- menu(c("Increase the maximum number of subjects", "Choose a different value of k", "Exit"),
                    title = cat("\n", "I don't have a design for just", maxsub, 
                                 "subjects.", "\n", "Please choose one of the following items.", 
                                 "\n"))
      if (dummy == 1) {
        cat("Please specify the maximum number of subjects.", "\n")
        maxsub <- menu(as.character(2:10000))
      }
      if (dummy == 2) {
        cat("Please specify the value of k.", "\n")
        k <- menu(as.character(2:trt))
      }
      if (dummy == 3) {
        stop("Exit function.", "\n")
      }
      choi <- choices(trt, k, maxsub)
    }
  }
  
  choichar <- c("all.combin", "williams", "des.MOLS", "williams.BIB", "balmin.RMD")[choi[[1]]]
  cat("Possible constructions and minimum numbers of subjects:", "\n")
  showchoices <- rbind(choichar, choi[[2]][choi[[1]]])
  rownames(showchoices) <- c("Method: ", "Number: ")
  colnames(showchoices) <- 1:sum(choi[[1]])
  print(showchoices, quote = FALSE)
  cat("\n")
  
  nextchoi <- menu(c(choichar, "Exit"), title = "Please choose one of the following constructions")
  # Choose one of the possible methods
  
  if (nextchoi == (length(choichar) + 1)) {
    stop("Exit function")
  }
  
  construct <- which(c("all.combin", "williams", "des.MOLS", "williams.BIB", "balmin.RMD") == choichar[nextchoi])
  maxsubchoicon <- maxsub%/%choi[[2]][construct]
  cat(choichar[nextchoi], "selected. How many 'replicates' do you wish (1 -", maxsubchoicon, ")?", "\n")
  replic <- menu(as.character(1:maxsubchoicon))
  if (replic > maxsubchoicon) {
    replic <- maxsubchoicon
  }
  cat(replic, "replicate(s) chosen", "\n")
  # Choose the number of 'replicates', determining the number of subjects which
  # is replic*choi[[2]][construct].  Additional subjects are assigned to the
  # replicates.
  
  # Generate the design
  
  # if( construct==3){ primefact<-matrix(
  # c(2,3,2,5,7,2,3,11,13,2,17,19,23,5,3,29,31,2,37,41,43,47,7,53,59,
  # 61,2,67,71,73,79,3,83,89,97,
  # 1,1,2,1,1,3,2,1,1,4,1,1,1,2,3,1,1,5,1,1,1,1,2,1,1,1,6,1,1,1,1,4,1,1,1),ncol=2)
  # trtpp <- primefact[primep100==trt,] } # The primepowers in primep100, in the
  # representation p^n # This was formerly used in des.MOLS.
  
  if (construct == 4) {
    bibdsub <- ifelse(!(k%%2), (choi[[2]][construct])/k, (choi[[2]][construct])/(2 *  k))
    # Note in Method 4: choi[[2]][con..] is the number of subjects for the
    # resulting design, the BIBD has only this number divided by k resp. 2k
    # subjects.
    lookforBIB <- find.BIB(trt, bibdsub, k)
    if (!all(isGYD(lookforBIB, FALSE, FALSE)[[1]][1:4])) {
      stop("Sorry. No BIBD found for these parameters. Please try again.")
    }
  }
  
  des <- switch(construct, allcombs(trt, k), williams(trt), des.MOLS(trt, k), 
                williams.BIB(lookforBIB), balmin.RMD(trt, choi[[2]][construct], k))
  
  # Now replicate the design as requested
  
  if (replic > 1) {
    des <- kronecker(rep(1, replic), des)
  }
  
  cat("Rows represent subjects, columns represent periods.", "\n", "\n")
  
  des
} 
