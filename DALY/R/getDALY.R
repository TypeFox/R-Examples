## Retrieve user inputs
## Calculate DALY
## Return object of S3 class 'DALY'

getDALY <-
function(button.call = FALSE, aw, dr){

  ## Retrieve & evaluate user inputs
  ## NOTE: all inputs must be sent as 1D vector !!

  # Population
  sendPop <- c(DALYget("pop"))
  sendPop[is.na(sendPop)] <- 0

  if (is.na(sum(sendPop))){
    if (button.call)
	  tkmessageBox(title = "Error", icon = "error", type = "ok",
                   message = "Population table contains\nnon-numeric value(s)")
    stop("Population table contains non-numeric value(s).", call. = FALSE)
  } else if (sum(sendPop) == 0){
    if (button.call)
      tkmessageBox(title = "Error", icon = "error", type = "ok",
                   message = "Population table is empty")
    stop("Population table is empty.", call. = FALSE)
  }

  # Age weighting
  aw_gui <- ifelse(DALYtclvalue(".aw") == "Yes", TRUE, FALSE)
  if (missing(aw)){
    sendAW  <- aw_gui
  } else {
    if (!is.logical(aw))
	  stop("Argument 'aw' must be a logical value", call. = FALSE)
	if (aw != aw_gui)
	  warning(paste("Age weighting constant specified in function call",
		            "differs from the one specified in GUI."),
              call. = FALSE)
    sendAW <- aw
  }

  # Discount rate
  if (length(DALYget(".dr")) == 0 & missing(dr)){
    if (button.call)
      tkmessageBox(title = "Error", icon = "error", type = "ok",
                   message = "Please enter a Discount Rate")
    stop("Please enter a Discount Rate.", call. = FALSE)
  } else {
    if (missing(dr)){  # dr only specified in GUI
      sendDR <- DALYtclvalue(".dr")
	} else {  # dr specified in both GUI and function call
	  sendDR <- dr * 100
	  if (length(DALYget(".dr")) != 0 && DALYtclvalue(".dr") != sendDR)
	    warning(paste("Discount rate specified in function call",
		              "differs from the one specified in GUI."),
                call. = FALSE)
	}
  }

  sendDR <- gsub("%", "", sendDR)  # remove % sign, if present
  if (!grepl("^[[:digit:]]*\\.?[[:digit:]]+$", sendDR)){
    if (button.call)
      tkmessageBox(title = "Error", icon = "error", type = "ok",
                   message = "Discount Rate contains non-numeric value")
    stop("Discount Rate contains non-numeric value.", call. = FALSE)
  }
  if (as.numeric(sendDR) < 0 | as.numeric(sendDR) > 100){
    if (button.call)
      tkmessageBox(title = "Error", icon = "error", type = "ok",
                   message = paste("Discount Rate must be",
                                   "a numeric value between 0 and 100%"))
    stop("Discount Rate must be a numeric value between 0 and 100%.",
         call. = FALSE)
  }

  # Data distributions; convert to number (cf C++ code)
  listDst <- NULL
  for (i in seq(8))
    for (j in seq(8))
      listDst[j+(i-1)*8] <-
        DALYget(paste("dist", DALYget("txtLbl")[j], i, sep = ""))
  sendDst <- numeric(length(listDst))
  for (i in seq_along(listDst))
    sendDst[i] <- which(DALYget("distributions") == listDst[i])

  # Data stratifications; convert to number (cf C++ code)
  listStr <- NULL
  for (i in seq(8))
    for (j in seq(8))
      listStr[j+(i-1)*8] <-
        DALYget(paste("strat", DALYget("txtLbl")[j], i, sep = ""))
  sendStr <- numeric(length(listStr))
  for (i in seq_along(listStr))
    sendStr[i] <- which(DALYget("stratifications") == listStr[i])

  sendStrAge <- numeric(length(sendStr))
  for (i in seq_along(sendStr))
    sendStrAge[i] <- c(1, 2, 0, 0)[sendStr[i]]

  sendStrSex <- numeric(length(sendStr))
  for (i in seq_along(sendStr))
    sendStrSex[i] <- ifelse(sendStr[i] == 1 || sendStr[i] == 3, 5, 0)

  # Data parameters
  sendInc <-
    c(getData(DALYget("inc1"), "data"), getData(DALYget("inc2"), "data"),
      getData(DALYget("inc3"), "data"), getData(DALYget("inc4"), "data"),
      getData(DALYget("inc5"), "data"), getData(DALYget("inc6"), "data"),
      getData(DALYget("inc7"), "data"), getData(DALYget("inc8"), "data"))
  sendTrt <-
    c(getData(DALYget("trt1"), "data"), getData(DALYget("trt2"), "data"),
      getData(DALYget("trt3"), "data"), getData(DALYget("trt4"), "data"),
      getData(DALYget("trt5"), "data"), getData(DALYget("trt6"), "data"),
      getData(DALYget("trt7"), "data"), getData(DALYget("trt8"), "data"))
  sendOns <-
    c(getData(DALYget("ons1"), "data"), getData(DALYget("ons2"), "data"),
      getData(DALYget("ons3"), "data"), getData(DALYget("ons4"), "data"),
      getData(DALYget("ons5"), "data"), getData(DALYget("ons6"), "data"),
      getData(DALYget("ons7"), "data"), getData(DALYget("ons8"), "data"))
  sendDur <-
    c(getData(DALYget("dur1"), "data"), getData(DALYget("dur2"), "data"),
      getData(DALYget("dur3"), "data"), getData(DALYget("dur4"), "data"),
      getData(DALYget("dur5"), "data"), getData(DALYget("dur6"), "data"),
      getData(DALYget("dur7"), "data"), getData(DALYget("dur8"), "data"))
  sendDWt <-
    c(getData(DALYget("DWt1"), "data"), getData(DALYget("DWt2"), "data"),
      getData(DALYget("DWt3"), "data"), getData(DALYget("DWt4"), "data"),
      getData(DALYget("DWt5"), "data"), getData(DALYget("DWt6"), "data"),
      getData(DALYget("DWt7"), "data"), getData(DALYget("DWt8"), "data"))
  sendDWn <-
    c(getData(DALYget("DWn1"), "data"), getData(DALYget("DWn2"), "data"),
      getData(DALYget("DWn3"), "data"), getData(DALYget("DWn4"), "data"),
      getData(DALYget("DWn5"), "data"), getData(DALYget("DWn6"), "data"),
      getData(DALYget("DWn7"), "data"), getData(DALYget("DWn8"), "data"))
  sendMrt <-
    c(getData(DALYget("mrt1"), "data"), getData(DALYget("mrt2"), "data"),
      getData(DALYget("mrt3"), "data"), getData(DALYget("mrt4"), "data"),
      getData(DALYget("mrt5"), "data"), getData(DALYget("mrt6"), "data"),
      getData(DALYget("mrt7"), "data"), getData(DALYget("mrt8"), "data"))
  sendDth <-
    c(getData(DALYget("lxp1"), "data"), getData(DALYget("lxp2"), "data"),
      getData(DALYget("lxp3"), "data"), getData(DALYget("lxp4"), "data"),
      getData(DALYget("lxp5"), "data"), getData(DALYget("lxp6"), "data"),
      getData(DALYget("lxp7"), "data"), getData(DALYget("lxp8"), "data"))

  # Life expectancy
  sendLxp <- getLifeExp()
  
  # Iterations
  sendIT  <- DALYget("it")

  # Derive which outcomes are to be used..
  outcomes <- numeric(8)
  for (i in seq(0, 7))
    outcomes[i+1] <-
      (sum(sendInc[((30 * i) + 1):(30 * (i + 1))]) != 0 ||
       sum(sendMrt[((30 * i) + 1):(30 * (i + 1))]) != 0) * (i + 1)
  outcomes <- outcomes[outcomes != 0]

  ## Send data to C++ function
  vector_length <- length(outcomes) * 5 * 2 * sendIT
  output <-
    .C("getMC",
       MRT = as.double(vector("numeric", vector_length)), 
       INC = as.double(vector("numeric", vector_length)), 
       YLD = as.double(vector("numeric", vector_length)), 
       YLL = as.double(vector("numeric", vector_length)),
	   samplesInc = as.double(vector("numeric", vector_length)),
	   samplesTrt = as.double(vector("numeric", vector_length)),
	   samplesOns = as.double(vector("numeric", vector_length)),
	   samplesDur = as.double(vector("numeric", vector_length)),
	   samplesDWt = as.double(vector("numeric", vector_length)),
	   samplesDWn = as.double(vector("numeric", vector_length)),
	   samplesMrt = as.double(vector("numeric", vector_length)),
	   samplesDth = as.double(vector("numeric", vector_length)),
       IT = as.integer(sendIT), 
       AW = as.integer(sendAW),
       DR = as.double(sendDR),
       OC = as.integer(outcomes),
       nOC = as.integer(length(outcomes)), 
       getDist = as.integer(sendDst),
       getStrat = as.integer(sendStr), 
       getStrAge = as.integer(sendStrAge),
       getStrSex = as.integer(sendStrSex),
       getPop = as.double(sendPop), 
       getDur = as.double(sendDur),
       getOns = as.double(sendOns),
       getInc = as.double(sendInc), 
       getTrt = as.double(sendTrt),
       getMrt = as.double(sendMrt),
       getDWt = as.double(sendDWt), 
       getDWn = as.double(sendDWn),
       getDth = as.double(sendDth),
       getLxp = as.double(sendLxp))

  ## Retrieve results
  # YLD, YLL, MRT, INC: list of length 'nOutcomes'
  nOutcomes <- length(outcomes)
  nSamples <- 5 * 2 * sendIT
  
  out <- vector("list", nOutcomes)
  as.DALY.matrix <-
  function(x){
    dim(x) <- c(sendIT, 5, 2)
	dimnames(x) <- list(NULL, DALYget("ageGroups"), c("Male", "Female"))
	return(x)
  }
  
  for (i in seq(nOutcomes)){
    range <- seq(nSamples) + nSamples * (i-1)
    YLD <- as.DALY.matrix(output$YLD[range])
    YLL <- as.DALY.matrix(output$YLL[range])
    DALY <- as.DALY.matrix(output$YLD[range] + output$YLL[range])
    MRT <- as.DALY.matrix(output$MRT[range])
    INC <- as.DALY.matrix(output$INC[range])
	out[[i]] <-
	  list(DALY = DALY, YLD = YLD, YLL = YLL,
           cases = INC, deaths = MRT,
           name = DALYtclvalue(paste("outcome", i, "Name", sep = "")),
		   input = list(inc = as.DALY.matrix(output$samplesInc[range]),
		                trt = as.DALY.matrix(output$samplesTrt[range]),
		                ons = as.DALY.matrix(output$samplesOns[range]),
		                dur = as.DALY.matrix(output$samplesDur[range]),
		                DWt = as.DALY.matrix(output$samplesDWt[range]),
		                DWn = as.DALY.matrix(output$samplesDWn[range]),
		                mrt = as.DALY.matrix(output$samplesMrt[range]),
                        dth = as.DALY.matrix(output$samplesDth[range])))
  }

  out[["pop"]] <- getData(DALYget("pop"), "pop")  # get '0' instead of 'NA'
  out[["name"]] <- DALYtclvalue("diseaseName")

  ## Define S3 class
  class(out) <- c("DALY", "list")

  if (button.call){
    show.relative <- DALYget("optRA") == "Relative (per 1000 pop)"
	show.outcomes <- DALYget("optOC") == "Per outcome"
	if (DALYget("optOP") == "Summed over age/sex classes")
	  print(out, relative = show.relative, outcomes = show.outcomes)
	if (DALYget("optOP") == "Per age/sex class")
	  summary(out, relative = show.relative, outcomes = show.outcomes)
	if (DALYget("optHist") == 1)
	  hist(out)
  }

  return(out)
}
