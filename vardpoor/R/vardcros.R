vardcros <- function(Y, H, PSU, w_final, id,
                     Dom = NULL,
                     Z = NULL, 
                     country, periods,
                     dataset = NULL,
                     linratio = FALSE,
                     percentratio=1,
                     use.estVar = FALSE,
                     household_level_max = TRUE,
                     withperiod = TRUE,
                     netchanges = TRUE,
                     confidence = .95) {
 
  ### Checking
  if (length(linratio) != 1 | !any(is.logical(linratio))) stop("'linratio' must be the logical value")
  if (length(percentratio) != 1 | !any(is.integer(percentratio) | percentratio > 0)) stop("'percentratio' must be the positive integer value")
  if (length(netchanges) != 1 | !any(is.logical(netchanges))) stop("'netchanges' must be the logical value")
  if (length(withperiod) != 1 | !any(is.logical(withperiod))) stop("'withperiod' must be the logical value")
  if (length(use.estVar) != 1 | !any(is.logical(use.estVar))) stop("'use.estVar' must be the logical value")
  if (length(household_level_max) != 1 | !any(is.logical(household_level_max))) stop("'household_level_max' must be the logical value")
  if(length(confidence) != 1 | any(!is.numeric(confidence) |  confidence < 0 | confidence > 1)) {
          stop("'confidence' must be a numeric value in [0,1]")  }
  if (is.null(Z)==linratio & linratio==TRUE) stop("'linratio' must be FALSE")

  if(!is.null(dataset)) {
      dataset <- data.frame(dataset)
      aY <- Y
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) {
                                Y <- data.frame(dataset[, Y], check.names=FALSE)
                                names(Y) <- aY }
      if(!is.null(H)) {
          aH <- H  
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) {
                                H <- as.data.frame(dataset[, aH], stringsAsFactors=FALSE)
                                names(H) <- aH }}
      if(!is.null(id)) {
          aid <- id  
          if (min(id %in% names(dataset))!=1) stop("'id' does not exist in 'dataset'!")
          if (min(id %in% names(dataset))==1) {
                                id <- as.data.frame(dataset[, aid], stringsAsFactors=FALSE)
                                names(id) <- aid }}
     if(!is.null(PSU)) {
          aPSU <- PSU  
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) {
                                PSU <- as.data.frame(dataset[, aPSU], stringsAsFactors=FALSE)
                                names(PSU) <- aPSU }}
      if(!is.null(w_final)) {
          aw_final <- w_final  
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) {
                                w_final <- data.frame(dataset[, aw_final])
                                names(w_final) <- aw_final }}
      if(!is.null(Z)) {
          aZ <- Z
          if (min(Z %in% names(dataset))!=1) stop("'Z' does not exist in 'dataset'!")
          if (min(Z %in% names(dataset))==1) {
                                Z <- data.frame(dataset[, aZ], check.names=FALSE, stringsAsFactors=FALSE)
                                names(Z) <- aZ }}
      if(!is.null(country)) {
          country2 <- country
          if (min(country %in% names(dataset))!=1) stop("'country' does not exist in 'dataset'!")
          if (min(country %in% names(dataset))==1) country <- as.data.frame(dataset[, country], stringsAsFactors=FALSE)
          names(country) <- country2  }

      if(!is.null(periods)) {
          periods2 <- periods
          if (min(periods %in% names(dataset))!=1) stop("'periods' does not exist in 'dataset'!")
          if (min(periods %in% names(dataset))==1) periods <- data.frame(dataset[, periods], stringsAsFactors=FALSE)
          names(periods) <- periods2  }
     
      if (!is.null(Dom)) {
          Dom2 <- Dom
          if (min(Dom %in% names(dataset))!=1) stop("'Dom' does not exist in 'data'!")
          if (min(Dom %in% names(dataset))==1) {  
                  Dom <- as.data.frame(dataset[, Dom2], stringsAsFactors=FALSE) 
                  names(Dom) <- Dom2 }    }
      }

  # Y
  Y <- data.table(Y, check.names=TRUE)
  n <- nrow(Y)
  m <- ncol(Y)
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numeric values")
  if (any(is.na(Y))) stop("'Y' has unknown values")
  if (is.null(names(Y))) stop("'Y' must be colnames")
  
  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' length must be equal with 'Y' row count")
  if (ncol(H) != 1) stop("'H' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colname")
  if (names(H)=="dataH_stratas") stop("'H' must be different colname")
  H[, (names(H)):=lapply(.SD, as.character)]

  # id
  id <- data.table(id)
  if (any(is.na(id))) stop("'id' has unknown values")
  if (nrow(id) != n) stop("'id' length must be equal with 'Y' row count")
  if (ncol(id) != 1) stop("'id' must be 1 column data.frame, matrix, data.table")
  if (is.null(names(id))||(names(id)=="id")) setnames(id, names(id), "h_ID")

  # PSU
  PSU <- data.table(PSU)
  if (any(is.na(PSU))) stop("'PSU' has unknown values")
  if (nrow(PSU) != n) stop("'PSU' length must be equal with 'Y' row count")
  if (ncol(PSU) != 1) stop("'PSU' has more than 1 column")
  
  # w_final 
  w_final <- data.frame(w_final)
  if (nrow(w_final) != n) stop("'w_final' must be equal with 'Y' row count")
  if (ncol(w_final) != 1) stop("'w_final' must be vector or 1 column data.frame, matrix, data.table")
  w_final <- w_final[,1]
  if (!is.numeric(w_final)) stop("'w_final' must be numerical")
  if (any(is.na(w_final))) stop("'w_final' has unknown values") 
  
  # country
  country <- data.table(country)
  if (any(is.na(country))) stop("'country' has unknown values")
  if (nrow(country) != n) stop("'country' length must be equal with 'Y' row count")
  if (ncol(country) != 1) stop("'country' has more than 1 column")
  if (!is.character(country[[names(country)]])) stop("'country' must be character")


  # periods
  if (withperiod) {
        periods <- data.table(periods)
        if (any(is.na(periods))) stop("'periods' has unknown values")
        if (nrow(periods) != n) stop("'periods' length must be equal with 'Y' row count")
    } else if (!is.null(periods)) stop("'periods' must be NULL for those data")

  # Dom
  namesDom <- NULL
  if (!is.null(Dom)) {
    Dom <- data.table(Dom)
    if (any(duplicated(names(Dom)))) 
           stop("'Dom' are duplicate column names: ", 
                 paste(names(Dom)[duplicated(names(Dom))], collapse = ","))
    if (nrow(Dom) != n) stop("'Dom' and 'Y' must be equal row count")
    if (any(is.na(Dom))) stop("'Dom' has unknown values")
    if (is.null(names(Dom))) stop("'Dom' must be colnames")
    namesDom <- names(Dom)
    Dom <- Dom[, lapply(.SD, as.character), .SDcols = namesDom]
    Dom_agg <- Dom[,.N, keyby=namesDom][,N:=NULL]
    Dom_agg1 <- Dom_agg[, lapply(namesDom, function(x) make.names(paste0(x,".", get(x))))]
    Dom_agg1[, Dom := Reduce(function(x, y) paste(x, y, sep="__"), .SD)]
    Dom_agg <- data.table(Dom_agg, Dom_agg1[, "Dom", with=FALSE])
  }

  # Z
  if (!is.null(Z)) {
    Z <- data.table(Z, check.names=TRUE)
    if (!all(sapply(Z, is.numeric))) stop("'Z' must be numeric values")
    if (nrow(Z) != n) stop("'Z' and 'Y' must be equal row count")
    if (ncol(Z) != m) stop("'Z' and 'Y' must be equal column count")
    if (any(is.na(Z))) stop("'Z' has unknown values")
    if (is.null(names(Z))) stop("'Z' must be colnames")
  }
      
  
  # Calculation
      
  # Domains
  N <- nameY <- nameZ <- variable <- NULL
  sample_size <- totalY <- totalZ <- Z1 <- NULL

  size <- data.table(size=rep(1, nrow(Y)))
  if (!is.null(Dom)) size1 <- domain(size, Dom) else size1 <- copy(size)

  if (!is.null(Dom)) Y1 <- domain(Y, Dom) else Y1 <- Y

  namesDom <- names(Dom)
  DTp <- copy(country)
  if (withperiod) DTp <- data.table(periods, country)
  namesperc <- names(DTp)
  namesperc2 <- c("period_country", names(DTp))
  period_country <- do.call("paste", c(as.list(DTp), sep="_")) 

  if (!is.null(Z)) {
       if (!is.null(Dom)) Z1 <- domain(Z, Dom) else Z1 <- Z               
       if (linratio){ 
                   sorts <- unlist(split(Y1[, .I], period_country))
                   lin1 <- lapply(split(Y1[, .I], period_country), 
                                  function(i) lin.ratio(Y1[i], Z1[i], w_final[i],
                                                        Dom=NULL, percentratio=percentratio))
                   Y2 <- rbindlist(lin1)[sorts]
                   if (any(is.na(Y2))) print("Results are calculated, but there are cases where Z = 0")
                  } else Y2 <- data.table(copy(Y1), copy(Z1))
    } else Y2 <- copy(Y1)

  namesY2 <- names(Y2)
  namesY2w <- paste0(namesY2, "w")
 
  Y2w <- Y2 * w_final
  setnames(Y2w, namesY2, namesY2w)

  size1w <- size1 * w_final
  size1wn <- paste0(names(size1), "w")
  setnames(size1w, names(size1), size1wn)
 
  DT <- data.table(period_country, DTp, H, PSU, id,
                    w_final, size1, size1w, Y2, Y2w) 

  namesY <- names(Y)
  namesZ <- names(Z)
  names_H <- names(H)
  namesY1 <- names(Y1)
  namesZ1 <- names(Z1)
  names_id <- names(id)
  names_PSU <- names(PSU)
  names_size1 <- names(size1)
  namesYZ <- c(namesY, namesZ)
  namesY1Z1 <- c(namesY1, namesZ1)
  names_country <- names(country)
  H <- PSU <- id <- DTp <- country <- NULL

  DTagg <- data.table(DT[, namesperc, with=FALSE], w_final)
  if (!is.null(Dom)) DTagg <- data.table(DTagg, Dom)
  DTagg <- data.table(DTagg, sample_size=1, 
                     pop_size=w_final, w_final*Y)
  if (!is.null(Z)) DTagg <- data.table(DTagg, w_final*Z)

  DTagg <- DTagg[, lapply(.SD, sum, na.rm=TRUE),
                    keyby=c(namesperc, namesDom),
                   .SDcols=c("sample_size", "pop_size", namesYZ)]
  DTaggs <- DTagg[, c(namesperc, namesDom,
                    "sample_size", "pop_size"), with=FALSE] 

  vars <- data.table(variable=namesY, namesY=namesY)
  if (!is.null(namesZ)) vars <- data.table(variable=as.character(1:length(namesY)), 
                                           namesY=namesY, namesZ=namesZ) 

  varsYZ <- list(namesY)
  if (!is.null(namesZ)) varsYZ <- list(namesY, namesZ)
  DTagg <- melt(DTagg, id=c(namesperc, namesDom),
                       measure=varsYZ,
                       variable.factor=FALSE)
  setnames(DTagg, "value1", "totalY")
  if (!is.null(Z)) setnames(DTagg, "value2", "totalZ")

  DTagg <- merge(DTagg, vars, by="variable")
  DTagg[, variable:=NULL]

  DTagg <- merge(DTagg, DTaggs, all.x=TRUE,
                 by=c(namesperc, namesDom))

  if (!is.null(namesDom)) DTagg[,(paste0(namesDom, "_new")):=lapply(namesDom,
                     function(x) make.names(paste0(x,".", get(x))))]

  varsYZ <- vars <- nameY1 <- nameZ1 <- valueY1 <- valueZ1 <- Dom <- NULL
  Z1 <- Y1 <- period_country <- Y2 <- Y2w <- total <- pop_size <- NULL
  stderr_nw <- nhcor <- num1 <- num <- den1 <- den <- num_den1 <- NULL
  grad1 <- grad2 <- estim <- sd_nw <- stderr_w <- sd_w <- se <- rse <- NULL
  cv <- CI_lower <- absolute_margin_of_error <- CI_upper <- totalZ <- NULL
  relative_margin_of_error <- NULL

  #--------------------------------------------------------*
  # AGGREGATION AT PSU LEVEL ("ULTIMATE CLUSTER" APPROACH) |
  #--------------------------------------------------------*

  DT1 <- DT[, lapply(.SD, sum, na.rm=TRUE), keyby=c(namesperc2,
               names_H, names_PSU), .SDcols = namesY2w]
  setnames(DT1, namesY2w, namesY2)

  DTnet <- copy(DT1)
  if (!netchanges)  DTnet <- NULL

  DT2 <- DT1[, lapply(.SD, sum, na.rm=TRUE),
            keyby=namesperc, .SDcols = namesY2]
  varsYZ <- list(namesY1)
  if (!is.null(namesZ1) & !linratio) varsYZ <- list(namesY1, namesZ1)

  DT2 <- melt(DT2, id=namesperc,
                      measure=varsYZ,
                      variable.factor=FALSE)
  if (!is.null(namesZ1) & !linratio){ setnames(DT2, c("value1", "value2"),
                                                    c("valueY1", "valueZ1"))
                     } else setnames(DT2, "value1", "valueY1")
  
  if (!is.null(namesZ1) & !linratio) {
                   vars <- data.table(variable=1:length(namesY1))
                       } else vars <- data.table(variable=namesY1)

  if (!is.null(namesDom)) { vars <- data.table(vars, nameY1=namesY1,
                                               t(data.frame(strsplit(namesY1, "__"))))
                            setnames(vars, names(vars)[3:length(vars)], 
                                        c("namesY", paste0(namesDom, "_new")))
                     } else {vars <- data.table(vars, nameY1=namesY1, namesY=namesY1) }
                           
  if (!is.null(namesZ1)) { vars <- data.table(vars, nameZ1=namesZ1)
                           if (!is.null(namesDom)) {
                                      varsZ <- data.table(nameZ1=namesZ1, 
                                                          t(data.frame(strsplit(namesZ1, "__"))))
                                      setnames(varsZ, names(varsZ)[2:length(varsZ)],
                                                    c("namesZ", paste0(namesDom, "_new"))) 
                                      varsZ[, (paste0(namesDom, "_new")):=NULL] 
                                      vars <- merge(vars, varsZ, by="nameZ1") 
                               } else vars[, namesZ:=nameZ1]  }

  vars <- vars[, lapply(vars, as.character)]

  DT2 <- merge(DT2, vars, by="variable")
  DT2[, variable:=NULL]
  vars <- varsZ <- NULL 

  vars <- c(namesperc, paste0(namesDom, "_new"), "namesY", "namesZ")
  vars <- names(DT2)[names(DT2) %in% vars]

  DTagg <- merge(DTagg, DT2, by=vars)
  DT2 <- vars <- NULL


  # VECTOR OF THE PARTIAL DERIVATIVES (GRADIENT FUNCTION)

  if (!is.null(namesZ1) & !linratio) {
      DTagg[, grad1:=1/valueZ1]
      DTagg[, grad2:=-valueY1/valueZ1^2]
  }

  # NUMBER OF PSUs PER STRATUM
  setkeyv(DT1, c(namesperc2, names_H))

  stratasf <- nh <- nhcor <- NULL
  DT1[, nh:=.N, by=c(namesperc2, names_H)]



 #--------------------------------------------------------------------------*
 # MULTIVARIATE REGRESSION APPROACH USING STRATUM DUMMIES AS REGRESSORS AND |
 # STANDARD ERROR ESTIMATION 						      |
 #--------------------------------------------------------------------------*


  DT1H <- DT1[[names_H]]
  DT1H <- factor(DT1H)
  if (length(levels(DT1H))==1) { DT1[, stratasf:=1]
                                 DT1H <- "stratasf"
                       }  else { DT1H <- data.table(model.matrix( ~ DT1H-1, DT1H,  contrasts = "contr.SAS"))
                                 DT1 <- cbind(DT1, DT1H)
                                 DT1H <- names(DT1H) }

  fits <-lapply(1:length(namesY1), function(i) {
           fitss <- lapply(split(DT1, DT1$period_country), function(DT1c) {

                        y <- namesY1[i]
                        if ((!is.null(namesZ1))&(!linratio)) z <- paste0(",", toString(namesZ1[i])) else z <- ""

                        funkc <- as.formula(paste("cbind(", trim(toString(y)), z, ")~ 0 + ",
                                       paste(c(0, DT1H), collapse= "+")))

                        res1 <- lm(funkc, data=DT1c)
                        if (use.estVar==TRUE) {res1 <- data.table(crossprod(res1$res))
                                } else res1 <- data.table(res1$res)
                        setnames(res1, names(res1)[1], "num") 
                        res1[, nameY1:=y]
                        if (!is.null(namesZ1) & !linratio) {
                              setnames(res1, names(res1)[2], "den")
                              res1[, nameZ1:=namesZ1[i]]
                            }


                        if (use.estVar==TRUE) {    
                              setnames(res1, "num", "num1") 
                              if (!is.null(namesZ1) & !linratio) {
                                       res1[, num_den1:=res1[["den"]][1]]
                                       res1[, den1:=res1[["den"]][2]] }
                              res1 <- data.table(res1[1], DT1c[1])
                          } else {
                              res1 <- data.table(res1, DT1c)
                              res1[, nhcor:=ifelse(nh>1, nh/(nh-1), 1)]
                              res1[, num1:=nhcor * num * num]
                              if (!is.null(namesZ1) & !linratio) {
                                   res1[, num_den1:=nhcor * num * den]
                                   res1[, den1:=nhcor * den * den]
                               }}
                         namep <- c("nameY1", "nameZ1")
                         namep <- namep[namep %in% names(res1)]
                         varsp <- c("num1", "den1", "num_den1")
                         varsp <- varsp[varsp %in% names(res1)]

                        fits <- res1[, lapply(.SD, sum), 
                                       keyby=c("period_country",
                                               namesperc, namep),
                                       .SDcols=varsp]
                        return(fits)
                    })
            return(rbindlist(fitss))
        })
  res <- rbindlist(fits)
  DT1 <- fits <- DT1H <- NULL

  vars <- c(namesperc, namesDom, "nameY1", "nameZ1")
  vars <- names(res)[names(res) %in% vars]

  res <- merge(DTagg, res, by=vars)

  DTagg <- total <- NULL
	
  res[, estim:=totalY]
  res[, var:=num1]
  if (!is.null(res$totalZ) & !linratio) { 
                    res[, estim:=totalY/totalZ * percentratio]
                    res[, var:= (grad1 * grad1 * num1) +
                                (grad2 * grad2 * den1) +
                              2*(grad1 * grad2 * num_den1)] 
                    res[, var:=var * (percentratio)^2] }
    
  main <- c(namesperc, namesDom, "namesY", "nameY1")
  if (!is.null(namesDom)) main <- c(main, paste0(namesDom, "_new"))
  if (!is.null(res$namesZ)) main <- c(main, "namesZ", "nameZ1") 
  main2 <- c(main, "estim", "totalY", "valueY1")
  if (!is.null(res$namesZ)) main <- c(main, "totalZ")
  if (!is.null(res$namesZ)) main2 <- c(main2, "totalZ")
  if (!is.null(namesZ1) & !linratio) main2 <- c(main2, "valueZ1")
  main2 <- c(main2, "num1")
  if (!is.null(namesZ1) & !linratio) main2 <- c(main2, "den1", "grad1", "grad2")

  if (netchanges) { res1 <- res[, main2[!(main2 %in% c("nameY1",
                                     paste0(namesDom, "_new"), "nameZ1"))], with=FALSE]                  
                  } else res1 <- NULL

  main <- c(main, "totalY", "sample_size", "pop_size", "estim", "var")
  res <- res[, main, with=FALSE]

  #-------------------------------------------------------------------------*
  # DESIGN EFFECT (DEFF) ESTIMATION - VARIANCE UNDER SIMPLE RANDOM SAMPLING |
  #-------------------------------------------------------------------------*

  # We aggregate the target variables at household level


  DTs <- DT[, lapply(.SD, sum, na.rm=TRUE), 
                          keyby=c(namesperc2, names_id, "w_final"),
                         .SDcols = c(names_size1, size1wn, namesY2)]
  if (household_level_max) {    
           DTm <- DT[, lapply(.SD, max, na.rm=TRUE), keyby=c(namesperc2, names_id), .SDcols = names_size1]
       } else {
          DTm <- DT[, lapply(.SD, sum, na.rm=TRUE), keyby=c(namesperc2, names_id), .SDcols = names_size1]
     }

  setnames(DTm, names_size1, paste0(names_size1, "m"))
  DTs <- merge(DTs, DTm, by=c(namesperc2, names_id))

  # Linearised variables

  if (!is.null(namesZ1) & !linratio) {
                   lin1 <- lapply(split(DTs[, .I], DTs$period_country), function(i) 
                                lin.ratio(Y=DTs[i, namesY1, with=FALSE],
                                          Z=DTs[i, namesZ1, with=FALSE],
                                          weight=DTs[["w_final"]][i], Dom=NULL,
                                          percentratio=percentratio))
                   Y2a <- rbindlist(lin1)
                   setnames(Y2a, names(Y2a), paste0("lin___", namesY1))
                   DTs <- data.table(DTs, Y2a)
                   Y2a <- paste0("lin___", namesY1)
               } else Y2a <- namesY1

  w_final <- DTs[["w_final"]]
  DTsd <- DTs[, lapply(.SD[, Y2a, with = FALSE], function(x) 
                       sum(w_final*((x-sum(w_final*x)/sum(w_final))^2))/(sum(w_final)-1)),
                       keyby="period_country"]

  setnames(DTsd, Y2a, paste0("sd_w__", namesY1))
  DTs <- merge(DTs, DTsd, by="period_country")

  DTm <- DTs[, lapply(.SD[, paste0(names_size1, "m"), with = FALSE], function(x) sum(w_final*x, na.rm=TRUE)),
                       keyby="period_country"]
  setnames(DTm, paste0(names_size1, "m"), paste0("pop_", names_size1))
  DTs <- merge(DTs, DTm, by="period_country")

  DTsd <- DTs[, lapply(.SD, sd, na.rm=TRUE), keyby="period_country", .SDcols = Y2a]
  setnames(DTsd, Y2a, paste0("sd_nw__", namesY1))
  DTs <- merge(DTs, DTsd, by="period_country")

  DTm <- DTs[, lapply(.SD, sum, na.rm=TRUE), keyby="period_country", .SDcols = names_size1]
  setnames(DTm, names_size1, paste0("samp_", names_size1))
  DTs <- merge(DTs, DTm, by="period_country")
    
  DTx <- DTs[, .N, keyby=c(namesperc, paste0("sd_w__", namesY1),
                               paste0("sd_nw__", namesY1),
                               paste0("pop_", names_size1),
                               paste0("samp_", names_size1))]
  DTx[, N:=NULL]

  main <- melt(DTx[, c(namesperc, paste0("sd_w__", namesY1)), with=FALSE],  id=namesperc)
  main[, nameY1:=substr(variable, 7, nchar(trim(as.character(variable))))] 
  main[, variable:=NULL]
  setnames(main, "value", "sd_w")
  res <- merge(res, main, all.x=TRUE, by=c(namesperc, "nameY1"))

  main <- melt(DTx[, c(namesperc, paste0("sd_nw__", namesY1)), with=FALSE],  id=namesperc)
  main[, nameY1:=substr(variable, 8, nchar(trim(as.character(variable))))] 
  main[, variable:=NULL]
  setnames(main, "value", "sd_nw")
  res <- merge(res, main, all=TRUE, by=c(namesperc, "nameY1"))

  main <- melt(DTx[, c(namesperc, paste0("pop_", names_size1)), with=FALSE],  id=namesperc)
  if (!is.null(namesDom)){ 
                   main[, Dom:=substr(variable, 11, nchar(trim(as.character(variable))))] 
                   vars <- unique(main[["Dom"]])
                   vars <- data.table(Dom=vars, t(data.frame(strsplit(vars, "__"))))
                   setnames(vars, names(vars)[2:length(vars)], paste0(namesDom, "_new"))
                   main <- merge(main, vars, all.x=TRUE, by="Dom")  }
  main[, variable:=NULL]
  setnames(main, "value", "pop")
  nds <- namesperc
  if (!is.null(namesDom)) nds <- c(namesperc, paste0(namesDom, "_new"))
  res <- merge(res, main, all.x=TRUE, by=nds)

  main <- melt(DTx[, c(namesperc, paste0("samp_", names_size1)), with=FALSE],  id=namesperc)
  if (!is.null(namesDom)) main[, Dom:=substr(variable, 12, nchar(trim(as.character(variable))))] 
  main[, variable:=NULL]
  setnames(main, "value", "sampl_siz")
  if (is.null(namesDom)) nds <- namesperc else nds <- c(namesperc, "Dom")
  res <- merge(res, main, all=TRUE, by=nds)

  res[sample_size<pop_size, stderr_nw:=100*sqrt((1-(sample_size/pop_size))/pop_size * sd_nw * sd_nw/sample_size)]
  res[sample_size<pop_size, stderr_w:=100*sqrt((1-(sample_size/pop_size))/pop_size * sd_w * sd_w/sample_size)]

  DT <- DTw <- DTx <- DTs <- DTsd <- sd1 <- nds <- NULL

  res[, se:=sqrt(var)]
  res[, rse:=se/estim]
  res[, cv:=rse*100]
  tsad <- qnorm(0.5*(1+confidence))
  res[, absolute_margin_of_error:=tsad*se]
  res[, relative_margin_of_error:=tsad*cv]
  res[, CI_lower:=estim - tsad*se]
  res[, CI_upper:=estim + tsad*se]

  main <- namesperc 
  if (!is.null(namesDom))  main <- c(main, namesDom)
  main <- c(main, "namesY")
  if (!is.null(res$namesZ)) main <- c(main, "namesZ")

  main <- c(main, "sample_size", "pop_size", "estim", "se", 
            "var", "rse", "cv", "absolute_margin_of_error",
            "relative_margin_of_error", "CI_lower", "CI_upper", 
            "sd_w", "sd_nw", "pop", "sampl_siz", "stderr_nw",
            "stderr_w")
  main <- main[main %in% names(res)]
  res <- res[, main, with=FALSE]
  list(data_net_changes=DTnet, var_grad=res1, results=res)
}   



