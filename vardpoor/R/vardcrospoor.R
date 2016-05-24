
vardcrospoor <- function(Y,  
                     age=NULL,
                     pl085=NULL,
                     month_at_work=NULL,
                     Y_den=NULL,
                     Y_thres = NULL, 
                     wght_thres = NULL, 
                     H, PSU, w_final, id,
                     Dom = NULL,
                     country, periods,
                     sort=NULL,
                     gender = NULL,
                     percentage=60,
                     order_quant=50,
                     alpha = 20,
                     dataset = NULL,
                     use.estVar = FALSE,
                     withperiod = TRUE,
                     netchanges = TRUE,
                     confidence = .95,
                     type="linrmpg") {
  ### Checking

  all_choices <- c("linarpr","linarpt","lingpg","linpoormed",
                   "linrmpg","lingini","lingini2","linqsr", "linrmir", "linarr")
  choices <- c("all_choices", all_choices)
  type <- tolower(type)

  type <- match.arg(type, choices, length(type)>1) 
  if (any(type == "all_choices"))  type <- all_choices

  # check 'p'
  p <- percentage
  if(length(p) != 1 | any(!is.numeric(p) | p < 0 | p > 100)) {
         stop("'percentage' must be a numeric value in [0, 100]")
     } else p <- percentage[1]

  # check 'order_quant'

  oq <- order_quant
   if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]")
      } else order_quant <- order_quant[1]

  if(length(alpha) != 1 | any(!is.numeric(alpha) | alpha < 0 | alpha > 100)) {
         stop("'alpha' must be a numeric value in [0, 100]")  }
 
  if (length(netchanges) != 1 | !any(is.logical(netchanges))) stop("'netchanges' must be the logical value")
  if (length(withperiod) != 1 | !any(is.logical(withperiod))) stop("'withperiod' must be the logical value")
  if (length(use.estVar) != 1 | !any(is.logical(use.estVar))) stop("'use.estVar' must be the logical value")

  if(length(confidence) != 1 | any(!is.numeric(confidence) | confidence < 0 | confidence > 1)) {
         stop("'confidence' must be a numeric value in [0, 1]")  }


  if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE] 

      if(!is.null(age)) {
          if (min(age %in% names(dataset))!=1) stop("'age' does not exist in 'dataset'!")
          if (min(age %in% names(dataset))==1) age <- dataset[, age, with=FALSE] }
      if(!is.null(pl085)) {
          if (min(pl085 %in% names(dataset))!=1) stop("'pl085' does not exist in 'dataset'!")
          if (min(pl085 %in% names(dataset))==1) pl085 <- dataset[, pl085, with=FALSE] }
      if(!is.null(month_at_work)) {
          if (min(month_at_work %in% names(dataset))!=1) stop("'month_at_work' does not exist in 'dataset'!")
          if (min(month_at_work %in% names(dataset))==1) month_at_work <- dataset[, month_at_work, with=FALSE] }
      if(!is.null(Y_den)) {
          if (min(Y_den %in% names(dataset))!=1) stop("'Y_den' does not exist in 'dataset'!")
          if (min(Y_den %in% names(dataset))==1) Y_den <- dataset[, Y_den, with=FALSE] }
      if(!is.null(Y_thres)) {
          if (min(Y_thres %in% names(dataset))!=1) stop("'Y_thres' does not exist in 'dataset'!")
          if (min(Y_thres %in% names(dataset))==1) Y_thres <- dataset[, Y_thres, with=FALSE] }    
      if(!is.null(wght_thres)) {
          if (min(wght_thres %in% names(dataset))!=1) stop("'wght_thres' does not exist in 'dataset'!")
          if (min(wght_thres %in% names(dataset))==1) wght_thres <- dataset[, wght_thres, with=FALSE] }
      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE]  }
      if(!is.null(id)) {
          if (min(id %in% names(dataset))!=1) stop("'id' does not exist in 'dataset'!")
          if (min(id %in% names(dataset))==1) id <- dataset[, id, with=FALSE] }
     if(!is.null(PSU)) {
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) PSU <- dataset[, PSU, with=FALSE] }
      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }
      if(!is.null(country)) {
          if (min(country %in% names(dataset))!=1) stop("'country' does not exist in 'dataset'!")
          if (min(country %in% names(dataset))==1) country <- dataset[, country, with=FALSE]  }

      if(!is.null(periods)) {
          if (min(periods %in% names(dataset))!=1) stop("'periods' does not exist in 'dataset'!")
          if (min(periods %in% names(dataset))==1) periods <- dataset[, periods, with=FALSE] }

      if(!is.null(gender)) {
          if (min(gender %in% names(dataset))!=1) stop("'gender' does not exist in 'dataset'!")
          if (min(gender %in% names(dataset))==1) gender <- dataset[, gender, with=FALSE] }

      if(!is.null(sort)) {
          if (min(sort %in% names(dataset))!=1) stop("'sort' does not exist in 'dataset'!")
          if (min(sort %in% names(dataset))==1) sort <- dataset[, sort, with=FALSE] }
     
      if (!is.null(Dom)) {
          if (min(Dom %in% names(dataset))!=1) stop("'Dom' does not exist in 'data'!")
          if (min(Dom %in% names(dataset))==1) Dom <- dataset[, Dom, with=FALSE] }
      }

  # Y
  Y <- data.frame(Y)
  n <- nrow(Y)
  if (ncol(Y) != 1) stop("'Y' must have vector or 1 column data.frame, matrix, data.table")
  Y <- Y[,1]
  if (!is.numeric(Y)) stop("'Y' must be numerical")
  if (any(is.na(Y))) stop("'Y' has unknown values")

   if (!is.null(Y_den)) {
          Y_den <- data.frame(Y_den)
          if (ncol(Y_den) != 1) stop("'Y_den' must be vector or 1 column data.frame, matrix, data.table")
          if (nrow(Y_den) != n) stop("'Y_den' must be the same length as 'Y'")
          Y_den <- Y_den[,1]
          if(!is.numeric(Y_den)) stop("'Y_den' must be numerical")
          if (any(is.na(Y_den))) stop("'Y_den' has unknown values")
  }

  # age
  if (!is.null(age)) {
       age <- data.frame(age)
       if (nrow(age) != n) stop("'age' must be the same length as 'Y'")
       if (ncol(age) != 1) stop("'age' must be vector or 1 column data.frame, matrix, data.table")
      age <- age[, 1]
      if (!is.numeric(age)) stop("'age' must be numerical")
      if (any(is.na(age))) stop("'age' has unknown values")
   }

   # pl085
   if (!is.null(pl085)) {
       pl085 <- data.frame(pl085)
       if (nrow(pl085) != n) stop("'pl085' must be the same length as 'Y'")
       if (ncol(pl085) != 1) stop("'pl085' must be vector or 1 column data.frame, matrix, data.table")
       pl085 <- pl085[, 1]
       if (!is.numeric(pl085)) stop("'pl085' must be numerical")
       if (any(is.na(pl085))) stop("'pl085' has unknown values")
   }

   # month_at_work
   if (!is.null(month_at_work)) {
        month_at_work <- data.frame(month_at_work)
        if (nrow(month_at_work) != n) stop("'month_at_work' must be the same length as 'Y'")
        if (ncol(month_at_work) != 1) stop("'month_at_work' must be vector or 1 column data.frame, matrix, data.table")
        month_at_work <- month_at_work[, 1]
        if (!is.numeric(pl085)) stop("'month_at_work' must be numerical")
        if (any(is.na(pl085))) stop("'month_at_work' has unknown values")
  }

  # Y_thres
  if (!is.null(Y_thres)) {
       Y_thres <- data.frame(Y_thres)
       if (nrow(Y_thres) != n) stop("'Y_thres' must have the same length as 'Y'")
       if (ncol(Y_thres) != 1) stop("'Y_thres' must have vector or 1 column data.frame, matrix, data.table")
       Y_thres <- Y_thres[,1]
       if (!is.numeric(Y_thres)) stop("'Y_thres' must be numerical")
       if (any(is.na(Y_thres))) stop("'Y_thres' has unknown values") 
     } else Y_thres <- Y

  # wght_thres
  if (is.null(wght_thres)) wght_thres <- w_final
  wght_thres <- data.frame(wght_thres)
  if (nrow(wght_thres) != n) stop("'wght_thres' must have the same length as 'Y'")
  if (ncol(wght_thres) != 1) stop("'wght_thres' must have vector or 1 column data.frame, matrix, data.table")
  wght_thres <- wght_thres[,1]
  if (!is.numeric(wght_thres)) stop("'wght_thres' must be a numeric vector")
 
  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' length must be equal with 'Y' row count")
  if (ncol(H) != 1) stop("'H' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colnames")
  H[, (names(H)):=lapply(.SD, as.character)]

  # id
  if (is.null(id)) id <- 1:n
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
  PSU[, (names(PSU)):=lapply(.SD, as.character)]
  
  # gender
  if (!is.null(gender)) {
      gender <- data.frame(gender)
      if (nrow(gender) != n) stop("'gender' must be the same length as 'Y'")
      if (ncol(gender) != 1) stop("'gender' must be vector or 1 column data.frame, matrix, data.table")
      gender <- gender[,1]
      if (!is.numeric(gender)) stop("'gender' must be numerical")
      if (length(unique(gender)) != 2) stop("'gender' must be exactly two values")
      if (!all.equal(unique(gender),c(1, 2))) stop("'gender' must be value 1 for male, 2 for females")
   }

  # sort
  if (!is.null(sort)) {
        sort <- data.frame(sort)
        if (length(sort) != n) stop("'sort' must have the same length as 'Y'")
        if (ncol(sort) != 1) stop("'sort' must be vector or 1 column data.frame, matrix, data.table")
        sort <- sort[,1]
   }

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
  country[, (names(country)):=lapply(.SD, as.character)]

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
    namesDom <- names(Dom)
    if (is.null(namesDom)) stop("'Dom' must be colnames")
    Dom[, (namesDom):=lapply(.SD, as.character)]
  }
    
  # Calculation
  Dom1 <- n_h <- stratasf <- name1 <- nhcor <- n_h <- var <- NULL
  num <- count_respondents <- value <- estim <- pop_size <- NULL
  period_country <- N <- se <- rse <- cv <- namesY <- H_sk <- NULL 

  estim <- c()
  countryper <- copy(country)
  if (!is.null(periods)) countryper <- data.table(periods, countryper)
  idper <- data.table(id, countryper)

  size <- copy(countryper)
  if (!is.null(namesDom)) size <- data.table(size, Dom)
  names_size <- names(size)
  size <- data.table(size, sk=1, w_final)
  size <- size[, .(count_respondents=.N,
                  pop_size=sum(w_final)), keyby=names_size]
 
  Y1 <- data.table(idper)
  Y1$period_country <- do.call("paste", c(as.list(Y1[,names(countryper),with=FALSE]), sep="_"))
  Y1 <- data.table(Y1, H, PSU, w_final, check.names=TRUE)
  namesY1 <- names(Y1)
  setkeyv(Y1, names(idper)) 

  if ("linarpt" %in% type) {
        varpt <- linarpt(Y=Y, id=id, weight=w_final,
                         sort=sort, Dom=Dom,
                         period=countryper,
                         dataset=NULL, percentage=percentage,
                         order_quant=order_quant, var_name="lin_arpt")
        Y1 <- merge(Y1, varpt$lin, all.x=TRUE)
        esti <- data.table("ARPT", varpt$value, NA)
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        varpt <- esti <- NULL
     }
  if ("linarpr" %in% type) {
        varpr <- linarpr(Y=Y, id=id, weight=w_final,
                         Y_thres=Y_thres,
                         wght_thres=wght_thres, sort=sort, 
                         Dom=Dom, period=countryper,
                         dataset=NULL, 
                         percentage=percentage,
                         order_quant=order_quant,
                         var_name="lin_arpr")
        Y1 <- merge(Y1, varpr$lin, all.x=TRUE)
        esti <- data.table("ARPR", varpr$value, NA)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        varpr <- esti <- NULL
      }
   if (("lingpg" %in% type)&&(!is.null(gender))) {
         vgpg <- lingpg(Y=Y, gender=gender, id=id,
                        weight=w_final, sort=sort,
                        Dom=Dom, period=countryper,
                        dataset=NULL, var_name="lin_gpg")
         Y1 <- merge(Y1, vgpg$lin, all.x=TRUE)
         esti <- data.table("GPG", vgpg$value, NA)  
         setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                    c("type", "value", "value_eu"))
         estim <- rbind(estim, esti)
         vgpg <- esti <- NULL
      }
   if ("linpoormed" %in% type) {
         vporm <- linpoormed(Y=Y, id=id, weight=w_final,
                             sort=sort, Dom=Dom, period=countryper, 
                             dataset=NULL, percentage=percentage,
                             order_quant=order_quant, var_name="lin_poormed")
         Y1 <- merge(Y1, vporm$lin, all.x=TRUE)
         esti <- data.table("POORMED", vporm$value, NA)  
         setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                    c("type", "value", "value_eu"))
         estim <- rbind(estim, esti)
         vporm <- esti <- NULL
      }
   if ("linrmpg" %in% type) {
         vrmpg <- linrmpg(Y=Y, id=id, weight=w_final,
                          sort=sort, Dom=Dom, period=countryper,
                          dataset=NULL, percentage=percentage,
                          order_quant=order_quant, var_name="lin_rmpg")
         Y1 <- merge(Y1, vrmpg$lin, all.x=TRUE)
         esti <- data.table("RMPG", vrmpg$value, NA)  
         setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                    c("type", "value", "value_eu")) 
         estim <- rbind(estim, esti)
         vrmpg <- esti <- NULL
      }
   if ("linqsr" %in% type) {
        vqsr <- linqsr(Y=Y, id=id, weight=w_final, 
                       sort=sort, Dom=Dom, period=countryper,
                       dataset=NULL, alpha=alpha, var_name="lin_qsr") 
        Y1 <- merge(Y1, vqsr$lin, all.x=TRUE)
        esti <- data.table("QSR", vqsr$value)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vqsr <- esti <- NULL
      }
   if ("lingini" %in% type) {
        vgini <- lingini(Y=Y, id=id, weight=w_final,
                         sort=sort, Dom=Dom, period=countryper,
                         dataset=NULL, var_name="lin_gini")
        Y1 <- merge(Y1, vgini$lin, all.x=TRUE)
        esti <- data.table("GINI", vgini$value)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vgini <- vginia <- esti <- NULL
      }
   if ("lingini2" %in% type) {
        vgini2 <- lingini2(Y=Y, id=id, weight=w_final,
                           sort=sort, Dom=Dom, period=countryper,
                           dataset=NULL, var_name="lin_gini2")
        Y1 <- merge(Y1, vgini2$lin, all.x=TRUE)
        esti <- data.table("GINI2", vgini2$value)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vgini2 <- esti <- NULL
      }
   if (("linrmir" %in% type)&&(!is.null(age))) {
        vrmir <- linrmir(Y=Y, id=id, age=age, weight=w_final, 
                       sort=sort, Dom=Dom, period=countryper,
                       dataset=NULL, order_quant=order_quant,
                       var_name="lin_rmir") 
        Y1 <- merge(Y1, vrmir$lin, all.x=TRUE)
 
        esti <- data.table("RMIR", vrmir$value, NA)  
        setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                   c("type", "value", "value_eu"))
        estim <- rbind(estim, esti)
        vrmir <-  esti <- NULL
      } 
   if (("linarr" %in% type)&&(!is.null(age))
                &&(!is.null(pl085))&&(!is.null(month_at_work))) {

       varr <- linarr(Y=Y, Y_den=Y_den, id=id, age=age, pl085=pl085, 
                             month_at_work=month_at_work, weight=w_final, 
                             sort=sort, Dom=Dom, period=countryper, dataset=NULL,
                             order_quant=order_quant,  var_name="lin_arr") 

       Y1 <- merge(Y1, varr$lin, all.x=TRUE)

       esti <- data.table("ARR", varr$value, NA)  
       setnames(esti, names(esti)[c(1, -1:0+ncol(esti))],
                                  c("type", "value", "value_eu"))
       estim <- rbind(estim, esti)
       varr <- esti <- NULL
     }


   setnames(estim, "value", "estim")
   estim$period_country <- do.call("paste", c(as.list(estim[,names(countryper),with=FALSE]), sep="_"))
   nams <- names(countryper)
   if (!is.null(namesDom)) nams <- c(nams, namesDom)
   estim <- merge(estim, size, all=TRUE, by=nams)

   namesY2 <- names(Y1)[!(names(Y1) %in% namesY1)]
   namesY2w <- paste0(namesY2, "w")
   Y1[, (namesY2w):=lapply(namesY2, function(x) get(x)*w_final)]

   DT1 <- copy(Y1)
   names_id <- names(id)
   names_H <- names(H)
   names_PSU <- names(PSU)

   namesperc <- c("period_country", names(countryper))
   namesDT1k <- c(namesperc, names_H, names_PSU)

   size <- id <- Dom <- country <-  NULL
   H <- PSU <- nh <- nh_cor <- NULL
  
   #--------------------------------------------------------*
   # AGGREGATION AT PSU LEVEL ("ULTIMATE CLUSTER" APPROACH) |
   #--------------------------------------------------------*

   DTY2 <- Y1[, lapply(.SD, sum, na.rm=TRUE), keyby=namesDT1k, .SDcols = namesY2w]
   setnames(DTY2, namesY2w, namesY2)
   DT1 <- copy(DTY2)
   DT1[, period_country:=NULL]
   if (!netchanges) DT1 <- NULL

   # NUMBER OF PSUs PER STRATUM
   setkeyv(DTY2, c(namesperc, names_H))
   DTY2[, nh:=.N, by=c(namesperc, names_H)]

   #--------------------------------------------------------------------------*
   # MULTIVARIATE REGRESSION APPROACH USING STRATUM DUMMIES AS REGRESSORS AND |
   # STANDARD ERROR ESTIMATION 						      |
   #--------------------------------------------------------------------------*

   DTY2H <- DTY2[[names_H]]
   DTY2H <- factor(DTY2H)
   if (length(levels(DTY2H))==1) { DTY2[, stratasf:=1]
                                   DTY2H <- "stratasf"
                          } else { DTY2H <- data.table(model.matrix( ~ DTY2H-1))
                                   DTY2 <- cbind(DTY2, DTY2H)
                                   DTY2H <- names(DTY2H) }
   namesY2m <-  make.names(namesY2)
   setnames(DTY2, namesY2, namesY2m)

   fits <-lapply(1:length(namesY2), function(i) {
             fitss <- lapply(split(DTY2, DTY2$period_country), function(DTY2c) {
                          y <- namesY2m[i]
                          funkc <- as.formula(paste("cbind(", trim(toString(y)), ")~",
                                         paste(c(-1, DTY2H), collapse= "+")))
                   	  res1 <- lm(funkc, data=DTY2c)
                            
           	          if (use.estVar==TRUE) {res1 <- data.table(crossprod(res1$res))
                                 } else res1 <- data.table(res1$res)
                          setnames(res1, names(res1)[1], "num") 
                          res1[, namesY:=y]
                           
                          if (use.estVar==TRUE) {
                                setnames(res1, "num", "var") 
                                res1 <- data.table(res1[1], DTY2c[1])
                            } else {
                                res1 <- data.table(res1, DTY2c)
                                res1[, nhcor:=ifelse(nh>1, nh/(nh-1), 1)]
                                res1[, var:=nhcor * num * num]
                              }
                          fits <- res1[, lapply(.SD, sum), 
                                         keyby=c(namesperc, "namesY"),
                                         .SDcols="var"]
                          return(fits)
                     })
            return(rbindlist(fitss))
      })
   res <- rbindlist(fits)
   
   estim[, namesY:=paste0("lin_", tolower(type))]
   if (!is.null(namesDom)) {
        Dom1 <- estim[, lapply(namesDom, function(x) make.names(paste0(x, ".", get(x))))]
        Dom1 <- Dom1[, Dom := Reduce(function(x, y) paste(x, y, sep="__"), .SD)]    
        estim <- data.table(estim, Dom1=Dom1[,Dom])
        estim[, namesY:=paste0(namesY, "__", Dom1)]
     }
  
   res <- merge(estim, res, all=TRUE, 
                 by=names(res)[!(names(res) %in% "var")])

   Dom1 <- estim <- DT3H <- NULL
   if (is.null(res$Dom1)) res[, Dom1:="1"]
   res[, (c("namesY", "Dom1", "period_country")):=NULL]

   res[, se:=sqrt(var)]
   res[, rse:=se/estim]
   res[, cv:=rse*100]
   
   res <- res[, c(names(countryper), namesDom, "type", "count_respondents",
                  "pop_size", "estim", "se", "var", "rse", "cv"), with=FALSE]

   list(data_net_changes=DT1, results=res)
 }   



