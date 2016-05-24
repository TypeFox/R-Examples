vardchangespoor <- function(Y,
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
                     period1, period2,
                     use.estVar = FALSE,
                     confidence=0.95,
                     type="linrmpg",
                     change_type="absolute") {
 
  ### Checking
  if (!change_type %in% c("absolute", "relative")) stop("'change_type' must be 'absolute' or 'relative'")
 
  all_choices <- c("linarpr","linarpt","lingpg",
                   "linpoormed",  "linrmpg","lingini",
                   "lingini2","linqsr", "linrmir", "linarr")
  choices <- c("all_choices", all_choices)
  type <- tolower(type)

  type <- match.arg(type, choices, length(type)>1) 
  if (any(type == "all_choices"))  type <- all_choices

  # check 'p'
  p <- percentage
   if(length(p) != 1 |  any(!is.numeric(p) | p < 0 | p > 100)) {
          stop("'percentage' must be a numeric value in [0, 100]")  }

  # check 'order_quant'
  oq <- order_quant
  if(length(oq) != 1 | any(!is.numeric(oq) | oq < 0 | oq > 100)) {
          stop("'order_quant' must be a numeric value in [0, 100]")  }

  if(length(alpha) != 1 | any(!is.numeric(alpha) | alpha < 0 | alpha > 100)) {
         stop("'alpha' must be a numeric value in [0, 100]")  }

  if (!is.logical(use.estVar)) stop("'use.estVar' must be the logical value")

  if(length(confidence) != 1 | any(!is.numeric(confidence) | confidence < 0 | confidence > 1)) {
         stop("'confidence' must be a numeric value in [0,1]")  }

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
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }
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
          if (min(country %in% names(dataset))==1) country <- dataset[, country, with=FALSE] }

      if(!is.null(periods)) {
          if (min(periods %in% names(dataset))!=1) stop("periods' does not exist in 'dataset'!")
          if (min(periods %in% names(dataset))==1) periods <- dataset[, periods, with=FALSE] }

      if(!is.null(gender)) {
          if (min(gender %in% names(dataset))!=1) stop("'gender' does not exist in 'dataset'!")
          if (min(gender %in% names(dataset))==1) gender <- dataset[, gender, with=FALSE] }

      if(!is.null(sort)) {
          if (min(sort %in% names(dataset))!=1) stop("'sort' does not exist in 'dataset'!")
          if (min(sort %in% names(dataset))==1) sort <- dataset[, sort, with=FALSE] }
     
      if (!is.null(Dom)) {
          if (min(Dom %in% names(dataset))!=1) stop("'Dom' does not exist in 'dataset'!")
          if (min(Dom %in% names(dataset))==1) Dom <- dataset[, Dom, with=FALSE]   }
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
       Y_thres <- Y_thres[, 1]
       if (!is.numeric(Y_thres)) stop("'Y_thres' must be numerical")
       if (any(is.na(Y_thres))) stop("'Y_thres' has unknown values") 
     } else Y_thres <- Y

  # wght_thres
  if (is.null(wght_thres)) wght_thres <- w_final
  wght_thres <- data.frame(wght_thres)
  if (nrow(wght_thres) != n) stop("'wght_thres' must have the same length as 'Y'")
  if (ncol(wght_thres) != 1) stop("'wght_thres' must have vector or 1 column data.frame, matrix, data.table")
  wght_thres <- wght_thres[, 1]
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
  if (is.null(names(id))||(names(id)=="id")) setnames(id, names(id), "ID")

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
        sort <- sort[, 1]
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
  if (ncol(country) != 1) stop("'country' has more than 1 column")
  if (!is.character(country[[names(country)]])) stop("'country' must be character")
 
  # periods
  periods <- data.table(periods, check.names=TRUE)
  if (any(is.na(periods))) stop("'periods' has unknown values")
  if (nrow(periods) != n) stop("'periods' length must be equal with 'Y' row count")

  # Dom
  if (!is.null(Dom)) {
    Dom <- data.table(Dom)
    if (any(duplicated(names(Dom)))) 
           stop("'Dom' are duplicate column names: ", 
                 paste(names(Dom)[duplicated(names(Dom))], collapse = ","))
    if (nrow(Dom) != n) stop("'Dom' and 'Y' must be equal row count")
    if (any(is.na(Dom))) stop("'Dom' has unknown values")
    if (is.null(names(Dom))) stop("'Dom' must be colnames")
    Dom[, (names(Dom)):=lapply(.SD, as.character)] 
  }
  
  # period1
  period1 <- data.table(period1, check.names=TRUE)
  if (ncol(period1) != 1) stop("'period1' must be 1 column")
  setnames(period1, names(period1), names(periods))
  if (any(is.na(period1))) stop("'period1' has unknown values")
  periodss <- copy(periods)
  periodss[, periodss:=1]
  if (any(is.na(merge(period1, periodss, all.x=TRUE, by=names(periods)))))
              stop("'period1' row must be exist in 'periods'")

  # period2
  period2 <- data.table(period2, check.names=TRUE)
  if (ncol(period2) != 1) stop("'period2' must be 1 column")
  setnames(period2, names(period2), names(periods))
  if (any(is.na(period2))) stop("'period2' has unknown values")
  if (any(is.na(merge(period2, periodss, all.x=TRUE, by=names(periods)))))
              stop("'period2' row must be exist in 'periods'")


  data <- vardcrospoor(Y=Y, age=age, pl085=pl085,
                        month_at_work=month_at_work,
                        Y_den=Y_den, Y_thres=Y_thres,
                        H=H, PSU=PSU, w_final=w_final,
                        id=id, Dom=Dom, country=country,
                        periods=periods, sort=sort, 
                        gender = gender,
                        percentage=percentage,
                        order_quant=order_quant,
                        alpha = alpha,
                        dataset = NULL,
                        use.estVar = use.estVar,
                        withperiod = TRUE,
                        netchanges = TRUE,
                        confidence=confidence,
                        type=type)

  crossectional_results <- data$results

  PSU <- names(PSU)
  Dom <- names(Dom)
  H <- names(H)
  Y <- type
  country <- names(country)
  per <- names(periods)
  sar <- c(country, Dom, "type")
  sarp <- c(country, H, PSU)

  N <- namesY <- w_final <- ind <- dataset <- rot <- NULL
  rot_1 <- rot_2 <- stratasf <- name1 <- num1 <- num1num1 <-NULL
  num2num2 <- num1num2 <- num2 <- num1_1 <- V12 <- NULL
  grad1 <- grad2 <- num1_2 <- estim <- estim_1 <- NULL
  nh <- nhcor <- var_1 <- var_2 <- typs <- estim_2 <- se  <- NULL
  significant <- rho <- q_1 <- q_2 <- sum1 <- sum2 <- NULL

  var_grad <- copy(crossectional_results)
  var_grad <- var_grad[, c(per, sar, "estim", "var"), with=FALSE]
  period1[, ind:=.I]
  period2[, ind:=.I]
  per1 <- paste0(per, "_1")
  per2 <- paste0(per, "_2")  
  setnames(period1, per, per1)
  setnames(period2, per, per2)
  period1 <- merge(period1, period2, by="ind")
  period2 <- NULL

  var_grad1 <- merge(period1, var_grad, all.x=TRUE,
                              by.x=per1, by.y=per,
                              allow.cartesian=TRUE)
  var_grad2 <- merge(period1, var_grad, all.x=TRUE,
                              by.x=per2, by.y=per,
                              allow.cartesian=TRUE)
  setnames(var_grad1, c("estim", "var"), paste0(c("estim", "var"), "_1"))
  setnames(var_grad2, c("estim", "var"), paste0(c("estim", "var"), "_2"))

  var_grad <- merge(var_grad1, var_grad2, all=TRUE, by=c("ind", per1, per2, sar))
  var_grad1 <- var_grad2 <- NULL

  data <- data.table(data$data_net_changes, check.names=TRUE)
  data[, rot:=1]
  data1 <- merge(period1, data, all.x=TRUE,
                    by.x=per1, by.y=per,
                    allow.cartesian=TRUE)
  data2 <- merge(period1, data, all.x=TRUE, 
                    by.x=per2, by.y=per,
                    allow.cartesian=TRUE)  

  sard <- names(data)[!(names(data) %in% c(sarp, per))]

  setnames(data1, sard, paste0(sard, "_1"))
  setnames(data2, sard, paste0(sard, "_2"))
  data <- merge(data1, data2, all=TRUE, by=c("ind", per1, per2, sarp))
  data1 <- data2 <- NULL

  recode.NA <- function(DT, cols = seq_len(ncol(DT))) {
     for (j in cols) if (is.numeric(DT[[j]]))
      set(DT, which(is.na(DT[[j]])), j, ifelse(is.integer(DT[[j]]), 0L, 0))
   }
  recode.NA(data, c(paste0(sard,"_1"), paste0(sard,"_2")))

  dataH <- data[[H]]
  dataH <- factor(dataH)
  if (length(levels(dataH))==1) { data[, stratasf:= 1]
                                  dataH <- "stratasf"
                         } else { dataH <- data.table(model.matrix( ~ dataH-1))
                                  data <- cbind(data, dataH)
                                  dataH <- names(dataH) }


  fit <- lapply(1:(length(sard)-1), function(i) {
         fitd <- lapply(split(data, data[["ind"]]), function(data1) {

                 fits <- lapply(split(data1, data1[[country]]), function(DT3c) {

                           y1 <- paste0(sard[i], "_1")
                           y2 <- paste0(sard[i], "_2")

                           funkc <- as.formula(paste0("cbind(", trim(toString(y1)), ", ", 
                                                                trim(toString(y2)), ")~-1+",
                                                                paste(t(unlist(lapply(dataH, function(x) 
                                                                         paste0("rot_1:", toString(x), "+",
                                                                                "rot_2:", toString(x), "+",
                                                                                "rot_1:rot_2:", toString(x))))),
                                                                                collapse= "+"))) 
                           res <- lm(funkc, data=DT3c)
                           ssumas <- DT3c[, .(sum1=sum(get(y1)), sum2=sum(get(y2)))]

                           if (use.estVar) { res <- data.table(estVar(res))
                                        } else res <- data.table(lm(funkc, data=DT3c)$res)
                           setnames(res, names(res), c("num1", "num2"))
                           res[, namesY:=sard[i]]
                          
                           if (use.estVar) { 
                               res[, num1num1:=res[["num1"]][1]]
                               res[, num2num2:=res[["num2"]][2]]
                               res[, num1num2:=res[["num1"]][2]]
                               res <- data.table(res[1], DT3c[1])
                             } else {
                               res[, num1num1:=num1 * num1]
                               res[, num2num2:=num2 * num2]
                               res[, num1num2:=num1 * num2]
                               res <- data.table(res, DT3c)}

                           keynames <- c(country, "ind", paste0(per, "_1"), paste0(per, "_2"), "namesY")
                           fits <- res[, lapply(.SD, sum), keyby=keynames,
                                      .SDcols=c("num1num1", "num2num2", "num1num2")]
                           fits <- data.table(fits, ssumas)
                          return(fits)
                      })
               rbindlist(fits)
            })
            rbindlist(fitd)      
        })
   res <- rbindlist(fit)

   set(res, j=country, value=as.character(res[[country]]))

   var_grad[, namesY:=paste0("lin_", tolower(type))]
   if (!is.null(Dom)) {
          var_grad[, paste0(Dom, "_ss"):=lapply(Dom, function(x) make.names(paste0(x,".", get(x))))]
          var_grad[, paste0(Dom[1], "_ss"):=paste0(get(paste0(Dom[1], "_ss")))]
          var_grad[, typs:=paste0("lin_", tolower(type))]
          var_grad[, namesY:=Reduce(function(x, y)
                                      paste(x, y, sep = "__"), .SD),
                                     .SDcols=c("typs", paste0(Dom, "_ss"))]
         var_grad[, (c("typs", paste0(Dom, "_ss"))):=NULL]
       }

   setkeyv(res, c(country, "ind", per1, per2, "namesY"))
   setkeyv(var_grad, c(country, "ind", per1, per2, "namesY"))
   data <- merge(res, var_grad, all=TRUE)
   res <- fit <- var_gr <- NULL
   data[, namesY:=NULL]

   data[, rho:=num1num2/sqrt(num1num1*num2num2)]
   data[, V12:=num1num2*sqrt(var_1*var_2/(num1num1*num2num2))]

   if (change_type=="relative") {
        data[, q_1:=-sum2/sum1^2]
        data[, q_2:=1/sum1]
      } else {
         data[, q_1:=-1]
         data[, q_2:=1]
       }

   if (change_type=="relative") {
        data[, estim:=estim_2/estim_1]
     } else data[, estim:=estim_2 - estim_1]
   data[, var:=q_1*q_1*var_1 + 2*q_1*q_2*V12 + q_2*q_2*var_2]

   data[var>=0, se:=sqrt(var)]
  
   CI_lower <- CI_upper <- NULL
   tsad <- qnorm(0.5*(1+confidence))
   data[, CI_lower:=estim - tsad*se]
   data[, CI_upper:=estim + tsad*se]

   var_grad <- data[, c(country, per1, per2, Dom, 
                               "type", "q_1", "q_2",
                               "rho", "var_1", "var_2"), with=FALSE]

   changes_results <- data[, c(country, per1, per2, Dom, 
                               "type", "estim_1", "estim_2",
                               "estim", "var", "se",
                               "CI_lower", "CI_upper"), with=FALSE]
   data <- NULL

   changes_results[, significant:=TRUE]
   boundss <- as.numeric(change_type=="relative")
   changes_results[CI_lower<=boundss & CI_upper>=boundss, significant:=FALSE]
   changes_results[significant==TRUE]


  list(crossectional_results=crossectional_results, changes_results=changes_results)
}   


