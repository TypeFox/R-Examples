
vardchangannual <- function(Y, H, PSU, w_final, id,
                        Dom = NULL, Z=NULL, 
                        country, years, subperiods,
                        dataset = NULL,
                        year1, year2,
                        percentratio = 1,
                        use.estVar = FALSE,
                        confidence=0.95) {
 
  ### Checking

  if (length(percentratio) != 1 | !any(is.integer(percentratio) | percentratio > 0)) stop("'percentratio' must be the positive integer value")
  if (length(use.estVar) != 1 | !any(is.logical(use.estVar))) stop("'use.estVar' must be the logical value")
  if(length(confidence) != 1 | any(!is.numeric(confidence) |  confidence < 0 | confidence > 1)) {
          stop("'confidence' must be a numeric value in [0,1]")  }

  if(!is.null(dataset)) {
      dataset <- data.table(dataset)
      if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
      if (min(Y %in% names(dataset))==1)  Y <- dataset[, Y, with=FALSE]

      if(!is.null(H)) {
          if (min(H %in% names(dataset))!=1) stop("'H' does not exist in 'dataset'!")
          if (min(H %in% names(dataset))==1) H <- dataset[, H, with=FALSE] }

      if(!is.null(id)) {
          if (min(id %in% names(dataset))!=1) stop("'id' does not exist in 'dataset'!")
          if (min(id %in% names(dataset))==1) id <- dataset[, id, with=FALSE]}

      if(!is.null(PSU)) {
          if (min(PSU %in% names(dataset))!=1) stop("'PSU' does not exist in 'dataset'!")
          if (min(PSU %in% names(dataset))==1) PSU <- dataset[, PSU, with=FALSE] }

      if(!is.null(w_final)) {
          if (min(w_final %in% names(dataset))!=1) stop("'w_final' does not exist in 'dataset'!")
          if (min(w_final %in% names(dataset))==1) w_final <- dataset[, w_final, with=FALSE] }

      if(!is.null(Z)) {
          if (min(Z %in% names(dataset))!=1) stop("'Z' does not exist in 'dataset'!")
          if (min(Z %in% names(dataset))==1) Z <- dataset[, Z, with=FALSE]}

      if(!is.null(country)) {
          if (min(country %in% names(dataset))!=1) stop("'country' does not exist in 'dataset'!")
          if (min(country %in% names(dataset))==1) country <- dataset[, country, with=FALSE] }

      if(!is.null(years)) {
          if (min(years %in% names(dataset))!=1) stop("years' does not exist in 'dataset'!")
          if (min(years %in% names(dataset))==1) years <- dataset[, years, with=FALSE] }

      if(!is.null(subperiods)) {
          if (min(subperiods %in% names(dataset))!=1) stop("subperiods' does not exist in 'dataset'!")
          if (min(subperiods %in% names(dataset))==1) subperiods <- dataset[, subperiods, with=FALSE] }
     
      if (!is.null(Dom)) {
          if (min(Dom %in% names(dataset))!=1) stop("'Dom1' does not exist in 'dataset'!")
          if (min(Dom %in% names(dataset))==1) Dom <- dataset[, Dom, with=FALSE]    }
   }

  # Y
  Y <- data.table(Y, check.names=TRUE)
  n <- nrow(Y)
  m <- ncol(Y)
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numerical")
  if (any(is.na(Y))) stop("'Y' has unknown values")
  if (is.null(names(Y))) stop("'Y' must be colnames")
  
  # H
  H <- data.table(H)
  if (nrow(H) != n) stop("'H' length must be equal with 'Y' row count")
  if (ncol(H) != 1) stop("'H' must be 1 column data.frame, matrix, data.table")
  if (any(is.na(H))) stop("'H' has unknown values")
  if (is.null(names(H))) stop("'H' must be colnames")
  
  # id
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
  if (ncol(country) != 1) stop("'country' must be 1 column")
  if (!is.character(country[[names(country)]])) stop("'country' must be character")

  # years
  years <- data.table(years, check.names=TRUE)
  if (any(is.na(years))) stop("'years' has unknown values")
  if (nrow(years) != n) stop("'years' length must be equal with 'Y' row count")
  if (ncol(years) != 1) stop("'years' must be 1 column")
  yearm <- names(years)

  # subperiods
  subperiods <- data.table(subperiods, check.names=TRUE)
  if (any(is.na(subperiods))) stop("'subperiods' has unknown values")
  if (nrow(subperiods) != n) stop("'subperiods' length must be equal with 'Y' row count")
  if (ncol(subperiods) != 1) stop("'subperiods' must be 1 column")
  subn <- data.table(years, subperiods)
  subn <- nrow(subn[,.N, by=names(subn)])/nrow(unique(years))
  subpm <- names(subperiods)

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
  
  namesZ <- NULL
  if (!is.null(Z)) {
    Z <- data.table(Z, check.names=TRUE)
    if (nrow(Z) != n) stop("'Z' and 'Y' must be equal row count")
    if (ncol(Z) != m) stop("'Z' and 'Y' must be equal column count")
    if (any(is.na(Z))) stop("'Z' has unknown values")
    if (is.null(names(Z))) stop("'Z' must be colnames")
    namesZ <- names(Z)
  }
 
   # year1
   year1 <- data.table(year1, check.names=TRUE)
   if (ncol(year1) != 1) stop("'year1' must be 1 column")
   setnames(year1, names(year1), names(years))
   if (any(is.na(year1))) stop("'year1' has unknown values")
   yearss <- copy(years)
   yearss[, yearss:=1]
   if (any(is.na(merge(year1, yearss, all.x=TRUE,
                       by=names(years), allow.cartesian=TRUE))))
                       stop("'year1' row must be exist in 'years'")

   # year2
   year2 <- data.table(year2, check.names=TRUE)
   if (ncol(year2) != 1) stop("'year2' must be 1 column")
   setnames(year2, names(year2), names(years))
   if (any(is.na(year2))) stop("'year2' has unknown values")
   if (any(is.na(merge(year2, yearss, all.x=TRUE,
                       by=names(years), allow.cartesian=TRUE))))
                       stop("'year2' row must be exist in 'years'")

   ids <- nams <- cros_se <- num1 <- totalY <- totalZ <- NULL
   estim_1 <- estim_2 <- avar <- N <- estim <- NULL
   var_est2 <- se  <- CI_lower <- CI_upper <- NULL

   pers <- data.table(years, subperiods, 
                      pers=paste0(years[[names(years)]], "__", subperiods[[names(subperiods)]]))
   sarak <- pers[,.N, keyby=names(pers)][, N:=NULL]
   
   namesDom <- names(Dom)

   apst <- lapply(1:nrow(year1), function(i) {

                 atsyear <- rbindlist(list(year1[i], year2[i]))
                 atsyear <- merge(atsyear, sarak, all.x=TRUE, by=yearm, sort = FALSE)
                 yr12 <- data.table(year1=year1[i][[1]], year2=year2[i][[1]])
                 setnames(yr12, paste0("year", c(1,2)), paste0(yearm, c(1,2)))
                 atsyrm <- names(atsyear)
                 atsyear[, ids:=.I]

                 nr1 <- nrow(atsyear)
                 yrs <- rbindlist(lapply(1:(nr1-1), function(j) { 
                           atsy1 <- atsyear[j]
                           atsy2 <- atsyear[(j+1):nr1]
                           setnames(atsy1, names(atsy1), paste0(names(atsy1), "_1"))
                           setnames(atsy2, names(atsy2), paste0(names(atsy2), "_2"))
                           data.table(atsy1, atsy2)                           
                         }))
                 yrs[, ids:=.I]

                 datas <- vardchanges(Y=Y, H=H, PSU=PSU, w_final=w_final,
                                      id=id, Dom=Dom, Z=Z, country=country,
                                      periods=pers[, "pers", with=FALSE], dataset=NULL,
                                      period1=yrs[["pers_1"]], period2=yrs[["pers_2"]],
                                      annual=TRUE, linratio=!is.null(Z),
                                      percentratio=percentratio,
                                      use.estVar = use.estVar,
                                      confidence=confidence,
                                      change_type="absolute")

                 crossectional_results <- datas$crossectional_results
                 crossectional_results <- merge(sarak, crossectional_results, all.y=TRUE, by="pers")

                 grad_var <- datas$grad_var
                 grad_var <- merge(yrs, grad_var, all.y=TRUE, by=c("pers_1", "pers_2"))

                 crossectional_var_grad <- datas$crossectional_var_grad
                 crossectional_var_grad <- merge(sarak, crossectional_var_grad, 
                                                 all.y=TRUE, by=c("pers"))

                 var_tau <- datas$var_tau
                 var_tau <- merge(yrs, var_tau, all.y=TRUE, by=c("pers_1", "pers_2"))
                 setkeyv(var_tau, "ids")

                 vardchanges_results <- datas$changes_results
                 vardchanges_results <- merge(yrs, vardchanges_results, all.y=TRUE, by=c("pers_1", "pers_2"))
                   
                 rho <- datas$rho
                 rho <- merge(yrs, rho, all.y=TRUE, by=c("pers_1", "pers_2"))
                 sar <- c("country", "namesY", "namesZ", namesDom)
                 sar <- sar[sar %in% names(rho)]
                 rhoj <- rho[,.N, keyby=sar][, N:=NULL]

                 apstr <- lapply(1:ncol(Y), function(j){                               

                               rho0 <- rhoj[j]
                               rho1 <- merge(rho0, rho, by=sar)[nams=="num2"]
                               A_matrix <- diag(1, nrow(atsyear), nrow(atsyear))

                               for (k in 1:nrow(rho1)) {

                                     at <- rho1[k==ids]
                                     A_matrix[at[["ids_1"]],at[["ids_2"]]] <- at[["rho_num1"]]
                                     A_matrix[at[["ids_2"]],at[["ids_1"]]] <- at[["rho_num1"]]
                                     if (at[["ids_2"]]> subn & at[["ids_1"]]< subn+1) {
                                                  A_matrix[at[["ids_1"]],at[["ids_2"]]] <- - 2 * at[["rho_num1"]]
                                                  A_matrix[at[["ids_2"]],at[["ids_1"]]] <- - 2 * at[["rho_num1"]]
                                           }
                                   }
                               rho1 <- merge(rho0, crossectional_var_grad, by=sar)
                               rho1 <- merge(atsyear[, c("pers", "ids"), with=FALSE], rho1, 
                                               by="pers", sort=FALSE, allow.cartesian=TRUE)
                               rho1[, cros_se:=sqrt(num1)]
                               X <- rho1[["cros_se"]]
                        
                               annual_var <- data.table(rho0, yr12, 1/(subn)^2 * (t(X)%*% A_matrix) %*% X)
                               setnames(annual_var, "V1", "var")
                         
                               A_matrix <- data.table(rho0, yr12, cols=paste0("V", 1:nrow(A_matrix)), A_matrix)
                               rho1[, ids:=paste0("V", ids)]
                               setnames(rho1, "ids", "cols")
                               rho1 <- data.table(yr12, rho1)

                               list(rho1, A_matrix, annual_var)})

                 rho1 <- rbindlist(lapply(apstr, function(x) x[[1]]))
                 A_matrix <- rbindlist(lapply(apstr, function(x) x[[2]]))
                 annual_var <- rbindlist(lapply(apstr, function(x) x[[3]]))

                 sars <- c(names(country), yearm, namesDom, "namesY", "namesZ")
                 sars <- sars[sars %in% names(crossectional_var_grad)]
                 sarsb <- sars[!(sars %in% yearm)]
                 sarc <- c("totalY", "totalZ")
                 sarc <- sarc[sarc %in% names(crossectional_var_grad)]
                 ysum <- crossectional_var_grad[,lapply(.SD, mean), by=sars, .SDcols=sarc]
                 if (!is.null(ysum$namesZ)) {
                              ysum[, estim:=totalY/totalZ * percentratio]
                       } else ysum[, estim:=totalY]
                 ysum1 <- ysum[get(yearm)==year1[i][[1]],c(sarsb, "estim"), with=FALSE]
                 ysum2 <- ysum[get(yearm)==year2[i][[1]],c(sarsb, "estim"), with=FALSE]
                 setnames(ysum1, "estim", "estim_1")
                 setnames(ysum2, "estim", "estim_2")
                 ysum1 <- data.table(yr12, merge(ysum1, ysum2, by=sarsb))
                 ysum1[, estim:=estim_2 - estim_1]
                 annual_changes <- merge(ysum1, annual_var, by=c(sarsb, names(yr12)))

                 list(crossectional_results,
                      crossectional_var_grad, grad_var,
                      rho, var_tau, vardchanges_results,
                      rho1, A_matrix, annual_changes, ysum)             
   })

  crossectional_results <- rbindlist(lapply(apst, function(x) x[[1]]))
  crossectional_var_grad <- rbindlist(lapply(apst, function(x) x[[2]]))
  grad_var <- rbindlist(lapply(apst, function(x) x[[3]]))
  rho <- rbindlist(lapply(apst, function(x) x[[4]]))
  var_tau <- rbindlist(lapply(apst, function(x) x[[5]]))
  vardchanges_results <- rbindlist(lapply(apst, function(x) x[[6]]))

  X_annual <- rbindlist(lapply(apst, function(x) x[[7]]))
  A_matrix <- rbindlist(lapply(apst, function(x) x[[8]]))
  annual_changes <- rbindlist(lapply(apst, function(x) x[[9]]))
  ysum <- rbindlist(lapply(apst, function(x) x[[10]]))

  crossectional_results[, pers:=NULL]
  crossectional_var_grad[, pers:=NULL]
  grad_var[, (c("pers_1", "pers_2", "ids_1", "ids_2", "ids")):=NULL]
  rho[, (c("pers_1", "pers_2", "ids_1", "ids_2", "ids")):=NULL]
  var_tau[, (c("pers_1", "pers_2", "ids_1", "ids_2", "ids")):=NULL]
  vardchanges_results[, (c("pers_1", "pers_2",
                           "ids_1", "ids_2", "ids")):=NULL]

  vars <- c(paste0(yearm, c(1,2)), yearm, names(country),
            namesDom, "namesY", "namesZ", "cols", "cros_se")
  X_annual <- X_annual[, vars[vars %in% names(X_annual)], with=FALSE]

  vars <- c(paste0(yearm, c(1,2)), names(country), namesDom, 
            "namesY", "namesZ", "cols", paste0("V", 1:8))
  A_matrix <- A_matrix[, vars[vars %in% names(A_matrix)], with=FALSE] 

  vars <- c(names(country), yearm, namesDom, "namesY", 
            "namesZ", "totalY", "totalZ", "estim")
  ysum <- ysum[, vars[vars %in% names(ysum)], with=FALSE] 

  vars <- c(paste0(yearm, c(1,2)), names(country), namesDom, "namesY", 
             "namesZ", paste0("estim_", c(1,2)), "estim", "var")       
  annual_changes <- annual_changes[, vars[vars %in% names(annual_changes)], with=FALSE] 

  annual_changes[, var_est2:=var]  
  annual_changes[xor(is.na(var_est2), var_est2 < 0), var_est2:=NA]  
  annual_changes[, se:=sqrt(var_est2)]
  annual_changes[, var_est2:=NULL]

  tsad <- qnorm(0.5*(1+confidence))
  annual_changes[, CI_lower:= estim - tsad * se]
  annual_changes[, CI_upper:= estim + tsad * se]

  significant <- NULL
  annual_changes[, significant:=TRUE]
  annual_changes[CI_lower<=0 & CI_upper>=0, significant:=FALSE]

  list(crossectional_results=crossectional_results,
       crossectional_var_grad=crossectional_var_grad,
       vardchanges_grad_var=grad_var,
       vardchanges_rho=rho,
       vardchanges_var_tau=var_tau,
       vardchanges_results=vardchanges_results,
       X_annual=X_annual, A_matrix=A_matrix,
       annual_sum=ysum,  
       annual_changes=annual_changes)

}
