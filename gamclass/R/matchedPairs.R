matchedPairs <-
function (years = 2007:2008, prefix = "fars", compareBYvar = c("airbagAvail",
                                                  "airbagDeploy", "Restraint"), bycat = list(airbagAvail = list(yes = c(1:9,
                                                                                                                20, 28:29, 31:32), no = 30, leaveout = c(0, 98, 99)), airbagDeploy = list(yes = c(1:9),
                                                                                                                                                                      no = c(20, 28, 30:32), leaveout = c(0, 29, 98, 99)), Restraint = list(yes = c(1:4,
                                                                                                                                                                                                                                            8, 10:12, 97), no = c(0, 5, 6, 7, 13:17), leaveout = (98:99))),
              restrict = "body%in%c(1:19,48,49,61,62)&!(mhevent%in%(2:5))",
              restrictvars = c("body", "mhevent", "seatpos", "injury"),
              retain = c("state", "age", "airbag", "injury", "restraint",
              "sex", "inimpact", "modelyr"), progress = TRUE)
{
    checkmiss <- c("injury", "age", "airbag", "restraint", "sex",
                   "any", "totnum", "nomatch", "dups")
    tabmiss <- array(0, c(51, length(checkmiss), length(years)))
    eachyear <- function(year, prefix, env = parent.frame()) {
        i <- year - min(years) + 1
        data <- get(paste(prefix, year, sep = ""))
        useINpassing <- c("casenum", "vnum", "ptype", "seatpos",
                          "sex", "body", "mhevent")
        usevar <- c(airbagAvail = "airbag", airbagDeploy = "airbag",
                    `deploy-0unk` = "airbag", `deploy-1unk` = "airbag",
                    Restraint = "restraint")
        notNA <- complete.cases(data[, restrictvars])
        condn <- with(data, notNA & seatpos %in% c(11, 13) &
                      (injury != 6)) & eval(parse(text = restrict), data)
        data <- data[condn, ]
        miss <- array(0, c(51, length(checkmiss)))
        rownames(miss) <- with(data, sort(unique(state)))
        colnames(miss) <- checkmiss
        data[, compareBYvar] <- data[, usevar[compareBYvar]]
        if (year <= 2008)
            data <- within(data, age[age %in% c(NA, 99)] <- 999)
        miss[, "injury"] <- table(data[,'state'], data[,'injury'] ==
                                             9)[, 2]
        miss[, "age"] <- table(data[,'state'], data[,'age'] == 999)[,
                                                 2]
        miss[, "airbag"] <- table(data[,'state'], data[,'airbag'] %in%
                                             c(NA, 98:99))[, 2]
        miss[, "restraint"] <- table(data[,'state'], data[,'restraint'] %in%
                                                c(NA, 98:99))[, 2]
        miss[, "sex"] <- table(data[,'state'], data[,'sex'] %in% c(NA,
                                                            9))[, 2]
        miss[, "totnum"] <- table(data[,'state'])
        data <- na.omit(data)
        idfields <- c("state", "casenum", "vnum")
        data[, 'caseid'] <- do.call("paste", c(data[, idfields], ":"))
        data[, 'anymiss'] <- (data[, 'airbag'] %in% c(98:99)) |
            (data[, 'sex'] %in% c(9)) |
             (data[, 'restraint'] %in% c(98:99)) |
             (data[, 'age'] %in% c(999))
        miss[, "any"] <- table(data[, 'state'], data[, 'anymiss'])[, 2]
        env[["tabmiss"]][, , i] <- miss
        levs <- factor(c("no", "yes", "NA-code"), levels = c("no",
                                                  "yes", "NA-code"))
        for (compareBY in compareBYvar) {
            var <- usevar[compareBY]
            toval <- val <- data[, var]
            bylev <- bycat[[compareBY]]
            toval[val %in% bylev[["no"]]] <- 1
            toval[val %in% bylev[["yes"]]] <- 2
            toval[val %in% bylev[["leaveout"]]] <- 3
            illegal <- val[!(toval %in% 1:3)]
            if (length(illegal) > 0)
                stop(paste("Illegal value(s) of", compareBY,
                           ":", paste(illegal, collapse = ",")))
            data[, compareBY] <- levs[toval]
        }
        ptype <- data[,'ptype']
        dfPass <- data[ptype==2,]
        dfDriver <- data[ptype==1,]
        dups <- dfPass$caseid[duplicated(dfPass$caseid)]
        indups1 <- dfPass$caseid %in% dups
        dups1 <- with(dfPass, table(state, indups1)[, 2])
        dfPass <- dfPass[!indups1, ]
        indups2 <- dfDriver$caseid %in% dups
        miss[, "dups"] <- dups1 + with(dfDriver, table(state,
                                                       indups2)[, 2])
        dfDriver <- dfDriver[!indups2, ]
        mapping <- match(dfPass[, "caseid"], dfDriver[, "caseid"],
                         nomatch = 0)
        nomatch1 <- with(dfPass, table(state, mapping == 0)[,
                                              2])
        dfPass <- dfPass[mapping > 0, c("caseid", retain, compareBYvar)]
        mapping <- mapping[mapping > 0]
        nomatch2 <- with(dfDriver, table(state, !(1:nrow(dfDriver)) %in%
                                         mapping)[, 2])
        miss[, "nomatch"] <- nomatch1 + nomatch2
        dfPass$D_injury <- dfDriver[mapping, "injury"]
        dvars <- paste("D_", compareBYvar, sep = "")
        dfPass[, dvars] <- dfDriver[mapping, compareBYvar]
        omit <- with(dfPass, injury == 9 | D_injury == 9)
        dfPass$year <- rep(year, nrow(dfPass))
        dfPass[!omit, ]
    }
    dflist <- vector("list", length(years))
    names(dflist) <- years
    i <- 0
    for (year in years) {
        if (progress)
            print(year)
        i <- i + 1
        dflist[[i]] <- eachyear(year, prefix)
    }
    df <- do.call("rbind", dflist)
    stateid <- sort(unique(df$state))
    dimnames(tabmiss) <- list(state = stateid, variables = checkmiss,
                              years = years)
    list(data = df, miss = tabmiss)
}
