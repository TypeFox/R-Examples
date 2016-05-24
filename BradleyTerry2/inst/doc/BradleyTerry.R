### R code from vignette source 'BradleyTerry.Rnw'

###################################################
### code chunk number 1: set_options
###################################################
options(prompt = "R> ", continue = "+  ", width = 70,
        useFancyQuotes = FALSE, digits = 7)


###################################################
### code chunk number 2: LoadBradleyTerry2
###################################################
library("BradleyTerry2")


###################################################
### code chunk number 3: CitationData
###################################################
data("citations", package = "BradleyTerry2")


###################################################
### code chunk number 4: CitationData2
###################################################
citations


###################################################
### code chunk number 5: countsToBinomial
###################################################
citations.sf <- countsToBinomial(citations)
names(citations.sf)[1:2] <- c("journal1", "journal2")
citations.sf


###################################################
### code chunk number 6: citeModel
###################################################
citeModel <- BTm(cbind(win1, win2), journal1, journal2, ~ journal,
    id = "journal", data = citations.sf)
citeModel


###################################################
### code chunk number 7: citeModelupdate
###################################################
update(citeModel, refcat = "JASA")


###################################################
### code chunk number 8: citeModelupdate2
###################################################
update(citeModel, br = TRUE)


###################################################
### code chunk number 9: lizModel
###################################################
options(show.signif.stars = FALSE)
data("flatlizards", package = "BradleyTerry2")
lizModel <- BTm(1, winner, loser, ~ SVL[..] + (1|..),
                data = flatlizards)


###################################################
### code chunk number 10: summarize_lizModel
###################################################
summary(lizModel)


###################################################
### code chunk number 11: lizModel2
###################################################
lizModel2 <- BTm(1, winner, loser,
            ~ throat.PC1[..] + throat.PC3[..] +
            head.length[..] + SVL[..] + (1|..),
            data = flatlizards)
summary(lizModel2)


###################################################
### code chunk number 12: baseball
###################################################
data("baseball", package = "BradleyTerry2")
head(baseball)


###################################################
### code chunk number 13: baseballModel
###################################################
baseballModel1 <- BTm(cbind(home.wins, away.wins), home.team, away.team,
                      data = baseball, id = "team")
summary(baseballModel1)


###################################################
### code chunk number 14: baseballDataUpdate
###################################################
baseball$home.team <- data.frame(team = baseball$home.team, at.home = 1)
baseball$away.team <- data.frame(team = baseball$away.team, at.home = 0)


###################################################
### code chunk number 15: baseballModelupdate
###################################################
baseballModel2 <- update(baseballModel1, formula = ~ team + at.home)
summary(baseballModel2)


###################################################
### code chunk number 16: CEMSmodel
###################################################
data("CEMS", package = "BradleyTerry2")
table8.model <-  BTm(outcome = cbind(win1.adj, win2.adj),
    player1 = school1, player2 = school2, formula = ~ .. +
    WOR[student] * LAT[..] +  DEG[student] * St.Gallen[..] +
    STUD[student] * Paris[..] + STUD[student] * St.Gallen[..] +
    ENG[student] * St.Gallen[..] + FRA[student] * London[..] +
    FRA[student] * Paris[..] + SPA[student] * Barcelona[..] +
    ITA[student] * London[..] + ITA[student] * Milano[..] +
    SEX[student] * Milano[..],
    refcat = "Stockholm", data = CEMS)


###################################################
### code chunk number 17: BTabilities
###################################################
BTabilities(baseballModel2)


###################################################
### code chunk number 18: BTabilities2
###################################################
head(BTabilities(lizModel2), 4)


###################################################
### code chunk number 19: residuals
###################################################
res.pearson <- round(residuals(lizModel2), 3)
head(cbind(flatlizards$contests, res.pearson), 4)


###################################################
### code chunk number 20: BTresiduals
###################################################
res <- residuals(lizModel2, type = "grouped")
#  with(flatlizards$predictors, plot(throat.PC2, res))
#  with(flatlizards$predictors, plot(head.width, res))


###################################################
### code chunk number 21: residualWLS
###################################################
lm(res ~ throat.PC1, weights = attr(res, "weights"),
   data = flatlizards$predictors)
lm(res ~ head.length, weights = attr(res, "weights"),
   data = flatlizards$predictors)


###################################################
### code chunk number 22: baseballModel2_call
###################################################
 baseballModel2$call


###################################################
### code chunk number 23: str_baseball
###################################################
str(baseball, vec.len = 2)


###################################################
### code chunk number 24: first_comparison
###################################################
baseball$home.team[1,]
baseball$away.team[1,]


###################################################
### code chunk number 25: first_outcome
###################################################
 baseball[1, c("home.wins", "away.wins")]


###################################################
### code chunk number 26: str_CEMS
###################################################
str(CEMS, vec.len = 2)


###################################################
### code chunk number 27: student-specific_data
###################################################
library("prefmod")
student <- cemspc[c("ENG", "SEX")]
student$ENG <- factor(student$ENG, levels = 1:2,
                      labels = c("good", "poor"))
student$SEX <- factor(student$SEX, levels = 1:2,
                      labels = c("female", "male"))


###################################################
### code chunk number 28: student_factor
###################################################
cems <- list(student = student)
student <- gl(303, 1, 303 * 15) #303 students, 15 comparisons
contest <- data.frame(student = student)


###################################################
### code chunk number 29: binomial_response
###################################################
win <- cemspc[, 1:15] == 0
lose <- cemspc[, 1:15] == 2
draw <- cemspc[, 1:15] == 1
contest$win.adj <- c(win + draw/2)
contest$lose.adj <- c(lose + draw/2)


###################################################
### code chunk number 30: school_factors
###################################################
lab <- c("London", "Paris", "Milano", "St. Gallen", "Barcelona",
         "Stockholm")
contest$school1 <- factor(sequence(1:5), levels = 1:6, labels = lab)
contest$school2 <- factor(rep(2:6, 1:5), levels = 1:6, labels = lab)


###################################################
### code chunk number 31: cems_data
###################################################
cems$contest <- contest


###################################################
### code chunk number 32: functions
###################################################
## cf. prompt
options(width = 55)
for (fn in getNamespaceExports("BradleyTerry2")) {
    name <- as.name(fn)
    args <- formals(fn)
    n <- length(args)
    arg.names <- arg.n <- names(args)
    arg.n[arg.n == "..."] <- "\\dots"
    is.missing.arg <- function(arg) typeof(arg) == "symbol" &&
        deparse(arg) == ""
    Call <- paste(name, "(", sep = "")
        for (i in seq_len(n)) {
        Call <- paste(Call, arg.names[i], if (!is.missing.arg(args[[i]]))
            paste(" = ", paste(deparse(args[[i]]),
                collapse = "\n"), sep = ""), sep = "")
        if (i != n)
            Call <- paste(Call, ", ", sep = "")
    }
    Call <- paste(Call, ")", sep = "")
    cat(deparse(parse(text = Call)[[1]], width.cutoff = 50), fill = TRUE)
}
options(width = 60)


