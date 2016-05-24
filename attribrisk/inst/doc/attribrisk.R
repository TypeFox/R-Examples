## ----echo=FALSE, include=FALSE, cache=FALSE--------------------
library(knitr)
options(width=65)
require(attribrisk)
# Make output look like output, input like input
opts_chunk$set(fig.with=7, fig.height=5.5, fig.path="figures/",
               comment=NA, prompt=TRUE, tidy=FALSE, cache=FALSE,
               warning=FALSE, error=FALSE,
               out.width="\\textwidth", out.height="!")

## ----results="as.is", echo=FALSE-------------------------------
temp <- with(chapter.dat, table(cases, hbp)) 
temp2 <- colSums(temp)
tpct <- round(temp[2,] / colSums(temp) ,2)
ar <- round(diff(tpct)/tpct[2], 2)

## ----results="as.is", echo=FALSE-------------------------------
cat("Controls &", temp[1,1], " &", temp[1,2], "\\\\ \n")
cat("Stroke   &", temp[2,1], " &", temp[2,2], "\\\\ \\hline \n")
cat("Total    &", temp2[1] , " &", temp2[2] , "\n")

## --------------------------------------------------------------
require(attribrisk)
data(chapter.dat)

#Show first and last row.
chapter.dat[c(1,2644),]

#Summarize the relationship between hbp and case/control status.
count <- table(chapter.dat$hbp, chapter.dat$cases)
count

## ----cache=TRUE, eval=TRUE-------------------------------------
example1 <- attribrisk(cases ~ expos(hbp), data=chapter.dat)

example1

## ----wt1-------------------------------------------------------
tdata <- data.frame(case=c(0, 1, 0, 1),
                    hbp =c(0, 0, 1, 1),
                    count = c(559, 384, 763, 938))
example1b <- attribrisk(case ~ expos(hbp), data=tdata, weight=count)
example1b

## ----boot------------------------------------------------------
example1boot <- attribrisk(cases ~ expos(hbp), data=chapter.dat, 
   varmethod = "boot")
example1boot

## ----cache=T, tidy=FALSE, eval=TRUE----------------------------
chapter.dat[chapter.dat$match.id==1,]

example2 <- attribrisk(cases ~ strata(match.id) + expos(hbp), 
   data=chapter.dat)

example2

## ----cache=T, tidy=FALSE, eval=TRUE----------------------------
# Build Targe
data(stroke.dat)

stroke.target <- data.frame(smoke = "Never",
                            diastolic = .9*stroke.dat$diastolic)

set.seed(21790)
example4a <- attribrisk(
   cases ~ age + expos(smoke) + expos(diastolic), 
   data=stroke.dat, varmethod="boot", baseline = stroke.target)

example4a

## ----cache=T, tidy=FALSE, eval=TRUE----------------------------
# Build baseline 
target <- cut(stroke.dat$diastolic, c(0, 85, 100, 120, 150, 500))
reduce <- c(0, .05, .1, .15, .25)[as.numeric(target)]
newbp <-with(stroke.dat, diastolic *(1-reduce))
newsm <- with(stroke.dat, ifelse(smoke=="Current", "Former", smoke))

stroke.target2 <- data.frame(diastolic = newbp,
                             smoke = newsm)
example4b <- attribrisk(
   cases ~ age + expos(smoke) + expos(diastolic), 
   data=stroke.dat, baseline = stroke.target2)

print(example4b, digits=3)

## ----bgfit, tidy=FALSE, echo=FALSE, results='asis'-------------
data(benichou)

bg_model <- 1:18
bg_formula <- c('$\\alpha + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot Ag$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S \\cdot Ag$',
                '$\\alpha + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot Ag$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S \\cdot Ag$',
                '$\\alpha + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot Ag$',
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S$', 
                '$\\alpha_1 Ag + \\alpha_2 S + \\alpha_3 Ag \\cdot S + \\beta Al + \\gamma Al \\cdot S \\cdot Ag$',
                '$\\alpha + \\beta X$',
                '$\\alpha Ag + \\beta X$',
                '$\\alpha Ag + \\beta X + \\gamma Ag X$')

bg_ar <- c(0.395, 0.382, 0.380, 0.381, 0.380,
           0.709, 0.719, 0.723, 0.703, 0.700,
           0.709, 0.721, 0.726, 0.703, 0.701,
           0.862, 0.866, 0.868)

bg_sd <- c(0.042, 0.044, 0.044, 0.044, 0.044,
           0.051, 0.050, 0.050, 0.054, 0.056,
           0.051, 0.050, 0.050, 0.054, NA,
           0.046, 0.045, 0.044)

ar_formula <- c(cases ~ expos(alcohol80),
                cases ~ age * smoke + expos(alcohol80),
                cases ~ age * (smoke + expos(alcohol80)),
                cases ~ smoke * (age + expos(alcohol80)),
                cases ~ age * smoke * expos(alcohol80),

                cases ~ expos(alcohol40),
                cases ~ age * smoke + expos(alcohol40),
                cases ~ age * (smoke + expos(alcohol40)),
                cases ~ smoke * (age + expos(alcohol40)),
                cases ~ age * smoke * expos(alcohol40),

                cases ~ expos(alcohol),
                cases ~ age * smoke + expos(alcohol),
                cases ~ age * (smoke + expos(alcohol)),
                cases ~ smoke * (age + expos(alcohol)),
                cases ~ age * smoke * expos(alcohol),

                cases ~ expos(fsmoke.alc),
                cases ~ age + expos(fsmoke.alc),
                cases ~ age * expos(fsmoke.alc))

ar_jk <- rep(0, 18)
sd_jk <- rep(0, 18)

benichou$fsmoke.alc <-  factor(benichou$smoke.alc,
                               levels=c('Unexposed', 'Exposed'))

for (i in 1:length(ar_formula)){
  t <- attribrisk(ar_formula[[i]], data=benichou, 
                  varmethod = "jackknife")

  ar_jk[i] <- t$attribrisk
  sd_jk[i] <- sqrt(t$var)
}

bg_df <- data.frame(model=bg_model, 
                    #formula = bg_formula,
                    formula = unlist(lapply(ar_formula, function(x) 
                                     deparse(x[[3]]))), 
                    ar = format(round(bg_ar,2)), 
                    sd = format(round(bg_sd,3)), 
                    ar_jk = format(round(ar_jk,2)),
                    sd_jk= format(round(sd_jk,3)),
                    stringsAsFactors = FALSE)


#bg_df <- bg_df[,-2] #TMT addition: take out the formula


for (i in 1:nrow(bg_df)) {
    cat (paste(bg_df[i,], collapse=" & "))
    if (i == nrow(bg_df)) cat("\n") else cat("\\\\ \n")
}

