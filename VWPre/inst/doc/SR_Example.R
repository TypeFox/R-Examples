## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, warning=FALSE)

## ---- echo=FALSE, eval=TRUE, message=FALSE-------------------------------
library(VWPre)
data(VWdat)

## ---- eval= FALSE, echo=TRUE, results='asis'-----------------------------
#  library(VWPre)
#  VWdat <- read.table("1000HzData.txt", header = T, sep = "\t", na.strings = c(".", "NA"))

## ---- eval= FALSE, echo=TRUE, results='asis'-----------------------------
#  data(VWdat)

## ---- eval=TRUE, echo=TRUE, results='asis'-------------------------------
dat0 <- prep_data(data = VWdat, Subject = "RECORDING_SESSION_LABEL", Item = "itemid")

## ---- eval= FALSE, echo=TRUE, results='asis'-----------------------------
#  dat0 <- select(dat0, -starts_with("AVERAGE"), -starts_with("DATA_"),
#                 -starts_with("HTARGET"), -starts_with("IP"),
#                 -starts_with("LEFT_ACCELLERATION"), -starts_with("LEFT_GAZE"),
#                 -starts_with("LEFT_IN_"), -starts_with("LEFT_PUPIL"),
#                 -starts_with("LEFT_VELOCITY"), -starts_with("RESOLUTION"),
#                 -starts_with("RIGHT_ACCELLERATION"), -starts_with("RIGHT_GAZE"),
#                 -starts_with("RIGHT_IN_"), -starts_with("RIGHT_PUPIL"),
#                 -starts_with("RIGHT_VELOCITY"), -starts_with("SAMPLE"),
#                 -starts_with("TARGET_"), -starts_with("TRIAL_START"),
#                 -starts_with("VIDEO"))

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat1 <- relabel_na(data = dat0, NoIA = 4)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat2 <- create_time_series(data = dat1, Offset = 100)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
check_time_series(data = dat2)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
check_eye_recording(data = dat2)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat3 <- select_recorded_eye(data = dat2, Recording = "R", WhenLandR = "Right")

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
check_samplingrate(dat3)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
ds_options(SamplingRate = 1000)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat4 <- bin_prop(dat3, NoIA = 4, BinSize = 20, SamplingRate = 1000)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
check_samplingrate(dat4)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
check_samples_per_bin(dat4)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat5 <- transform_to_elogit(dat4, NoIA = 4, SamplesPerBin = 20)

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat5a <- create_binomial(data = dat4, NoIA = 4)

## ---- eval = TRUE, echo=TRUE, results='asis'-----------------------------
dat5b <- fasttrack(data = VWdat, Subject = "RECORDING_SESSION_LABEL", Item = "itemid", 
	EventColumns = c("Subject", "TRIAL_INDEX"), NoIA = 4, Offset = 100, Recording = "LandR", 
  WhenLandR = "Right", BinSize = 20, SamplingRate = 1000,
  SamplesPerBin = 20, Constant = 0.5, Output = "ELogit")

## ---- eval= TRUE, echo=TRUE, results='asis'------------------------------
dat6 <- rename_columns(dat5, Labels = c(IA1="Target", IA2="Rhyme", 
                                       IA3="OnsetComp", IA4="Distractor")) 

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"),
    Condition1 = NA, Condition2 = NA, Cond1Labels = NA, Cond2Labels = NA,
    ErrorBar = TRUE, VWPreTheme = TRUE) 

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"),
    Condition1 = NA, Condition2 = NA, Cond1Labels = NA, Cond2Labels = NA,
    ErrorBar = TRUE, VWPreTheme = TRUE) + ggtitle("Grand Average Plot")

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"),
    Condition1 = NA, Condition2 = NA, Cond1Labels = NA, Cond2Labels = NA,
    ErrorBar = TRUE, VWPreTheme = FALSE) + theme(axis.text = element_text(size = 15))

## ---- eval= TRUE, fig.show='hold', fig.height=5, results='asis', message=FALSE----
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"), Condition1 = "talker", 
    Condition2 = NA, Cond1Labels = c(CH1 = "Chinese 1", CH10 = "Chinese 3", 
                                     CH9 = "Chinese 2", EN3 = "English 1"),
    Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE)

## ---- eval= TRUE, fig.show='hold', fig.width=10, fig.height=5, results='asis', message=FALSE----
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"), Condition1 = NA, 
    Condition2 = "talker", Cond1Labels = NA, Cond2Labels = c(CH1 = "Chinese 1", 
                                                             CH10 = "Chinese 3", 
                                                             CH9 = "Chinese 2", 
                                                             EN3 = "English 1"), 
    ErrorBar = TRUE, VWPreTheme = TRUE)

## ---- eval= TRUE, fig.show='hold', fig.width=10, fig.height=8, results='asis', message=FALSE----
plot_avg(data = dat4, type = "proportion", xlim = c(0, 1000), 
    IAColumns = c(IA_1_P = "Target", IA_2_P = "Rhyme", IA_3_P = "OnsetComp", 
                  IA_4_P = "Distractor"), Condition1 = "talker", 
    Condition2 = "Exp", Cond1Labels = c(CH1 = "Chinese 1", CH10 = "Chinese 3", 
                                     CH9 = "Chinese 2", EN3 = "English 1"),
    Cond2Labels = c(High = "High Exp", Low = "Low Exp"), ErrorBar = TRUE, 
    VWPreTheme = TRUE)

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg_diff(data = dat4, xlim = c(0, 1000), DiffCols = c(IA_1_P = "Target", IA_2_P = "Rhyme"), 
            Condition1 = NA, Condition2 = NA, Cond1Labels = NA,
            Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE)

## ---- eval= TRUE, fig.show='hold', fig.height=5, results='asis', message=FALSE----
plot_avg_diff(data = dat4, xlim = c(0, 1000), DiffCols = c(IA_1_P = "Target", IA_2_P = "Rhyme"), 
            Condition1 = "talker", Condition2 = NA, Cond1Labels = c(CH1 = "Chinese 1", 
            CH10 = "Chinese 3", CH9 = "Chinese 2", EN3 = "English 1"),
            Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE)

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg_contour(data = dat4, IA = "IA_1_P", type = "proportion", Var = "Rating", 
VarLabel = "Accent Rating", xlim = c(0,1000), Theme = FALSE, 
Color = c("gray20", "gray90"))

## ---- eval= TRUE, fig.show='hold', results='asis', message=FALSE---------
plot_avg_contour(data = dat4, IA = "IA_1_P", type = "proportion", Var = "Rating", 
VarLabel = "Accent Rating", xlim = c(0,1000), Theme = FALSE, 
Color = c("red", "green")) + ggtitle("Looks to target")

## ---- eval=FALSE, echo=TRUE, results='asis'------------------------------
#  plot_var_app(dat4)

## ---- eval=FALSE, echo=TRUE, results='asis'------------------------------
#  plot_indiv_app(dat4)

## ---- eval=FALSE, echo=TRUE, results='asis'------------------------------
#  FinalDat <- dat5 %>%
#    # un-do any previous groupings
#    ungroup() %>%
#    # Select just the columns you want
#    select(., Subject, Item, Time, starts_with("IA"), Event, TRIAL_INDEX, Rating,
#           InteractChinese, Exp, target, rhymecomp, onsetcomp, distractor) %>%
#    # Order the data by Subject, Trial, and Time
#    arrange(., Subject, TRIAL_INDEX, Time)

## ---- eval=FALSE, echo=TRUE, results='asis'------------------------------
#  save(FinalDat, file = "FinalDat.rda", compress = "xz")

