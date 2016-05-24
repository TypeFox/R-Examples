## ----message=FALSE, warning=FALSE, echo=FALSE----------------------------
knitr::opts_chunk$set(message=FALSE, warning=FALSE)
library(magrittr)

## ----Coxph_tt_example, eval=FALSE----------------------------------------
#  library(survival)
#  coxph(Surv(time, event) ~ age + sex +
#          type_of_surgery + tt(type_of_surgery) +
#          tt(surgery_length),
#        data = my_surgical_data,
#        tt = list(
#          function(type_of_surgery, time, ...){
#            # As type_of_surgery is a categorical variable
#            # we must transform it into a design matrix that
#            # we can use for multiplication with time
#            # Note: the [,-1] is for dropping the intercept
#            mtrx <- model.matrix(~type_of_surgery)[,-1]
#            mtrx * time
#          },
#          function(surgery_length, time, ...){
#            # Note that the surgery_length and the time
#            # should probably have similar ranges. If you
#            # report operating time in minutes, follow-up
#            # in years the t will be dwarfed by the
#            pspline(surgery_length + time)
#          }
#        ))

## ----eval=FALSE----------------------------------------------------------
#  library(plyr)
#  library(magrittr)
#  models <-
#    timeSplitter(your_data,
#                 by = 4,
#                 event_var = "status",
#                 time_var = "years",
#                 event_start_status = "alive",
#                 time_related_vars = c("age", "date")) %>%
#    dlply("Start_time",
#          function(data){
#            coxph(Surv(Start_time, End_time, status) ~ age + sex + treatment, data = data)
#          })

## ------------------------------------------------------------------------
library(Greg)
test_data <- data.frame(
  id = 1:4,
  time = c(4, 3.5, 1, 5),
  event = c("censored", "dead", "alive", "dead"),
  age = c(62.2, 55.3, 73.7, 46.3),
  date = as.Date(
    c("2003-01-01", 
      "2010-04-01", 
      "2013-09-20",
      "2002-02-23"))
)

test_data$time_str <- sprintf("0 to %.1f", test_data$time)

## ----Display_data, echo=FALSE, fig.height=4, fig.width=7-----------------
library(grid)
library(magrittr)
getMaxWidth <- function(vars){
  sapply(vars, 
         USE.NAMES = FALSE,
         function(v){
           grobWidth(x = textGrob(v)) %>% 
             convertX(unitTo = "mm")
         }) %>% 
    max %>% 
    unit("mm")
}
plotTitleAndPushVPs <- function(title_txt){
  pushViewport(viewport(width = unit(.9, "npc"),
                        height = unit(.9, "npc")))
  
  title <- textGrob(title_txt, gp = gpar(cex = 2))
  title_height <- grobHeight(title) %>% 
    convertY(unitTo = "mm", valueOnly = TRUE) * 2 %>% 
    unit("mm")
  grid.layout(nrow = 3,
              heights = unit.c(title_height,
                               unit(.1, "npc"),
                               unit(1, "npc") - 
                                 title_height - 
                                 unit(.1, "npc") -
                                 unit(2, "line"),
                               unit(2, "line"))) %>% 
    viewport(layout = .) %>% 
    pushViewport()
  viewport(layout.pos.row = 1) %>% 
    pushViewport()
  grid.draw(title)
  upViewport()
  
  viewport(layout.pos.row = 3) %>% 
    pushViewport()
}

plotLine <- function (row_no,
                      start_time,
                      stop_time,
                      event,
                      data_range = c(0, max(test_data$time)),
                      print_axis = FALSE) {
  viewport(layout.pos.row = row_no,
           layout.pos.col = 6,
           xscale = data_range) %>% 
    pushViewport()
  on.exit(upViewport())
  
  if (event){
    grid.lines(x = unit(c(start_time, 
                          stop_time), "native"), 
               y = rep(0.5, 2))
    grid.points(x = unit(stop_time, "native"), y = 0.5, 
                pch = "*", 
                gp = gpar(cex = 2))
  }else{
    grid.lines(x = unit(c(start_time, 
                          stop_time), "native"), 
               y = rep(0.5, 2),
               arrow = arrow(length = unit(3, "mm"),
                             type = "closed"),
               gp = gpar(fill = "#000000"))
  }
  grid.points(x = unit(start_time, "native"), y = 0.5, pch = 20)
  if (print_axis)
    grid.xaxis()
}

plotIDcell <- function(row_no, id){
  viewport(layout.pos.row = row_no,
           layout.pos.col = 2) %>% 
    pushViewport()
  grid.text(id)
  upViewport()
}
plotTimeStrcell <- function(row_no, time_str){
  viewport(layout.pos.row = row_no,
           layout.pos.col = 4) %>% 
    pushViewport()
  grid.text(time_str)
  upViewport()
}

plotRowColor <- function(row_no, clr = "#F6F6FF"){
  viewport(layout.pos.row = row_no) %>% 
    pushViewport()
  grid.rect(gp = gpar(col = clr, fill = clr))
  upViewport()
}


# Do the actual plot
grid.newpage()
plotTitleAndPushVPs("Time spans")
widths <- 
  unit.c(unit(.1, "npc"),
         getMaxWidth(test_data$id),
         unit(.1, "npc"),
         getMaxWidth(test_data$time_str),
         unit(.1, "npc")) %>% 
  unit.c(., 
         unit(1, "npc") - sum(.) - unit(.1, "npc"),
         unit(.1, "npc"))

grid.layout(nrow = nrow(test_data),
            ncol = length(widths),
            widths = widths) %>% 
  viewport(layout = .) %>% 
  pushViewport()


for (i in 1:nrow(test_data)){
  if (i %% 2 == 0)
    plotRowColor(i)
  plotIDcell(i, test_data$id[i])
  plotTimeStrcell(i, test_data$time_str[i])

  plotLine(row_no = i, 
           start_time = 0,
           stop_time = test_data$time[i],
           event = test_data$event[i] == "dead",
           print_axis = i == nrow(test_data))
}
upViewport(2)

## ----Split_data----------------------------------------------------------
library(dplyr)
split_data <- 
  test_data %>% 
  select(id, event, time, age, date) %>% 
  timeSplitter(by = 2, # The time that we want to split by
               event_var = "event",
               time_var = "time",
               event_start_status = "alive",
               time_related_vars = c("age", "date"))

knitr::kable(head(split_data, 10))

## ----Complex_split_plot, fig.height=6, fig.width=7, echo=FALSE-----------
# Do the actual plot
plotTitleAndPushVPs("Time spans with split")

grid.layout(nrow = nrow(test_data) + nrow(split_data),
            ncol = length(widths),
            widths = widths) %>% 
  viewport(layout = .) %>% 
  pushViewport()

current_id <- NULL
no_ids <- 0
for (i in 1:nrow(split_data)){
  if (is.null(current_id) ||
      split_data$id[i] != current_id){
    current_id <- split_data$id[i]
    subjects_splits <- subset(split_data, id == current_id)
    rowspan <- (i + no_ids):(i + no_ids + nrow(subjects_splits))
    if (no_ids %% 2 == 1)
      plotRowColor(rowspan)
    plotIDcell(row_no = rowspan, id = current_id)
    plotTimeStrcell(row_no = rowspan, time_str = subset(test_data,
                                                        id == current_id,
                                                        "time_str"))
    with(subset(test_data,
                id == current_id),
         plotLine(row_no = i + no_ids, 
                  start_time = 0,
                  stop_time = time,
                  event = event == "dead"))
    no_ids = no_ids + 1
  }
  
  plotLine(row_no = i + no_ids, 
           start_time = split_data$Start_time[i],
           stop_time = split_data$Stop_time[i],
           event = split_data$event[i] == "dead",
           print_axis = i == nrow(split_data))
}
upViewport(2)


## ------------------------------------------------------------------------
# First we start with loading the dataset
data("melanoma", package = "boot")

# Then we munge it according to ?boot::melanoma
library(dplyr)
library(magrittr)
melanoma %<>% 
  mutate(status = factor(status,
                        levels = 1:3,
                        labels = c("Died from melanoma", 
                                   "Alive", 
                                   "Died from other causes")),
        ulcer = factor(ulcer,
                       levels = 0:1,
                       labels = c("Absent", "Present")),
        time = time/365.25, # All variables should be in the same time unit
        sex = factor(sex,
                     levels = 0:1,
                     labels = c("Female", "Male")))


## ------------------------------------------------------------------------
library(survival)
regular_model <- coxph(Surv(time, status == "Died from melanoma") ~
                         age + sex + year + thickness + ulcer,
                       data = melanoma,
                       x = TRUE, y = TRUE)
summary(regular_model)

## ------------------------------------------------------------------------
spl_melanoma <-
  melanoma %>% 
  timeSplitter(by = .5,
               event_var = "status",
               event_start_status = "Alive",
               time_var = "time",
               time_related_vars = c("age", "year"))

interval_model <-
  update(regular_model, 
         Surv(Start_time, Stop_time, status == "Died from melanoma") ~ .,
         data = spl_melanoma)

summary(interval_model)

## ------------------------------------------------------------------------
library(htmlTable)
cbind(Regular = coef(regular_model),
      Interval = coef(interval_model),
      Difference = coef(regular_model) - coef(interval_model)) %>% 
  txtRound(digits = 5) %>% 
  knitr::kable(align = "r")

## ------------------------------------------------------------------------
cox.zph(regular_model) %>% 
  extract2("table") %>% 
  txtRound(digits = 2) %>% 
  knitr::kable(align = "r")

## ------------------------------------------------------------------------
time_int_model <- 
  update(interval_model,
         .~.+thickness:Start_time)
summary(time_int_model)

## ------------------------------------------------------------------------
# First we need to manually add an interaction term
spl_melanoma %<>% 
  mutate(thickness_start = thickness * Start_time) 
anova(time_int_model,
      update(time_int_model, .~.+I(thickness_start^2)))

## ------------------------------------------------------------------------
update(time_int_model, .~.-thickness:Start_time+pspline(thickness_start))

## ------------------------------------------------------------------------
# Lets create an evenly distributed categorical thickness variable
# and interactions
spl_melanoma %<>% 
  mutate(thickness_cat = cut(thickness,
                             breaks = c(0, 1, 5, Inf),
                             labels = c("less than 1.0",
                                        "1.0 to 4.9",
                                        "at least 5.0")))
# Now create interaction variables
for (l in levels(spl_melanoma$thickness_cat)[-1]){
  spl_melanoma[[sprintf("thickness_%s_time", gsub(" ", "_", l))]] <- 
    (spl_melanoma$thickness_cat == l)*spl_melanoma$Start_time
}

# Now for the model specification where we use a 
# pspline for the two interaction variables
adv_int_model <-
  coxph(Surv(Start_time, Stop_time, status == "Died from melanoma") ~ 
          age + sex + year + ulcer + 
          thickness_cat + pspline(thickness_1.0_to_4.9_time) + pspline(thickness_at_least_5.0_time), 
        data = spl_melanoma, 
        x = TRUE, y = TRUE,
        iter.max = 1000)

# To get the estimates we use the predict function
new_data <- data.frame(thickness_cat = rep(levels(spl_melanoma$thickness_cat)[-1],
                                           each = 100),
                       Start_time = 2^seq(-3, 3, length.out = 100)) %>% 
  mutate(thickness_1.0_to_4.9_time = (thickness_cat == levels(spl_melanoma$thickness_cat)[2]) *
           Start_time,
         thickness_at_least_5.0_time = (thickness_cat == levels(spl_melanoma$thickness_cat)[3]) *
           Start_time)
new_data$sex = "Female"
new_data$age = median(melanoma$age)
new_data$year = median(melanoma$year)
new_data$ulcer = "Absent"

adv_pred <- predict(adv_int_model,
                    newdata = new_data,
                    type = "terms",
                    terms = c("thickness_cat", 
                              "pspline(thickness_1.0_to_4.9_time)",
                              "pspline(thickness_at_least_5.0_time)"),
                    se.fit = TRUE)

new_data$fit <- rowSums(adv_pred$fit)
new_data$se.fit <- apply(adv_pred$se.fit, 1, function(x) x^2) %>% 
  colSums %>% 
  sqrt
new_data %<>% 
  mutate(risk = exp(fit),
         upper = exp(fit + 1.96*se.fit),
         lower = exp(fit - 1.96*se.fit))


## ----fig.width=8, fig.height=6-------------------------------------------
library(ggplot2)
new_data %>% 
  mutate(adapted_risk = sapply(risk, function(x) min(max(2^-4, x), 2^5)),
         adapted_upper = sapply(upper, function(x) min(max(2^-4, x), 2^5)),
         adapted_lower = sapply(lower, function(x) min(max(2^-4, x), 2^5))) %>% 
  ggplot(aes(y = adapted_risk, 
             x = new_data$Start_time, 
             group = thickness_cat,
             col = thickness_cat,
             fill = thickness_cat)) + 
  # The confidence intervals are too wide to display in this case
  # geom_ribbon(aes(ymax = adapted_upper, ymin = adapted_lower), fill = "red") + 
  geom_line() +
  scale_x_log10(breaks = 2^(-3:4),
                limits = c(2^-3, 8),
                expand = c(0, 0)) + 
  scale_y_log10(breaks = 2^(-4:5), 
                labels = txtRound(2^(-4:5), 2),
                expand= c(0,0)) +
  scale_color_brewer(type = "qual", guide = guide_legend(title = "Thickness (mm)")) + 
  ylab("Hazard ratio") +
  xlab("Time (years)") +
  theme_bw()

