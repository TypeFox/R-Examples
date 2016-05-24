#' @name fishers.exact
#' @title Fisher's Exact Test
#' @param targetcolumn  Character string, name of response column to be tested.
#' @param alternative  Character string(s) specifying the direction of the alternative hypothesis.
#' Must be one or more of "greater", "two.sided", or "less".
#' @param alpha  Significance level (numeric) to be used.
#' @param control  Level of dose to be used as the control for dichotomous data.
#' @param tot.obs  Total number of observations for each dose level for dichotomous data.
#' @param label  Output label.
#' @param data  Input dataframe.
#' @keywords internal

fishers.exact <- function (targetcolumn, alternative, alpha, control, tot.obs, label, data) {
  
  dose_fac <- data$dose_fac
  level_list <- levels(dose_fac)
  
  if (control == "") {control <- min(level_list)}

    for (c in 1:length(level_list)) {
      
      if (level_list[c] != control) {
        
        Events <- data[data$dose_fac == level_list[c],]
        Sum_events <- sum(Events[get('targetcolumn')])
        Sum_noevents <- sum(Events[get('tot.obs')])
      
        Control_events <- data[data$dose_fac == get('control'),]
        Control_sum_events <- sum(Control_events[get('targetcolumn')])
        Control_sum_noevents <- sum(Control_events[get('tot.obs')])
      
        Matrix <-matrix(c(Sum_events, Control_sum_events, Sum_noevents, Control_sum_noevents),
                        nrow = 2, dimnames = list(c(paste("Dose", level_list[c], sep = " "), paste("Control", control, sep = " ")), 
                                                  c(targetcolumn, "Total Obs")))
        conf.level <- 1-alpha
        Label <- paste("Tested Contrast,", label)
        writeLines(Label)
        print(Matrix)
        print(stats::fisher.test(Matrix, alternative = alternative, conf.level = conf.level))
      } else {next}
    }
  }

