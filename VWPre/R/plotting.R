# Custom theme for ggplot2
theme_mybw <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
    )
}


#' Plots average looks to interest areas.
#' 
#' \code{plot_avg} calculates the grand or conditional averages of 
#' looks to each interest area along with standard error. It then plots the results.
#' N.B.: This function will work for data with a maximum of 8 interest areas
#' and two conditions.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' @import ggplot2
#' @import mgcv
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @param type A character string indicating "proportion" or "elogit".
#' @param xlim A vector of two integers specifying the limits of the x-axis.
#' @param IAColumns A named character vector specifying the desired interest 
#' area columns with custom strings for the legend.
#' @param Condition1 A string containing the column name corresponding to the 
#' first condition, if available. 
#' @param Condition2 A string containing the column name corresponding to the 
#' second condition, if available.
#' @param Cond1Labels A named character vector specifying the desired custom 
#' labels of the levels of the first condition. 
#' @param Cond2Labels A named character vector specifying the desired custom 
#' labels of the levels of the second condition. 
#' @param ErrorBar A logical indicating whether standard error bars should 
#' included in the plot.
#' @param VWPreTheme A logical indicating whether the theme included with the 
#' function should be applied, or ggplot2's base theme (to which any other 
#' custom theme could be added).
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For plotting the grand average with the included theme
#' plot_avg(data = dat, type = "elogit", xlim = c(0, 1000), 
#'    IAColumns = c(IA_1_ELogit = "Target", IA_2_ELogit = "Rhyme", 
#'    IA_3_ELogit = "OnsetComp", IA_4_ELogit = "Distractor"),
#'    Condition1 = NA, Condition2 = NA, Cond1Labels = NA, Cond2Labels = NA,
#'    ErrorBar = TRUE, VWPreTheme = TRUE) 
#'
#' # For plotting conditional averages (one condition) with the included theme.
#' # This produces plots arranged vertically.
#' plot_avg(data = dat, type = "elogit", xlim = c(0, 1000), 
#'    IAColumns = IAColumns = c(IA_1_ELogit = "Target", IA_2_ELogit = "Rhyme", 
#'    IA_3_ELogit = "OnsetComp", IA_4_ELogit = "Distractor"),
#'    Condition1 = "talker", Condition2 = NA, 
#'    Cond1Labels = c(CH1 = "Chinese 1", CH10 = "Chinese 3", CH9 = "Chinese 2", 
#'    EN3 = "English 1"), Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE)
#' 
#' # For plotting conditional averages (one condition) with the included theme.
#' # This produces plots arranged horizontally
#' plot_avg(data = dat, type = "elogit", xlim = c(0, 1000), 
#'    IAColumns = IAColumns = c(IA_1_ELogit = "Target", IA_2_ELogit = "Rhyme", 
#'    IA_3_ELogit = "OnsetComp", IA_4_ELogit = "Distractor"),
#'    Condition1 = NA, Condition2 = "talker", Cond1Labels = NA, 
#'    Cond2Labels = c(CH1 = "Chinese 1", CH10 = "Chinese 3", CH9 = "Chinese 2", 
#'    EN3 = "English 1"), ErrorBar = TRUE, VWPreTheme = TRUE)
#' 
#' # For plotting conditional averages (two conditions) with the included theme.
#' # This produces plots arranged in grid format.
#' plot_avg(data = dat, type = "elogit", xlim = c(0, 1000),
#'    IAColumns = IAColumns = c(IA_1_ELogit = "Target", IA_2_ELogit = "Rhyme", 
#'    IA_3_ELogit = "OnsetComp", IA_4_ELogit = "Distractor"),
#'    Condition1 = "talker", Condition2 = "Exp",
#'    Cond1Labels = c(CH1 = "Chinese 1", CH10 = "Chinese 3", CH9 = "Chinese 2", 
#'    EN3 = "English 1"), Cond2Labels = c(High = "H Exp", Low = "L Exp"),
#'    ErrorBar = TRUE, VWPreTheme = TRUE)
#' }
plot_avg <- function(data = data, type = NA, xlim = NA, IAColumns = IAColumns, 
                     Condition1 = NA, Condition2 = NA, Cond1Labels = NA, 
                     Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE) {
  
  NoIA <- length(names(IAColumns))
  IAColumns <- IAColumns
  Condition1 <- Condition1
  Condition2 <- Condition2
  Cond1Labels <- Cond1Labels
  Cond2Labels <- Cond2Labels
  ErrorBar <- ErrorBar
  Theme <- VWPreTheme
  
  if (is.na(xlim)[1]) {
    xlim <- c(range(data$Time)[1], range(data$Time)[2])
  } else {
    xlim <- xlim
  }
  
  if (type == "proportion") {
    ylabel = "Proportion Looks"
	ylim = c(0,1)
  } else if (type == "elogit") {
    ylabel = "Empirical Logit Looks"
	ylim = c(-4,4)
    } else {
	stop("You must specify either 'proportion' or 'elogit'.")
	}
  
  if (NoIA>8) {
    stop("You have more than 8 interest areas; you must modify this function.")
  } else if (NoIA == 1) {
    sel_names = c("Subject", "Time", names(IAColumns)[1])
    leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1]),
                          labels = c(IAColumns[[1]]))
      )
    } 
  } else if (NoIA == 2) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2]),
                          labels = c(IAColumns[[1]], IAColumns[[2]]))
      )
    } 
  } else if (NoIA == 3) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]]))
      )
    } 
  } else if (NoIA == 4) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], 
                  names(IAColumns)[3], names(IAColumns)[4])
    gath_col = "IA_"
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], 
                                     names(IAColumns)[3], names(IAColumns)[4]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]], 
                                     IAColumns[[4]]))
      )
    } 
  } else if (NoIA == 5) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                  names(IAColumns)[4], names(IAColumns)[5])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                                     names(IAColumns)[4], names(IAColumns)[5]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]], 
                                     IAColumns[[4]], IAColumns[[5]]))
      )
    } 
  } else if (NoIA == 6) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                  names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                                     names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]], 
                                     IAColumns[[4]], IAColumns[[5]], IAColumns[[6]]))
      )
    } 
  } else if (NoIA == 7) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                  names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6], names(IAColumns)[7])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                                     names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6], 
                                     names(IAColumns)[7]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]], 
                                     IAColumns[[4]], IAColumns[[5]], IAColumns[[6]], 
                                     IAColumns[[7]]))
      )
    } 
  } else if (NoIA == 8) {
    sel_names = c("Subject", "Time", names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                  names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6], 
                  names(IAColumns)[7], names(IAColumns)[8])
    my_leg <- function() {
      eval(
        scale_colour_grey(name = "Interest Area",
                          breaks = c(names(IAColumns)[1], names(IAColumns)[2], names(IAColumns)[3], 
                                     names(IAColumns)[4], names(IAColumns)[5], names(IAColumns)[6], 
                                     names(IAColumns)[7], names(IAColumns)[8]),
                          labels = c(IAColumns[[1]], IAColumns[[2]], IAColumns[[3]], 
                                     IAColumns[[4]], IAColumns[[5]], IAColumns[[6]], 
                                     IAColumns[[7]], IAColumns[[8]]))
      )
    } 
  }
  

  my_labeller <- labeller(
    CustCond1 = Cond1Labels,
    CustCond2 = Cond2Labels
  )
  
  
  if (!is.na(Condition1) & !is.na(Condition2)) {
  
  	GrandAvg <- data %>% select(match(sel_names,names(.)), get(Condition1), get(Condition2)) %>% 
      tidyr::gather_("IA", "VALUE", unique(names(IAColumns)), na.rm = FALSE, convert = FALSE) %>%
      group_by_("IA", "Time", Condition1, Condition2) %>%
      summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE))) %>%
    	rename_(CustCond1 = as.name(eval(Condition1)), CustCond2 = as.name(eval(Condition2)))
  
    if (type == "elogit") {
	ylim[1] = min(GrandAvg$mean) - 0.5
	ylim[2] = max(GrandAvg$mean) + 0.5
	}
	
    ggplot(GrandAvg, aes(x = Time, y = mean, colour = IA)) + 
      geom_point() +
      geom_line() +
      ylab(ylabel) +
      facet_grid(CustCond1 ~ CustCond2, labeller = my_labeller) +
      my_leg() +
      scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
      scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
        if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
      } + {
        if (Theme == TRUE) theme_mybw()
      }
    
  }
  
  else if (!is.na(Condition1) & is.na(Condition2)) {
    
    GrandAvg <- data %>% select(match(sel_names,names(.)), get(Condition1)) %>% 
      tidyr::gather_("IA", "VALUE", unique(names(IAColumns)), na.rm = FALSE, convert = FALSE) %>%
      group_by_("IA", "Time", Condition1) %>%
      summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE))) %>%
      rename_(CustCond1 = as.name(eval(Condition1)))
    
	if (type == "elogit") {
	ylim[1] = min(GrandAvg$mean) - 0.5
	ylim[2] = max(GrandAvg$mean) + 0.5
	}
	
      ggplot(GrandAvg, aes(x = Time, y = mean, colour = IA)) + 
        geom_point() +
        geom_line() +
        ylab(ylabel) +
        facet_grid(CustCond1 ~ ., labeller = my_labeller) +
        my_leg() +
        scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
        scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
          if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
        } + {
          if (Theme == TRUE) theme_mybw()
        }
    
  }
  
  else if (is.na(Condition1) & !is.na(Condition2)) {
    
    GrandAvg <- data %>% select(match(sel_names,names(.)), get(Condition2)) %>% 
      tidyr::gather_("IA", "VALUE", unique(names(IAColumns)), na.rm = FALSE, convert = FALSE) %>%
      group_by_("IA", "Time", Condition2) %>%
      summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE))) %>%
      rename_(CustCond2 = as.name(eval(Condition2)))
    

	if (type == "elogit") {
	ylim[1] = min(GrandAvg$mean) - 0.5
	ylim[2] = max(GrandAvg$mean) + 0.5
	}
	
      ggplot(GrandAvg, aes(x = Time, y = mean, colour = IA)) + 
        geom_point() +
        geom_line() +
        ylab(ylabel) +
        facet_grid(. ~ CustCond2, labeller = my_labeller) +
        my_leg() +
        scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
        scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
          if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
        } + {
          if (Theme == TRUE) theme_mybw()
        }
    
  }
  
  else if (is.na(Condition1) & is.na(Condition2)) {
    
    GrandAvg <- data %>% select(match(sel_names,names(.))) %>% 
      tidyr::gather_("IA", "VALUE", unique(names(IAColumns)), na.rm = FALSE, convert = FALSE) %>%
      group_by_("IA", "Time") %>%
      summarise(mean = mean(VALUE, na.rm = T), se = sd(VALUE) / sqrt(length(VALUE)))

	if (type == "elogit") {
	ylim[1] = min(GrandAvg$mean) - 0.5
	ylim[2] = max(GrandAvg$mean) + 0.5
	}
	
        ggplot(GrandAvg, aes(x = Time, y = mean, colour = IA)) + 
        geom_point() +
        geom_line() +
        ylab(ylabel) +
        my_leg() +
        scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
        scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
          if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
        } + {
          if (Theme == TRUE) theme_mybw()
        }
    
  }
  
}


#' Plots average contour surface of looks to a given interest area.
#' 
#' \code{plot_avg_contour} calculates the conditional average of empirical logit 
#' looks to a given interest area by Time and a specified continuous variable. 
#' It then applies a 3D smooth (derived using \code{\link[mgcv]{gam}}) over 
#' the surface and plots the results as a contour plot. 
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' @import ggplot2
#' @import mgcv
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @param IA A string specifying the column name of the IA to use. 
#' @param type A character string indicating "proportion" or "elogit".
#' @param xlim A vector of two integers specifying the limits of the x-axis.
#' @param Var A string containing the column name corresponding to the continuous
#' variable.
#' @param VarLabel A string specifying the axis label to use for \code{Var}.
#' @param Theme A logical indicating whether the theme included with the 
#' function, or ggplot2's base theme (which any other custom theme could be added).
#' @param Colors A vector of two strings specifying the colrs of the contour shading - The default values represent grayscale.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For plotting a conditional contour surface...
#' plot_avg_contour(data = dat, IA = "IA_1_ELogit", type = "elogit", 
#'                Var = "Rating", VarLabel = "Accent Rating", xlim = c(0,1000), 
#'                Theme = FALSE, Color = c("red", "white"))
#' }
plot_avg_contour <- function(data = data, IA = "IA_1_P", type = NA, Var = Var, 
                             VarLabel = VarLabel, xlim=NA, Theme=TRUE, 
                             Colors=c("gray20", "gray90")) {
  data <- data
  IA <- IA
  Var <- Var
  zlim <- c(-4,4)
  sel_names <- c("Time", IA, Var)
  Colors <- Colors
  
  if (is.na(xlim)[1]) {
    xlim <- c(range(data$Time)[1], range(data$Time)[2])
  } else {
    xlim <- xlim
  }
  
  Avg <- data %>% select(match(sel_names,names(.))) %>% 
    group_by_("Time", Var) %>%
    summarise_(mean = interp(~mean(IA, na.rm = T), IA = as.name(IA)))
  
  if (type == "proportion") {
    Avg$meanon <- round(Avg$mean*100)
    Avg$meanoff <- 100 - Avg$meanon
    Avg$meanbinom <- cbind(Avg$meanon, Avg$meanoff)
    var<-list(Var)
    model <- lapply(var, function(x) {
      mgcv::gam(substitute(meanbinom ~ te(Time, i), list(i = as.name(x))), data = Avg, family = binomial)
    })
  } else if (type == "elogit") {
    var<-list(Var)
    model <- lapply(var, function(x) {
      mgcv::gam(substitute(mean ~ te(Time, i), list(i = as.name(x))), data = Avg)
    })
  }
  
  model<-model[[1]]
  
  names(Avg)[names(Avg)==Var] <- "Variable"
  Varcol<-Avg$Variable
  Varcol<-unique(Varcol)
  
  # create a vector for variable A of length np
  data<-data[order(data$Event, data$Time),]
  rate <- 1000 / (data$Time[2] - data$Time[1])
  nptime<-round((xlim[2]-xlim[1])/100*(100/(1000/rate)))
  time<-seq(xlim[1],xlim[2],length.out=nptime) 
  # create a vector for variable B of length np
  npvar<-length(unique(Varcol))
  variable<-seq(range(Varcol)[1],range(Varcol)[2],length.out=npvar)
  tmp0<-expand.grid(time,variable)
  colnames(tmp0)<-c("Time", Var)
  newdat<-tmp0
  
  if (type == "proportion") {
    X1<-predict(model,newdat,type="lpmatrix")
  } else if (type == "elogit") {
    X1<-predict(model,newdat,type="lpmatrix")
  }
  
  # compute lpmatrix
  #X1<-predict(model,newdat,type="lpmatrix")
  mean<-X1%*%coef(model)
  predsmooth<-cbind(newdat, mean)
  names(predsmooth)[names(predsmooth)==Var] <- "Variable"
  
  if (type == "elogit") {
    zlim[1] = min(predsmooth$mean) - 0.25
    zlim[2] = max(predsmooth$mean) + 0.25
  } else if (type == "proportion") {
    predsmooth$mean <- plogis(predsmooth$mean)
    zlim[1] = 0
    zlim[2] = 1
  }
  
  if (Theme == TRUE) {
    ggplot(predsmooth) + 
      aes(x = Time, y = Variable, z = mean, fill = mean) + 
      geom_tile(aes(fill = mean)) +
      stat_contour() +
      geom_contour(color = "white") + 
      scale_fill_gradient(limit=c(zlim[1], zlim[2]), low = Colors[1], high = Colors[2]) + 
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      ylab(VarLabel) +
      theme_bw()
  } else {
    ggplot(predsmooth) + 
      aes(x = Time, y = Variable, z = mean, fill = mean) + 
      geom_tile(aes(fill = mean)) +
      stat_contour() +
      geom_contour(color = "white") + 
      scale_fill_gradient(limit=c(zlim[1], zlim[2]), low = Colors[1], high = Colors[2]) + 
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      ylab(VarLabel) 
  }
  
}


#' Plots average difference between looks to interest areas.
#' 
#' \code{plot_avg_diff} calculates the grand or conditional averages of 
#' differences between looks to two interest area along with standard error. 
#' It then plots the results.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' @import ggplot2
#' @import mgcv
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @param xlim A vector of two integers specifying the limits of the x-axis.
#' @param DiffCols A named character vector specifying the desired columns 
#' corresponding to the interest areas. 
#' @param Condition1 A string containing the column name corresponding to the 
#' first condition, if available. 
#' @param Condition2 A string containing the column name corresponding to the 
#' second condition, if available.
#' @param Cond1Labels A named character vector specifying the desired labels 
#' of the levels of the first condition. 
#' @param Cond2Labels A named character vector specifying the desired labels 
#' of the levels of the second condition. 
#' @param ErrorBar A logical indicating whether standard error bars should 
#' included in the plot.
#' @param VWPreTheme A logical indicating whether the theme included with the 
#' function should be applied, or ggplot2's base theme (which any other 
#' custom theme could be added).
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For plotting grand average differences...
#' plot_avg_diff(data = dat, xlim = c(0, 1000), DiffCols = c(IA_1_P = "Target", IA_2_P = "Rhyme"),
#'              Condition1 = NA, Condition2 = NA, Cond1Labels = NA, Cond2Labels = NA,
#'              ErrorBar = TRUE, VWPreTheme = TRUE)
#'              
#' # For plotting conditional average differences (one condition) with the included theme.
#' plot_avg_diff(data = dat, xlim = c(0, 1000), DiffCols = c(IA_1_P = "Target", IA_2_P = "Rhyme"), 
#'            Condition1 = "talker", Condition2 = NA, Cond1Labels = c(CH1 = "Chinese 1", 
#'            CH10 = "Chinese 3", CH9 = "Chinese 2", EN3 = "English 1"),
#'            Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE)
#'            
#' # For plotting conditional average differences (two conditions) with the included theme.
#' plot_avg_diff(data = dat, xlim = c(0, 1000), DiffCols = c(IA_1_P = "Target", IA_2_P = "Rhyme"), 
#'            Condition1 = "talker", Condition2 = "Exp", Cond1Labels = c(CH1 = "Chinese 1", 
#'            CH10 = "Chinese 3", CH9 = "Chinese 2", EN3 = "English 1"),
#'            Cond2Labels = c(High = "H Exp", Low = "L Exp"), ErrorBar = TRUE, 
#'            VWPreTheme = TRUE)
#' }
plot_avg_diff <- function(data = data, DiffCols = DiffCols, xlim = NA,
                          Condition1 = NA, Condition2 = NA, Cond1Labels = NA, 
                          Cond2Labels = NA, ErrorBar = TRUE, VWPreTheme = TRUE) {
  dat <- data
  
  DiffCols <- DiffCols
  DiffCol1 <- names(DiffCols)[1]
  DiffCol2 <- names(DiffCols)[2]
  ylim <- c(0,0)
  dat <- dat %>% mutate_(., Diff = interp(~DiffCol1 - DiffCol2, DiffCol1 = as.name(DiffCol1), DiffCol2 = as.name(DiffCol2)))
  Condition1 <- Condition1
  Condition2 <- Condition2
  Cond1Labels <- Cond1Labels
  Cond2Labels <- Cond2Labels
  ErrorBar <- ErrorBar
  Theme <- VWPreTheme
  sel_names <- c("Time", "Diff")
  
  if (is.na(xlim)[1]) {
    xlim <- c(range(dat$Time)[1], range(dat$Time)[2])
  } else {
    xlim <- xlim
  }

  if (is.na(Condition1) & is.na(Condition2)) {
    
    GrandAvg <- dat %>% select(match(sel_names,names(.))) %>% 
      group_by(Time) %>%
      summarise(mean = mean(Diff, na.rm = T), se = sd(Diff) / sqrt(length(Diff)))
    
    ylim[1] = min(GrandAvg$mean) - 0.15
    ylim[2] = max(GrandAvg$mean) + 0.15
    
    ggplot(GrandAvg, aes(x = Time, y = mean)) + 
      geom_point() +
      geom_line() +
      ylab("Difference") +
      geom_hline(yintercept=0) + 
      scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
      scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
        if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
      } + {
        if (Theme == TRUE) theme_mybw()
      } +
      ggtitle(paste("Difference between", DiffCols[[1]], "and", DiffCols[[2]], sep=" "))
    
  } else if (!is.na(Condition1) & is.na(Condition2)) {
    
    GrandAvg <- dat %>% select(match(sel_names,names(.)), get(Condition1)) %>% 
      group_by_("Time", Condition1) %>%
      summarise(mean = mean(Diff, na.rm = T), se = sd(Diff) / sqrt(length(Diff))) %>%
      mutate_(., Cond = interp(~Condition1, Condition1 = as.name(Condition1))) %>%
      mutate(., Cond = as.factor(Cond))
    
    lev1 <- unique(levels(GrandAvg$Cond))
    for (x in 1:length(names(Cond1Labels))) {
      for(i in 1:length(lev1)) {
        if (lev1[i] == names(Cond1Labels)[x]) {
          lev1[i] <- Cond1Labels[[x]]
        }
      }
    }
    levels(GrandAvg$Cond) <- lev1
    
    ylim[1] = min(GrandAvg$mean) - 0.15
    ylim[2] = max(GrandAvg$mean) + 0.15
    
    
    ggplot(GrandAvg, aes(x = Time, y = mean, colour = Cond)) + 
      geom_point() +
      geom_line() +
      ylab("Difference") +
      scale_colour_grey(name = "Condition") +
      geom_hline(yintercept=0) + 
      scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
      scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
        if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
      } + {
        if (Theme == TRUE) theme_mybw()
      } +
      ggtitle(paste("Difference between", DiffCols[[1]], "and", DiffCols[[2]], sep=" "))
    
  }  else if (!is.na(Condition1) & !is.na(Condition2)) {
    
    GrandAvg <- dat %>% select(match(sel_names,names(.)), get(Condition1), get(Condition2)) %>% 
      group_by_("Time", Condition1, Condition2) %>%
      summarise(mean = mean(Diff, na.rm = T), se = sd(Diff) / sqrt(length(Diff))) %>%
      mutate_(., Cond1 = interp(~Condition1, Condition1 = as.name(Condition1)),
              Cond2 = interp(~Condition2, Condition2 = as.name(Condition2))) %>%
      mutate(., Cond1 = as.factor(Cond1), Cond2 = as.factor(Cond2))
    
    
    lev1 <- unique(levels(GrandAvg$Cond1))
    for (x in 1:length(names(Cond1Labels))) {
      for(i in 1:length(lev1)) {
        if (lev1[i] == names(Cond1Labels)[x]) {
          lev1[i] <- Cond1Labels[[x]]
        }
      }
    }
    levels(GrandAvg$Cond1) <- lev1
    
    lev2 <- unique(levels(GrandAvg$Cond2))
    for (x in 1:length(names(Cond2Labels))) {
      for(i in 1:length(lev2)) {
        if (lev2[i] == names(Cond2Labels)[x]) {
          lev2[i] <- Cond2Labels[[x]]
        }
      }
    }
    levels(GrandAvg$Cond2) <- lev2
    
    GrandAvg$Cond <- interaction(GrandAvg$Cond1, GrandAvg$Cond2, sep="_")
    
    ylim[1] = min(GrandAvg$mean) - 0.15
    ylim[2] = max(GrandAvg$mean) + 0.15
    
    
    ggplot(GrandAvg, aes(x = Time, y = mean, colour = Cond)) + 
      geom_point() +
      geom_line() +
      ylab("Difference") +
      scale_colour_grey(name = "Conditions") +
      geom_hline(yintercept=0) + 
      scale_x_continuous(limits = c(xlim[1], xlim[2])) + 
      scale_y_continuous(limits = c(ylim[1], ylim[2])) + {
        if (ErrorBar == TRUE) geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .3)
      } + {
        if (Theme == TRUE) theme_mybw()
      } +
      ggtitle(paste("Difference between", DiffCols[[1]], "and", DiffCols[[2]], sep=" "))
  }  
  
}

