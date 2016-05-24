# Part of the migui package for multiple imputation of missing data
# Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

GUI_ENV <- new.env()

# Invoke mi GUI
migui <- function(tk = "RGtk2") {
  options("guiToolkit" = tk)
  initialize();
  nextScreen();
}

# TOP-LEVEL FUNCTIONS
# =========================================

# initialize data and GUI state (once per GUI invocation)
initialize <- function() {

  GUI_ENV$window <- gwindow("mi", visible = FALSE);
  GUI_ENV$state <- "get_data";
  GUI_ENV$group <- ggroup(container=GUI_ENV$window);
  GUI_ENV$n.chains = 4;
  GUI_ENV$n.iter = 30;
  GUI_ENV$mdf_args <- list();
  GUI_ENV$mdf <- NULL
  GUI_ENV$user_df <- NULL
  nextScreen();
  addSpring(GUI_ENV$group);
  visible(GUI_ENV$window) <- TRUE;
  visible(GUI_ENV$group) <- TRUE;
  return(invisible(NULL))
}

# top-level loop to show next screen
nextScreen <- function() {
  delete(GUI_ENV$window,GUI_ENV$group);
  #   visible(GUI_ENV$window) <- FALSE;
  GUI_ENV$group <- ggroup(container=GUI_ENV$window,
                          horizontal=FALSE,
                          expand=TRUE, 
                          fill="both");
  #   addSpring(GUI_ENV$group);
  functionName <- paste("do_", GUI_ENV$state, sep = "");
  do.call(functionName, args = list());
  #   visible(GUI_ENV$window) <- TRUE;
  return(invisible(NULL))
}

# SCREENS
# ========================================

# STEP 1
do_get_data <- function() {
  innerGroup <- setupScreen(step = "1", 
                            title = "Get Data",
                            nextHandler = function(h,...) {    
                              selected = svalue(dataSourceRadio);
                              if (selected == dataSourceOptions["disk"]) {
                                GUI_ENV$state = "get_data_from_disk";
                              } else if (selected == dataSourceOptions["r"]) {
                                GUI_ENV$state = "get_data_from_r";
                              } else if (selected == dataSourceOptions["mi"]) {
                                GUI_ENV$state = "get_data_from_mi";
                              }
                              nextScreen(); 
                            },
                            backHandler = NULL );
  choiceLabel("Select location of data", container=innerGroup);
  dataSourceOptions <- c(disk="disk drive [csv, dta, por, sav, RData]",
                         r="R's workspace",
                         mi="mi package example [built in]");
  dataSourceRadio <- gradio(dataSourceOptions, selected=1,
                            container = indent(innerGroup));
}

# STEP 1.5.a
do_get_data_from_mi <- function() {
  innerGroup <- setupScreen(step = "1.5",
                            title = "Get Data from mi Package",
                            backHandler = navHandler("get_data"),
                            nextHandler = function(h,...) {
                              do.call("data", args = list(svalue(selector), 
                                                          package = "mi"))
                              GUI_ENV$user_df <- get(svalue(selector), 
                                                     envir= .GlobalEnv);
                              GUI_ENV$state <- "describe_data";
                              nextScreen();
                            });
  
  choiceLabel("Select an example", container = innerGroup);
  data.frames <- data(package = "mi")$results[,3];
  selector <- gcombobox(data.frames, container=indent(innerGroup),
                        editable = TRUE, expand=FALSE, anchor=c(-1,-1));
}

# STEP 1.5.b (pops up own dialog for file using gfile(), then calls next screen)
do_get_data_from_disk <- function() {
  extensions <- paste("*", c("csv", "dta", "por", "sav", "RData", "Rdata"), sep = ".");
  file_name <- gfile(text="Choosing Data File for mi GUI", type="open", 
                     filter = list("All files" = list(patterns = extensions)));
  if (is.na(file_name)) {
    GUI_ENV$state = "get_data";
  } else {
    if (grepl("\\.csv$", file_name)) {
      user_df <- read.csv(file_name);
      GUI_ENV$user_df <- user_df;
      GUI_ENV$state = "describe_data";
    } else if (grepl("\\.dta$", file_name)) {
      user_df <- foreign::read.dta(file_name); # FIXME: untested
      GUI_ENV$user_df <- user_df;
      GUI_ENV$state = "describe_data";
    } else if (grepl("\\.\\(por,sav\\)$", file_name)) {
      user_df <- foreign::read.spss(file_name, to.data.frame = TRUE);
      GUI_ENV$user_df <- user_df;
      GUI_ENV$state = "describe_data";
    } else {
      dfs <- load(file_name, envir = GUI_ENV);
      GUI_ENV$user_df <- GUI_ENV[[dfs[1]]]
      GUI_ENV$state = "describe_data";
    }
  }
  nextScreen();
  return(invisible(NULL));
}

# STEP 1.5.c
do_get_data_from_r <- function() {
  data.frames <- sapply(ls(envir = .GlobalEnv), FUN = function(x) {
    is.data.frame(get(x, envir = .GlobalEnv))
  });

  if(length(data.frames) == 0) {
    GUI_ENV$state <- "get_data";
    gmessage("No data.frames found in R's global environment, starting over", 
             icon = "error", parent = GUI_ENV$group);
    nextScreen();
    return(invisible(NULL))
  }
  data.frames <- data.frames[data.frames];
  innerGroup <- setupScreen(step = "1.5",
                            title = "Get Data Frame from R's Environment",
                            nextHandler = function(h,...) {
                              user_df <- get(svalue(selector), envir= .GlobalEnv);
                              if(is(user_df, "missing_data.frame")) {
                                GUI_ENV$user_df <- as.data.frame(user_df);
                                GUI_ENV$mdf <- user_df;
                                GUI_ENV$state <- "var_options";
                              }
                              else {
                                GUI_ENV$user_df <- user_df
                                GUI_ENV$state <- "describe_data";
                              }
                              nextScreen();
                            },
                            backHandler = navHandler("get_data"));
  
  choiceLabel("Select a data frame in R from the list", container = innerGroup);
  selector <- gcombobox(names(data.frames), container=indent(innerGroup),
                        editable = TRUE);
}

# STEP 2
do_describe_data <- function() {
  innerGroup <- setupScreen(step = "2",
                            title = "Describe the Structure of Your Data",
                            backHandler = navHandler("get_data"),
                            nextHandler = function(h,...) {
                              subclass <- svalue(data_desc);
                              if(!any(grepl(paste("^", subclass, "_", sep = ""), 
                                            mi::.prune("missing_data.frame")))) subclass <- NULL
                              GUI_ENV$mdf_args$subclass <- subclass;
                              GUI_ENV$state <- "get_var_classes";
                              nextScreen();
                            });
  
  choiceLabel("Select the best description of your data", container = innerGroup);
  options <- mi::.prune("missing_data.frame")
  # FIXME: Update this every time you define a new class from missing_data.frame
  human_options <- sapply(strsplit(options, "_"), FUN = function(x) x[1]);
  human_options[1] <- "no special structure";
  names(human_options) <- options;
  
  data_desc <- gradio(human_options, container=indent(innerGroup));
}

# STEP 3
do_get_var_classes <- function() {
  innerGroup <- setupScreen(step = "3",
                            title = "Set Class of Variables",
                            backHandler = navHandler("describe_data"),
                            nextHandler = navHandler("var_options"));
  
  label_group <- ggroup(horizontal=TRUE, container=innerGroup);
  glabel("<b>Select default parameters for variable types:</b>", 
                    markup=TRUE, container= label_group);
  addSpring(label_group);
  
  defaults <- formals(.guess_type);
  threshold_group <- ggroup(horizontal=TRUE, container = innerGroup);
  contLabel <- glabel("Continuity Threshold", container=threshold_group);
  thresholdButton <- gspinbutton(from=3, to=100, by=1, value=defaults$threshold, 
                                            container=threshold_group, handler = function(h, ...) {
                                   GUI_ENV$mdf <- NULL;
                                   types <- sapply(colnames(GUI_ENV$user_df), FUN = function(y) {
                                     .guess_type(GUI_ENV$user_df[,y],
                                                   favor_ordered = svalue(favor_ord),
                                                   favor_positive = svalue(favor_pos),
                                                   threshold = svalue(thresholdButton),
                                                   variable_name = y)
                                   })
                                   names(types) <- colnames(GUI_ENV$user_df);
                                   GUI_ENV$mdf_args$types <- types;
                                   for(i in seq_along(types)) {
                                     tmp <- selectionListGroup[i+1,2]
                                     svalue(tmp) <- types[i];
                                   }
                                 });
#   tooltip(contLabel) <- tooltip(thresholdButton) <- 
#     "Specified value indicates minimum number of unique values required to favor continuous over categorical."
  addSpace(threshold_group,15);
  favor_ord <- gcheckbox("favor ordered", checked=defaults$favor_ordered, container=threshold_group,
                         handler = function(h, ...) {
                           GUI_ENV$mdf <- NULL;
                           types <- sapply(colnames(GUI_ENV$user_df), FUN = function(y){
                             .guess_type(GUI_ENV$user_df[,y],
                                           favor_ordered = svalue(favor_ord),
                                           favor_positive = svalue(favor_pos),
                                           threshold = svalue(thresholdButton),
                                           variable_name = y)
                           })
                           names(types) <- colnames(GUI_ENV$user_df);
                           GUI_ENV$mdf_args$types <- types;
                           categoricals <- c("unordered-categorical", "ordered-categorical")
                           for(i in seq_along(types)) if(types[i] %in% categoricals) {
                             tmp <- selectionListGroup[i+1,2]
                             svalue(tmp) <- types[i];
                           }
                         });
  # tooltip(favor_ord) <- "When checked, favor ordered-categorical over unordered-categorical for initial values."
  addSpace(threshold_group,15);
  favor_pos <- gcheckbox("favor positive", checked=defaults$favor_positive, container=threshold_group,
                         handler = function(h, ...) {
                           GUI_ENV$mdf <- NULL;
                           types <- sapply(colnames(GUI_ENV$user_df), FUN = function(y) {
                             .guess_type(GUI_ENV$user_df[,y],
                                           favor_ordered = svalue(favor_ord),
                                           favor_positive = svalue(favor_pos),
                                           threshold = svalue(thresholdButton),
                                           variable_name = y)
                           })
                           names(types) <- colnames(GUI_ENV$user_df);
                           GUI_ENV$mdf_args$types <- types;
                           cons <- mi::.prune("continuous")
                           for(i in seq_along(types)) if(types[i] %in% cons) {
                             tmp <- selectionListGroup[i+1,2]
                             svalue(tmp) <- types[i];
                           }
                         });
  # tooltip(favor_pos) <- "When checked, favor positive-continuous to continuous for initial values."
  
  gseparator(horizontal=TRUE,container=innerGroup);
  
  addSpring(innerGroup);
  
  glabel("<b>Variable Class Selector:</b>", markup = TRUE, container=innerGroup);
  
  # FIXME: this sucks on Linux with tcltk; also better to do use.scrollwindow = TRUE
  #   size(innerGroup) <- c(500,500);
  selectionGroup <- ggroup(container=innerGroup, use.scrollwindow= TRUE, #fill = "both",
                           horizontal=FALSE, expand=TRUE, visible = TRUE, 
                           width = 300, height = 300);
  #   size(selectionGroup) <- c(300,300);
  selectionListGroup <- glayout(container=selectionGroup, expand=TRUE);
  
  #  FIXME: maintain tooltips for all types
  tt <- mi::.prune("missing_variable");
  names(tt) <- tt
  tt_names <- names(tt)
  tt[] <- "stub"
  tt["irrelevant"] <- "EXCLUDE the variable from all models"
  tt["count"] <- "non-negative integers"
  tt["continuous"] <- "continuous but NOT special like positive-continuous, etc"
  tt["fixed"] <- "one unique value"
  tt["group"] <- "grouping factors for multilevel models"
  tt["grouped-binary"] <- "binary variable that is grouped into strata"
  tt["unordered-categorical"] <- "discrete variable with >= 3 values but WITHOUT ordering, sometimes called nominal or factor"
  tt["ordered-categorical"] <- "discrete variable with >= 3 values WITH ordering by increasing value"
  tt["interval"] <- "interval-censored from the left and / or right"
  tt["binary"] <- "discrete variable with exactly 2 values"
  tt["bounded-continuous"] <- "continuous variable with lower and / or upper bounds (but not a proportion)"
  tt["positive-continuous"] <- "continuous variable with no values of 0 or negative (but not a proportion)"
  tt["proportion"] <- "proportion on the (0,1) interval"
  tt["semi-continuous"] <- "continuous variable with a point mass at the lower and / or upper end of the distribution"
  tt["nonnegative-continuous"] <- "continuous variable with a point mass at zero and no negative values"
  tt["SC_proportion"] <- "proportion on the [0,1] interval"
  if(any(tt == "stub")) print(tt)
  tt <- paste(names(tt), tt, sep = ": ")
  names(tt) <- tt_names
  
  types <- sapply(colnames(GUI_ENV$user_df), FUN = function(y) {
    .guess_type(GUI_ENV$user_df[,y], favor_ordered = svalue(favor_ord),
                  favor_positive = svalue(favor_pos), threshold = svalue(thresholdButton),
                  variable_name = y)
  })
  names(types) <- colnames(GUI_ENV$user_df);
  GUI_ENV$mdf_args$types <- types;
  
  selectionListGroup[1,1] = glabel("<b>Variable</b>", container=selectionListGroup, markup=TRUE, anchor = c(1,0));
  selectionListGroup[1,2] = glabel("<b>Class</b>", container=selectionListGroup, markup=TRUE);
  
  if(!is.null(GUI_ENV$mdf_args$subclass) && GUI_ENV$mdf_args$subclass == "experiment") {
    concept <- c("covariate", "outcome", "treatment")
    GUI_ENV$mdf_args$concept <- factor(rep("covariate", length(types)), levels = concept)
  }
  for(i in seq_along(types)) {
    selectionListGroup[i+1,1] = glabel(names(types)[i], container = selectionListGroup,
                                       anchor=c(1,0));
    likelyTypeTF <- mi::.possible_missing_variable(GUI_ENV$user_df[,i])
    likelyTypes <- names(likelyTypeTF[likelyTypeTF])
    likelyTypes <- data.frame(likelyTypes, NA, tt[likelyTypes], stringsAsFactors = FALSE)
    selectedId <- which(likelyTypes == types[i]);
    selectionListGroup[i+1,2] = gcombobox(items = likelyTypes, anchor=c(-1,0),
                                          action=i, editable = TRUE,
                                          selected = selectedId,
                                          container = selectionListGroup,
                                          handler = function(h, ...) {
                                            GUI_ENV$mdf <- NULL;
                                            GUI_ENV$mdf_args$types[h$action] <- svalue(h$obj)
#                                             tooltip(selectionListGroup[h$action + 1,1]) <- 
#                                               tooltip(selectionListGroup[h$action + 1,2]) <- tt[svalue(h$obj)];
                                          });
#     tooltip(selectionListGroup[i+1,1]) <- 
#       tooltip(selectionListGroup[i+1,2]) <- likelyTypes[selectedId,3];
    if(!is.null(GUI_ENV$mdf_args$subclass) && GUI_ENV$mdf_args$subclass == "experiment") {
      selectionListGroup[i+1,3] = gcombobox(items = concept, anchor=c(-1,0), 
                                            editable = TRUE, action = i, selected = 1,
                                            container = selectionListGroup, handler = function(h, ...) {
                                              GUI_ENV$mdf_args$concept[h$action] <- svalue(h$obj)
                                            })
    }
  }
  #   visible(selectionList) <- TRUE;
}

# STEP 4
do_var_options <- function() {
  if (is.null(GUI_ENV$mdf)) {
    GUI_ENV$mdf_args$y <- GUI_ENV$user_df
    GUI_ENV$mdf <- do.call("missing_data.frame", args = GUI_ENV$mdf_args)
    GUI_ENV$imputations <- NULL
  }
  innerGroup <- setupScreen(step = "4",
                            title = "Options for Variables",
                            backHandler = navHandler("get_var_classes"),
                            nextHandler = navHandler("pre_mi"));
  topGroup <- ggroup(horizontal = FALSE, container = innerGroup);
  glabel("Choose a variable and then its associated options:", container = topGroup);
  #   addSpring(topGroup)
  
  varnames <- colnames(GUI_ENV$mdf)

  varHandler <- function(h, ...) {
    mv <- GUI_ENV$mdf@variables[[svalue(h$obj)]];
    varHandler2(mv);
    return(invisible(NULL));
  }
  
  varHandler2 <- function(mv) {
    classes <- mi::.possible_missing_variable(GUI_ENV$user_df[,mv@variable_name])
    classes <- names(classes[classes])
    tmp <- classCombo
    if(!identical(svalue(tmp), class(mv)[1])) {
      svalue(tmp) <- class(mv)[1]
      tmp[] <- classes;
      svalue(tmp) <- class(mv)[1]
    }
    tmp <- miss
    svalue(tmp) <- paste(" # missing: ", mv@n_miss, sep = "");
    varHandler3(mv);
    
    codes <- unique(mv@raw_data);
    codes <- codes[!is.na(codes)]
    if(is(mv, "continuous")) codes <- sort(codes)
    tmp <- table.group[5,2]
    svalue(tmp) <- "";
    tmp[] <- codes
    table.group[5,4] <- gcombobox(items = c("NA", "unpossible", codes), 
                                  selected = 0, editable = TRUE, 
                                  container = table.group, action = "recode", 
                                  handler = changeHandler);
    
#     tmp <- 
#     svalue(tmp) <- "";
#     if(is(mv, "continuous")) tmp[] <- c("NA", "unpossible")
#     else tmp[] <- c("NA", "unpossible", codes)

    impModel <- .default_model(mv, GUI_ENV$mdf)[1]
    #   FIXME: do something about the options for the embedded ordered-categorical
    tmp <- model
    if(is.na(impModel)) {
      svalue(tmp) <- "No model needed because no missing values"
    }
    else svalue(tmp) <- paste("Model implied by the above options:", impModel)
    
    return(invisible(NULL));
  }
  
  varHandler3 <- function(mv) {
    if (!mv@all_obs) {
      for(i in 1:3) {
        tmp <- table.group[i,2]
        enabled(tmp) <- TRUE
      }
      im <- getClass(class(mv))@prototype@imputation_method      
      im <- im[!is.na(im)]
      im <- sapply(im, rosettaStone)
      tmp <- table.group[1,2]
      tmp[] <- im # FIXME: consider better names, but if so change the changeHandler()
      svalue(tmp) <- rosettaStone(mv@imputation_method)
#       tooltip(table.group[1,1]) <-
#         tooltip(table.group[1,2]) <- "Way of imputing missing values given a model"
      tmp <- table.group[2,2]
      svalue(tmp) <- mv@family$family
      tmp[] <- getClass(class(mv))@prototype@known_families
      svalue(tmp) <- mv@family$family
#       tooltip(table.group[2,1]) <- 
#         tooltip(table.group[2,2]) <- "Likelihood function used to model variable" 
      tmp <- table.group[3,2]
      svalue(tmp) <- mv@family$link
      tmp[] <- getClass(class(mv))@prototype@known_links
      svalue(tmp) <- mv@family$link
#       tooltip(table.group[3,1]) <- 
#         tooltip(table.group[3,2]) <- "Link function used to map between linear predictor and support of outcome"
    }
    else for(i in 1:3) {
      tmp <- table.group[i,2]
      svalue(tmp) <- ""
      enabled(tmp) <- FALSE
    }
    
    if (is(mv, "continuous")) {
      tmp <- table.group[4,2]
      enabled(tmp) <- TRUE
      svalue(tmp) <- .parse_trans(mv@transformation)
      tmp[] <- mv@known_transformations
      svalue(tmp) <- .parse_trans(mv@transformation)
#       tooltip(table.group[4,1]) <- 
#         tooltip(table.group[4,2]) <- paste("Transformation applied to variable both when it is being modeled and when",
#                                            "it is used to model other variables")
    }
    else {
      tmp <- table.group[4,2]
      svalue(tmp) <- ""
      enabled(tmp) <- FALSE
    }
    return(invisible(NULL))
  }
  
  rosettaStone <- function(im, inverse = FALSE) {
    if(inverse) switch(im,
                       "posterior predictive distribution" = "ppd",
                       "predictive mean matching" = "pmm",
                       "conditional mean of observed" = "mean",
                       "conditional median of observed" = "median",
                       "conditional mean" = "expectation",
                       "highest conditional probability" = "mode",
                       NA_character_)
    else switch(im,
                ppd = "posterior predictive distribution",
                pmm = "predictive mean matching",
                mean = "conditional mean of observed",
                median = "conditional median of observed",
                expectation = "conditional mean",
                mode = "highest conditional probability",
                "not recognized")
  }
  
  changeHandler <- function(h, ...) {
    to <- svalue(h$obj)
    if(is.null(to)) return(invisible(NULL))
    if(to == "") return(invisible(NULL))
    what <- h$action
    if(what == "imputation_method") to <- rosettaStone(to, inverse = TRUE)
    
    # try to avoid unnecessary changes
    mv <- GUI_ENV$mdf@variables[[svalue(selectedVar)]];
    if(what == "imputation_method" && mv@imputation_method == to) return(invisible(NULL))
    if(what == "family" && mv@family$family == to) return(invisible(NULL))
    if(what == "link" && mv@family$link == to) return(invisible(NULL))
    if(is(mv, "continuous") && what == "transformation" && 
       .parse_trans(mv@transformation) == to) return(invisible(NULL))
    if(what == "class" && class(mv) == to) return(invisible(NULL))
    
    if(what == "recode") {
      what <- svalue(table.group[5,2]);
      if(to == "NA") to <- NA;
    }
    
    GUI_ENV$mdf <- change(data = GUI_ENV$mdf, y = svalue(selectedVar), 
                          what = what, to = to);
    GUI_ENV$imputations <- NULL;
    mv <- GUI_ENV$mdf@variables[[svalue(selectedVar)]];
    
    varHandler3(mv);
    return(invisible(NULL));
  }
  varGroup <- ggroup(horizontal = TRUE, container = topGroup)
  selectedVar <- gcombobox(items = varnames, selected = 1, container = varGroup, 
                           editable = TRUE, handler = varHandler);
  mv <- GUI_ENV$mdf@variables[[1]]
  glabel(" Class: ", container = varGroup);
  classes <- mi::.prune(class(mv));
  classCombo <- gcombobox(classes, selected = which(class(mv) == classes), 
                          action = "class", editable = TRUE,
                          container = varGroup, handler = changeHandler);
  miss <- glabel(paste(" # missing: ", mv@n_miss, sep = ""), container = varGroup);        
  gbutton("histogram ...", container = varGroup,
          handler = function(h, ...) {
            mv <- GUI_ENV$mdf@variables[[svalue(selectedVar)]]
            hist(mv);
            title(main = paste("Observed", mv@variable_name));
          });
  addSpring(innerGroup)
  gseparator(horizontal = TRUE, container = innerGroup);
  addSpace(innerGroup, value = 30);
  table.group <- glayout(container = innerGroup, expand = TRUE, fill = "both");
  
  # imputation_method
  ims <- getClass(class(mv))@prototype@imputation_method;
  ims <- ims[!is.na(ims)]
  table.group[1,1] <- glabel("imputation method", container = table.group);
  table.group[1,2] <- gcombobox(items = sapply(ims, rosettaStone),
                                selected = 1, editable = TRUE,
                                handler = changeHandler,
                                action = "imputation_method",
                                container = table.group);
  
  # family
  table.group[2,1] <- glabel("family (likelihood)", container = table.group);
  table.group[2,2] <- gcombobox(items = ifelse(is.na(mv@family), "none", mv@family$family),
                                selected = 1, editable = TRUE,
                                handler = changeHandler,
                                action = "family",
                                container = table.group);
  
  # link
  table.group[3,1] <- glabel("link function", container = table.group);
  table.group[3,2] <- gcombobox(items = ifelse(is.na(mv@family), "none", mv@family$link),
                                selected = 1, editable = TRUE,
                                handler = changeHandler,
                                action = "link",
                                container = table.group);
  
  # transformation (continuous variables only)
  table.group[4,1] <- glabel("transformation function", container = table.group);
  table.group[4,2] <- gcombobox(items = if(is(mv, "continuous")) .parse_trans(mv@transformation) else "none",
                                selected = 1, editable = TRUE,
                                handler = changeHandler,
                                action = "transformation",
                                container = table.group);
  for(i in 1:4) {
    tmp <- table.group[i,2]
    size(tmp) <- c(250,25);
  }
  #   visible(table.group) <- TRUE;
  #   addSpring(innerGroup);
  #   recodeGroup <- ggroup(horizontal = TRUE, container = innerGroup);
  table.group[5,1] <- glabel("Recode", container = table.group);
  codes <- unique(mv@raw_data);
  codes <- codes[!is.na(codes)];
  if(is(mv, "continuous")) codes <- sort(codes);
  table.group[5,2] <- gcombobox(items = codes, selected = 0, 
                                editable = TRUE, container = table.group);
  table.group[5,3] <- glabel("as", container = table.group);
  table.group[5,4] <- gcombobox(items = c("NA", "unpossible", codes), 
                                selected = 0, editable = TRUE, 
                                container = table.group, action = "recode", handler = changeHandler);
  gseparator(horizontal = TRUE, container = innerGroup);
  modelGroup <- ggroup(horizontal = TRUE, container = innerGroup);
  model <- glabel("", container = modelGroup);
  #   addSpring(modelGroup);
  #   addSpring(innerGroup);
  varHandler2(mv);
}

# STEP 5
do_pre_mi <- function() {
  innerGroup <- setupScreen(step = "5",
                            title = "Pre-imputation Diagnostics",
                            backHandler = navHandler("var_options"),
                            nextHandler = navHandler("impute"));
  
  choiceLabel("Click one or more", container = innerGroup);
  button.group <- ggroup(horizontal = TRUE, container = indent(innerGroup));
  table.group <- glayout(horizontal = TRUE, container = button.group, visible = TRUE);
  table.group[1,1] <- gbutton("image plot ...", container = table.group,
                              handler = function(h, ...) {
                                image(GUI_ENV$mdf, ask = FALSE);
                                return(invisible(NULL));
                                # The following works for RGtk2
                                w <- gwindow(visible = FALSE);
                                ggraphics(container = w);
                                visible(w) <- TRUE
                                image(GUI_ENV$mdf, ask = FALSE);
                                # but right-click to save throws an error
                              });
  table.group[2,1] <- gbutton("histograms ...", container = table.group,
                              handler = function(h, ...) hist(GUI_ENV$mdf, ask = FALSE));
  table.group[3,1] <- gbutton("numerical summary ...", container = table.group,
                              handler = function(h, ...) print(summary(as.data.frame(GUI_ENV$mdf))));
  #   visible(table.group) <- TRUE;
  addSpring(button.group);
}

# STEP 6
do_impute <- function() {
  innerGroup <- setupScreen(step = "6",
                            title = "Options for Imputation",
                            backHandler = navHandler("pre_mi"),
                            nextHandler = function(h,...) { runImputationInDialog(); });
  
  choiceLabel("Choose options for imputation",container = innerGroup);
  table.group <- glayout(horizontal = TRUE, container = indent(innerGroup), visible = TRUE);
  table.group[1,1] <- glabel("Number of chains:", container=table.group);
  table.group[1,2] <- gspinbutton(from = 1, to=100, by=1, value=GUI_ENV$n.chains, container=table.group,
                                  handler= function(h,...) {
                                    GUI_ENV$n.chains = svalue(h$obj);
                                  });
  table.group[2,1] <- glabel("Number of iterations:", container=table.group);
  table.group[2,2] <- gspinbutton(from = 1, to=1000, by=1, value=GUI_ENV$n.iter, container=table.group,
                                  handler = function(h,...) {
                                    GUI_ENV$n.iter = svalue(h$obj);
                                  });
  #   visible(table.group) <- TRUE;
}

# STEP 7 (in its own popup dialog)
runImputationInDialog <- function() {
  exitImputationDialog <- function() {
    dispose(dialogWindow);
    nextScreen();
    return(invisible(NULL))
  }
  foo <- function() {
    on.exit(exitImputationDialog());
    if (is.null(GUI_ENV$imputations)) {
      GUI_ENV$state <- "impute"
      GUI_ENV$imputations <- mi(GUI_ENV$mdf, n.chains = GUI_ENV$n.chains, 
                                n.iter = GUI_ENV$n.iter, parallel = svalue(parallel));
    }
    else {
      GUI_ENV$state <- "convergence_eval";
      GUI_ENV$imputations <- mi(GUI_ENV$imputations, n.iter = GUI_ENV$n.iter);
    }
    GUI_ENV$state <- "convergence_eval";
    return(invisible(NULL))
  }
  dialogWindow <- gwindow("Running mi");
  group <- ggroup(container=dialogWindow, horizontal=FALSE, expand=TRUE);
  addTitle("Step 7: Running mi", group);
  addSpring(group);
  optionsGroup <- ggroup(container=group, horizontal=TRUE, expand=TRUE);
  parallel <- gcheckbox("Run in parallel", checked=TRUE, container=optionsGroup);
#   launch   <- gcheckbox("Launch browser to view progress", checked=TRUE, container=optionsGroup);
  txt <- if(.Platform$OS.type == "windows") "ESC" else "<CTRL>+c"
  exec <- gbutton("Run mi",container=group, handler = function(h, ...) foo())
  focus(exec) <- TRUE
  addSpring(group);
}

# STEP 8
do_convergence_eval <- function() {
  innerGroup <- setupScreen(step = "8",
                            title = "Evaluate Convergence",
                            backHandler=navHandler("impute"),
                            nextHandler=navHandler("parametric_eval"));
  
  choiceLabel("Click to show convergence evaluations", container=innerGroup);
  BUGS <- NULL;
  gbutton("R-hat statistics...",
          container=innerGroup,
          handler = function(h,...) {
            if (is.null(BUGS)) BUGS <- mi::Rhats(GUI_ENV$imputations)
            print(round(BUGS, 2));
          });

  addSpring(innerGroup);
  gseparator(horizontal=TRUE, container=innerGroup);
  addSpring(innerGroup);
  
  #   choiceLabel("Continue for more iterations", container = innerGroup);
  cont.group <- ggroup(horizontal=TRUE, expand=TRUE, container=innerGroup);
  gbutton("Continue ...",
          handler = function(h, ...) {
            runImputationInDialog();
          },
          container= cont.group);
  glabel("mi for", container = cont.group);
  gspinbutton(from = 1, to = 1000, by = 1,
              value = GUI_ENV$n.iter, container = cont.group,
              handler = function(h, ...) {
                GUI_ENV$n.iter = svalue(h$obj);
              });
  glabel("more iterations", markup = FALSE, container = cont.group);
}

# STEP 9
do_parametric_eval <- function() {
  innerGroup <- setupScreen(step = "9",
                            title = "Parametric Diagnostics",
                            backHandler = navHandler("convergence_eval"),
                            nextHandler = navHandler("analysis"));
  
  choiceLabel("Run one or more of the following diagnostics", container = innerGroup);
  button.group <- ggroup(horizontal = TRUE, container = indent(innerGroup));
  variable_names <- colnames(GUI_ENV$mdf);
  variable_names <- c("Choose variable", variable_names);
  table.group <- glayout(horizontal = TRUE, container = button.group, visible = TRUE);
  
  table.group[1,1] <- glabel("Plot model fits:", container = table.group);
  table.group[1,2] <- gcombobox(items = variable_names, selected = 1, 
                                editable = TRUE, container = table.group,
                                handler = function(h, ...) {
                                  y <- svalue(h$obj)
                                  if (y != "Choose variable") plot(x = GUI_ENV$imputations, y = y, ask = FALSE)
                                  else dev.off();
                                });
  table.group[2,2] <- gbutton("image plot ...", container = table.group,
                              handler = function(h, ...) image(GUI_ENV$imputations, ask = FALSE));
  table.group[3,2] <- gbutton("histograms for chain ...", container = table.group,
                              handler = function(h, ...) {
                                m <- svalue(table.group[3,3])
                                hist(GUI_ENV$imputations, m = m, ask = FALSE);
                                title(sub = paste("Chain", m), outer = TRUE)
                              });
  table.group[3,3] <- gspinbutton(from = 1, to = length(GUI_ENV$imputations), by = 1, value = 1,
                                  container = table.group)
  table.group[4,2] <- gbutton("numerical summary ...", container = table.group,
                              handler = function(h, ...) 
                                print(summary(do.call(rbind, complete(GUI_ENV$imputations)))));
  #   visible(table.group) <- TRUE;
  addSpring(button.group);
}

# STEP 10
# do_structural_eval <- function() {
#   stop("this function should not be called")
#   innerGroup <- setupScreen(step = "10",
#                             title = "Structural Diagnostics",
#                             backHandler = navHandler("parametric_eval"),
#                             nextHandler = navHandler("analysis"));
#   
#   choiceLabel("Estimate switching regression", container = innerGroup);
#   button.group <- ggroup(horizontal = TRUE, container = indent(innerGroup));
#   variable_names <- colnames(GUI_ENV$mdf);
#   variable_names <- c("Choose variable", variable_names);
#   table.group <- glayout(horizontal = TRUE, container = button.group, visible = TRUE);
#   table.group[1,1] <- glabel("Dependent variable:", container = table.group);
#   table.group[1,2] <- gcombobox(items = variable_names, selected = 1, 
#                                 editable = TRUE, container = table.group,
#                                 handler = function(h, ...) {
#                                   y <- svalue(h$obj)
#                                   if (y != "Choose variable") {
#                                     GUI_ENV$SR <- mi:::tobin5(GUI_ENV$imputations@data[[1]], y = y) # FIXME
#                                     print(summary(GUI_ENV$SR))
#                                   }
#                                 });
#   #   visible(table.group) <- TRUE;
# }

# STEP 10
do_analysis <- function() {
  innerGroup <- setupFinalScreen(step = "10",
                                 title = "Analyze Completed Datasets",
                                 backHandler = navHandler("parametric_eval"),
                                 nextHandler = function(h,...) { exit_gui(h); });
  
  gbutton("Estimate pooled model ...", container = innerGroup, handler = function(h, ...) {
    runPooledAnalysisDialog();
  });
  gbutton("Create completed data.frames in R's environment", container = innerGroup, handler = function(h, ...) {
    .GlobalEnv$dfs <- mi::complete(GUI_ENV$imputations);
    msg <- "data.frames are in R's environment with object name 'dfs'"
    gmessage(msg, , icon = "info")
  });
  gbutton("Save imputations in R binary format", container = innerGroup, handler = function(h, ...) {
    file_name <- gfile(text="Choosing File Name to Save As", type="save");
    if (is.na(file_name))
      return();
    if (!grepl("\\.R[dD]ata$", x = file_name, perl = TRUE)) {
      file_name <- paste(file_name, "RData", sep = ".");
    }
    saveMe <- GUI_ENV$imputations
    save(saveMe, file = file_name)
  });
  gbutton("Export imputations for Stata", container = innerGroup, handler = function(h, ...) {
    file_name <- gfile(text="Choosing File Name to Save As", type="save");
    if (is.na(file_name))
      return();
    if (!grepl("\\.dta$", x = file_name, perl = TRUE)) {
      file_name <- paste(file_name, "dta", sep = ".");
    }
    mi::mi2stata(GUI_ENV$imputations, m = length(GUI_ENV$imputations), file = file_name)
  })
  addSpring(innerGroup)
  
}

# STEP 11  (pops up something to run, then closes)
runPooledAnalysisDialog <- function() {
  dialogWindow <- gwindow("Pooled analysis");
  dialog.group <- ggroup(horizontal = FALSE, container = dialogWindow, expand = TRUE);  
  addTitle("Step 12: Specify analysis model", dialog.group);
  addSpring(dialog.group);
  
  missing <- colnames(GUI_ENV$mdf@X);
  missing <- missing[grepl("^missing_", missing)];
  var.names <- c(colnames(GUI_ENV$imputations), missing);
  var.group <- glabel(strwrap(paste("Variables:", paste(var.names, collapse = ", "))), container = dialog.group);
  #   var.group <- ggroup(horizontal = TRUE, container = dialog.group, use.scrollwindow = TRUE, expand = TRUE);
  #   var.table.group <- glayout(horizontal = FALSE, container = var.group, visible = TRUE);
  #   row.index <- 1
  #   col.index <- 1
  #   for(i in seq_along(var.names)) {
  #     var.table.group[row.index,col.index] <- glabel(var.names[i], container = var.table.group);
  #     if(col.index == 5) {
  #       row.index <- row.index + 1;
  #       col.index <- 1;
  #     }
  #     else col.index <- col.index + 1;
  #   }
  #   visible(var.table.group) <- TRUE;
  addSpring(dialog.group);
  
  possibleFunctions <- c("bayesglm", "glm", "lm");
  estimatingFunction <- possibleFunctions[1]
  estimatingFormula <- NA_character_;
  m <- length(GUI_ENV$imputations)
  dotdotdot <- ""
  table.group <- glayout(horizontal = FALSE, container = dialog.group, visible = TRUE);
  table.group[1,1] <- glabel("estimating function", container = table.group);
  table.group[1,2] <- gcombobox(possibleFunctions, selected = 1, 
                                editable = TRUE, container = table.group);
  table.group[2,1] <- glabel("formula", container = table.group);
  table.group[2,2] <- gedit("", container = table.group, initial.msg = "formula like: y ~ x1 + x2 ...");
  table.group[3,1] <- glabel("number of chains", container = table.group);
  table.group[3,2] <- gspinbutton(from = 1, to = 2 * length(GUI_ENV$imputations), value = length(GUI_ENV$imputations),
                                  by = 1, container = table.group);
  table.group[4,1] <- glabel("...", container = table.group);
  table.group[4,2] <- gedit("", container = table.group, width = 30,
                            initial.msg = "any additional arguments like: foo = bar",
                            handler = function(h, ...) {
                              dotdotdot <- svalue(h$obj);
                            });
  #   visible(table.group) <- TRUE;
  
  gbutton("Estimate", container = dialog.group, handler = function(h, ...) {
    poolCall <- paste("pool(formula = ", svalue(table.group[2,2]),
                      ", data = GUI_ENV$imputations",
                      ", m = ", svalue(table.group[3,2]),
                      ", FUN = ", svalue(table.group[1,2]), sep = "")
    if (svalue(table.group[4,2]) != "") poolCall <- paste(poolCall, svalue(table.group[4,2]), sep = ", ")
    poolCall <- paste(poolCall, ")")
    model <- eval(parse(text = poolCall))
    print(model)
  });
  addSpring(dialog.group)
  gbutton("Return to previous screen", container = dialog.group, handler = function(h, ...) {
    dispose(dialogWindow);
  });
}

# GENERAL UTILITIES
# ========================================

helpHandler <- function(h,...) {
  URL <- "http://www.google.com";
  # FIXME: add help pages to migui package under migui/inst/htmlHelp/g
  # URL <- file.path(dirname(system.file(package = "migui")),
  #                          "migui","inst","htmlHelp",
  #                          paste(GUI_ENV$state,".html",sep=""));
  browseURL(URL);
};


setupScreen <- function(step, title, backHandler, nextHandler) {
  setupScreenHelper(step,title,backHandler,nextHandler,"NEXT >");
}  

setupFinalScreen <- function(step, title, backHandler, nextHandler) {
  setupScreenHelper(step,title,backHandler,nextHandler,"FINISH >");
}

setupScreenHelper <- function(step, title, backHandler, nextHandler, nextButtonText) {
  topGroup <- ggroup(horizontal=FALSE, expand=TRUE, fill="both",
                     container = GUI_ENV$group, visible = FALSE);
  glabel(paste("<big><b>","Step ", step, ".  <i>", title,"</i></b></big>",
               sep = ""), markup=TRUE,
         container=topGroup, anchor=c(0,-1));
  
  addSpring(topGroup);
  
  innerGroup <- ggroup(horizontal=FALSE, expand=TRUE, fill="both",
                       container = topGroup);
  
  addSpring(topGroup);
  
  gseparator(horizontal=TRUE, container = topGroup);
  button.group <- ggroup(container = topGroup, horizontal=TRUE);
  gbutton(" QUIT! ", container=button.group,
          handler = function(h,...) { exit_gui(h) });
  addSpring(button.group);
  if (!is.null(backHandler)) {
    gbutton(" < BACK ", container=button.group, handler=backHandler);
    addSpace(button.group,10);
  }
  gbutton(" HELP... ", container=button.group, handler=helpHandler);
  addSpace(button.group,10);
  nextButton <- gbutton(nextButtonText, container=button.group, handler=nextHandler);
  # focus(nextButton) <- TRUE
  return(innerGroup);
}

navHandler <- function(stateName) {
  return(function(h, ...) {
    GUI_ENV$state = stateName;
    nextScreen();
  })
}

choiceLabel <- function(title, container) {
  glabel(paste("<b>",title,":</b>"), markup=TRUE,
         expand=FALSE, fill="x", anchor=c(-1,0),
         container=container);
  addSpace(container,3);
}

indent <- function(container) {
  group <- ggroup(horizontal=TRUE, container=container, fill=FALSE);
  addSpace(group,5);
  return(group);
}

exit_gui <- function(h) {
  #   dispose(h$obj);
  reallyQuit <- gconfirm(paste("Do you really want to quit?\nCancel means 'no' while 'OK' means 'yes'"), icon = "warning")
  if(reallyQuit) dispose(GUI_ENV$window);
  return(invisible(NULL)); 
}

addTitle <- function(title,container) {
  title_group <- ggroup(horizontal=TRUE, expand=TRUE, fill="both",
                        container=container);
  glabel(paste("<b>",title,"</b>"),
         container=title_group,
         markup=TRUE);
}

.default_model <-
  function(y, data) {
    if(is(data, "allcategorical_missing_data.frame")) return("Gibbs")
    if(y@all_obs) {
      if(is(y, "semi-continuous")) return(rep(NA_character_, 2))
      else return(NA_character_)
    }
    if(is(y, "irrelevant")) return(NA_character_)
    if(y@imputation_method == "mcar") return(NA_character_)
    if(!is.method_in_mi("fit_model", y = class(y), data = class(data))) {
      if(is(y, "semi-continuous")) return(rep("user-defined", 2))
      else return("user-defined")
    }
    fam <- y@family$family
    link <- y@family$link
    if(is(y, "count")) {
      if(fam == "quasipoisson" && link == "log") return("qpoisson")
      else if(fam == "poisson" && link == "log") return("poisson")
      else return("****")
    }
    else if(is(y, "binary")) {
      if(is(y, "grouped-binary")) return("clogit")
      if(fam == "quasibinomial")  return(paste("q", link, sep = ""))
      else if(fam == "binomial")  return(link)
      else return("****")
    }
    else if(is(y, "interval")) return("survreg")
    else if(is(y, "ordered-categorical"))   return(paste("o", link, sep = ""))
    else if(is(y, "unordered-categorical")) {
      if(fam == "binomial") out <- "RN"
      else out <- "m"
      return(paste(out, link, sep = ""))
    }
    else if(is(y, "proportion")) return(if(fam == "gaussian") "linear" else "betareg")
    else if(is(y, "SC_proportion")) {
      out <- .default_model(y@indicator, data)
      return(c("betareg", out))
    }
    else if(is(y, "semi-continuous")) {
      out <- .default_model(y@indicator, data)
      if(fam == "gaussian") { 
        if(link == "identity") return(c("linear", out))
        else if(link == "log") return(c("loglinear", out))
        else if(link == "inverse") return(c("inverselinear", out))
        else return(c("****", out))
      }
      else if(fam == "Gamma") return(c("****", out))
      else if(fam == "inverse.gaussian") return(c("****", out))
      else if(fam == "quasi") return(c("quasi", out))
      else return(c("****", out))
    }
    else if(is(y, "continuous")) {
      if(fam == "gaussian") { 
        if(link == "identity") return("linear")
        else if(link == "log") return("loglinear")
        else if(link == "inverse") return("inverselinear")
        else return("****")
      }
      else if(fam == "Gamma") return("****")
      else if(fam == "inverse.gaussian") return("****")
      else if(fam == "quasi") return("quasi")
      else return("****")
    }
    else return("user-defined")
  }

.guess_type <-
  function(y, favor_ordered = TRUE, favor_positive = FALSE, threshold = 5,
           variable_name = deparse(substitute(y))) {
    
    if(!is.null(dim(y))) stop(paste(variable_name, ": must be a vector"))
    if(is.factor(y)) y <- factor(y) # to drop unused levels    
    values <- unique(y)
    values <- sort(values[!is.na(values)])
    len <- length(values)
    if(len == 0) {
      warning(paste(variable_name, ": cannot infer variable type when all values are NA, guessing 'irrelevant'"))
      type <- "irrelevant"
    }
    else if(len == 1)             type <- "fixed"
    else if(grepl("^[[:punct:]]", 
                  variable_name)) type <- "irrelevant"
    else if(identical("id",
                      tolower(variable_name))) type <- "irrelevant"
    else if(len == 2) {
      if(!is.numeric(values))     type <- "binary"
      else if(all(values == 
                  as.integer(values))) type <- "binary"
      else if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else if(is.ts(y)) {
      if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else if(is.ordered(y))        type <-   "ordered-categorical"
    else if(is.factor(y))         type <- "unordered-categorical"
    else if(is.character(y))      type <- "unordered-categorical"
    else if(is.numeric(y)) {
      if(all(values >= 0) && 
         all(values <= 1)) {
        
        if(any(values %in% 0:1))  type <- "SC_proportion"
        else                      type <- "proportion"
      }
      else if(len <= threshold && 
              all(values == as.integer(values)))
        type <- if(favor_ordered) "ordered-categorical" else "unordered-categorical"
      else if(favor_positive) {
        if(all(values > 0))       type <- "positive-continuous"
        else if(all(values >= 0)) type <- "nonnegative-continuous"
        else                      type <- "continuous"
      }
      else                        type <- "continuous"
    }
    else stop(paste("cannot infer variable type for", variable_name))
    
    return(type)
  }

.parse_trans <-
  function(trans) {
    if(identical(names(formals(trans)), c("y", "mean", "sd", "inverse"))) return("standardize")
    if(identical(names(formals(trans)), c("y", "a", "inverse"))) return("logshift")
    if(identical(body(trans), body(.squeeze_transform))) return("squeeze")
    if(identical(body(trans), body(.identity_transform))) return("identity")
    if(identical(body(trans), body(log))) return("log")
    if(identical(body(trans), body(sqrt))) return("sqrt")
    if(identical(body(trans), body(.cuberoot))) return("cuberoot")
    if(identical(body(trans), body(qnorm))) return("qnorm")
    return("user-defined")
  }

.squeeze_transform <- function(y, inverse = FALSE) {
  n <- length(y)
  if(inverse) (y * n - .5) / (n - 1)
  else (y * (n - 1) + .5) / n
}

.cuberoot <-
  function(y, inverse = FALSE) {
    if(inverse) y^3
    else        y^(1/3)
  }

.identity_transform <- function(y, ...) return(y)

is.method_in_mi <-
  function(generic, ...) {
    method <- selectMethod(generic, signature(...))
    return(environmentName(environment(method@.Data)) == "mi")
  }
