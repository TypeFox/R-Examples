#' Fuzzy Forest Object
#'
#' Fuzzy forests returns an object of type
#' fuzzyforest.
#' @export
#' @param feature_list      List of selected features along with variable
#'                          importance measures.
#' @param final_rf          A final random forest fit using the features
#'                          selected by fuzzy forests.
#' @param module_membership Module membership of each feature.
#' @param WGCNA_object      If applicable, output of WGCNA analysis.
#' @param survivor_list     List of features that have survived screening step.
#' @param selection_list    List of features retained at each iteration of
#'                          selection step.
#' @return An object of type fuzzy_forest.
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
fuzzy_forest <- function(feature_list, final_rf, module_membership,
                         WGCNA_object=NULL, survivor_list, selection_list) {
  out <- list()
  out[[1]] <- feature_list
  out[[2]] <- final_rf
  out[[3]] <- module_membership
  out[[4]] <- module_membership
  out[[5]] <- survivor_list
  out[[6]] <- selection_list
  names(out) <- c("feature_list", "final_rf", "module_membership",
                  "WGCNA_object", "survivor_list", "selection_list")
  class(out) <- "fuzzy_forest"
  return(out)
}


#' Print fuzzy_forest object.
#' Prints output from fuzzy forests algorithm.
#' @export
#' @param x   A fuzzy_forest object.
#' @param ... Additional arguments not in use.
#' @return data.frame with list of selected features and variable
#'          importance measures.
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
print.fuzzy_forest <- function(x, ...) {
  print(x$feature_list)
  if(!is.null(x$final_rf$test)) {
    if(!is.null(x$final_rf$test$mse)) {
      cat(c("test set error: ", x$final_rf$test$mse[x$final_rf$ntree]))
    }
    if(!is.null(x$final_rf$test$err.rate)) {
      cat(c("test set error: ", x$final_rf$test$err.rate[x$final_rf$ntree]))
    }
  }
}

#' Predict method for fuzzy_forest object.
#' Obtains predictions from fuzzy forest algorithm.
#' @export
#' @param object   A fuzzy_forest object.
#' @param new_data A matrix or data.frame containing new_data.
#'                 Pay close attention to ensure feature names
#'                 match between training set and test set
#'                 data.frame.
#' @param ...      Additional arguments not in use.
#' @return A vector of predictions
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
predict.fuzzy_forest <- function(object, new_data, ...) {
  out <- predict(object$final_rf, new_data)
  return(out)
}

#' Plots relative importance of modules.
#'
#' The plot is designed
#' to depict the size of each module and what percentage of selected
#' features fall into each module.  In particular, it is easy to
#' determine which module is over-represented in the group of selected
#' features.
#' @export
#' @param object   A fuzzy_forest object.
#' @param main Title of plot.
#' @param xlab Title for the x axis.
#' @param ylab Title for the y axis.
#' @param module_labels Labels for the modules.  A data.frame
#'                      or character matrix with first column giving
#'                      the current name of module and second column giving
#'                      the assigned name of each module.
#'
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
modplot <- function(object, main=NULL, xlab=NULL, ylab=NULL,
                              module_labels=NULL) {
  if(is.null(main)) {
    main <- "Module Membership Distribution"
  }
  if(is.null(xlab)) {
   xlab <- "Module"
  }
  if(is.null(ylab)) {
    ylab <- "Percentage of features in module"
  }

  #allows user to supply new names for modules
  if(!is.null(module_labels)) {
    old_labels <- object$module_membership$module
    #module_labels should be re-ordered so that the old labels are in
    #alphabetical order.  This is because factor(old_labels) has levels in
    #alphabetical order. Note that the `labels` below is contains new labels.
    module_labels <- module_labels[order(module_labels[, 1]), ]
    new_labels <- as.character(factor(old_labels, labels=module_labels[, 2]))
    object$module_membership$module <- new_labels

    #Now module labels need to be changed for the table of variable importances.
    select_mods <- as.factor(object$feature_list$module_membership)
    select_module_table <- module_labels[which(module_labels[, 1] %in%
                                        levels(select_mods)), ,drop=FALSE]

    #This line of code may be slightly dangerous depending on where "." is.
    #It should be ok because after removing "." the remaining levels are in
    #alphabetical order.
    if( "." %in% levels(select_mods)) {
      dot_index <- which(levels(select_mods) == ".")
      levels(select_mods)[-dot_index] <- select_module_table[, 2]
      }
      else {
        levels(select_mods) <- select_module_table[, 2]
      }
      object$feature_list$module_membership <- as.character(select_mods)
  }
  mods <- object$module_membership$module
  mod_length <- length(mods)
  mod_tab <- table(mods)
  mod_name <- names(mod_tab)
  feature_list <- object$feature_list$module_membership
  #This line is here in the case that some covariates are not in a module
  mod_feature_list <- feature_list[feature_list != "."]
  #Table showing how many important features are in each selected module.
  imp_feature_tab <- table(mod_feature_list)
  imp_names <- names(imp_feature_tab)
  feature_tab <- rep(0, length(mod_tab))
  names(feature_tab) <- mod_name
  for(i in 1:length(feature_tab)) {
    if(mod_name[i] %in% names(imp_feature_tab)) {
      feature_tab[i] <- imp_feature_tab[which(imp_names == mod_name[i])]
      }
    }
  unimportant_pct <- (mod_tab - feature_tab)/mod_length
  important_pct <- feature_tab/mod_length
  mod_name <- rep(mod_name, 2)
  pct <- c(unimportant_pct, important_pct)
  pct_type <- rep(c("% Unimportant", "% Important"), each=length(mod_tab))
  importance_pct <- data.frame(Module=mod_name, Status=pct_type,
                               Percentage=pct)
  #test whether labels are numeric
  #reorder the labels if they are numeric
  num_mods <- suppressWarnings(as.numeric(object$module_membership[, 2]))
  num_test <- sum(is.na(num_mods))
  if(num_test == 0) {
    importance_pct[, 1] <- as.factor(importance_pct[, 1])
    levels(importance_pct[, 1]) <- sort(unique(num_mods))
  }

  #this is a work-around to get rid of notes in R CMD Check
  #Probably a better way to address the issue
  Module <- NULL
  Percentage <- NULL
  Status <- NULL
  ###############
  imp_plot <- ggplot(importance_pct, aes(x=Module, y=Percentage, fill=Status)) +
              geom_bar(stat="identity") +
              ggtitle(main) + labs(x = xlab, y = ylab) +
              theme(plot.title = element_text(lineheight=.8, face="bold"),
                    legend.title = element_blank())
  plot(imp_plot)
}
