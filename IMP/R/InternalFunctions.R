
error_handler <- function(list_models,arg,method) {

  if (!is.numeric(arg)) {
    if (method == "Performance_Measure") {
      stop("Error: The argument g should be numeric")
    } else {
      stop("Error: The argument t should be numeric")
    }
  }

  count_cols <- sapply(list_models,function(x) ncol(x))
  if (sum(count_cols != 2) > 0) {

    stop("Error: Each dataframe in the list should consist of only 2 columns.
          The first column with class labels (0/1) and the second column indicating the predicted probability")

  }

  check_col_type <- sum(unlist(lapply(list_models,function(x)
    lapply(x, function(y) !is.numeric(y)))))

  if (check_col_type > 0) {

    stop("Error: All columns in each dataframe should be numeric")

  }

  check_class_lables <- unlist(lapply(list_models, function(x) unique(x[,1])))

  if (sum (!(check_class_lables == 0 | check_class_lables == 1)) > 0) {

    stop("Error: Class labels should be either 0 or 1")

  }

  check_pred_prob <- unlist(lapply(list_models,function(x) range(x[,2])))

  if (sum (!(check_pred_prob >= 0 & check_pred_prob <= 1)) > 0) {

    stop("Error: Predicted probability should be between zero and one")

  }

}

obs_exp_summary <- function(x) {
    cut_points <- obs <- pred <- NULL
    colnames(x) <- c('obs','pred','cut_points')
    obs_exp <- x %>% group_by(cut_points) %>% summarise(obs_zero = length(obs[obs == 0]), obs_one = length(obs[obs ==
        1]), exp_zero = (1 - mean(pred)) * n(), exp_one = mean(pred) * n())
    obs_exp

}


hl_function <- function(x, g, sample_size_concord = NULL) {

    cut_points <- obs <- pred <- NULL
    x <- x %>% mutate(cut_points = cut(x[, 2], breaks = unique(quantile(x[, 2], seq(0, 1, 1/g))), include.lowest = T))

    obs_exp <- obs_exp_summary(x)
    return(list(obs_exp))
}

hl_test <- function(x, g, sample_size_concord = NULL) {

    hl_df <- hl_function(x, g)
    hl_df <- as.data.frame(hl_df)
    colnames(hl_df) <- gsub("hl_df.", "", colnames(hl_df))

    chi_square <- sum(with(hl_df, (obs_zero - exp_zero)^2/exp_zero + (obs_one - exp_one)^2/exp_one))
    p_value <- 1 - pchisq(chi_square, g - 2)
    hl_out <- data.frame(chi_square, p_value)
    return(list(hl_out))

}

calib_function <- function(x, g, sample_size_concord = NULL) {

    cut_points <- obs <-  NULL

    colnames(x) <- c('obs','pred')
    cut_points_mid <- data.frame(bin_mid = round(((seq(0, 1, 1/g) + lag(seq(0, 1, 1/g)))/2)[-1], 2))
    obs_prob <- x %>% mutate(cut_points = cut(x[, 2], breaks = seq(0, 1, 1/g), include.lowest = T)) %>%
        group_by(cut_points) %>% summarise(obs_rate = sum(obs)/n())
    mid_point <- function(x) {
        round(mean(as.numeric(unlist(regmatches(x, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", x))))), 2)
    }


    calib_df <- data.frame(bin_mid = sapply(1:nrow(obs_prob), function(i) mid_point(obs_prob$cut_points[i])),
        obs_prob) %>% select(-cut_points)
    calib_df <- merge(cut_points_mid, calib_df[c("bin_mid", "obs_rate")], by.x = "bin_mid", by.y = "bin_mid",
        all.x = T)
    calib_df[is.na(calib_df)] <- 0
    rownames(calib_df) <- NULL

    return(list(calib_df))

}

lift_function <- function(x, g, sample_size_concord = NULL) {

    obs_zero <- obs_one <- NULL

    if (g < length(unique(x[, 2]))) {

        x <- x %>% mutate(cut_points = cut(x[, 2], breaks = unique(quantile(x[, 2], seq(0, 1, 1/g))), include.lowest = T))

    } else {

        x <- x %>% mutate(cut_points = as.factor(round(x[, 2], 2)))

    }

    obs_exp <- as.data.frame(obs_exp_summary(x))
    obs_exp <- obs_exp[nrow(obs_exp):1, ] %>% mutate(Total = obs_zero + obs_one)
    cum_capture_1 <- with(obs_exp, round((cumsum(obs_one)/(sum(obs_one))) * 100, 2))
    fpr <- with(obs_exp, round((cumsum(obs_zero)/(sum(obs_zero))) * 100, 2))
    lift_index <- with(obs_exp, cum_capture_1/(cumsum(Total)/sum(Total)))
    lift_df <- data.frame(obs_exp["cut_points"], cum_capture_1, fpr, lift_index, KS = cum_capture_1 - fpr)
    rownames(lift_df) = NULL
    return(list(lift_df))

}

conc_disc <- function(x, g, sample_size_concord) {

    obs <- pred <- p_one <- p_zero <- compare <- Count <- NULL

    colnames(x) <- c('obs','pred')

    if (nrow(x) > 5000) {
        x <- x[sample(nrow(x), sample_size_concord), ]
    }

    df_one <- unlist(x %>% filter(obs == 1) %>% select(pred))
    df_zero <- unlist(x %>% filter(obs == 0) %>% select(pred))

    con_dis <- expand.grid(df_one, df_zero)
    colnames(con_dis)[1:2] <- c("p_one", "p_zero")
    con_dis <- con_dis %>% mutate(compare = c("Lesser than", "tied", "Greater than")[sign(p_one - p_zero) +
        2])

    con_dis <- con_dis %>% group_by(compare) %>% summarize(Count = n())
    con_dis <- con_dis %>% mutate(Perct = paste(round((Count/sum(Count)) * 100, 2), "%", sep = ""))
    count <- as.numeric(con_dis$Count)
    c_stat <- paste(round(((count[1] + 0.5 * count[3])/sum(count)) * 100, 2), "%", "")
    c_stat_df <- data.frame(label = "C-statistic", val = c_stat)
    con_dis <- con_dis[-2]

    colnames(c_stat_df) <- colnames(con_dis)

    con_dis <- rbind(con_dis, c_stat_df)
    return(list(con_dis))

}


combine_df <- function(list_df, index) {

    comb_df <- do.call(rbind, sapply(list_df, "[[", index))
    rep_nos <- sapply(sapply(list_df, "[[", index), function(x) nrow(x))
    comb_df <- data.frame(Model = unlist(lapply(seq_along(rep_nos), function(i) rep(paste("Model", i), rep_nos[i]))),
        comb_df)

    comb_df
}


"Plots"

plot_HL <- function(df) {

    Model <- Expected <- Value <- obs_one <- exp_one <- bins <- NULL

    df$bins <- unlist(lapply(1:length(unique(df$Model)), function(i) paste("Bin", seq(1, nrow(filter(df,
        Model == paste("Model", i)))))))
    df$bins <- factor(df$bins, levels = unique(df$bins))
    df <- df %>% gather(Expected, Value, obs_one, exp_one)
    g <- ggplot(df, aes(x = bins, y = Value, fill = Expected))
    g <- g + geom_bar(stat = "identity", position = "dodge") + facet_wrap(~Model)
    g
}

plot_calib <- function(df) {

    bin_mid <- obs_rate <- Model <- NULL
    calib_plot <- ggplot(df, aes(bin_mid, obs_rate, colour = Model)) + geom_line(size = 1) + geom_point(size = 3)
    calib_plot <- calib_plot + geom_abline(intercept = 0, slope = 1)
    calib_plot
}

plot_lift <- function(df) {

    Model <- bins <- lift_index <- NULL

    df$bins <- unlist(lapply(1:length(unique(df$Model)), function(i) paste("Bin", seq(1, nrow(filter(df,
        Model == paste("Model", i)))))))
    df$bins <- factor(df$bins, levels = unique(df$bins))
    g <- ggplot(df, aes(x = bins, y = lift_index, group = Model, colour = Model))
    g <- g + geom_line(size=1) + geom_point(size=3)
    g
}

plot_condis <- function(df) {

    obs <- pred <- NULL

    for (i in 1:length(df)) {
      colnames(df[[i]]) <- c('obs','pred')
    }

    comb_df <- do.call(rbind, df)
    colnames(comb_df) <- c('obs','pred')
    reps_no <- sapply(df, function(x) nrow(x))
    comb_df <- data.frame(Model = unlist(lapply(seq_along(reps_no), function(i) rep(paste("Model", i), reps_no[i]))),
        comb_df)

    g <- ggplot(comb_df, aes(x = pred, fill = as.factor(obs))) + geom_density(alpha = 0.5) + facet_wrap(~Model)
    g

}


############## Confusion ####################################

conf_mat <- function(x, t) {


    x <- x %>% mutate(pred_prob = as.factor(ifelse(x[, 2] >= t, "Pred-1", "Pred-0")))
    x$pred_prob <- factor(x = x$pred_prob, levels = c("Pred-0", "Pred-1"))
    conf_mat <- table(x[, 3], x[, 1])
    conf_mat

}

conf_mat_metrics <- function(x, t) {

    matrix <- conf_mat(x, t)
    Acc <- (matrix[1, 1] + matrix[2, 2])/sum(matrix)
    Acc <- round(Acc * 100, 2)
    TPR <- matrix[2, 2]/sum(matrix[, 2])
    TPR <- round(TPR * 100, 2)
    FPR <- matrix[2, 1]/sum(matrix[, 1])
    FPR <- round(FPR * 100, 2)
    Prec <- matrix[2, 2]/sum(matrix[2, ])
    Prec <- round(Prec * 100, 2)
    output <- data.frame(Threshold = t, Acc = Acc, TPR = TPR, FPR = FPR, Prec = Prec)
    output

}

conf_range <- function(x, reps, all.unique = F) {

    if (all.unique == T) {
      prob_values <- sort(unique(x[,2]))
    } else {
      prob_values <- seq(0, 1, 1/reps)
    }

    out <- lapply(seq_along(prob_values), function(i) conf_mat_metrics(x, prob_values[i]))
    out <- do.call(rbind, out)
    out

}

