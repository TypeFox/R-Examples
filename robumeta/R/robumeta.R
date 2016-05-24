#-------------------------------------------------------------------------------
forest.robu <- function(x, es.lab, study.lab, ...){
  
  if (paste(x$ml[3]) != 1){
    stop("Requires an intercept-only model.")
  }  
  ellipsis        <- lapply(substitute(list(...))[-1L], deparse)
  n_user_cols     <- length(ellipsis) # num. of additional columns
  reg_table       <- x$reg_table    
  N               <- x$N # num. of studies
  M               <- x$M # num. of clusters 
  n_rows          <- M + (2 * N) + 4
  data            <- as.data.frame(x$data) 
  data.full       <- as.data.frame(x$data.full) 
  data$orig.study <- as.factor(x$study_orig_id)
  data            <- data[order(data$orig.study),]            
  data$r.weights  <- data.full$r.weights
  data$effect.size  <- data.full$effect.size
  data$var.eff.size <- data.full$var.eff.size
  data$study.num  <- data.full$study 
  add_col_titles  <- as.list(names(ellipsis)) # user supplied titles
  add_col_values  <- as.list(data[, unlist(ellipsis, use.names = FALSE)]) 
  id_col_title    <- "Studies" 
  id_col_study_values <- unique(data[,study.lab]) 
  id_col_es_values    <- as.character(data[,es.lab]) 
  data$obs_num    <- seq(1, M)
  data$study_num  <- data$study.num 
  data$es_rows    <- as.numeric(data$obs_num + (2 * data$study_num) + 1) 
  data$study_rows <- as.numeric(ave(data$es_rows, data$study_num, FUN = min)- 1)
  es_rows         <- data$es_rows 
  study_rows      <- unique(data$study_rows) 
  total_row       <- max(n_rows)
  title_row       <- min(n_rows)
  data_col_values <- data[, c("r.weights", "effect.size", "var.eff.size")]
  data_col_values <- cbind(data_col_values, es_rows)
  grand.ES        <- reg_table$b.r
  grand.CI.L      <- reg_table$CI.L 
  grand.CI.U      <- reg_table$CI.U

  is.numeric_df   <- function(x) all(sapply(x, is.numeric))
  specify_decimal <- function(x, k) format(round(x, k), nsmall = k)
  
  makeTextGrob <- function(values, rows, just = "left", bold = FALSE){  
    if (is.numeric_df(values)) 
      values <- lapply(values, function (x) specify_decimal(x, 3))
    if (!bold){
      t <- lapply(values, function(x) textGrob(paste(x), x = 0, just = just))
    } else {
      t <- lapply(values, function(x) textGrob(paste(x), x = 0, just = just,
                                                 gp = gpar(fontface = "bold")))
    }
    return(list(values = t, rows = rows)) 
  }
  
  addTitleToGrob <- function(col, title){
    titleGrob  <- makeTextGrob(values = title, rows = 1, bold = TRUE)
    values     <- c(col$values, titleGrob$values)
    rows       <- c(col$rows, titleGrob$rows)
    return(list(values = values, rows = rows)) 
  }
  
  addGrobToGrob <- function(col1, col2){
    values <- c(col1$values, col2$values)
    rows   <- c(col1$rows, col2$rows)
    return(list(values = values, rows = rows)) 
  } 

  makeDataGrob <- function(x){  
    ES      <- x$effect.size
    size    <- x$r.weights / max(x$r.weights)
    CI.U    <- x$effect.size + (1.96 * sqrt(x$var.eff.size))
    CI.L    <- x$effect.size - (1.96 * sqrt(x$var.eff.size))
    type    <- rep("n", M)   
    rows    <- x$es_rows
    return(list(type = type, rows = rows, size = size, CI.L = CI.L, CI.U = CI.U, 
                ES = ES)) 
  }

  addSummaryToDataGrob <- function(x){ 
    type  <- c(x$type, "s")
    rows  <- c(x$rows, max(x$rows) + 2)
    size  <- as.numeric(x$size)
    size  <- x$size / max(x$size)
    ES    <- c(x$ES, grand.ES)
    CI.L  <- c(x$CI.L, grand.CI.L)
    CI.U  <- c(x$CI.U, grand.CI.U)
    min   <- floor(as.numeric(min(CI.L)))
    max   <- ceiling(as.numeric(max(CI.U)))
    range <- c(min, max)
    return(list(type = type, rows = rows, size = size, CI.L = CI.L, CI.U = CI.U, 
                ES = ES, min = min, max = max, range = range)) 
  }

  if (n_user_cols > 1){
    add_col <- lapply(add_col_values, function(x) makeTextGrob(x, es_rows))
    add_col <- Map(function(x, y) addTitleToGrob(x, y), add_col, add_col_titles)
  }
  
  if (n_user_cols == 1){
    add_col <- makeTextGrob(add_col_values, es_rows)
    add_col <- addTitleToGrob(add_col, add_col_titles)
  }

  id_col_study_grob <- makeTextGrob(id_col_study_values, study_rows, bold =TRUE)
  id_col_es_grob    <- makeTextGrob(id_col_es_values, es_rows)
  id_col            <- addGrobToGrob(id_col_study_grob, id_col_es_grob)
  id_col            <- addTitleToGrob(id_col, id_col_title) 
  data_col          <- makeDataGrob(data_col_values)
  data_col          <- addSummaryToDataGrob(data_col)

  drawLabelCol <- function(col, j) {
    for (i in 1:length(col$rows)) {
      pushViewport(viewport(layout.pos.row = col$rows[i], layout.pos.col = j))
      grid.draw(col$values[[i]])
      popViewport()
    }
  }

  drawNormalCI <- function(CI.L, ES, CI.U, size) {
    grid.rect(x = unit(ES, "native"), 
              width = unit(size, "snpc"), 
              height = unit(size, "snpc"),
              gp = gpar(fill = "black"))
    
    if (convertX(unit(CI.U, "native"), "npc", valueOnly = TRUE) > 1)
      grid.lines(x = unit(c(CI.L, 1), c("native", "npc")), 
                 y = .5,
                 arrow = arrow(length = unit(0.05, "inches")))
    else { 
      lineCol <- "black"
      grid.lines(x = unit(c(CI.L, CI.U), "native"), 
                 y = 0.5,
                 gp = gpar(col = lineCol))
    }
  }

  drawSummaryCI <- function(CI.L, ES, CI.U) {
    grid.polygon(x=unit(c(CI.L, ES, CI.U, ES), "native"),
                 y=unit(0.5 + c(0, 0.25, 0, -0.25), "npc"))
  }

  drawDataCol <- function(col, j) { # j = col_place
    pushViewport(viewport(layout.pos.col = j, xscale = col$range))
    grid.lines(x = unit(col$ES[length(col$ES)], "native"),
               y = unit(0:(n_rows - 2), "lines"), 
               gp = gpar(lty = "dashed"))
    grid.xaxis(gp=gpar(cex = 1))
    grid.text("Effect Size",                  
              y = unit(-3, "lines"), 
              x = unit(0.5, "npc"), 
              just = "centre", 
              gp = gpar(fontface = "bold"))
    popViewport()
    x = unit(0.5, "npc")
    
    for (i in 1:length(col$rows)) {
      pushViewport(viewport(layout.pos.row = col$rows[i], 
                            layout.pos.col = j,
                            xscale = col$range))
      if (col$type[i] == "n")
        drawNormalCI(col$CI.L[i], col$ES[i], col$CI.U[i], col$size[i])
      else
        drawSummaryCI(col$CI.L[i], col$ES[i], col$CI.U[i])
      popViewport()
    }
  }

  id_col_width   <- max(unit(rep(1, length(id_col$values)), "grobwidth", 
                           id_col$values))
  data_col_width <- unit(3, "inches")
  gap_col        <- unit(10, "mm")
  cols           <- unit.c(id_col_width, gap_col, data_col_width, gap_col)
  
  add_col_widths <- c()
  if (n_user_cols > 1){
    for (i in 1:n_user_cols) {
      add_col_widths[[i]] <-  max(unit(rep(1, length(add_col[[i]]$values)), 
                                     "grobwidth", add_col[[i]]$values))
      cols <- unit.c(cols, add_col_widths[[i]])
      cols <- unit.c(cols, gap_col)
    }
  }
  if (n_user_cols == 1){
      add_col_widths <-  max(unit(rep(1, length(add_col[1]$values)), 
                                     "grobwidth", add_col[1]$values))
      cols <- unit.c(cols, add_col_widths[1])
      cols <- unit.c(cols, gap_col)
  }

  pushViewport(viewport(layout = grid.layout(n_rows, (4 + (2 * n_user_cols)),
                        widths = cols,
                        heights = unit(c(1, rep(1, n_rows)), "lines"))))

  pushViewport(viewport(layout.pos.row = 1))
  
  grid.text("Forest Plot", 
            y = unit(+3, "lines"),
            just = "center", 
            gp = gpar(fontface = "bold"))
  
  grid.text(paste(x$model.lab1), 
            y = unit(+2, "lines"),
            just = "center", 
            gp = gpar(fontface = "italic"))
  
  popViewport()

  drawLabelCol(id_col, 1)
  
  if (n_user_cols > 1){
    for (i in 1:n_user_cols) {
      drawLabelCol(add_col[[i]], ((i * 2) + 3))
    }
  }
  if (n_user_cols == 1){
    for (i in 1:n_user_cols) {
      drawLabelCol(add_col, 5)
    }
  }

  drawDataCol(data_col, 3)
  popViewport()
}
#-------------------------------------------------------------------------------
sensitivity <- function(x) UseMethod("sensitivity")

  sensitivity.robu <- function(x){
    
    modelweights   <- x$modelweights
    user_weighting <- x$user_weighting
    
    if(modelweights == "HIER")  
      stop("Sensitivity analysis is not available for hierarchical effects.")
    
    if(user_weighting == TRUE)  
      stop("Sensitivity analysis is not available for user specified weights.")
        
        mod_info  <- x$mod_info
        p         <- x$p
        N         <- x$N
        Xreg      <- x$Xreg
        y         <- x$y
        X         <- x$X
        data.full <- x$data.full
        X.full    <- x$X.full
        k         <- data.full$k
        k_list    <- x$k_list
        ml        <- x$ml
        term1     <- mod_info$term1
        term2     <- mod_info$term2
        small     <- x$small
        labels    <- x$labels
        mod_label <- x$mod_label
        rho.test  <- seq(0, 1, .2)
    
        rho_labels   <- c(paste("Rho = ", seq(0, 1, .2), sep=""))
        var_labels   <- rep("", 2 * (p + 1))
        var_labels[seq(1, length(var_labels), by = 2)] <- labels
        var_labels   <- c(var_labels, "Tau.sq")
        
        col2_labels  <- rep("", 2 * (p + 1))
        col2_labels[seq(1, length(col2_labels), by = 2)] <- "Coefficient"
        col2_labels[seq(2, length(col2_labels), by = 2)] <- "Std. Error"
        col2_labels  <- c(col2_labels, "Estimate")
        sen          <- data.frame(cbind(var_labels, col2_labels))
        
        for (i in (1: length(rho.test))){
          
           tau.sq1             <- term1 + rho.test[i] * term2 
           tau.sq              <- ifelse(tau.sq1 < 0, 0, tau.sq1)
           data.full$r.weights <- 1 / (data.full$k * 
                                  (data.full$avg.var.eff.size + tau.sq))
           W.r.big             <- diag(data.full$r.weights)  # W
           W.r                 <- by(data.full$r.weights, data.full$study, 
                                     function(x) diag(x, nrow = length(x)))
           sumXWX.r            <- Reduce("+", Map(function(X, W) 
                                                  t(X) %*% W %*% X, 
                                                  X, W.r))
           sumXWy.r            <- Reduce("+", Map(function(X, W, y) 
                                                  t(X) %*% W %*% y, 
                                                  X, W.r, y))
           b.r                 <- solve(sumXWX.r) %*% sumXWy.r 
           data.full$pred.r    <- Xreg %*% b.r
           data.full$e.r       <- cbind(data.full$effect.size) - 
                                  data.full$pred.r
           data.full$e.r       <- as.numeric(data.full$e.r)
           sigma.hat.r         <- by(data.full$e.r, data.full$study, 
                                    function(x) tcrossprod(x))
           
              if (!small) { # Begin small = FALSE 
                
                 sumXWeeWX.r <- Reduce("+", Map(function(X, W, V) 
                                                t(X) %*% W %*% V %*% W %*% X, 
                                                X, W.r, sigma.hat.r))
                 VR.r        <- solve(sumXWX.r) %*% sumXWeeWX.r %*% 
                                solve(sumXWX.r)  
                 SE          <- sqrt(diag(VR.r)) * sqrt(N / (N - (p + 1)))

              } else { 
                
                Q             <- solve(sumXWX.r) #
                Q.list        <- rep(list(Q), N)
                H             <- Xreg %*% Q %*% t(Xreg) %*% W.r.big 
                ImH           <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) - H
                data.full$ImH <- cbind(ImH)
                ImHj          <- by(data.full$ImH, data.full$study, 
                                    function(x) as.matrix(x))
                dfS           <- c(rep(0, p + 1))
                diag_one      <- by(rep(1, nrow(X.full)), X.full$study, 
                                    function(x) diag(x, nrow = length(x)))
                ImHii         <- Map(function(X, Q, W, D) 
                                     D - X %*% Q %*% t(X) %*% W,
                                     X, Q.list, W.r, diag_one)
                eigenvec <- lapply(ImHii, function(x) eigen(x)$vectors) 
                eigenval <- lapply(ImHii, function(x) eigen(x)$values)
                I        <- ImHii
                A.MBB    <- Map(function (eigenvec, eigenval, k_list) 
                                eigenvec %*% diag(1/sqrt(eigenval), 
                                                  k_list, k_list) 
                                %*% t(eigenvec),
                                eigenvec, eigenval, k_list)
                A.MBB1    <- Map(function(K, A, I) 
                                 if (K > 1) A else matrix(sqrt(solve(I))), 
                                 k_list, A.MBB, I)
                A.MBB2    <- A.MBB 
                
                sumXWA.MBBeeA.MBBWX.r <- Reduce("+", Map(function(X,W,A,S) 
                                                         t(X) %*% W %*% A %*% 
                                                           S %*% A %*% W %*% X, 
                                                         X, W.r, A.MBB2, 
                                                           sigma.hat.r))
                data.full$ImH         <- ImH
                ImH                   <- lapply(split(data.full$ImH, 
                                                      data.full$study), 
                                                matrix, ncol=nrow(data.full))
                giTemp                <- Map(function(I, A, W, X, Q)
                                             t(I) %*% A %*% W %*% X %*% Q, 
                                             ImHj, A.MBB2, W.r, X, Q.list) 
                dfs                   <- c(rep(0, p + 1))
                
                for (i in 1:(p + 1)) { 
                   L      <- c(rep(0,p + 1))
                   L[i]   <- 1
                   Ll     <- rep(list(L), N)
                   gi     <- Map(function(G, L) G %*% cbind(L), giTemp, Ll)
                   G      <- Reduce("+", lapply(gi, function(x) tcrossprod(x)))
                   B      <- solve(sqrt(W.r.big) )%*% G %*% solve(sqrt(W.r.big))
                   e.val2 <- eigen(B)
                   dfs[i] <- sum(e.val2$values)^2 / sum(e.val2$values^2)
                }
                
               VR.MBB1 <- solve(sumXWX.r) %*% sumXWA.MBBeeA.MBBWX.r %*% 
                          solve(sumXWX.r)
               VR.r    <- VR.MBB1
               SE      <- sqrt(diag(VR.r))         
             } 
           
           vals <- c()
           temp_vals <- c()
           for (i in 1:(p + 1)) { 
              temp_vals <- c(b.r[i], SE[i])
              vals <- c(vals, temp_vals)
           }
           vals <- c(vals, tau.sq)
           vals <- format(vals, digits=3, justify="centre")
           sen <- cbind(sen, vals)
        }
    colnames(sen)     <- c(" ", " ", rho_labels)
    format(sen[,1], justify = "left")
    cat(mod_label, "\n")
    cat("Model:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
    cat(paste("Sensitivity Analysis"), "\n\n")
    print.data.frame(sen, quote = FALSE, row.names = FALSE, right = FALSE)
  }
# ------------------------------------------------------------------------------
group.mean <- function(var, grp) {
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(tapply(var, grp, mean, na.rm = TRUE)[grp])
}
#-------------------------------------------------------------------------------
group.center <- function(var, grp) {
  grp <- as.factor(grp)
  grp <- as.numeric(grp)
  var <- as.numeric(var)
  return(var - tapply(var, grp, mean, na.rm = TRUE)[grp])
}
#-------------------------------------------------------------------------------
print.robu <- function(x, digits = 3,...){

user_weighting <- x$user_weighting
modelweights   <- x$modelweights
mod_info       <- x$mod_info

output              <- x$reg_table
output              <- format(output, trim = TRUE, digits = digits, 
                              scientific = FALSE)
colnames(output)    <- c("", "Estimate","StdErr", "t-value", "dfs", "P(|t|>)", 
                             "95% CI.L","95% CI.U", "Sig")

if(!user_weighting){
  
  switch(modelweights,
         
    HIER = { # Begin HIER 
          
      cat(x$mod_label, "\n")
      cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
      cat(paste("Number of clusters ="), x$N, "\n")
      cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
          paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
          median(x$k), paste(", max ="), max(x$k),")\n")
      cat(paste("Omega.sq ="), mod_info$omega.sq, "\n")
      cat(paste("Tau.sq ="), mod_info$tau.sq, "\n\n")
      print(output, digits = 3)
      cat("---\n")
      cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
      cat("---\n")
      cat(x$mod_notice)
           
    }, # End HIER
         
    CORR = { # Begin CORR
           
      cat(x$mod_label, "\n")
      cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
      cat(paste("Number of studies ="), x$N, "\n")
      cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
          paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
          median(x$k), paste(", max ="), max(x$k),")\n")
      cat(paste("Rho ="), mod_info$rho, "\n")
      cat(paste("I.sq ="), mod_info$I.2, "\n")
      cat(paste("Tau.sq ="), mod_info$tau.sq, "\n\n")
      print(output)
      cat("---\n")
      cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
      cat("---\n")
      cat(x$mod_notice)      
             
    } # End CORR
    
  ) 
  
  } else { # Begin userweights
    
      cat(x$mod_label, "\n")
      cat("\nModel:",paste(x$ml[2]), paste(x$ml[[1]]), paste(x$ml[3]),"\n\n")
      cat(paste("Number of studies ="), x$N, "\n")
      cat(paste("Number of outcomes ="), sum(x$k), paste("(min ="), min(x$k),
          paste(", mean ="), format(mean(x$k), digits = 3), paste(", median ="), 
          median(x$k), paste(", max ="), max(x$k),")\n")
      print(output)
      cat("---\n")
      cat("Signif. codes: < .01 *** < .05 ** < .10 *\n")
      cat("---\n")
      cat(x$mod_notice) 
    
  } 
}
#-------------------------------------------------------------------------------
robu     <- function(formula, data, studynum,var.eff.size, userweights,
                     modelweights = c("CORR", "HIER"), rho = 0.8, 
                     small = TRUE, ...) {
  
  # Evaluate model weighting scheme.
  modelweights <- match.arg(modelweights)

  if (modelweights == "CORR" && rho > 1 | rho < 0)  
      stop ("Rho must be a value between 0 and 1.")
  
  if (missing(userweights)){
     user_weighting = FALSE 
  } else { 
     user_weighting = TRUE
  }

  cl                       <- match.call() # Full model call
  mf                       <- match.call(expand.dots = FALSE)
  ml                       <- mf[[2]] # Model formula 
  m                        <- match(c("formula", "data", "studynum", 
                                      "var.eff.size", "userweights"), names(mf))
  mf                       <- mf[c(1L, m)] 
  mf$drop.unused.levels    <- TRUE
  mf[[1L]]                 <- as.name("model.frame")
  mf                       <- eval(mf, parent.frame()) 
  
  if(!user_weighting){ 
    
    dframe                 <- data.frame(effect.size = mf[,1],
                                        model.matrix(formula, mf),
                                        studynum = mf[["(studynum)"]],
                                        var.eff.size = mf[["(var.eff.size)"]])
    
    X.full.names           <- names(dframe)[-match(c("effect.size", 
                                                     "studynum", 
                                                     "var.eff.size"), 
                                                   names(dframe))] 
    
  } else { # Begin userweights
    
    dframe                 <- data.frame(effect.size = mf[,1],
                                        model.matrix(formula, mf), 
                                        studynum = mf[["(studynum)"]], 
                                        var.eff.size = mf[["(var.eff.size)"]],
                                        userweights = mf[["(userweights)"]])
    
    X.full.names           <- names(dframe)[-match(c("effect.size", 
                                                     "studynum", 
                                                     "userweights", 
                                                     "var.eff.size"), 
                                                   names(dframe))] 
  } # End userweights
  study_orig_id            <- dframe$studynum
  dframe$study             <- as.factor(dframe$studynum)
  dframe$study             <- as.numeric(dframe$study)
  dframe                   <- dframe[order(dframe$study),]
  k_temp                   <- as.data.frame(unclass(rle(sort(dframe$study))))
  dframe$k                 <- k_temp[[1]][ match(dframe$study, k_temp[[2]])]
  dframe$avg.var.eff.size  <- ave(dframe$var.eff.size, dframe$study)
  dframe$sd.eff.size       <- sqrt(dframe$var.eff.size)
  
  switch(modelweights, 
         
    HIER = { # Begin HIER
      
        dframe$weights <- 1 / dframe$var.eff.size
      
    }, # End HIER
    
    CORR = { # Begin CORR
      
        dframe$weights <- 1 / (dframe$k * dframe$avg.var.eff.size)
      
    } # End CORR
    
  ) 
  
  X.full           <- dframe[c("study", X.full.names)]
  
  data.full.names  <- names(dframe)[-match(c("studynum",X.full.names), 
                                          names(dframe))] 
  data.full        <- dframe[c(data.full.names)]
  k                <- data.full[ !duplicated(data.full$study), ]$k
  k_list           <- as.list(k) 
  M                <- nrow(data.full) # Number of units in analysis
  p                <- ncol(X.full) - 2 # Number of (non-intercept) covariates 
  N                <- max(data.full$study) # Number of studies
  W                <- as.matrix(by(data.full$weights, data.full$study, 
                                   function(x) diag(x, nrow = length(x))))
  X                <- data.matrix(X.full)
  X                <- lapply(split(X[,2:(p + 2)], X[,1]), matrix, ncol = p + 1)
  y                <- by(data.full$effect.size, data.full$study, 
                         function(x) matrix(x))
  J                <- by(rep(1, nrow(X.full)), X.full$study, 
                         function(x) matrix(x, nrow = length(x), 
                                               ncol = length(x)))
  sigma            <- by(data.full$sd.eff.size, data.full$study, 
                         function(x) tcrossprod(x))
  vee              <- by(data.full$var.eff.size, data.full$study, 
                         function(x) diag(x, nrow = length(x)))
  SigmV            <- Map(function(sigma, V) 
                          sigma - V, sigma, vee)
  sumXWX           <- Reduce("+", Map(function(X, W) 
                                      t(X) %*% W %*% X, 
                                      X, W))
  sumXWy           <- Reduce("+", Map(function(X, W, y) 
                                      t(X) %*% W %*% y, 
                                      X, W, y))
  sumXWJWX         <- Reduce("+", Map(function(X, W, J) 
                                      t(X) %*% W %*% J %*% W %*% X, 
                                      X, W, J))
  sumXWVWX         <- Reduce("+", Map(function(X, W, V) 
                                      t(X) %*% W %*% V %*% W %*% X, 
                                      X, W, vee))
  sumXW.sig.m.v.WX <- Reduce("+", Map(function(X, W, V) 
                                      t(X) %*% W %*% V %*% W %*% X, 
                                      X, W, SigmV))
  
  switch(modelweights, 
    
    HIER = { # Begin HIER
        
        tr.sumJJ <- Reduce("+", Map(function(J) 
                                    sum(diag(J %*% J)), 
                                    J)) 
        sumXJX   <- Reduce("+", Map(function(X, J) 
                                    t(X) %*% J %*% X, 
                                    X, J))
        sumXWJJX <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% W %*% J %*% J %*% X, 
                                    X, W, J))
        sumXJJWX <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% J %*% J %*% W %*% X, 
                                    X, W, J))
        sumXWWX  <- Reduce("+", Map(function(X, W) 
                                    t(X) %*% W %*% W %*% X, 
                                    X, W))
        sumXJWX  <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% J %*% W %*% X, 
                                    X , W, J))
        sumXWJX  <- Reduce("+", Map(function(X, W, J) 
                                    t(X) %*% W %*% J %*% X, 
                                    X, W, J))
    } # End HIER
    
  ) 
  
  b              <- solve(sumXWX) %*% sumXWy 
  Xreg           <- as.matrix(X.full[-c(1)], dimnames = NULL)
  data.full$pred <- Xreg %*% b 
  data.full$e    <- data.full$effect.size - data.full$pred 
  
  if (!user_weighting) { 
    
  switch(modelweights, 
    
    HIER = { # Begin HIER
      
      # Sigma_aj = tau.sq * J_j + omega.sq * I_j + V_j 
      # Qe is sum of squares 1
      # Qe = Sigma(T'WT)-(Sigma(T'WX)(Sigma(X'WX))^-1(Sigma(X'WT)
      # where W = V^(-1) and V = data.full$var.eff.size
      # Also, Qe = (y-xb)' W (y-xb)
      sumV <- sum(data.full$var.eff.size)
      W    <- diag(1 / data.full$var.eff.size) 
      sumW <- sum(W)
      Qe   <- t(data.full$e) %*% W %*% data.full$e
     
      # Qa is sum of squares 2
      # Qa = sum(T-XB.hat)'J(T-XB.hat)
      # where B.hat = (X'WX)^-1(X'WT)
      # Also, Qa = (y-xb)'A (y-xb), A=diag(J)
      e      <- by(data.full$e, data.full$study, function(x) matrix(x))
      sumEJE <- Reduce("+", Map(function(e, J) t(e) %*% J %*% e, e, J))
      Qa     <- sumEJE 
      
      # MoM estimators for tau.sq and omega.sq can be written as
      # omega.sq.h = A2(Qa-C1)-A1(Qe-C2) / B1A2-B2A1
      # tau.sq.h = Qe-C2/A2 - omega.sq.h(B2/A2) where
      
      # Vi = (t(X)WX)^-1
      V.i    <- solve(sumXWX)
      
      # A1 = Sigma(kj^2) - tr(V*Sigma(kj*t(Xj)*Jj*Wj*Xj)) - 
      #                    tr(V*Sigma(kj*t(Xj)*Jj*Wj*Xj)) +
      #                    tr(V*[Sigma(t(Xj)*Jj*Xj)]*V*Sigma(t(Xj)*Wj*Jj*Wj*Xj))
      # B1 = Sigma(kj)   - tr(V Sigma(t(Xj)*Jj*Wj*Xj)) -
      #                    tr(V Sigma(t(Xj)*Wj*Jj*Xj)) +
      #                    tr(V*[Sigma(t(Xj)*Jj*Xj)]*V*Sigma(t(Xj)*Wj^2*Xj)) 
      # C1 = tr(W^-1)    - tr(V*Sigma(t(X)*Jj*Xj))
      
      A1    <- tr.sumJJ  - sum(diag(V.i %*% sumXJJWX)) - 
                           sum(diag(V.i %*% sumXWJJX)) + 
                           sum(diag(V.i %*% sumXJX %*% V.i %*% sumXWJWX))
      
      B1    <- length(data.full$study) -
                           sum(diag(V.i %*% sumXWJX)) -
                           sum(diag(V.i %*% sumXJWX)) +
                           sum(diag(V.i %*% sumXJX%*%V.i %*% sumXWWX))
      C1    <- sumV - sum(diag(V.i %*% sumXJX))
     
      # A2 = tr(W) - tr(V*Sigma(t(X)*Wj*Jj*Wj*Xj))
      # B2 = tr(W) - tr(V*Sigma(t(X)*Wj^2*Xj))
      # C2 = Sigma(kj-p)
      
      A2   <- sumW - sum(diag(V.i %*% sumXWJWX)) 
      B2   <- sumW - sum(diag(V.i %*% sumXWWX)) 
      C2   <- length(data.full$study) - (p + 1) 
      
      # MoM estimator for omega.sq.h = A2(Qa-C1)-A1(Qe-C2) / B1A2-B2A1
      # Estimate of between-studies-wthin-cluster variance component
      omega.sq1  <- ((Qa - C1) * A2 - (Qe - C2) * A1) / (B1 * A2 - B2 * A1)
      omega.sq   <- ifelse(omega.sq1 < 0, 0, omega.sq1)
      
      # MoM estimators for tau.sq: Qe-C2/A2 - omega.sq.h(B2/A2)
      # Estimate of between-clusters variance component 
      tau.sq1  <- ((Qe - C2) / A2) - omega.sq  * (B2 / A2)
      tau.sq   <- ifelse(tau.sq1 < 0, 0, tau.sq1) 

      # Approximate inverse variance weights
      data.full$r.weights <- (1 / (data.full$var.eff.size + tau.sq + omega.sq))
  
      # Model info list for hierarchical effects
      mod_info            <- list(omega.sq = omega.sq, tau.sq = tau.sq)

    }, # End HIER
         
    CORR = { # Begin CORR
      
      W       <- diag (data.full$weights) 
      sumW    <- sum(data.full$weights) # Sum (k.j*w.j)
      Qe      <- t(data.full$e) %*% W %*% data.full$e 
      
      # The following components (denom, termA, termB, term1, term2)
      # are used in the calculation of the estimate of the residual 
      # variance component tau.sq.hat. 
      # Note: The effect of correlation on the estimates occurs entirely 
      # through the rho*term2 component.
      
      denom   <- sumW - sum(diag(solve(sumXWX) %*% sumXWJWX)) 
      termA   <- sum(diag(solve(sumXWX) %*% sumXWVWX))
      termB   <- sum(diag(solve(sumXWX) %*% sumXW.sig.m.v.WX ))
      term1   <- (Qe - N + termA) / denom 
      term2   <- termB / denom 
      tau.sq1 <- term1 + rho * term2 
      tau.sq  <- ifelse(tau.sq1 < 0, 0, tau.sq1)
      df      <- N - termA - rho * (termB) 
      I.2.1   <- ((Qe - df) / Qe) * 100
      I.2     <- ifelse(I.2.1 < 0, 0, I.2.1)
      
      # Approximate inverse variance weights
      data.full$r.weights <- 1 / (data.full$k * 
                             (data.full$avg.var.eff.size + tau.sq))
      
      # Model info list for correlated effects
      mod_info            <- list(rho = rho, I.2 = I.2, tau.sq = tau.sq,
                                  term1 = term1, term2 = term2)
      
    } # End CORR 
    
  ) 
  
  } else { # Begin userweights
  
    data.full$r.weights <- data.full$userweights
    
    # Model info list for userweights
    mod_info            <- list(k = k, N = N, p = p, M = M)
    
  } # End userweights
  
  W.r.big          <- diag(data.full$r.weights)  # W
  W.r              <- by(data.full$r.weights, data.full$study, # Wj
                         function(x) diag(x, nrow = length(x)))
  sumXWX.r         <- Reduce("+", Map(function(X, W) 
                                      t(X) %*% W %*% X, 
                                      X, W.r))
  sumXWy.r         <- Reduce("+", Map(function(X, W, y) 
                                      t(X) %*% W %*% y, 
                                      X, W.r, y))
  b.r              <- solve(sumXWX.r) %*% sumXWy.r 
  data.full$pred.r <- Xreg %*% b.r
  data.full$e.r    <- cbind(data.full$effect.size) - data.full$pred.r
  data.full$e.r    <- as.numeric(data.full$e.r)
  sigma.hat.r      <- by(data.full$e.r, data.full$study, 
                         function(x) tcrossprod(x))
  
  if (!small) { # Begin small = FALSE
    
    sumXWeeWX.r  <- Reduce("+", Map(function(X, W, V) 
                                    t(X) %*% W %*% V %*% W %*% X, 
                                    X, W.r, sigma.hat.r))
    
    VR.r         <- solve(sumXWX.r) %*% sumXWeeWX.r %*% solve(sumXWX.r)  
    SE           <- sqrt(diag(VR.r)) * sqrt(N / (N - (p + 1)))
    t            <- b.r / SE
    dfs          <- N - (p + 1)
    prob         <- 2 * (1 - pt(abs(t), dfs))
    CI.L         <- b.r - qt(.975, dfs) * SE
    CI.U         <- b.r + qt(.975, dfs) * SE
    
  } else { # Begin small = TRUE

    Q             <- solve(sumXWX.r) # Q = (X'WX)^(-1)
    Q.list        <- rep(list(Q), N)
    H             <- Xreg %*% Q %*% t(Xreg) %*% W.r.big # H = X * Q * X' * W
    ImH           <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) - H
    data.full$ImH <- cbind(ImH)
    ImHj          <- by(data.full$ImH, data.full$study, 
                        function(x) as.matrix(x))
    diag_one      <- by(rep(1, M), X.full$study, 
                                   function(x) diag(x, nrow = length(x)))
    ImHii         <- Map(function(X, Q, W, D) 
                         D - X %*% Q %*% t(X) %*% W,
                         X, Q.list, W.r, diag_one)
    
    if (!user_weighting){
    
    switch(modelweights, 
      
      HIER = { # Begin HIER
        
         # inside = Wj^(-1/2) * (I-Hjj) * Wj^(-3/2)
         inside   <- Map(function(W, I) 
                         solve(sqrt(W)) %*% I %*% solve(sqrt(W)^3),
                         W.r, ImHii)
         I        <- inside
         eigenvec <- lapply(inside, function(x) eigen(x)$vectors) 
         eigenval <- lapply(inside, function(x) eigen(x)$values)
        
      }, # End HIER
      
      CORR = { # Begin CORR
      
         eigenvec <- lapply(ImHii, function(x) eigen(x)$vectors) 
         eigenval <- lapply(ImHii, function(x) eigen(x)$values)
         I        <- ImHii
        
      } # End CORR
           
    ) 
    
    } else { # Begin userweights

        V.big        <- diag(c(1), dim(Xreg)[1], dim(Xreg)[1]) %*% 
                        diag(data.full$avg.var.eff.size)
        V.big.list   <- rep(list(V.big), N)
        v.j          <- by(data.full$avg.var.eff.size, data.full$study, 
                           function(x) diag(x, nrow = length(x)))
        v.j.sqrt     <- lapply(v.j, function (x) sqrt(x))
        inside       <- Map(function(V, I) 
                            I  %*% V %*% t(I),
                            V.big.list, ImHj)
        eigenvec     <- lapply(inside, function(x) eigen(x)$vectors)
        eigenval     <- lapply(inside, function(x) eigen(x)$values)
        I            <- inside
        
      } # End userweights
    
    A.MBB  <- Map(function (eigenvec, eigenval, k_list) 
                  eigenvec %*% 
                    diag(1/sqrt(eigenval), k_list, k_list) %*% t(eigenvec),
                  eigenvec, eigenval, k_list)
    A.MBB1 <- Map(function(K, A, I) 
                  if (K > 1) A else matrix(sqrt(solve(I))), 
                  k_list, A.MBB, I)
    
    if (!user_weighting){
      
    switch(modelweights, 
           
     HIER = { # Begin HIER
        
        A.MBB2                <- Map(function(W, A) 
                                     solve(sqrt(W)) %*% A %*% solve(sqrt(W)),
                                     W.r, A.MBB1) 
        sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
                                     t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
                                     X, W.r, A.MBB2, sigma.hat.r)
      }, # End HIER
      
      CORR = { # Begin CORR
        
        A.MBB2                <- A.MBB1
        sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
                                     t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
                                     X, W.r, A.MBB2, sigma.hat.r)
      } # End CORR
           
    ) 
    
    } else { # Begin userweights
        
        A.MBB2                <- Map(function(V, A) 
                                     V %*% A,
                                     v.j.sqrt, A.MBB1) 
        sumXWA.MBBeeA.MBBWX.r <- Map(function(X,W,A,S) 
                                     t(X) %*% W %*% A %*% S %*% A %*% W %*%X, 
                                     X, W.r, A.MBB2, sigma.hat.r)
      } # End userweights
    
    sumXWA.MBBeeA.MBBWX.r <- Reduce("+", sumXWA.MBBeeA.MBBWX.r) 
    giTemp                <- Map(function(I, A, W, X, Q)
                                 t(I) %*% A %*% W %*% X %*% Q, 
                                 ImHj, A.MBB2, W.r, X, Q.list)
   
    dfs <- c(rep(0, p + 1))
    
      for (i in 1:(p+1)) { 
        
        L      <- c(rep(0,p+1))
        L[i]   <- 1
        Ll     <- rep(list(L), N)
        gi     <- Map(function(G, L) G %*% cbind(L), giTemp, Ll)
        G      <- Reduce("+", lapply(gi, function(x) tcrossprod(x)))
        
        if (!user_weighting){
          
          switch(modelweights, 
                 
            HIER = { # Begin HIER
              
              B <- solve(sqrt(W.r.big) )%*% G %*% solve(sqrt(W.r.big))
              
            }, # End HIER
            
            CORR = { # Begin CORR
              
              B <- solve(sqrt(W.r.big) )%*% G %*% solve(sqrt(W.r.big))
              
            } # End CORR
            
          ) 
          
        } else { # Begin userweights
          
            B <- solve(sqrt(V.big) )%*% G %*% solve(sqrt(V.big)) 
            
        } # End userweights
        
        e.val2 <- eigen(B)
        dfs[i] <- sum(e.val2$values)^2/sum(e.val2$values^2)
        
      } # End loop
    
    VR.MBB1 <- solve(sumXWX.r) %*% sumXWA.MBBeeA.MBBWX.r %*% solve(sumXWX.r)
    VR.r    <- VR.MBB1
    SE      <- sqrt(diag(VR.r))
    t       <- b.r / SE
    prob    <- 2 * (1 - pt(abs(t), df = dfs)) 
    CI.L    <- b.r - qt(.975, dfs) * SE
    CI.U    <- b.r + qt(.975, dfs) * SE
    
  } # End small = TRUE
        
    reg_table           <- data.frame(cbind(b.r, SE, t, dfs, prob, CI.L, CI.U))
    names(X.full)[2]    <- "intercept"
    labels              <- c(colnames(X.full[2:length(X.full)]))
    sig                 <- ifelse(prob < .01, "***", 
                           ifelse(prob > .01 & prob < .05, "**",
                           ifelse(prob > .05 & prob < .10, "*", "")))
    reg_table           <- cbind(labels, reg_table, sig)
    colnames(reg_table) <- c("labels", "b.r", "SE", "t", "dfs", "prob", "CI.U", 
                             "CI.L", "sig")
  
   if (!small) { # Begin small = FALSE
 
      mod_label_sm   <- ""
      mod_notice     <- ""
      
   } else { # Begin small = TRUE
     
      mod_label_sm  <- "with Small-Sample Corrections"            
      mod_notice    <- "Note: If df < 4, do not trust the results"
    
   } # End small = TRUE
       
  
   if (!user_weighting) {
  
     switch(modelweights,
           
       HIER = { # Begin HIER
             
          mod_label <- c("RVE: Hierarchical Effects Model", mod_label_sm)

        }, # End HIER
         
        CORR = { # Begin CORR
             
          mod_label <- c("RVE: Correlated Effects Model", mod_label_sm)
             
        } # End CORR
           
     ) 
         
   } else { # Begin userweights

     mod_label <- c("RVE: User Specified Weights", mod_label_sm)
         
   } # End userweights

  res <- list(data.full = data.full, X.full = X.full, reg_table = reg_table, 
              mod_label = mod_label, mod_notice = mod_notice, modelweights = 
              modelweights, mod_info = mod_info, user_weighting = 
              user_weighting, ml = ml, cl = cl, N = N, M = M, k = k, 
              k_list = k_list, p = p, X = X, y = y, Xreg = Xreg, b.r = b.r, 
              VR.r = VR.r, dfs = dfs, small = small, data = data, labels = 
              labels, study_orig_id = study_orig_id)
               
  class(res) <- "robu"
  res
}
