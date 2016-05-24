
#' Metier Editing GUI
#' 
#' The \code{gui_lb_met_edi} function implements the graphical user interface for the
#'  editing of a LogBook Metier file
#' 
#' In this gui, with a LogBook Metier file, the user can change the simple discovery
#'  ordering of the clusters to a corresponding international metier code.
#'
#' @return This function does not return a value. 
#' After the execution, the edited names will be added to the LogBook Metier file.
#' 
#' @usage gui_lb_met_edi()
#' 
#' @export gui_lb_met_edi
#'

gui_lb_met_edi <- function(){
  
  lb_CLA <- log_Cla$new()
  
  cla_file <- gfile(text = "Select Metier file",
                    type = "open",
                    filter = list("Metier data" = list(patterns = c("*.rData"))))
  
  lb_CLA <- readRDS(cla_file)
  
  # temp_med <- as.data.frame(lb_CLA$data$medoids[, -which(apply(lb_CLA$data$medoids ,2,sum)==0)])
  temp_med <- as.data.frame(lb_CLA$data$medoids)
  num_clu <- nrow(temp_med)
  
  if(length(lb_CLA$options) == 1)
  {
    temp_med <- cbind(1:num_clu, 1:num_clu, temp_med)
  }else{
    temp_med <- cbind(1:num_clu, lb_CLA$options, temp_med)
  }
  colnames(temp_med)[1] <- "Cluster"
  colnames(temp_med)[2] <- "Metier"
  col_zer <- as.numeric(which(apply(temp_med[,3:ncol(temp_med)], 2, sum) == 0))+2
  temp_med <- temp_med[,-col_zer]
  fao_spe <- sub("FAO_", "", colnames(temp_med))
  
  sta_met <- read.table(file = system.file("extdata/EU_CODES_MET.csv", package="vmsbase"),
                        header = FALSE)
  
  met_edi_win <- gwindow("Metier Editing Tool", 
                         width = 700, height = 500,
                         horizontal= FALSE,
                         visible = FALSE)
  
  big_g <- ggroup(horizontal = FALSE,
                  use.scrollwindow = TRUE,
                  container = met_edi_win,
                  expand = TRUE)
  
  
  
  g_top <- ggroup(horizontal = FALSE, container = big_g, expand = TRUE)
  g_bot <- ggroup(horizontal = FALSE, container = big_g)
  g_bot_t <- ggroup(horizontal = TRUE, container = g_bot)
  g_bot_u <- ggroup(horizontal = FALSE, container = g_bot)
  addSpring(g_bot_t)
  
  tab_dat <- gtable(items = temp_med,
                    name = "Metier data",
                    filter.column = NULL,
                    expand = TRUE,
                    container = g_top)
  
  #font(tab_dat) <- list(size = 14)
  colnames(tab_dat) <- fao_spe
  
  e1 <- environment()
  gbutton(text = "\n\tLoad Metier List\t\t\n", container = g_bot_t, handler = function(h,...){
    
    cus_met_file <- gfile(text = "Select Custom Metier List file",
                          type = "open",
                          filter = list("Metier Names List" = list(patterns = c("*.csv"))))
    #new_dro <- vector(mode = "integer", length = num_clu)
    sta_met <- read.table(file = cus_met_file,
                          header = FALSE)
    delete(g_bot, g_bot_u)
    g_bot_u <<- ggroup(horizontal = FALSE, container = g_bot)
    for(i in 1:num_clu)
    {
      if(((i-1) %% 4) == 0 )
      {
        new_gr <- paste("g_met_", i, sep = "")
        assign(new_gr, ggroup(horizontal = TRUE, container = g_bot_u), envir = e1)
        addSpring(get(new_gr))
      }
      glabel(text = paste("Group ", i, sep = ""), container = get(new_gr))
      new_dro <- paste("d_met_", i, sep = "")
      #     gr_met_nam[i] <- new_dro
      assign(new_dro, 
             gdroplist(as.character(sta_met[,1]), 
                       horizontal = TRUE, 
                       container = get(new_gr)), envir = e1)
      
      if(length(lb_CLA$options) != 1)
      {
        aho <- get(new_dro, envir = e1)
        svalue(aho) <- lb_CLA$options[i]
      }
      addSpring(get(new_gr))
    }
  })
  
  addSpring(g_bot_t)  
  
  new_name <- array(1:nrow(temp_med))
  
  gbutton(text = "\n\tAssign New Names\t\t\n", container = g_bot_t, handler = function(h,...){
    #     for(k in 1:num_clu)
    #     {
    #       new_name[k] <<- as.character(svalue(get(paste("d_met_", k, sep = ""), envir = e1)))
    #       temp_med[k,2] <<- new_name[k]
    #       tab_dat[] <- temp_med
    #     }
    th_nene <- mget(paste("d_met_", 1:num_clu, sep = ""), envir = e1)
    new_name <<- as.character(lapply(th_nene, FUN = svalue))
    temp_med[,2] <<- new_name
    tab_dat[] <- temp_med
  })
  addSpring(g_bot_t)  
  
  gbutton(text = "\n\tSave Editing\t\t\n", container = g_bot_t, handler = function(h,...){
    #     new_name <<- as.character(lapply(mget(paste("d_met_", 1:num_clu, sep = ""), envir = e1), FUN = svalue))
    #     
    lb_CLA$options <- new_name
    
    res_file <- gfile(text = "Save Metier Names Editing",
                      type = "save")
    
    saveRDS(lb_CLA, file = paste(res_file, "_ed.rData", sep = ""))
    
    gconfirm("Metier Name Editing saved",
             title = "Confirm",
             parent = met_edi_win,
             handler = function(h,...){dispose(met_edi_win)})
    
  })
  addSpring(g_bot_t)  
  
  #new_dro <- vector(mode = "integer", length = num_clu)
  #   gr_met_nam <- character(length = num_clu)
  if(is.numeric(temp_med[,2]) == TRUE){
    for(i in 1:num_clu)
    {
      if(((i-1) %% 4) == 0 )
      {
        new_gr <- paste("g_met_", i, sep = "")
        assign(new_gr, ggroup(horizontal = TRUE, container = g_bot_u), envir = e1)
        addSpring(get(new_gr))
      }
      glabel(text = paste("Group ", i, sep = ""), container = get(new_gr))
      new_dro <- paste("d_met_", i, sep = "")
      #     gr_met_nam[i] <- new_dro
      assign(new_dro, 
             gdroplist(as.character(sta_met[,1]), 
                       horizontal = TRUE, 
                       container = get(new_gr)), envir = e1)
      if(length(lb_CLA$options) != 1)
      {
        aho <- get(new_dro, envir = e1)
        svalue(aho) <- lb_CLA$options[i]
      }
      addSpring(get(new_gr))
    }
  }else{
    for(i in 1:num_clu)
    {
      if(((i-1) %% 4) == 0 )
      {
        new_gr <- paste("g_met_", i, sep = "")
        assign(new_gr, ggroup(horizontal = TRUE, container = g_bot_u), envir = e1)
        addSpring(get(new_gr))
      }
      glabel(text = paste("Group ", i, sep = ""), container = get(new_gr))
      new_dro <- paste("d_met_", i, sep = "")
      #     gr_met_nam[i] <- new_dro
      assign(new_dro, 
             gdroplist(as.character(unique(temp_med[,2])), 
                       horizontal = TRUE, 
                       container = get(new_gr)), envir = e1)
      if(length(lb_CLA$options) != 1)
      {
        aho <- get(new_dro, envir = e1)
        svalue(aho) <- lb_CLA$options[i]
      }
      addSpring(get(new_gr))
    }
  }
  visible(met_edi_win) <- TRUE
}