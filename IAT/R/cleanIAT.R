#' Clean IAT data using the updated D-Scoring algorithm
#' 
#' Transform a dataframe with trial latencies (stored as one line per trial) 
#' for a standard format IAT (7 blocks) into a one line summary per subject of the IAT effect using
#' the standard scoring algorithm recommended in Greenwald, Nosek, & Banaji (2003). The goal
#' is to prepare IAT data for subsequent analysis. However, this does not relieve the researcher 
#' from making conceptual decisions about how best to analyze IAT data. There are decisions to 
#' make about how the function is applied, and the function does not remove participants.
#' All subject exclusions must be made deliberately by the researcher. Note that the output of this function
#' is identical to that of the standard SAS macro (link in reference) for all meaningful columns.
#' 
#' @param my_data The raw dataframe to be used
#' @param block_name A string of the variable name for the blocks
#' @param trial_blocks A vector of the four essential blocks in the seven-block IAT (i.e., B3, B4, B6, and B7).
#' @param session_id A string of the variable name identifying each unique participant.
#' @param trial_latency A string of the variable name for the latency of each trial.
#' @param trial_error A string of the variable name identifying whether a trial was an error or not, where 1 indicates an error.
#' @param v_error If 1 (current standard), then means are calculated for the entire set of latencies. If 2, error latencies will be replaced by the block mean + 600ms
#' @param v_extreme If 1, then no extreme value treatment. If 2 (current standard), delete trial latencies < 400ms
#' @param v_std If 1 (current standard), block SD is performed including error trials (corrected or not). If 2, block SD is performed on correct responses only
#' @return Outputs a dataframe that must be saved to an object. The variable IAT is the calculated D-Score for each individual. SUBEXCL notes
#' any exclusion criteria, with 0 being inclusion data, 1 for exclusion due to fast response, and 2 for exclusion due to missing blocks. C indicates
#' standard deviation for combined blocks (correct trial only), while A indicates standard deviations for combined blocks (all trials). M (mean), E (percent error),
#' N (number of trials used), and F (percent fast responses), are reported for each block included in the original dataframe. 
#' @references \href{http://faculty.washington.edu/agg/pdf/GB&N.JPSP.2003.pdf}{Understanding and Using the Implicit Association Test: I. An Improved Scoring Algorithm}
#' @references \href{http://projectimplicit.net/nosek/papers/scoringalgorithm.sas.txt}{IAT SAS macro}
#' @examples
#' # Get Ps who recieve Math-Male sorting task in first blocks
#' 
#' cong_first <- IATData[IATData$isCongruentFirst == 1, ]
#' 
#' dscore_first <- cleanIAT(my_data = cong_first,
#'                          block_name = "BLOCK_NAME_S",
#'                          trial_blocks = c("BLOCK2", "BLOCK3", "BLOCK5", "BLOCK6"),
#'                          session_id = "SESSION_ID",
#'                          trial_latency = "TRIAL_LATENCY",
#'                          trial_error = "TRIAL_ERROR",
#'                          v_error = 1, v_extreme = 2, v_std = 1)
#' 
#' # Get Ps who recieve Math-Female sorting task in first blocks
#' 
#' cong_second <- IATData[IATData$isCongruentFirst == 0, ]
#' 
#' dscore_second <- cleanIAT(my_data = cong_second,
#'                           block_name = "BLOCK_NAME_S",
#'                           trial_blocks = c("BLOCK2", "BLOCK3", "BLOCK5", "BLOCK6"),
#'                           session_id = "SESSION_ID",
#'                           trial_latency = "TRIAL_LATENCY",
#'                           trial_error = "TRIAL_ERROR",
#'                           v_error = 1, v_extreme = 2, v_std = 1)
#' 
#' d_score <- rbind(dscore_first, dscore_second)
#'                                                       
#' @import dplyr lazyeval
#' @importFrom stats reshape sd
#' @export

cleanIAT <- function(my_data, block_name, trial_blocks, session_id, trial_latency, trial_error, v_error, v_extreme, v_std){
  
  # Step 1: Include data from B3, B4, B6, B7
  # Step 1 has been removed so all data is at least partially analyzed
  
  # Step 2a: Eliminate trial latencies > 10,000ms and < 0ms
  
  my_tbl <- my_data %>%
    tbl_df() %>%
    group_by_(session_id) %>%
    mutate(SUBEXCL = 0) %>%
    filter_(interp(~ x < 10000 & x >= 0, x = as.name(trial_latency)))
  
  # Step 2b: Mark subjects for whom more than 10% of trials have latencies < 300ms with a 1
  
  fast_tbl <- my_tbl %>%
    filter_(interp(~ x %in% trial_blocks, x = as.name(block_name))) %>%
    group_by_(session_id, block_name) %>%
    summarise_(FAST = interp(~ sum(x < 300)/length(x), x = as.name(trial_latency))) %>%
    summarise_(FASTM = interp(~ mean(y), y = as.name("FAST")))
              
  my_tbl$SUBEXCL[data.frame(my_tbl)[, session_id] %in% 
                   data.frame(filter_(fast_tbl, interp(~ z > 0.10, z = as.name("FASTM"))))[, session_id]] <- 1
  
  # Step 2c: Mark subjects with any missing blocks with a 2
  
  num_blocks <- my_tbl %>%
    summarise_(total = interp(~ sum(unique(x) %in% trial_blocks), x = as.name(block_name))) %>%
    filter_(interp(~ x < 4, x = as.name("total")))
  
  my_tbl$SUBEXCL[data.frame(my_tbl)[, session_id] %in% data.frame(num_blocks)[, session_id]] <- 2

  # Step 3: Use all trials (in the conventional algorithm) the first two trials of each
  # block would be dropped here
  
  # Step 4: No extreme value treatment <or> delete trial with latencies < 400ms
  
  if(v_extreme == 1){
    
    my_tbl_sub <- my_tbl %>%
      group_by_(session_id, block_name)
  
  } else if(v_extreme == 2){
    
    my_tbl_sub <- my_tbl %>%
      filter_(interp(~ x >= 400, x = as.name(trial_latency))) %>%
      group_by_(session_id, block_name)
    
  } else stop("Please enter a v_extreme value of 1 or 2.")
  
  # Step 5: Compute mean of correct latencies for each block
  # If v_error = 1 then means and SDs will be calculated for the entire set of latencies
  # IF v_error = 2 then the algorithm will replace error trial latencies with blockmean + 600
  # (where blockmean is mean of correct responses only)
  
  block_speed <- my_tbl %>%
    group_by_(session_id, block_name) %>%
    summarise_(F = interp(~ sum(y < 300)/length(y), y = as.name(trial_latency)))
  
  block_errors <- my_tbl_sub %>%
    group_by_(session_id, block_name) %>%
    summarise_(E = interp(~ sum(x)/length(x), x = as.name(trial_error)))
  
  block_extras <- full_join(block_errors, block_speed, by = c(session_id, block_name))
  
  if(v_error == 1){
    
    block_means <- my_tbl_sub %>%
      summarise_(SUBEXCL = interp(~ unique(z), z = as.name("SUBEXCL")),
                 M = interp(~ mean(x), x = as.name(trial_latency)),
                 N = interp(~ length(x), x = as.name(trial_latency))) %>%
      full_join(block_extras, by = c(session_id, block_name)) %>%
      full_join(fast_tbl, by = session_id)
    
  } else if(v_error == 2){
    
    my_tbl_sub <- my_tbl_sub %>%
      summarise_(new_latency = interp(~ mean(x[y == 0]) + 600,
                                      x = as.name(trial_latency),
                                      y = as.name(trial_error))) %>%
      full_join(my_tbl_sub, by = c(session_id, block_name))
    
    my_tbl_sub[, trial_latency] <- ifelse(data.frame(my_tbl_sub)[, trial_error] == 0,
                                          data.frame(my_tbl_sub)[, trial_latency],
                                          data.frame(my_tbl_sub)[, "new_latency"])
    
    my_tbl_sub$new_latency <- NULL
    
    block_means <- my_tbl_sub %>%
      group_by_(session_id, block_name) %>%
      summarise_(SUBEXCL = interp(~ unique(z), z = as.name("SUBEXCL")),
                 M = interp(~ mean(x), x = as.name(trial_latency)),
                 N = interp(~ length(x), x = as.name(trial_latency))) %>%
      full_join(block_extras, by = c(session_id, block_name)) %>%
      full_join(fast_tbl, by = session_id)
    
  } else stop("Please pick a value of 1 or 2 for v_error.")
  
  # Step 6: Compute pooled SD for B3 & B6, and separately for B4 & B7.
  # If v_std is 1, the block SD is performed including error trials (corrected or not)
  # If v_std is 2, the block SD is performed on correct responses only
  
  my_tbl_sds <- my_tbl_sub %>%
    filter_(interp(~ x %in% trial_blocks, x = as.name(block_name))) %>%
    mutate_(block_pairs = interp(~ ifelse(x %in% trial_blocks[c(1, 3)], "S1", "S2"),
                                 x = as.name(block_name))) %>%
    group_by_(session_id, "block_pairs") %>%
    summarise_(A = interp(~ sd(x[y == 0 | y == 1]),
                          x = as.name(trial_latency),
                          y = as.name(trial_error)),
               C = interp(~ sd(x[y == 0]),
                          x = as.name(trial_latency),
                          y = as.name(trial_error))) %>%
    data.frame() %>%
    reshape(v.names = c("A", "C"),
            idvar = session_id,
            timevar = "block_pairs",
            direction = "wide",
            sep = "") %>%
    full_join(suppressWarnings(reshape(data.frame(block_means),
                      v.names = c("M", "N", "E", "F"),
                      idvar = session_id,
                      timevar = block_name,
                      direction = "wide",
                      sep = "")), by = session_id)

  # Step 7: Replace error latencies with block mean + 600ms
  # Already done in step above
  
  # Step 8: No transformation of latencies
  
  # Step 9: Average latencies for each of the four blocks
  # This was already done above
  
  # Step 10: Compute the two differences B6 - B3 and B7 - B4
  # Step 11: Divide each difference by associated pooled SD from step 6
  
  tbl_diff <- my_tbl_sds %>%
    mutate_(DIFF1 = interp(~ c - a,
                           c = as.name(paste0("M", trial_blocks[3])),
                           a = as.name(paste0("M", trial_blocks[1]))),
            DIFF2 = interp(~ d - b,
                           d = as.name(paste0("M", trial_blocks[4])),
                           b = as.name(paste0("M", trial_blocks[2]))))
  
  if(v_std == 1){
    
    tbl_result <- tbl_diff %>%
      mutate_(IAT1 = interp(~ DIFF1/AS1, DIFF1 = as.name("DIFF1"), AS1 = as.name("AS1")),
              IAT2 = interp(~ DIFF2/AS2, DIFF2 = as.name("DIFF2"), AS2 = as.name("AS2"))) %>%
      mutate_(IAT =  interp(~ (IAT1 + IAT2)/2, IAT1 = as.name("IAT1"), IAT2 = as.name("IAT2")))
      
  } else if(v_std == 2){
    
    tbl_result <- tbl_diff %>%
      mutate_(IAT1 = interp(~ DIFF1/CS1, DIFF1 = as.name("DIFF1"), CS1 = as.name("CS1")),
              IAT2 = interp(~ DIFF2/CS2, DIFF2 = as.name("DIFF2"), CS2 = as.name("CS2"))) %>%
      mutate_(IAT =  interp(~ (IAT1 + IAT2)/2, IAT1 = as.name("IAT1"), IAT2 = as.name("IAT2")))
    
  } else stop("Please enter a v_std value of 1 or 2.")

  return(tbl_result)
  
}

