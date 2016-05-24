
globalVariables(c("id", "messages", "postchat_survey", "events", "visitor", "queue"))

#' Create a list that contains 4 data frames containing chat sessions and raw chat text.
#'
#' @param account defined by the livechatCreateAccount function
#' @param date_from pull all chats from the date in the format YYYY-MM-DD
#'
#' @return A list of 4. The being a data.table of chat sessions, messages, postchat surveys and events
#' @export
#'
#' @examples
#' # account <- livechatCreateAccount("email_here", "api_key_here")
#' # livechat_data <- livechatGetChats(account, date_from = "2016-02-23")
#' # str(livechat_data)
#' # livechat_data[[1]] %>% str # chat_sessions
#' # livechat_data[[2]] %>% str # messages
#' # livechat_data[[3]] %>% str # postchat_survey
#' # livechat_data[[4]] %>% str # events

livechatGetChats <- function(
  account,
  date_from
) {

  email <- account$email
  api_key <- account$api_key
  email_api_key <- paste0(email, ":", api_key)

  # check if there are multiple pages
  q <- paste0("curl 'https://api.livechatinc.com/chats?",
            "date_from=", date_from,
            "&'",
            " -u ",
            email_api_key,
            " -H X-API-Version:2")
  raw_json <- system(q, intern = TRUE)
  raw_json_char <- paste(raw_json, collapse = "")
  chats <- fromJSON(raw_json_char)
  pages <- chats$pages
  pages

  chats_list = list() # create empty chat list
  for(i in 1:pages){ # increment through pages
    q <- paste0("curl 'https://api.livechatinc.com/chats?",
                "date_from=", date_from,
                "&",
                "page=", i,
                "&'",
                " -u ",
                email_api_key,
                " -H X-API-Version:2")
    raw_json <- system(q, intern = TRUE)
    raw_json_char <- paste(raw_json, collapse = "")
    chats <- fromJSON(raw_json_char)
    tmp <- chats$chats
    chats_list[[i]] <- tmp
  }

  # select columns
  exclude_visitor_col <- function(df){
    # exclude visitor if it exists, otherwise select all, visitor wasn't working in rbind.fill
    if("visitor" %in% names(df)) {
      select(df, -visitor)
    } else {
      select(df, 1:ncol(df))
    }
  }

  exclude_queue_col <- function(df){
    if("queue" %in% names(df)) {
      select(df, -queue)
    } else {
      select(df, 1:ncol(df))
    }
  }

  chats_list_exclude_visitor <- map(chats_list, exclude_visitor_col)
  chats_list_excluded_cols <- map(chats_list_exclude_visitor, exclude_queue_col)
  all_data <- rbindlist(chats_list_excluded_cols, fill = TRUE)



  # explore all chat data ---------------------------------------------------
  # null columns that are blank
  all_data$supervisors <- NULL
  all_data$tags <- NULL
  all_data$tickets <- NULL # no tickets data in live chat, tickets are created in zendesk

  # unlist groups and add them back over the columns that were still lists of 1
  groups_unlisted <- unlist(all_data$group) # just 1 column
  all_data$group <- groups_unlisted

  # NULLs cannot be in a vector for some reason
  agents_no_null <- map(all_data$agents, function(x) ifelse(is.null(x), NA, x))
  agents_df <- purrr::flatten(agents_no_null) %>% as.matrix() %>% as.character()

  d <- cbind(all_data, agents_df)
  d$agents <- NULL # get rid of old agents column

  # treat actually list columns separately
  chat_data_lists <- d %>% select(id, messages, postchat_survey, events)
  chat_data <- d %>% select(-messages, -postchat_survey, -events)

  setkey(chat_data_lists, id)

  # NULL message column, database doesn't like multi line strings
  chat_data$message <- NULL

  # parse the parts that are still lists ---------------------------

  # MESSAGES
  # vars1 <- c("id", "messages")
  id_messages <- subset(chat_data_lists, select = c(id, messages))
  # for some reason, some messages data frames are null
  null_messages <- map(id_messages$messages, is.null) %>% unlist()

  if(sum(null_messages) == 0) {
    id_messages_no_null <- id_messages
  } else {
    id_messages_no_null <- id_messages[!null_messages, ] # filter out rows with null messages
  }

  # id_messages_no_null <- id_messages[!null_messages] # filter out rows with null messages

  df_lengths <- map(id_messages_no_null$messages, nrow) %>% unlist() # first create a vector with the length of every messages df
  # replicate id with df_lengths
  id_repped <- rep(id_messages_no_null$id, df_lengths)
  messages_data_frame <- rbindlist(id_messages_no_null$messages, fill = TRUE)
  messages_data_frame$id <- id_repped

  # POST CHAT SURVEY
  id_pcs <- subset(chat_data_lists, select = c(id, postchat_survey))
  # for some reason, some messages data frames are null
  null_pcs <- map(id_pcs$postchat_survey, is.null) %>% unlist()
  id_pcs_no_null <- id_pcs[!null_pcs, ] # filter out rows with null messages
  df_lengths <- map(id_pcs_no_null$postchat_survey, nrow) %>% unlist() # first create a vector with the length of every messages df
  # replicate id with df_lengths
  id_repped <- rep(id_pcs_no_null$id, df_lengths)
  pcs_df <- rbindlist(id_pcs_no_null$postchat_survey, fill = TRUE)
  pcs_df$id <- id_repped

  # EVENTS
  id_events <- subset(chat_data_lists, select = c(id, events))
  # for some reason, some messages data frames are null
  null_events <- map(id_events$events, is.null) %>% unlist()
  id_events_no_null <- id_events[!null_events, ] # filter out rows with null messages
  df_lengths <- map(id_events_no_null$events, nrow) %>% unlist() # first create a vector with the length of every messages df
  # replicate id with df_lengths
  id_repped <- rep(id_events_no_null$id, df_lengths)
  events_df <- rbindlist(id_events_no_null$events, fill = TRUE)
  events_df$id <- id_repped

  # combine into one list for export
  livechat_data <- list(chat_data, messages_data_frame, pcs_df, events_df)
  names(livechat_data) <- c("chat_sessions", "messages", "postchat_survey", "events")
  livechat_data
}

# livechatGetChats(account, date_from = "2016-02-23")
