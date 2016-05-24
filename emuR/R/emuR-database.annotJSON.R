#######################################################
# annotJSON representation to annotDFs conversion functions

# convert annotJSON to list of data.frames including 
# meta information (name, annotates, samplerate)
annotJSONcharToBundleAnnotDFs <- function(annotJSONchar){
  
  json = tidyjson::as.tbl_json(annotJSONchar)
  
  # get top level data
  tlData = json %>%
    tidyjson::spread_values(name = tidyjson::jstring("name"), annotates = tidyjson::jstring("annotates"), sampleRate = tidyjson::jstring("sampleRate"))
  
  # gen. links data.frame
  links = json %>%
    tidyjson::enter_object("links") %>%
    tidyjson::gather_array()  %>%
    tidyjson::spread_values(fromID = tidyjson::jstring("fromID"), toID = tidyjson::jstring("toID")) %>%
    dplyr::select_(~fromID, ~toID) %>%
    dplyr::rename_("from_id" = "fromID", "to_id" = "toID")
  
  # gen. items list of data.frame
  items = json %>%
    tidyjson::spread_values(sampleRate = tidyjson::jstring("sampleRate")) %>%
    tidyjson::enter_object("levels") %>%
    tidyjson::gather_array()  %>%
    tidyjson::spread_values(level = tidyjson::jstring("name"), type = tidyjson::jstring("type")) %>%
    tidyjson::enter_object("items") %>%
    tidyjson::gather_array(column.name = "seq_idx") %>%
    tidyjson::spread_values(itemID = tidyjson::jstring("id"), samplePoint = tidyjson::jstring("samplePoint"), sampleStart = tidyjson::jstring("sampleStart"), sampleDur = tidyjson::jstring("sampleDur")) %>%
    dplyr::select_(~itemID, ~level, ~type, ~seq_idx, ~sampleRate, ~samplePoint, ~sampleStart, ~sampleDur) %>%
    dplyr::rename_("item_id" = "itemID", "sample_rate" = "sampleRate", "sample_point" = "samplePoint", "sample_start" = "sampleStart", "sample_dur" = "sampleDur")
    
  # gen. label list of data.frame
  labels = json %>%
    tidyjson::enter_object("levels") %>%
    tidyjson::gather_array()  %>%
    tidyjson::spread_values(level = tidyjson::jstring("name")) %>%
    tidyjson::enter_object("items") %>%
    tidyjson::gather_array() %>%
    tidyjson::spread_values(itemID = tidyjson::jstring("id")) %>%
    tidyjson::enter_object("labels") %>%
    tidyjson::gather_array(column.name = "label_idx") %>%
    tidyjson::spread_values(name = tidyjson::jstring("name"), label = tidyjson::jstring("value")) %>%
    dplyr::select_(~itemID, ~label_idx, ~name, ~label)
  
  return(list(name = tlData$name, annotates = tlData$annotates, sampleRate = tlData$sampleRate, items = items, links = links, labels = labels))
  
}

# convert annotDFs (annotation list of data.frame representation) to annotJSON
bundleAnnotDFsToAnnotJSONchar <- function(emuDBhandle, annotDFs){
  # load DBconfig to generate levelNames vector (although levels are not ordered per say)
  levelDefs = list_levelDefinitions(emuDBhandle)
  
  levels = list()
  
  for(l in levelDefs$name){
    levelItems = dplyr::filter_(annotDFs$items, ~(level == l))
    
    levels[[length(levels) + 1]] = list(
      items = apply(levelItems, 1, function(r) {
      
      labels = apply(dplyr::filter_(annotDFs$labels, ~(item_id == as.numeric(r[1]))), 1, function(r2) list(name = as.character(r2[3]), value = as.character(r2[4])))
      res = NULL
      if(r[3] == "ITEM"){
        res = list(id = as.numeric(r[1]),
                   labels = labels)
      }else if(r[3] == "SEGMENT"){
        res = list(id = as.numeric(r[1]),
                   sampleStart = as.numeric(r[7]),
                   sampleDur = as.numeric(r[8]),
                   labels = labels)
      }else if(r[3] == "EVENT"){
        res = list(id = as.numeric(r[1]),
                   samplePoint = as.numeric(r[6]),
                   labels = labels)
      }
      return(res)
    }),
    name = l,
    type = levelDefs$type[levelDefs$name == l]
    )
  }
  
  if(nrow(annotDFs$links) >0){
    links = apply(annotDFs$links, 1, function(r) list(fromID = as.numeric(r[1]), toID = as.numeric(r[2])))
  }else{
    links = list()
  }
  
  annotJSON = list(name = annotDFs$name,
                   annotates = annotDFs$annotates,
                   sampleRate = annotDFs$sampleRate,
                   levels = levels, links = links)
  
  return(jsonlite::toJSON(annotJSON, auto_unbox = T, force = T, pretty = T))
}

