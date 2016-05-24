requireNamespace("RSQLite", quietly = T)
requireNamespace("stringr", quietly = T)

## EmuDB Parser for Bas Partitur Files
## 
## @param bpfPath
## @param samplerate
## @param encoding
## @param dbName
## @param bundle
## @param session
## @param dbUUID
## @param refLevel
## @param segmentToEventLevels
## @param levelClasses
## @return list(levelInfo, linkInfo, warningsInfo)
## @import stringr RSQLite
## @keywords emuR BPF Emu

parse_BPF <- function(emuDBhandle,
                      bpfPath,
                      encoding = "UTF-8",
                      bundle,
                      session,
                      refLevel,
                      extractLevels,
                      samplerate,
                      segmentToEventLevels,
                      unifyLevels,
                      levelClasses)
{
  # ---------------------------------------------------------------------------
  # --- Containers for info to be passed out to caller (converter) function ---
  # ---------------------------------------------------------------------------
  
  levelInfo = list()
  linkInfo = list()
  semicolonFound = FALSE
  
  # ---------------------------------------------------------------------------
  # ------------------------ Read BPF file from disk --------------------------
  # ---------------------------------------------------------------------------
  
  bpfLines = try(readLines(bpfPath))
  if(class(bpfLines) == "try-error") 
  {
    stop("Cannot read from file ", bpfPath)
  }
  
  if(length(bpfLines) == 0)
  {
    stop("File ", bpfPath, " has length 0. This does not conform to BPF specifications.")
  }
  
  # ---------------------------------------------------------------------------
  # -------------------------- Parse header -----------------------------------
  # ---------------------------------------------------------------------------
  
  returnContainer = parse_bpfHeader(
    bpfLines = bpfLines, 
    bpfPath = bpfPath,
    samplerate = samplerate
  )
  
  header = returnContainer$header
  bsKeyPosition = returnContainer$bsKeyPosition
  samplerate = returnContainer$samplerate
  
  # ---------------------------------------------------------------------------
  # ---------- Write 'Utterance' item to items and lables tables --------------
  # ---------------------------------------------------------------------------
  levelInfo = write_bpfUtteranceToDb(emuDBhandle,
                                     header = header,
                                     session = session,
                                     bundle = bundle,
                                     samplerate = samplerate)
  
  # Utterance item will be written even if the BPF body is empty!
  
  # ---------------------------------------------------------------------------
  # -------------------------- Parse body -------------------------------------
  # ---------------------------------------------------------------------------
  
  if(bsKeyPosition < length(bpfLines))
  {
    returnContainer = parse_bpfBody(bpfLines = bpfLines,
                                    bpfPath = bpfPath,
                                    bsKeyPosition = bsKeyPosition,
                                    extractLevels = extractLevels,
                                    levelClasses = levelClasses,
                                    unifyLevels = unifyLevels,
                                    refLevel = refLevel,
                                    segmentToEventLevels = segmentToEventLevels)
    
    levels = returnContainer$levels
    currentItemID = returnContainer$currentItemID
    semicolonFound = returnContainer$semicolonFound
    
    # -------------------------------------------------------------------------
    # --- Change classes of levels in segmentToEventLevels (2->3 and 4->5) ----
    # -------------------------------------------------------------------------
    
    # (done after parsing because parser needs original classes to make sense of BPF lines)
    
    for(key in segmentToEventLevels)
    {
      levelClasses[[key]] = levelClasses[[key]] + 1
    }
    
    # -------------------------------------------------------------------------
    # ------ Check for temporal overlap within levels with time information ---
    # -------------------------------------------------------------------------
    
    check_bpfOverlap(levels = levels,
                     bpfPath = bpfPath,
                     segmentToEventLevels = segmentToEventLevels,
                     levelClasses = levelClasses)
    
    # -------------------------------------------------------------------------
    # ---------- Pad segment tiers between segments with empty items ----------
    # -------------------------------------------------------------------------
    
    levels = pad_bpfSegments(levels = levels,
                             currentItemID = currentItemID,
                             levelClasses = levelClasses)
    
    # -------------------------------------------------------------------------
    # ------------------------------ Assign seqIdx ----------------------------
    # -------------------------------------------------------------------------
    
    levels = assign_bpfSeqIdx(levels = levels,
                              levelClasses = levelClasses)
    
    # -------------------------------------------------------------------------
    # ----------- Write item and label information to database ----------------
    # -------------------------------------------------------------------------
    
    levelInfo = write_bpfItemsLabelsToDb(emuDBhandle,
                                         levels = levels,
                                         session = session,
                                         bundle = bundle,
                                         samplerate = samplerate,
                                         unifyLevels = unifyLevels,
                                         levelInfo = levelInfo)
    
    # -------------------------------------------------------------------------
    # --------------- Write link information to database ----------------------
    # -------------------------------------------------------------------------
    
    if(!is.null(refLevel))
    {
      linkIdxMap = get_bpfLinkIdxMap(levels = levels,
                                     refLevel = refLevel)
      
      # -----------------------------------------------------------------------
      # --------------- Write link information to database --------------------
      # -----------------------------------------------------------------------
      
      linkInfo = write_bpfLinksToDb(emuDBhandle,
                                    levels = levels,
                                    levelClasses = levelClasses,
                                    linkIdxMap = linkIdxMap,
                                    refLevel = refLevel,
                                    session = session,
                                    bundle = bundle,
                                    unifyLevels = unifyLevels,
                                    bpfPath = bpfPath)
      
      # -----------------------------------------------------------------------
      # -------- Unify levels in unifyLevels with the reference level ---------
      # -----------------------------------------------------------------------
      
      if(!is.null(unifyLevels))
      {
        levelInfo = unify_bpfLevels(emuDBhandle,
                                    levels = levels,
                                    linkIdxMap = linkIdxMap,
                                    refLevel = refLevel,
                                    bpfPath = bpfPath,
                                    levelInfo = levelInfo,
                                    unifyLevels = unifyLevels,
                                    session = session,
                                    bundle = bundle)
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # --------- Return info containers to caller (converter) function -----------
  # ---------------------------------------------------------------------------
  
  returnContainer = list(levelInfo = levelInfo, linkInfo = linkInfo, semicolonFound = semicolonFound)
  return(returnContainer)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# -----------------------------------------------------------------------------
# ------------------------- HELPER FUNCTIONS ----------------------------------
# -----------------------------------------------------------------------------

## Parser for Bas Partitur header
## 
## @param bpfLines
## @param bpfPath
## @param samplerate
## @keywords emuR BPF Emu
## @return list(header, bsKeyPosition, missingHeaderKeys, samplerate)

parse_bpfHeader <- function(bpfLines,
                            bpfPath,
                            samplerate)
{
  # ---------------------------------------------------------------------------
  # --------------------------- Necessary constants ---------------------------
  # ---------------------------------------------------------------------------
  
  # Supplements for ranges in regular expressions.
  UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  DIGITS = "0123456789"
  
  # Obligatory header keys. If these are not present in the BPF header -> warningTracker
  OBLIGATORY_HEADER_KEYS = list("LHD", "REP", "SNB", "SAM", "SBF", "SSB", "NCH", "SPN")
  
  # All lines of the BPF must conform to the following regular expression:
  GLOBAL_REGEX = paste0("^[", UPPER, DIGITS, "]{3}:.*")
  
  # Key that marks the beginning of the BPF body / end of the BPF header:
  BODY_START_KEY = "LBD"
  
  # ---------------------------------------------------------------------------
  # --------------------------- Initialize containers -------------------------
  # ---------------------------------------------------------------------------
  
  # Container for key value pairs from the BPF header.
  header = list()
  
  # Container for found keys (to check for duplicates).
  foundKeys = c()
  
  # Line index of the body start key (needed for parse_bpfBody to know where to start).
  bsKeyPosition = NULL
  
  # ---------------------------------------------------------------------------
  # ---------------- Parse until body start key is found ----------------------
  # ---------------------------------------------------------------------------
  
  for (idx in 1:length(bpfLines))
  {
    # Skip empty lines.
    if(str_length(bpfLines[idx]) == 0)
    {
      next
    }
    
    # Check line's format.
    if (!str_detect(bpfLines[idx], GLOBAL_REGEX))
    {
      stop("Line ", idx, " of the following BPF does not conform to BPF specifications: ", bpfPath)
    }
    
    # Get key value pair.
    splitline = str_split_fixed(bpfLines[idx], ":", 2)
    key = splitline[1]
    
    # Remove trailing white space and escape single quotes (compatibility with SQL).
    value = str_replace(str_replace_all(splitline[2], "'", "''"), "^\\s+", "")
    
    # Once the body start key is found, remember its position and break.
    if(key == BODY_START_KEY)
    {
      bsKeyPosition = idx
      break
    }
    
    header[[key]] = value
    foundKeys = c(foundKeys, key)
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------- Some checks  ---------------------------------
  # ---------------------------------------------------------------------------
  
  # Throw exception if the body start key has not been found.
  if(is.null(bsKeyPosition))
  {
    stop("The following BPF does not contain the body start key 'LBD': ", bpfPath)
  }
  
  # Throw exception if a key was found more than once in the header.
  if(length(unique(foundKeys)) < length(foundKeys))
  {
    stop("There is a duplicate header key in the following BPF: ", bpfPath)
  }
  
  # ---------------------------------------------------------------------------
  # ----- Compare samplerate of audio with the one declared in BPF header -----
  # ---------------------------------------------------------------------------
  
  samplerate = compare_bpfSamplerate(samplerate = samplerate,
                                     header = header,
                                     bpfPath = bpfPath)
  
  # ---------------------------------------------------------------------------
  # ---------------------------------- Return  --------------------------------
  # ---------------------------------------------------------------------------
  
  returnContainer = list(header = header, bsKeyPosition = bsKeyPosition, samplerate = samplerate)
  
  return(returnContainer)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Compare samplerate of audio file with the one declared in the BPF header
## 
## @param header
## @param samplerate
## @keywords emuR BPF Emu
## @return samplerate

compare_bpfSamplerate <- function(header,
                                  samplerate,
                                  bpfPath)
{
  # Throw an exception if we can't get a sample rate from the BPF or the audio file.
  if(is.null(header$SAM) && is.null(samplerate))
  {
    stop("Sample rate has not been read from audio and is therefore needed in the following BPF: ", bpfPath)
  }
  
  # If we don't have a sample rate from the audio, get samle rate from BPF.
  else if(!is.null(header$SAM) && is.null(samplerate))
  {
    samplerate = as.integer(header$SAM)
  }
  
  # If we have one sample rate from the audio and one from the BPF, check if they match.
  else if(!is.null(header$SAM) && !is.null(samplerate))
  {
    if(as.integer(header$SAM) != samplerate)
    {
      stop("Declared sample rate in the following BPF does not match the sample rate of the audio: ", bpfPath)
    }
  }
  
  return(samplerate)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Write item and label information for Utterance Item
## 
## @param emuDBhandle
## @param session
## @param bundle
## @param samplerate
## @param header
## @keywords emuR BPF Emu
## @return levelInfo

write_bpfUtteranceToDb <- function(emuDBhandle,
                                   session,
                                   bundle,
                                   samplerate,
                                   header)
{
  # Utterance gets itemID 1 (other items will start at ID 2)
  utteranceItemID = 1
  
  # Collect label keys ("Utterance" + all header keys found).
  labelTracker = list("Utterance")
  
  queryTxt = paste0("INSERT INTO items VALUES"," ('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', ",
                    utteranceItemID, ", 'Utterance', 'ITEM', 1, ", samplerate, ", NULL, NULL, NULL)")
  
  DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
  
  labelIdxCounter = 1
  
  # First label: 'Utterance' -> name of bundle.
  queryTxt = paste0("INSERT INTO labels VALUES","('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', ",
                    utteranceItemID, ", ", labelIdxCounter, ", 'Utterance', '", bundle, "')")
  
  DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
  
  labelIdxCounter = labelIdxCounter + 1
  
  
  # Subsequent labels: Key -> value pairs found in BPF header.
  for(key in names(header))
  {
    queryTxt = paste0("INSERT INTO labels VALUES","('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', ",
                      utteranceItemID, ", ", labelIdxCounter, ", '", key,"', '", header[[key]], "')")
    
    DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
    
    labelTracker[[length(labelTracker) + 1L]] = key
    labelIdxCounter = labelIdxCounter + 1
  }
  
  
  levelInfo = list(list(key = "Utterance", type = "ITEM", labels = labelTracker))
  
  return(levelInfo)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Parser for Bas Partitur Body
## 
## @param bpfLines
## @param bpfPath
## @param bsKeyPosition
## @param extractLevels
## @param levelClasses
## @param unifyLevels
## @param currentItemID
## @param refLevel
## @keywords emuR BPF Emu
## @return list(levels, currentItemID, semicolonFound)

parse_bpfBody <- function(bpfLines,
                          bpfPath,
                          bsKeyPosition,
                          extractLevels, 
                          levelClasses,
                          unifyLevels,
                          refLevel,
                          segmentToEventLevels)
{
  # ---------------------------------------------------------------------------
  # -------------------------- NECESSARY CONSTANTS ----------------------------
  # ---------------------------------------------------------------------------
  
  # Supplements for ranges in regular expressions.
  UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  DIGITS = "0123456789"
  
  # Regular expressions for BPF lines of classes 1-5.
  CLASS_REGEXES = c(
    paste0("^[", UPPER, DIGITS, "]{3}:\\s+-?[", DIGITS, "][", DIGITS, ",;]*\\s+.*"),
    paste0("^[", UPPER, DIGITS, "]{3}:\\s+[", DIGITS, "]+\\s+[", DIGITS, "]+\\s+.*"),
    paste0("^[", UPPER, DIGITS, "]{3}:\\s+[", DIGITS, "]+\\s+.*"),
    paste0("^[", UPPER, DIGITS, "]{3}:\\s+[", DIGITS, "]+\\s+[", DIGITS, "]+\\s+-?[", DIGITS, "][", DIGITS, ",;]*\\s+.*"),
    paste0("^[", UPPER, DIGITS, "]{3}:\\s+[", DIGITS, "]+\\s+-?[", DIGITS, "][", DIGITS, ",;]*\\s+.*")
  )
  
  # The number of pieces a BPF line should be split into according to its class.
  CLASS_SPLITNUMS = c(3, 4, 3, 5, 4)
  
  # All lines of the BPF must conform to the following regular expression:
  GLOBAL_REGEX = paste0("^[", UPPER, DIGITS, "]{3}:.*")
  
  # Item type according to class.
  CLASS_TO_TYPE = c("ITEM", "SEGMENT", "EVENT", "SEGMENT", "EVENT")
  
  # ---------------------------------------------------------------------------
  # ------------------------------- Containers --------------------------------
  # ---------------------------------------------------------------------------
  
  # Boolean indicating whether a semicolon link operator (not supported) was found in the present BPF.
  semicolonFound = FALSE
  
  # Container for levels.
  levels = list()
  
  # ---------------------------------------------------------------------------
  # --------------------------- Parsing ---------------------------------------
  # ---------------------------------------------------------------------------
  
  # Initialize current itemID at 2 (since 'Utterance' is 1).
  currentItemID = 2
  
  # Start parsing from body start key onwards (body start key position + 1).
  for(idx in (bsKeyPosition+1):length(bpfLines))
  {
    # Skip empty lines.
    if(str_length(bpfLines[idx]) == 0)
    {
      next
    }
    
    # Throw an exception if a line does not match the global regular expression.
    if (!str_detect(bpfLines[idx], GLOBAL_REGEX))
    {
      stop("Line ", idx, " of the following BPF does not conform to BPF specification: ", bpfPath)
    }
    
    # Get level name (first three characters)
    key = str_sub(bpfLines[idx], start = 1, end = 3)
    
    # If only a subset of levels should be extracted, and this level is not one of them, next.
    if(!is.null(extractLevels))
    {
      if(!key %in% extractLevels)
      {
        next
      }
    }
    
    if(!key %in% names(levelClasses))
    {
      stop("Unknown level name in line ", idx, " of the following BPF: ", bpfPath, ". If this level is not one of the standard BPF tiers, you have to declare it using the newLevels argument.")
    }
    
    if(!key %in% names(levels))
    {
      levels[[key]] = list()
    }
    
    # Throw an exception if the line does not conform to the regular expression of its class.
    # WARNING: Cannot detect all errors!
    if(!str_detect(bpfLines[idx], CLASS_REGEXES[levelClasses[[key]]]))
    {
      stop("Line ", idx, " in the following BPF does not match the Bas Partitur File Specifications: ", bpfPath,
           ". Level '", key, "' should be of class ", levelClasses[[key]], ".")
    }
    
    # Split the line according to its class.
    splitline = str_split_fixed(bpfLines[idx], "\\s+", CLASS_SPLITNUMS[levelClasses[[key]]])
    
    # Assign and increment global index.
    if(!key %in% unifyLevels)
    {
      itemID = currentItemID
      currentItemID = currentItemID + 1
    }
    
    # If the key in unifyLevels, assign no index (since this won't become an independent item but a label).
    else
    {
      itemID = NA
    }
    
    # Assign type, based on key class.
    type = CLASS_TO_TYPE[levelClasses[[key]]]
    
    # Initialize seq index as NA (assigned later).
    seqIdx = NA
    
    # -------------------------------------------------------------------------
    # --------------- Parse BPF line accrding to level class ------------------
    # -------------------------------------------------------------------------
    
    returnContainer = parse_bpfLine(levelClass = levelClasses[[key]],
                                    splitline = splitline)
    
    start = returnContainer$start
    duration = returnContainer$duration
    point = returnContainer$point
    labelString = returnContainer$labelString
    linksString = returnContainer$linksString
    
    # -------------------------------------------------------------------------
    # ---------------- Evaluate information in labelString --------------------
    # -------------------------------------------------------------------------
    
    labels = evaluate_bpfLabelString(labelString = labelString,
                                     key = key)
    
    # ---------------------------------------------------------------------------
    # -------------------- Evaluate information in linksString ------------------
    # ---------------------------------------------------------------------------
    
    returnContainer = evaluate_bpfLinksString(linksString = linksString,
                                              bpfPath = bpfPath,
                                              refLevel = refLevel,
                                              key = key)
    
    links = returnContainer$links
    
    if(returnContainer$semicolon)
    {
      semicolonFound = TRUE
    }
    
    # -------------------------------------------------------------------------
    # ------- Turn segment into event if level in segmentToEventLevels --------
    # -------------------------------------------------------------------------
    
    if(key %in% segmentToEventLevels)
    {
      itemID_start = itemID
      itemID_end = currentItemID
      currentItemID = currentItemID + 1
      
      point_start = start
      point_end = start + duration
      
      start = "NULL"
      duration = "NULL"
      
      type = "EVENT"
      
      labels_start = list()
      labels_end = list()
      
      
      for(key in names(labels))
      {
        labels_start[[key]] = paste0(labels[[key]], "_start")
        labels_end[[key]] = paste0(labels[[key]], "_end")
      }
      
      levels[[key]][[length(levels[[key]]) + 1L]] = list(itemID = itemID_start, start = start, duration = duration, point = point_start,
                                                         labels = labels_start, links = links, seqIdx = seqIdx, type = type)
      
      levels[[key]][[length(levels[[key]]) + 1L]] = list(itemID = itemID_end, start = start, duration = duration, point = point_end,
                                                         labels = labels_end, links = links, seqIdx = seqIdx, type = type)
    }
    
    else
    {
      levels[[key]][[length(levels[[key]]) + 1L]] = list(itemID = itemID, start = start, duration = duration, point = point,
                                                         labels = labels, links = links, seqIdx = seqIdx, type = type)
    }
  }
  
  # ---------------------------------------------------------------------------
  # -------------------------------- Return -----------------------------------
  # ---------------------------------------------------------------------------
  
  returnContainer = list(levels = levels, currentItemID = currentItemID, semicolonFound = semicolonFound)
  return(returnContainer)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Parse a single BPF line according to the line's level class
## 
## @param levelClass
## @param splitline
## @import stringr
## @keywords emuR BPF Emu
## @return list(start, duration, point, labelString, linksString)

parse_bpfLine <- function(levelClass,
                          splitline)
{
  if(levelClass == 1)
  {
    start = "NULL"
    duration = "NULL"
    point = "NULL"
    
    # Escape single quotes with double quotes (conformity with SQL).
    labelString = str_replace_all(splitline[3], "'", "''")
    linksString = splitline[2]
  }
  
  else if(levelClass == 2)
  {
    start = as.integer(splitline[2])
    duration = as.integer(splitline[3])
    point = "NULL"
    labelString = str_replace_all(splitline[4], "'", "''")
    linksString = NA
  }
  
  else if(levelClass == 3)
  {
    start = "NULL"
    duration = "NULL"
    point = as.integer(splitline[2])
    labelString = str_replace_all(splitline[3], "'", "''")
    linksString = NA
  }
  
  else if(levelClass == 4)
  {
    start = as.integer(splitline[2])
    duration = as.integer(splitline[3])
    point = "NULL"
    labelString = str_replace_all(splitline[5], "'", "''")
    linksString = splitline[4]
  }
  
  else if(levelClass == 5)
  {
    start = "NULL"
    duration = "NULL"
    point = as.integer(splitline[2])
    labelString = str_replace_all(splitline[4], "'", "''")
    linksString = splitline[3]
  }
  
  return(list(start = start, duration = duration, point = point, 
              labelString = labelString, linksString = linksString))
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Turn raw label string into one or several value - key pairs
## 
## @param labelString
## @param key
## @import stringr
## @keywords emuR BPF Emu
## @return labels

evaluate_bpfLabelString <- function(labelString,
                                    key)
{
  # ---------------------------------------------------------------------------
  # -------------------------- NECESSARY CONSTANTS ----------------------------
  # ---------------------------------------------------------------------------
  
  # Supplements for ranges in regular expressions.
  UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  DIGITS = "0123456789"
  
  # If the label string contains more than one label, expand into separate key -> value pairs.
  # Syntax for label string with multiple labels: "ABC: value; DEF: value; GHI: value".
  # If not, the level name becomes the label key, and the full label string becomes the value.
  
  labels = list()
  
  if(str_detect(labelString, paste0("^[", UPPER, DIGITS, "]{3}:\\s+.*;")) &&
     str_detect(labelString, paste0(";\\s*[", UPPER, DIGITS, "]{3}:\\s+.*$")))
  {
    extractedLabels = str_split(labelString, "\\s*;\\s*")[[1]]
    for(extractedLabel in extractedLabels)
    {
      splitLabel = str_split(extractedLabel, ":\\s+", n=2)[[1]]
      
      labels[[splitLabel[1]]] = splitLabel[2]
    }
  }
  
  else
  {
    labels[[key]] = labelString
  }
  
  return(labels)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Turn raw links string into NA (no links) or a vector of integers
## 
## @param linksString
## @param bpfPath
## @param refLevel
## @param key
## @import stringr
## @keywords emuR BPF Emu
## @return list(links, semicolon)

evaluate_bpfLinksString <- function(linksString,
                                    bpfPath,
                                    refLevel,
                                    key)
{
  # Variable to be returned (TRUE if a semicolon has been found in this BPF -> warningTracker).
  semicolon = FALSE
  
  # If there was no link entry in the first place, or the link was '-1' -> no link information.
  if(is.na(linksString) || str_detect(linksString, "-1"))
  {
    links = NA
  }
  
  # Ignore links containing the ';' operator.
  else if(str_detect(linksString, ";"))
  {
    semicolon = TRUE
    links = NA
  }
  
  # Store links as a vector of integers. 
  else
  {
    links = as.integer(unlist(str_split(linksString, ",")))
    
    # Throw an exception if an item links to the same item more than once.
    for(link in links)
    {
      if(sum(links == link) > 1)
      {
        stop("An item cannot link to the same item more than once. BPF: ", bpfPath)
      }
    }
  }
  
  # If the current level is the reference level, check whether all links are valid and atomic.
  if(!is.null(refLevel))
  {
    if(key == refLevel)
    {
      if(length(links) > 1)
      {
        stop("The reference level must contain atomic links. Not the case in the following BPF: ", bpfPath)
      }
      if(is.na(links[[1]]))
      {
        stop("The reference level must contain valid symbolic links. Valid symbolic links are neither '-1', nor do they contain the ';' operator. BPF: ", bpfPath)
      }
    }
  }
  
  return(list(links = links, semicolon = semicolon))
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Check for temporal overlap on event and segment tiers
## 
## @param levels
## @param bpfPath
## @param segmentToEventLevels
## @param levelClasses
## @import stringr
## @keywords emuR BPF Emu
## @return

check_bpfOverlap <- function(levels, 
                             bpfPath, 
                             levelClasses,
                             segmentToEventLevels)
{
  for(key in names(levels))
  {
    # If the level is time consuming, check whether there is segmental overlap.
    if(levelClasses[[key]] %in% c(2, 4))
    {
      for(idx in 1:length(levels[[key]]))
      {
        for(jdx in 1:length(levels[[key]]))
        {
          if(
            idx!=jdx && 
            levels[[key]][[idx]][["start"]] >= levels[[key]][[jdx]][["start"]] && 
            levels[[key]][[idx]][["start"]] <= levels[[key]][[jdx]][["start"]] + levels[[key]][[jdx]][["duration"]]
          )
          {
            stop("The following BPF contains overlapping segments on level '", key, "': ", bpfPath)
          }
        }
      }
    }
    
    # If the level is not time consuming, check whether there are two events pointing to the same sample.
    if(levelClasses[[key]] %in% c(3, 5))  
    {
      for(idx in 1:length(levels[[key]]))
      {
        for(jdx in idx:length(levels[[key]]))
        {
          if(
            idx!=jdx &&
            levels[[key]][[idx]][["point"]] == levels[[key]][[jdx]][["point"]]
          )
          {
            if(key %in% segmentToEventLevels)
            {
              stop("The following BPF contains simultaneous events on level '", key, "' after segment overlap resolution: ", bpfPath, ". Check whether there are any segments with simultaneous starting and/or end points in this BPF.")
            }
            else
            {
              stop("The following BPF contains simultaneous events on level '", key, "': ", bpfPath)
            }
          }
        }
      }
    }
  }
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Pad space between segment items with empty strings
## 
## @param levels
## @param currentItemID
## @param levelClasses
## @keywords emuR BPF Emu
## @return levels

pad_bpfSegments <- function(levels,
                            currentItemID,
                            levelClasses)
{
  # Pad segment tiers with empty segments.
  # No padding before the first segment and after the last segment.
  
  for(key in names(levels))
  {
    # If there is only one item on this level, there is no padding required: Jump to next level.
    
    if(length(levels[[key]]) == 1 || levelClasses[[key]] %in% c(1, 3, 5))
    {
      next
    }
    
    start_order = sapply(levels[[key]], "[[", "start")
    levels[[key]] = levels[[key]][order(start_order)]
    
    for(idx in 1:(length(levels[[key]])-1))
    {
      if((levels[[key]][[idx]][["start"]] + levels[[key]][[idx]][["duration"]] + 1) < levels[[key]][[idx+1]][["start"]])
      {
        start = levels[[key]][[idx]][["start"]] + levels[[key]][[idx]][["duration"]] + 1
        duration = levels[[key]][[idx+1]][["start"]] - start -1
        
        # Create new item for the pad.
        labels = list()
        labels[[key]] = ""
        levels[[key]][[length(levels[[key]]) + 1L]] = 
          list(itemID = currentItemID, start = start, duration = duration, point = "NULL",
               labels = labels, links = NA, seqIdx = NA, type = "SEGMENT")
        
        currentItemID = currentItemID + 1
      }
    }
  }
  
  return(levels)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Assign seqIdx to items
## 
## @param levels
## @param levelClasses
## @keywords emuR BPF Emu
## @return levels

assign_bpfSeqIdx <- function(levels,
                             levelClasses)
{
  for(key in names(levels))
  {
    # If there is only one item on this level (and the level is thus already ordered chronologically):
    # Assign seqIdx as 1 and jump to next level.
    if(length(levels[[key]]) == 1)
    {
      levels[[key]][[1]][["seqIdx"]] = 1
      next
    }
    
    
    # Counter for assigned indices (starting at 1 for each level).
    currentSeqIdx = 1
    
    # Order items on levels with temporal information chronologically.
    if(levelClasses[[key]] %in% c(2, 4))
    {
      startOrder = sapply(levels[[key]], "[[", "start")
      levels[[key]] = levels[[key]][order(startOrder)]
    }
    
    else if(levelClasses[[key]] %in% c(3, 5))
    {
      pointOrder = sapply(levels[[key]], "[[", "point")
      levels[[key]] = levels[[key]][order(pointOrder)]
    }
    
    # Assign seqIdx from top to bottom.
    for(idx in 1:length(levels[[key]]))
    {
      levels[[key]][[idx]][["seqIdx"]] = currentSeqIdx
      currentSeqIdx = currentSeqIdx + 1
    }
  }
  
  return(levels)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Write item and label info to EmuDB
##
## @param emuDBhandle
## @param levels
## @param session
## @param bundle
## @param samplerate 
## @param unifyLevels
## @keywords emuR BPF Emu
## @return levelInfo

write_bpfItemsLabelsToDb <- function(emuDBhandle,
                                     levels,
                                     session,
                                     bundle,
                                     samplerate,
                                     unifyLevels,
                                     levelInfo)
{  
  for(key in names(levels))
  {
    # Skip current level if it is to be unified with the reference level.
    if(key %in% unifyLevels)
    {
      next
    }
    
    labelTracker = list()
    
    for(idx in 1:length(levels[[key]]))
    {
      
      # Write item information.
      queryTxt = paste0("INSERT INTO items VALUES"," ('", emuDBhandle$UUID, "', '", session, "', '", bundle, "', ",
                        levels[[key]][[idx]][["itemID"]], ", '", key, "', '", levels[[key]][[idx]][["type"]], "', ",
                        levels[[key]][[idx]][["seqIdx"]], ", ", samplerate, ", ", levels[[key]][[idx]][["point"]], ", ",
                        levels[[key]][[idx]][["start"]], ", ", levels[[key]][[idx]][["duration"]], ")")
      
      DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
      
      labelIdxCounter = 1
      
      for(labelKey in names(levels[[key]][[idx]][["labels"]]))
      {
        queryTxt = paste0("INSERT INTO labels VALUES","('", emuDBhandle$UUID, "', '", session, "', '", bundle, 
                          "', ", levels[[key]][[idx]][["itemID"]], ", ", labelIdxCounter,", '", 
                          labelKey, "', '", levels[[key]][[idx]][["labels"]][[labelKey]], "')")
        
        DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
        
        if(!labelKey %in% labelTracker)
        {
          labelTracker[[length(labelTracker) + 1L]] = labelKey
        }
        labelIdxCounter = labelIdxCounter + 1
      }
    }
    levelInfo[[length(levelInfo) + 1L]] = list(key = key, type = levels[[key]][[1]][["type"]], labels = labelTracker)
  }
  
  return(levelInfo)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Get a list mapping from link strings in the BPF to itemIDs on the reference level
## 
## @param levels
## @param refLevel
## @keywords emuR BPF Emu
## @return linkIdxMap

get_bpfLinkIdxMap <- function(
  levels,
  refLevel
)
{
  # Map from link name to indices on reference level.
  linkIdxMap = list()
  for(idx in 1:length(levels[[refLevel]]))
  {
    linkIdxMap[[toString(levels[[refLevel]][[idx]][["links"]][1])]] = levels[[refLevel]][[idx]][["itemID"]]
  }
  
  return(linkIdxMap)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Locally determine link directions and write link info to EmuDB
## 
## @param emuDBhandle
## @param levels
## @param levelClasses
## @param refLevel
## @param session
## @param bundle
## @param unifyLevels
## @param bpfPath
## @param linkIdxMap
## @keywords emuR BPF Emu
## @return linkInfo

write_bpfLinksToDb <- function(emuDBhandle,
                               levels,
                               levelClasses,
                               refLevel,
                               session,
                               bundle,
                               unifyLevels,
                               bpfPath,
                               linkIdxMap)
{
  # Container for information on levels found in this BPF. Will be returned to conversion function.
  linkInfo = list()
  
  for(key in names(levels))
  {
    if(key == refLevel || levelClasses[[key]] %in% c(2, 3) || key %in% unifyLevels)
    {
      next
    }
    
    # -------------------------------------------------------------------------
    # - Get direction and type of links between refLevel and current level ----
    # -------------------------------------------------------------------------
    
    returnContainer = get_bpfLinkCounts(levels,
                                        key)
    
    # If we haven't seen any links, skip and don't make entries to linkInfo or the temp DB
    if(is.null(returnContainer$seenLinks))
    {
      next
    }
    
    oneToMany = returnContainer$oneToMany
    manyToOne = returnContainer$manyToOne
    
    linkInfoEntry = bpf_get_link_info_entry(key = key,
                                            refLevel = refLevel,
                                            oneToMany = oneToMany, 
                                            manyToOne = manyToOne)
    
    upper = linkInfoEntry$fromkey
    lower = linkInfoEntry$tokey
    
    linkInfo[[length(linkInfo) + 1L]] = linkInfoEntry
    
    # -------------------------------------------------------------------------
    # ----------------------- Insert links into temp DB -----------------------
    # -------------------------------------------------------------------------
    
    for(idx in 1:length(levels[[key]]))
    {
      if(is.na(levels[[key]][[idx]][["links"]][1]))
      {
        next
      }
      
      for(link in levels[[key]][[idx]][["links"]])
      {
        if(!(link %in% names(linkIdxMap)))
        {
          stop("There is a symbolic link on level ", key, " in the following BPF that does not point to any item on the reference level: ", bpfPath)
        }
        
        if(upper == refLevel)
        {
          queryTxt = paste0("INSERT INTO links VALUES","('", emuDBhandle$UUID, "', '", session, "', '", 
                            bundle, "', ", linkIdxMap[[toString(link)]], ", ", 
                            levels[[key]][[idx]][["itemID"]],", NULL)")
          
          DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
        }
        else if(lower == refLevel)
        {
          queryTxt =  paste0("INSERT INTO links VALUES","('", emuDBhandle$UUID, "', '", session, "', '", 
                             bundle, "', ", levels[[key]][[idx]][["itemID"]], ", ", 
                             linkIdxMap[[toString(link)]],", NULL)")
          
          
          DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
        }
      }
    }
  }
  
  return(linkInfo)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## 
## 
## @param levels
## @param key
## @keywords emuR BPF Emu
## @return list(oneToMany, manyToOne, seenLinks)

get_bpfLinkCounts <- function(levels,
                              key)
{
  seenLinks = NULL
  oneToMany = 0   
  # oneToMany increments for items on current level linking to more than one item on refLevel 
  # (+1 for each extra item on refLevel)
  manyToOne = 0   
  # increments for items on refLevel linking to more than one item on current level
  # (+1 for each extra item on current level)
  
  for(idx in 1:length(levels[[key]]))
  {
    # Skip if current item does not have any links.
    if(is.na(levels[[key]][[idx]][["links"]])[1])
    {
      next
    }
    
    oneToMany = oneToMany + (length(levels[[key]][[idx]][["links"]]) - 1)
    
    for(link in levels[[key]][[idx]][["links"]])
    {
      if(link %in% seenLinks)
      {
        manyToOne = manyToOne + 1
      }
      else
      {
        seenLinks = c(seenLinks, link)
      }
    }
  }
  
  return(list(oneToMany = oneToMany, manyToOne = manyToOne, seenLinks = seenLinks))
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Evaluate link counts to determine link direction and type
## 
## @param key
## @param refLevel
## @param oneToMany
## @param manyToOne
## @keywords emuR BPF Emu
## @return list(linkInfoEntry)

bpf_get_link_info_entry <- function(key,
                                    refLevel,
                                    oneToMany,
                                    manyToOne)
{
  # ---------------------------------------------------------------------------
  # ------------------------ Determine link type ------------------------------
  # ---------------------------------------------------------------------------
  
  if(oneToMany == 0 && manyToOne == 0)
  {
    linkType = "ONE_TO_ONE"
  }
  
  else if(oneToMany == 0 || manyToOne == 0)
  {
    linkType = "ONE_TO_MANY"
  }
  
  else if(oneToMany > 0 && manyToOne > 0)
  {
    linkType = "MANY_TO_MANY"
  }
  
  # ---------------------------------------------------------------------------
  # ------------------------ Determine link direction -------------------------
  # ---------------------------------------------------------------------------
  
  if(oneToMany <= manyToOne)
  {
    upper = refLevel
    lower = key
    countRight = manyToOne
    countWrong = oneToMany
  }
  
  else
  {
    upper = key
    lower = refLevel
    countRight = oneToMany
    countWrong = manyToOne
  }
  
  # ---------------------------------------------------------------------------
  # ------------------------------------ Return -------------------------------
  # ---------------------------------------------------------------------------
  
  linkInfoEntry = list(fromkey = upper, tokey = lower, type = linkType, 
                       countRight = countRight, countWrong = countWrong)
  
  return(linkInfoEntry)
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

## Unify levels from unifyLevels with the reference level
## 
## @param emuDBhandle
## @param levels
## @param unifyLevels
## @param linkIdxMap
## @param levelInfo
## @param refLevel
## @param bpfPath
## @param session
## @param bundle
## @keywords emuR BPF Emu
## @return levelInfo


unify_bpfLevels <- function(emuDBhandle,
                            levels,
                            unifyLevels,
                            linkIdxMap,
                            levelInfo,
                            refLevel,
                            bpfPath,
                            session,
                            bundle)
{
  # Start currentLabelIdx at 2, since refLevel already has one label (namely the refLevel's name).
  currentLabelIdx = 2
  newLabelsforRefLevel = NULL
  
  for(key in unifyLevels)
  {
    if(!key %in% names(levels))
    {
      next
    }
    
    seenLinks = list()
    
    for(idx in 1:length(levels[[key]]))
    {
      
      if(levels[[key]][[idx]][["links"]][1] %in% seenLinks)
      {
        stop("If you want to unify level ", key, " with the reference level, you cannot have more than one item on ", key, " pointing to one item on the reference level. BPF: ", bpfPath)
      }
      
      if(is.na(levels[[key]][[idx]][["links"]][1]))
      {
        stop("If you want to unify level, ", key, " with the reference level, it must not contain any link-less items. BPF: ", bpfPath)
      }
      
      if(!toString(levels[[key]][[idx]][["links"]][1]) %in% names(linkIdxMap))
      {
        stop("There is a symbolic link on level ", key, " in the following BPF that does not point to any item on the reference level: ", bpfPath)
      }
      
      for(labelKey in names(levels[[key]][[idx]][["labels"]]))
      {
        for(link in 1:length(levels[[key]][[idx]][["links"]]))
        {
          queryTxt = paste0("INSERT INTO labels VALUES('", emuDBhandle$UUID, "', '", session, "', '", 
                            bundle, "', ", linkIdxMap[[toString(levels[[key]][[idx]][["links"]][link])]], 
                            ", ", currentLabelIdx, ", '", labelKey, "', '", 
                            levels[[key]][[idx]][["labels"]][[labelKey]], "')")
            
          DBI::dbGetQuery(emuDBhandle$connection, queryTxt)
          seenLinks[[length(seenLinks) + 1L]] = levels[[key]][[idx]][["links"]][link]
        }
        
        if(!labelKey %in% newLabelsforRefLevel)
        {
          newLabelsforRefLevel = c(newLabelsforRefLevel, labelKey)
        }
        
        currentLabelIdx = currentLabelIdx + 1
      }
    }
  }
  
  for(idx in 1:length(levelInfo))
  {
    if(levelInfo[[idx]][["key"]] == refLevel)
    {
      levelInfo[[idx]][["labels"]] = c(levelInfo[[idx]][["labels"]], newLabelsforRefLevel)
    }
  }
  
  return(levelInfo)
}


# TODO: Find a better solution for the ";"-case (links to space in between items)
# TODO: Build syntax tree
# TODO: unify levels with levels other than the reference level
# TODO: unify levels that are not class 1
# TODO: OOP-Implementation to avoid passing/returning so many variables
