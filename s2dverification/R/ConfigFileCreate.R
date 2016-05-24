ConfigFileCreate <- function(file_path, confirm = TRUE) {
  success <- ConfigFileSave(list(definitions = list(
    DEFAULT_EXP_MAIN_PATH = "$EXP_NAME$", 
    DEFAULT_EXP_FILE_PATH = "$STORE_FREQ$/$VAR_NAME$_$START_DATE$.nc",
    DEFAULT_NC_VAR_NAME = "$VAR_NAME$", 
    DEFAULT_SUFFIX = "", DEFAULT_VAR_MIN = "", 
    DEFAULT_VAR_MAX = "", DEFAULT_OBS_MAIN_PATH = "$OBS_NAME$", 
    DEFAULT_OBS_FILE_PATH = "$STORE_FREQ$/$VAR_NAME$_$YEAR$$MONTH$.nc", 
    DEFAULT_DIM_NAME_LONGITUDES = "longitude", DEFAULT_DIM_NAME_LATITUDES = "latitude", 
    DEFAULT_DIM_NAME_MEMBERS = "ensemble")), file_path, confirm = confirm)
  if (success) {
    cat("WARNING: You have just created an empty configuration file. You can edit it with ConfigAddEntry(). You can edit the defaults according to your needs with the functions ConfigFileOpen(), ConfigEditDefinition() and ConfigFileSave() or edit the file manually as specified in ?ConfigFileOpen.\n")
  }
}
