###############################################################################
# Mango Solutions, Chippenham SN14 0SQ 2008
# BLCOPOptions
# Author: Francisco
###############################################################################
# DESCRIPTION: Sets or retrieves the package's global settings .  See the online help for these.
# KEYWORDS: environment
###############################################################################

BLCOPOptions <- function
(
  opt,             # string with the option to be retrieved or changed
  setting          # new setting for the value.  Note: currently not error checked!
)
{
  if(missing(opt)&missing(setting))
    return(.BLEnv$settings)
  if(missing(setting))
    return(.BLEnv$settings[[opt]])
  .BLEnv$settings[[opt]] <- setting   
}