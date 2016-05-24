# I don't think this function is needed.  It is called from the odfTable
# methods to set the type of the cells in a table.  However, I think
# this is only really needed for spreadsheets.  And if anything but
# "string" is specified as the data type, then you have to include the
# value of the cell with the "office:value" attribute, which we could
# do, but I don't know if there's any point.  In fact, I don't think
# that Open Office Writer even allows you to set the types of the
# table cells.  That only seems to be supported by Calc.
#
# For now, I'm leaving the old code in, but adding in a final line
# that always returns the value "string".

"odfDataType" <-
function(x)
{
   internalMode <- typeof(x)
   dataMode <- internalMode
   dataMode[dataMode == "character"] <- "string"
   dataMode[dataMode == "integer"] <- "float"
   dataMode[dataMode == "double"] <- "float"
   # possible office:value-type values are:
   #float, percentage, currency, date, time, boolean, string
   if (any(internalMode == "integer" & is.factor(x)))
      dataMode[internalMode == "integer" & is.factor(x)] <- "string"
   dataMode

   # Forget everything we just did: always return the type "string".
   # See the notes at the top of the file for my reasons.
   "string"
}
