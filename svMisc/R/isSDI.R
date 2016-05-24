isSDI <- function ()
{
	# This function is specific to Windows, but it is defined everywhere
	# so that we don't have to test the platform before use!
	# Check if Rgui was started in SDI mode (needed by some GUI clients)

	# First, is it Rgui?
	if (!isRgui()) return(FALSE)
    # RGui SDI mode: returns "R Console", in MDI mode: returns "RGui"
    if (getIdentification() == "R Console") TRUE else FALSE
}
