.onAttach <- function(lib, pkg)
{
	pkg.info <- utils::packageDescription('autopls')
	packageStartupMessage ('autopls ', utils::packageDescription ('autopls',
                          field= "Version"), appendLF = TRUE)
}
