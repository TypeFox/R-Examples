
# == title
# Global options for GetoptLong() 
#
# == param
# -... options, see 'details' section
# -RESET Whether to reset options to their default values
# -READ.ONLY only return read-only options?
# -LOCAL switch local mode
#
# == detail
# Supported options are following:
#
# -startingMsg message that will be printed before the helping message when running ``Rscript foo.R --help``
# -endingMsg message that will be printed after the helping message when running ``Rscript foo.R --help``
# -config configuration of ``Getopt::Long``, check http://perldoc.perl.org/Getopt/Long.html#Configuring-Getopt\%3a\%3aLong
#
# ``GetoptLong.options(...)`` should be put before calling `GetoptLong` function.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
GetoptLong.options = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
GetoptLong.options = setGlobalOptions(
	startingMsg = list(.value = "",
					   .length = 1),
	endingMsg = list(.value = "",
					 .length = 1),
	config = list(.value = NULL,
		          .class = "character"),
	"__argv_str__" = list(.value = NULL,
						.length = c(0, 1),
						.private = TRUE,
						.visible = FALSE),
	"__script_name__" = list(.value = NULL,
		                    .private = TRUE,
		                    .visible = FALSE)
)

# == title
# Global options for qq() related functions
#
# == param
# -... options, see 'details' section
# -RESET Whether to reset options to their default values
# -READ.ONLY only return read-only options?
# -LOCAL switch local mode
#
# == detail
# Supported options are following:
#
# -cat_prefix prefix of the string which is printed by `qqcat`
# -cat_verbose whether to print text by `qqcat`
# -code.pattern code pattern for variable interpolation
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# a = 1
# qq.options(cat_prefix = "[INFO] ")
# qqcat("a = @{a}\n")
# qq.options(cat_verbose = FALSE)
# qqcat("a = @{a}\n")
# qq.options(RESET = TRUE)
# qq.options(code.pattern = "`CODE`")
# qqcat("a = `a`\n")
# qq.options(RESET = TRUE)
qq.options = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
qq.options = setGlobalOptions(
	cat_prefix = list(.value = "",
					  .length = c(0, 1),
					  .class = c("character", "numeric", "NULL"),
					  .filter = function(x) {
						if(is.null(x)) {
							return('')
						} else {
							return(x)
						}
					  }),
	cat_verbose = list(.value = TRUE,
					   .class = "logical"),
	code.pattern = list(.value = "@\\{CODE\\}",
						.length = 1,
						.class = "character")
)
