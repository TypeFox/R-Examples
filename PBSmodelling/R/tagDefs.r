#Defined in tagDefs.r:
#
#In .tagDefs list:
#-The names of each valid tag (names(.tagDefs))
#   -tag names given should only contain lowercase letters
#-The valid option names for each tag (names(.tagDefs$<tagName>))
#-The class of each tag option (.tagDefs$<tagName>$<optionName>$class
#   -currently supported classes: character, logical, integer (all scalars)
#-A vector of names of defined validation checks for each option
# (.tagDefs$<tagName>$<optionName>$checks)
#   -these are defined in .tagOptionChecks
#   -validation of correct class type is handled separately
#-The default value for each option (.tagDefs$<tagName>$<optionName>$default)
#   -the absence of a default implies the option is required
#
#In .tagOptionChecks:
#-Names to identify each option validation check, which are used in .tagDefs
# checks vectors (names(.tagOptionChecks))
#-perl regular expression for the check (.tagOptionChecks$<checkName>$regexp)
#-error message to display if the check is failed
# (.tagOptionChecks$<checkName>$error)

.tagDefs=list()
.tagDefs$talk=list(
	"name"=list(class="character", checks=c("name")),
	"button"=list(class="logical", default=FALSE),
	"col"=list(class="integer", checks="positiveInteger", default=1)
)
.tagDefs$section=list(
	"name"=list(class="character", checks=c("name")),
	"button"=list(class="logical", default=FALSE),
	"col"=list(class="integer", checks="positiveInteger", default=2)
)
.tagDefs$code=list(
	"show"=list(class="logical", default=TRUE),
	"print"=list(class="logical", default=TRUE),
	"break"=list(class="character", checks="codeBreakOpt", default="print")
)
.tagDefs$text=list(
	"break"=list(class="logical", default=TRUE)
)
.tagDefs$file=list(
	"name"=list(class="character", checks=c("name")),
	"button"=list(class="logical", default=FALSE),
	"col"=list(class="integer", checks="positiveInteger", default=3),
	"break"=list(class="logical", default=TRUE)
)

#checks are all case insensitive
.tagOptionChecks=list(
	name=list(regexp="^[a-z][a-z0-9_]*$", error=paste("Name values must contain",
		"only alphanumeric and underscore ('_') characters and start with a",
		"letter.")),
	codeBreakOpt=list(regexp="^(show)|(print)|(all)|(none)$",
		error=paste("Invalid choice for break option; must be \"show\", \"print\",",
		"\"all\", or \"none\".")),
	positiveInteger=list(regexp="^[1-9]\\d*$", error=paste("Option value must be",
		"positive integer."))
)

#tag hierarchy: a tree-- talk is the root, each parent is a list containing its
#children tags, and leaves are character scalars.
#a tag must be nested within its parent tag
#a tag can have more than one type of parent tag
#(the parent tag is given as an argument to handleTag.tagName and so the
#tag can be handled differently depending on its parent)
#NOT CURRENTLY IMPLEMENTED
#.tagHierarchy=list(
#	talk=list(
#		section=list(c("code", "text", "file")),
#		"file"
# )
#)
	