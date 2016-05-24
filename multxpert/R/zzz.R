# Multxpert package implements commonly used p-value-based and parametric
# multiple testing procedures and parallel gatekeeping procedures

# For more information about the Multxpert package, visit the Multiplicity Expert web site
# http://multxpert.com/wiki/MultXpert_package

.onLoad <- function(lib, pkg) {
	if (interactive()) {
		cat('multxpert: Implementation of commonly used p-value based and parametric\n',
				'multiple testing procedures and parallel gatekeeping procedures.\n',
				'For more information visit http://multxpert.com/wiki/MultXpert_package\n', sep='')
	}
}