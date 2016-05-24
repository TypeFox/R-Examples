{###############################################################################
# test-styles.assignment.R
# (c)2012 Andrew Redd 
# This is file part of the lint R package, a code style check package for R.
# 
# DESCRIPTION
# -----------
# Tests for assignment styles defined in R/styles.assignment.R
# 
}###############################################################################

context("Styles/Assignment")
autotest_style('styles.assignment.noeq')
autotest_style('styles.assignment.norightassign')
autotest_style('styles.assignment.notinfcall')


