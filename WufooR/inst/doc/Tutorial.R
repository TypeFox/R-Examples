## ------------------------------------------------------------------------
library(WufooR)

options(Wufoo_Name = "johnmalc", Wufoo_API = "F1QH-Q64B-BSBI-JASJ")

## ------------------------------------------------------------------------
auth_name(NULL)
auth_key(NULL)

## ------------------------------------------------------------------------
t(user_info())

## ------------------------------------------------------------------------
t(form_info())

# Show responses to the form
fe_1 <- form_entries(formIdentifier = "z5kqx7h1gtvg4g")
t(fe_1)

sapply(fe_1, class)

## ------------------------------------------------------------------------
# How many responses did you get ?
form_entriesCount(formIdentifier = "z5kqx7h1gtvg4g", domain = "wufoo.eu")

## ------------------------------------------------------------------------
fields_info(formIdentifier = "z5kqx7h1gtvg4g", showRequestURL = TRUE)

## ------------------------------------------------------------------------
t(reports_info())

