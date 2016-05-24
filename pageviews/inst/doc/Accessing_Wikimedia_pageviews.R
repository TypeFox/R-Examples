## ---- eval = FALSE-------------------------------------------------------
#  str(article_pageviews(project = "de.wikipedia", article = "R_(Programmiersprache)",
#                        start = "2015110100", end = "2015110300"))
#  'data.frame':	3 obs. of  6 variables:
#   $ project  : chr  "de.wikipedia" "de.wikipedia" "de.wikipedia"
#   $ article  : chr  "R_(Programmiersprache)" "R_(Programmiersprache)" "R_(Programmiersprache)"
#   $ timestamp: chr  "2015110100" "2015110200" "2015110300"
#   $ access   : chr  "all-access" "all-access" "all-access"
#   $ agent    : chr  "all-agents" "all-agents" "all-agents"
#   $ views    : num  308 536 537

## ---- eval = FALSE-------------------------------------------------------
#  str(project_pageviews())
#  'data.frame':	1 obs. of  6 variables:
#   $ project    : chr "en.wikipedia"
#   $ access     : chr "all-access"
#   $ agent      : chr "all-agents"
#   $ granularity: chr "daily"
#   $ timestamp  : chr "2015100100"
#   $ views      : num 2.72e+08

## ---- eval = FALSE-------------------------------------------------------
#  > str(top_articles())
#  'data.frame':	1000 obs. of  8 variables:
#   $ article: chr  "Main_Page" "Special:Search" "Special:BlankPage" "-" ...
#   $ views  : int  18840697 3191975 1862191 1660878 293537 289710 271152 195670 163707 124751 ...
#   $ rank   : int  1 2 3 4 5 6 7 8 9 10 ...
#   $ project: chr  "en.wikipedia" "en.wikipedia" "en.wikipedia" "en.wikipedia" ...
#   $ access : chr  "all-access" "all-access" "all-access" "all-access" ...
#   $ year   : chr  "2015" "2015" "2015" "2015" ...
#   $ month  : chr  "10" "10" "10" "10" ...
#   $ day    : chr  "01" "01" "01" "01" ...

