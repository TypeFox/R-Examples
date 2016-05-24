## ---- echo = FALSE, message = FALSE--------------------------------------
library(knitr)
library(ggvis)
opts_chunk$set(comment = "#>", error = FALSE, tidy = FALSE)

## ------------------------------------------------------------------------
merge_props <- ggvis:::merge_props
merge_props(props(x = ~x), props(y = ~y))
merge_props(props(x = ~a), props(x = ~b))
merge_props(props(x = ~a, y = ~a), props(x = ~b, inherit = FALSE))

## ---- eval = FALSE-------------------------------------------------------
#  ggplot(Minard.cities, aes(x = long, y = lat)) +
#    geom_path(
#      aes(size = survivors, colour = direction, group = group),
#      data = Minard.troops
#    ) +
#    geom_point() +
#    geom_text(aes(label = city), hjust=0, vjust=1, size=4)

## ---- eval = FALSE-------------------------------------------------------
#  ggvis(data = NULL, x = ~long, y = ~lat) %>%
#    layer_points(size = ~survivors, stroke = ~direction, data = Minard.troops) %>%
#    layer_text(text := ~city, dx := 5, dy := -5, data = Minard.cities)

