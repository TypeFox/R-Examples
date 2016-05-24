## ----basic---------------------------------------------------------------
library(scatterD3)
scatterD3(x = mtcars$wt, y = mtcars$mpg)

## ----basic_cust----------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, point_size = 15, point_opacity = 0.5, fixed = TRUE)

## ----labels--------------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, lab = rownames(mtcars), labels_size = 9)

## ----mapping-------------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, col_var = mtcars$cyl, symbol_var = mtcars$gear)

## ----map_size------------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, col_var = mtcars$cyl, size_var = mtcars$hp, 
          size_range = c(10,1000), point_opacity = 0.7)

## ----axis_limits---------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, xlim=c(0,10), ylim=c(10,35))

## ----cust_labels---------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, col_var = mtcars$cyl, symbol_var = mtcars$gear,
          xlab = "Weight", ylab = "Mpg", col_lab = "Cylinders", symbol_lab = "Gears")

## ----cust_tooltips-------------------------------------------------------
tooltips <- paste("This is an incredible <strong>", rownames(mtcars),"</strong><br />with ", 
                  mtcars$cyl, "cylinders !")
scatterD3(x = mtcars$wt, y = mtcars$mpg, tooltip_text = tooltips)

## ----ellipses------------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, ellipses = TRUE)

## ----ellipses_col--------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, col_var = mtcars$cyl, ellipses = TRUE)

## ----lasso---------------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, lab = rownames(mtcars), lasso = TRUE)

## ----lasso_callback------------------------------------------------------
scatterD3(x = mtcars$wt, y = mtcars$mpg, lab = rownames(mtcars), 
          lasso = TRUE,
          lasso_callback = "function(sel) {alert(sel.data().map(function(d) {return d.lab}).join('\\n'));}")

## ----cust_arrows---------------------------------------------------------
scatterD3(x = c(1, 0.9, 0.7, 0.2, -0.4, -0.5), xlab = "x",
          y = c(1, 0.1, -0.5, 0.5, -0.6, 0.7), ylab = "y",
          lab = LETTERS[1:6], type_var = c("point", rep("arrow", 5)),
          unit_circle = TRUE, fixed = TRUE, xlim = c(-1.2, 1.2))

