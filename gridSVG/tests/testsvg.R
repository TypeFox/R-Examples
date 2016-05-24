library(grid)
library(gridSVG)
svgdev <- svgOpen(width=400, height=400)
svgRect(0, 0, 400, 400, 
        style=svgStyle(fill="none", stroke="black"),
        svgdev=svgdev)
svgStartGroup(svgdev=svgdev)
svgRect(2, 2, 240, 320,
        style=svgStyle(fill="none", stroke="green"),
        svgdev=svgdev)
svgLines(c(23.81, 45.63, 67.45, 89.27, 111.09, 132.90,
           154.72, 176.54, 198.36, 220.18),
         c(292.90, 263.81, 234.72, 205.63, 176.54,
           147.45, 118.36, 89.27, 60.18, 31.09),
         style=svgStyle(stroke="green"),
         svgdev=svgdev)
svgPolygon(c(23.81, 67.45, 89.27, 23.81),
           c(31.09, 31.09, 147.45, 118.36),
           style=svgStyle(fill="grey"),
           svgdev=svgdev)
svgRect(132.90, c(89.27, 205.63), 43.63, 29.09,
        style=svgStyle(fill="cyan"),
        svgdev=svgdev)
svgText(45.63, c(234.728, 205.63),
        c("some text", "some more text!"),
        style=svgStyle(fill="red"),
        svgdev=svgdev)
svgCircle(176.54, 89.27, c(2.18, 43.63),
          style=svgStyle(stroke="blue", fill="none"),
          svgdev=svgdev)
svgText(89.27, 147.45, "centred text", hjust="centre", vjust="centre", rot=20,
        style=svgStyle(fill="yellow", stroke="black"),
        svgdev=svgdev)


  svgStartGroup(svgdev=svgdev)
  svgRect(132.90, 147.45, 65.45, 29.09,
          style=svgStyle(fill="none", stroke="black"),
          svgdev=svgdev)
  svgText(139.45, 162, "text in a box",
          svgdev=svgdev)
  svgEndGroup(svgdev=svgdev)


svgRect(111.09, 60.18, 43.63, 203.63,
        style=svgStyle(fill="green", opacity=.5),
        svgdev=svgdev)
svgEndGroup(svgdev=svgdev)
svgClose(svgdev)

