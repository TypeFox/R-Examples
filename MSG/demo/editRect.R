library(grid)
grid.rect(x = 0, y = 0, width = 0.1, height = 0.1,
    gp = gpar(col = NA, fill = "red"), name = "rect0")
grid.rect(x = 0.1, y = 0.9, width = 0.1, height = 0.1,
    gp = gpar(col = NA, fill = "green"), name = "rect1")
for (i in 1:100) {
    grid.edit("rect0", x = unit(i/100, "npc"), y = unit(i/100,
        "npc"), gp = gpar(fill = rainbow(100)[i]))
    Sys.sleep(0.05)
}
