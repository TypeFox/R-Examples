portfolio <-
function (x, y, z, label_x="X", label_y="Y", heading="Portfolio", pcol="given", colsp=0, leg=FALSE, leg_vec=0, leg_fsize=1, leg_x=-max_val, leg_y=-max_val/2) {

x_count <- length(x)
max_val <- max(x,y)

point_col <- character()
if (pcol=="random") { point_col <- sample(colours(), x_count) }

if (pcol=="given") {
for (i in 1:x_count) {
point_col[i] <- as.character(colsp[[i]])
}
}

point_size <- (z/max(z))*10
max_val_inc <- 0
plot(x,y, cex=point_size, xlim=c(-max_val,max_val), ylim=c(-max_val-max_val_inc,max_val), pch=16, col=point_col, xlab=label_x, ylab=label_y, main=heading)
abline (h=0, v=0)

if (leg==TRUE) {
legend(leg_x, leg_y, leg_vec, point_col, cex=leg_fsize, ncol=2, bg="white")
}
}
