mmds.plot <- function (x,project = NULL, new.plot = TRUE, pdf.file = NULL) {

  if (!inherits(x, "mmds"))
    stop("object of class 'mmds' expected")

  if (project != NULL && !inherits(project, "project"))
    stop("object of class 'project' expected")

  if(is.null(pdf.file)){
	if (new.plot)
		dev.new()
	}
        else {
              	pdf(pdf.file,height=10,width=10)
        }

  if(length(x$active.coord) == 2) {
    layout(matrix(1:2))

    scree.plot(x$eigen.perc, lab = TRUE, title = "Scree plot of metric MDS", new.plot = FALSE)
    mmds.2D.plot(x,project, axis = c(1, 2), legend = FALSE, new.plot = FALSE)

  }
  else {
    layout(matrix(1:4, 2, 2))

    scree.plot(x$eigen.perc, lab = TRUE, title = "Scree plot of metric MDS", new.plot = FALSE)
    mmds.2D.plot(x,project, axis = c(1, 3), legend = FALSE, new.plot = FALSE)
    mmds.2D.plot(x,project, axis = c(1, 2), legend = FALSE, new.plot = FALSE)
    mmds.2D.plot(x,project, axis = c(2, 3), legend = FALSE, new.plot = FALSE)

  }
  if(!is.null(pdf.file)){
    dev.off()
    print("File created")
  }
}
