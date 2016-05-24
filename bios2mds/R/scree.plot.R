scree.plot <- function (x, lab = FALSE, title = NULL, xlim = NULL,
  ylim = NULL, new.plot = TRUE, pdf.file = NULL) {

if(is.null(pdf.file)){
	if (new.plot)
		dev.new()
	}
        else {
              	pdf(pdf.file,height=10,width=10)
        }

  #default names of axis
  x.lab <- paste("Component number")
  y.lab <- paste("Eigenvalue %")
  
  names(x) <- seq_len(length(x))

  #create space to display
  #the first label
  if (lab) {
    if (is.null(ylim))
      ylim <- c(0, max(x) * 1.2)
  }
    
  if (is.null(xlim) && is.null(ylim)) {
    #get back positions of bars
    pos <- barplot(x, xlab = x.lab, ylab = y.lab, xlim = xlim, ylim = ylim)
  }
  else {
     pos <- barplot(x, xlab = x.lab, ylab = y.lab, xlim = xlim, ylim = ylim, xpd = FALSE)
  }
  if (is.null(title))
    title ("Scree plot")
  else
    title (title)

  if (lab)
    text(pos, x, labels = x, pos = 3)

        if(!is.null(pdf.file)){
                dev.off()
                print("File created")
        }

}