Venn <- function (A, B, AandB = 0, labA = "A", labB = "B", cex.text = .95, ...) 
{
    if(A > 1) stop("Violation of probability rules! \nP(A) cannot be greater than 1")
    if(B > 1) stop("Violation of probability rules! \nP(B) cannot be greater than 1")   
    if(AandB > 1) stop("Violation of probability rules! \nP(AandB) cannot be greater than 1") 
    if(A + B - AandB > 1)stop("Violation of probability rules! \nP(A) + P(B) -P(AandB) cannot be greater than 1")  
    # if(A + B + AandB > 1) stop("Violation of probability rules! \nP(A) + P(B) + P(AandB) cannot be greater than 1")   
    if(A < 0) stop("Violation of probability rules! \nP(A) cannot be less than 0")  
    if(B < 0) stop("Violation of probability rules! \nP(B) cannot be less than 0")
    if(AandB < 0) stop("Violation of probability rules! \nP(AandB) cannot be less than 0")    
    if(AandB > A | AandB > B) stop("Violation of probability rules! \nP(AandB) cannot be greater than P(A) or P(B)")

    if (A == 0 & B == 0) {
        S <- plot(seq(0, 1), seq(0, 1), type = "n", xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "",...)
        text(0.5, 0.6, bquote(paste(italic(P), "(", italic(.(labA)), 
            ") = 0", sep = "")), cex = cex.text)
        text(0.5, 0.4, bquote(paste(italic(P), "(", italic(.(labB)), 
            ") = 0", sep = "")), cex = cex.text)
    }
    if (A + B == 1 & AandB == 0) {
        S <- plot(seq(0, 1), seq(0, 1), type = "n", xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", bty = "n",...)
        rect(c(0, A), c(0, 0), c(A, 1), c(1, 1), col = c(rgb(0.7, 
            0.7, 0.7, 0.8), rgb(0.3, 0.3, 0.3, 0.8)))
        text(A/2, 0.5, bquote(paste(italic(P), "(", italic(.(labA)), 
            ") = ", .(A), sep = "")), cex = cex.text)
        text(A + ((1 - A)/2), 0.5, bquote(paste(italic(P), "(", 
            italic(.(labB)), ") = ", .(B), sep = "")), cex = cex.text)
    }
    if (A + B == 1 & AandB > 0) {
        S <- plot(seq(0, 1), seq(0, 1), type = "n", xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", bty = "n",...)
        rect(c(0, 0.5, A-(AandB/2)), c(0, 0, 0), c(A, 1, A+(AandB/2)), c(1, 1, 1), col = c(rgb(0.7, 
            0.7, 0.7, 0.8), rgb(0.3, 0.3, 0.3, 0.8), rgb(0.3, 0.3, 0.3, 0.8)))
        text(A/2, 0.5, bquote(paste(italic(P), "(", italic(.(labA)), 
            ") = ", .(A), sep = "")), cex = cex.text)
        text(A + ((1 - A)/2), 0.5, bquote(paste(italic(P), "(", 
            italic(.(labB)), ") = ", .(B), sep = "")), cex = cex.text)
        mtext(bquote(paste(italic(P), "(", italic(.(labA)), intersect(""), 
            italic(.(labB)), ") = ", .(AandB), sep = "")), side = 3, 
            at = (A - (AandB/2) + A + (AandB/2))/2, cex = cex.text)
        arrows((A - (AandB/2) + A + (AandB/2))/2, 1.01, (A - (AandB/2) + A + (AandB/2))/2, 0.75, length = 0.05)
    }
    if ((AandB == A | AandB == B) & (A != 0 | B != 0)) {
        S <- plot(seq(-0.55, 0.55, length = 3), seq(-0.55, 0.55, 
            length = 3), type = "n", xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "",...)
        r.A <- sqrt(A/pi)
        r.B <- sqrt(B/pi)
        draw.circle(0, 0, radius = r.A, col = rgb(blue = 0.7, 
            red = 0.7, green = 0.7, alpha = 0.8))
        draw.circle(0, 0, radius = r.B, col = rgb(blue = 0.3, 
            red = 0.3, green = 0.3, alpha = 0.8))
        if (A >= B) {
            text(-0.5, 0.5, bquote(paste(italic(P), "(", italic(.(labA)), 
                ") = ", .(A), sep = "")), cex = cex.text)
            text(0, 0, bquote(paste(italic(P), "(", italic(.(labB)), 
                ") = ", italic(P), "(", italic(A), intersect(""), 
                italic(B), ") = ", .(B), sep = "")), cex = cex.text)
        }
        if (B > A) {
            text(-0.5, 0.5, bquote(paste(italic(P), "(", italic(.(labB)), 
                ") = ", .(B), sep = "")), cex = cex.text)
            text(0, 0, bquote(paste(italic(P), "(", italic(.(labA)), 
                ") = ", italic(P), "(", italic(A), intersect(""), 
                italic(B), ") = ", .(A), sep = "")), cex = cex.text)
        }
        arrows(-0.5, 0.46, -mean(c(r.A, r.B)), 0, length = 0.05)
    }
    if (A + B != 1 & AandB == 0 & (A != 0 | B != 0)) {
        S <- plot(seq(0, 2), seq(0, 2), type = "n", xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", bty = "n",...)
        r.A <- sqrt(A/pi)
        r.B <- sqrt(B/pi)
        draw.circle(0.5, 1, radius = r.A, col = rgb(blue = 0.7, 
            red = 0.7, green = 0.7, alpha = 0.8))
        draw.circle(1.25, 1, radius = r.B, col = rgb(blue = 0.3, 
            red = 0.3, green = 0.3, alpha = 0.8))
        text(0.5, 1, bquote(paste(italic(P), "(", italic(.(labA)), 
            ") = ", .(A), sep = "")), cex = cex.text)
        text(1.25, 1, bquote(paste(italic(P), "(", italic(.(labB)), 
            ") = ", .(B), sep = "")), cex = cex.text)
        rect(0.5 - sqrt(1/pi), 1 - sqrt(1/pi), 1.25 + sqrt(1/pi), 
            1 + sqrt(1/pi))
    }
    if ((A + B != 1 & AandB != 0) & (A != AandB) & (B != AandB)) {
        S <- plot(seq(-1, 1), seq(-1, 1), type = "n", xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "", bty = "n",...)
        r.A <- sqrt(A/pi)
        r.B <- sqrt(B/pi)
        r.diff <- r.A - r.B
        d <- 2 * ((1 - AandB) * r.A) - r.diff
        draw.circle(-0.5, 0, radius = r.A, col = rgb(blue = 0.7, 
            red = 0.7, green = 0.7, alpha = 0.8))
        draw.circle(d - 0.5, 0, radius = r.B, col = rgb(blue = 0.3, 
            red = 0.3, green = 0.3, alpha = 0.8))
        text(-0.5, 0, bquote(paste(italic(P), "(", italic(.(labA)), 
            ") = ", .(A), sep = "")), cex = cex.text)
        text(d - 0.5, 0, bquote(paste(italic(P), "(", italic(.(labB)), 
            ") = ", .(B), sep = "")), cex = cex.text)
        x <- ((d + r.diff)/2) - 0.5
        arrows(x, 0.6, x, 0, length = 0.05)
        text(x, 0.65, bquote(paste(italic(P), "(", italic(.(labA)), 
            intersect(""), italic(.(labB)), ") = ", .(AandB), 
            sep = "")), cex = 0.8)
        rect(-0.5 - sqrt(1/pi), sqrt(1/pi), (d - 0.5) + sqrt(1/pi), 
            -sqrt(1/pi))
    }
}

Venn.tck<-function()
{

    local({
        have_ttk <- as.character(tcl("info", "tclversion")) >= 
            "8.5"
        if (have_ttk) {
            tkbutton <- ttkbutton
            tkcheckbutton <- ttkcheckbutton
            tkentry <- ttkentry
            tkframe <- ttkframe
            tklabel <- ttklabel
            tkradiobutton <- ttkradiobutton
        }
        tclServiceMode(FALSE) 
        dialog.sd <- function() {
            tt <- tktoplevel()
            tkwm.title(tt, "Venn diagrams")
            A.entry <- tkentry(tt, textvariable = Ae, width =10)
            B.entry <- tkentry(tt, textvariable = Be, width =10)
            AandB.entry <- tkentry(tt, textvariable = AandBe, width =10)
            A.label.entry <- tkentry(tt, textvariable = Alabel, width =10)
            B.label.entry <- tkentry(tt, textvariable = Blabel, width =10)
            
            done <- tclVar(0)
 
            reset <- function() {
                tclvalue(Ae) <- ""
                tclvalue(Be) <- ""
                tclvalue(AandBe) <- ""
                tclvalue(Alabel)<-"A"
                tclvalue(Blabel)<-"B"
            }
            reset.but <- tkbutton(tt, text = "Reset", command = reset)
            submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)
            build <- function() {
                A <-tclvalue(Ae)
                B <- tclvalue(Be)
                AandB <- tclvalue(AandBe)
                labA <- tclvalue(Alabel)
                labB <- tclvalue(Blabel)
                substitute(Venn(A=as.numeric(A),B=as.numeric(B),AandB=as.numeric(AandB),labA=labA,labB=labB))
            }
            
            tkgrid(tklabel(tt, text = "Venn probability diagrams"), 
                columnspan = 2)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = "    P(A)"), A.entry, sticky = "w")
            tkgrid(tklabel(tt, text = "    P(B)"), B.entry, sticky = "w")
            tkgrid(tklabel(tt, text = "    P(A\u2229B)"), AandB.entry, sticky = "w")
            tkgrid(tklabel(tt, text = "    A label"), A.label.entry,  sticky = "w")
            tkgrid(tklabel(tt, text = "    B label"), B.label.entry, sticky = "w")
            tkgrid(tklabel(tt, text = ""))
            tkgrid(submit.but, reset.but,sticky="w")
            tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
            tkwait.variable(done)
            if (tclvalue(done) == "2") 
                stop("aborted")
            tkdestroy(tt)
            cmd <- build()
            eval.parent(cmd)
        }
        Ae <- tclVar("0.4")
        Be <- tclVar("0.3")
        AandBe <- tclVar("0.1")
        Alabel <- tclVar("A")
        Blabel <-tclVar("B")
        dialog.sd()
    })
    }
