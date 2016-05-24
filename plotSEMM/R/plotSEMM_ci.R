plotSEMM_ci <- function(SEMLIdatapks, linesearch, lnty = 3, lncol = 1, deltaci=TRUE,
                        deltace=TRUE, bsci=TRUE, ninty_five = TRUE, include_title = TRUE,
                        use_fixed_value, cex = 1.5, ...) {
    nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), c(4, 4), c(4, 1), TRUE)
    par(mar=c(5,5,5,2))
    if(use_fixed_value){
        fixed_values <- SEMLIdatapks[nrow(SEMLIdatapks), ]
        SEMLIdatapks <- SEMLIdatapks[-nrow(SEMLIdatapks), ]
    }

    dots <- list(...)
    input <- dots$input
    if(!is.null(input$xlab)) xlab <- input$xlab
    if(!is.null(input$ylab)) ylab <- input$ylab
    legend_location <- if(!is.null(input$legend_location)) input$legend_location else 'bottomleft'
    if(legend_location == 'default') legend_location <- 'bottomleft'

    #requires setup from plotSEMM_setup2
    if(!SEMLIdatapks$setup2[1L])
        stop('plotSEMM_ci requires a setup model that included the parameter ACOV matrix')
    def.par <- par(no.readonly = TRUE)
    pick <- SEMLIdatapks$agg_denEta1 > 0.02
    SEMLIdatapks <- SEMLIdatapks[pick, ]
    legendkeep <- c(TRUE, deltaci, deltace, bsci)

    # plot(SEMLIdatapks$Eta1,SEMLIdatapks$y,type='n',xlab='Latent Predictor', ylab='Latent Outcome')
    plot(SEMLIdatapks$Eta1, SEMLIdatapks$Eta2, type = "n", xlab = xlab, ylab = ylab,
         cex.lab=cex, cex.axis=cex)
    lines(SEMLIdatapks$Eta1, SEMLIdatapks$agg_pred, col = 1, lwd = 2)
    if(deltaci){
        points(SEMLIdatapks$Eta1, SEMLIdatapks$delta_CIlo, col = 2, lwd = 1, lty = 2)
        points(SEMLIdatapks$Eta1, SEMLIdatapks$delta_CIhi, col = 2, lwd = 1, lty = 2)
    }
    if(deltace){
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$delta_CElo, col = 2, lwd = 2, lty = 2)
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$delta_CEhi, col = 2, lwd = 2, lty = 2)
    }
    if(bsci){
        points(SEMLIdatapks$Eta1, SEMLIdatapks$bs_CIlo, col = 4, lwd = 1, lty = 3, pch = 4)
        points(SEMLIdatapks$Eta1, SEMLIdatapks$bs_CIhi, col = 4, lwd = 1, lty = 3, pch = 4)
    }
    if(ninty_five){
        legend = c("Aggregate Function", "Delta Method 95% Confidence Interval",
                   "Delta Method 95% Confidence Envelope",
                   "Bootstrap 95% Confidence Interval")[legendkeep]
    } else {
        legend = c("Aggregate Function", "Delta Method 90% Confidence Interval",
                   "Delta Method 90% Confidence Envelope",
                   "Bootstrap 90% Confidence Interval")[legendkeep]
    }

    lwd = c(2, 1, 2, 1)[legendkeep]
    lty = c(1, 0, 2, 0) [legendkeep]
    pch = c(NA, 1, NA, 4)[legendkeep]
    col = c(1, 2, 2, 4)[legendkeep]

    if(SEMLIdatapks$boot[1]){
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$bs_CElo, col = 4, lwd = 2, lty = 4, pch = 4)
        lines(SEMLIdatapks$Eta1, SEMLIdatapks$bs_CEhi, col = 4, lwd = 2, lty = 4, pch = 4)
        if(ninty_five) legend <- c(legend, 'Bootstrap 95% Confidence Envelope')
        else legend <- c(legend, 'Bootstrap 90% Confidence Envelope')
        lwd <- c(lwd, 2); lty <- c(lty, 4); pch = c(pch, NA); col = c(col, 4)
    }
    if(legend_location != 'none')
        legend(legend_location, legend = legend, lwd = lwd, lty = lty,
               pch = pch, col = col, bty = "n", cex=cex)
    if(linesearch){
        found <- FALSE
        search <- attr(SEMLIdatapks, 'search')
        if(length(search) && deltace){
            found <- TRUE
            nc <- ncol(search)
            lines(c(search[1,1], search[1,nc]), c(search[2,1], search[2,nc]),
                  col = 'pink', lwd=3, lty=1)
        }
        if(SEMLIdatapks$boot[1]){
            search <- attr(SEMLIdatapks, 'search.bs')
            if(length(search)){
                found <- TRUE
                nc <- ncol(search)
                lines(c(search[1,1], search[1,nc]), c(search[2,1], search[2,nc]),
                      col = 'lightblue', lwd=3, lty=1)
            }
        }
        if(!found && include_title)
            title('No Line was Found within the Confidence Envelope(s)', cex.main=cex*(4/3))
    }
    if(use_fixed_value){
        #plot some text
        percent <- '90%'
        if(ninty_five) percent <- '95%'
#         points(fixed_values$Eta1, (fixed_values$delta_CIlo + fixed_values$delta_CIhi)/2,
#                col='purple', cex=1.5, pch=17)
	    par(mar = c(0, 0, 0, 0))
        plot(c(0,5), c(0,5), axes=FALSE, frame.plot=FALSE, type='n', xlab='', ylab='')
        txt <- paste0('Latent Predictor Value: ', fixed_values$Eta1)
        txt <- c(txt, paste0('Predicted Latent Outcome Value: ',
                       round((fixed_values$delta_CIlo + fixed_values$delta_CIhi)/2, 4)))
        if(deltaci)
            txt <- c(txt, paste0(percent, ' Delta Method Wald-type Confidence Interval: (',
                                 round(fixed_values$delta_CIlo, 3), ', ',
                     round(fixed_values$delta_CIhi, 3), ')'))
        if(bsci)
            txt <- c(txt, paste0(percent, ' Parametric Bootstrap Confidence Interval: (',
                                   round(fixed_values$bs_CIlo, 3), ', ',
                            round(fixed_values$bs_CIhi, 3), ')'))
        if(!deltaci && !bsci)
            txt <- c(txt, paste0(percent, ' Confidence interval: NA NA (No method selected)'))
        legend('topleft', legend=txt, cex=cex, bty='n')
    }
}
