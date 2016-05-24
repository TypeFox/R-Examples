##' Provides a GUI panel offering dynamic control of the adjust-to settings
##' of the partial effects plot of an rms fit object.
##'
##' The parameters are identified automatically by inspection of the fit object,
##' and widgets appropriate to their data types are chosen automatically.
##' @param fit The \code{rms} fit object to be visualized
##' @param datadist Optionally, a \code{datadist} for the fitted object may be
##' provided, conveniently enabling use of this function with models in which
##' logical covariates are not typed as such, but represented in {0,1}, or in
##' which covariates more properly coded as ordered factors are represented
##' instead as a finite number of integer values.
##' @return A controller panel is popped up on the screen, along with a partial
##' effects plot for the \code{fit} object. No value is returned.
##' @author David C. Norris
##' @keywords iplot dynamic
##' @export controller
## TODO: Gain control over LAYOUT and SIZING of the control panel, if possible.
##       Horizontally oriented radiogroups might be a good start, but seem
##       unsupported by the 'rpanel' package.
## TODO: Extend to controlling inputs to (precomputed) simulation exercises.
## TODO: Generalize to 'ols' models, including a confidence interval type selector:
##       rp.radiogroup(summary.panel, conf.type, c("mean","individual"),
##                     title="Conf. interval", action=summary.draw)
## TODO: Speed the creation of the controller, perhaps by delaying invocation
##       of the 'action' as it is built -- if rpanel supports this.
controller <- function(fit, datadist=NULL){
  summary.draw <- function(panel){
    vars <- intersect(fit$Design$name, names(panel))
    fit$Design$limits['Adjust to', vars] <- panel[vars]
    print(plot(Predict(fit, fun=plogis)))
    panel
  }
  limits <- fit$Design$limits
  adjust.to <- as.list(limits['Adjust to',])
  args <- c(title="Adjust-to values", adjust.to)
  summary.panel <- do.call("rp.control", args)
  ## Unfortunately, it seems that an 'lrm' fit, at least, has a limited $values
  ## component, which does not acknowledge numeric vars with finitely many levels.
  ## This is why I must employ an optional 'datadist' argument, in order to
  ## support use of unrefined models in which ordered factors and logical
  ## variables and represented as integer values.
  all.wholenumbers <- function(x){
    !is.null(x) && all(abs(x - round(x)) < .Machine$double.eps)
  }
  for(var in names(limits)){
    varlabel <- fit$Design$label[match(var, names(limits))]
    values <- datadist$values[[var]] # NB: this will be NULL if datadist==NULL
    if(length(values)==2 && all(values %in% 0:1)
       || all(limits[[var]] %in% 0:1)) # 1. Recognize logical var, and implement checkbox control
      do.call('rp.checkbox', list(panel=summary.panel,
                                  var=as.name(var),
                                  action=summary.draw,
                                  title=varlabel))
    else if(is.factor(limits[[var]])) # 2. Recognize factors, and treat as radiogroups
      do.call('rp.radiogroup', list(panel=summary.panel,
                                    var=as.name(var),
                                    values=fit$Design$values[[var]],
                                    ## TODO: Why does radiogroup not show its
                                    ##       initval as selected?  Sliders and
                                    ##       checkboxes DO, by contrast, show
                                    ##       the initialized values.  The
                                    ##       problem occurs even when I
                                    ##       explicitly set 'initval' here.
                                    title=varlabel,
                                    action=summary.draw))
    else # 3. Default is numeric variable, treated via slider
      do.call('rp.slider', list(panel=summary.panel,
                                var=as.name(var),
                                from=limits['Low', var],
                                to=limits['High', var],
                                showvalue=TRUE,
                                action=summary.draw,
                                title=varlabel,
                                ## use discrete slider if there are finitely many, entirely integral values,
                                resolution=as.integer(all.wholenumbers(values))))
  }
}
