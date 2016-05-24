setGeneric('losses', function(object){standardGeneric('losses')})

setMethod('losses',
          signature=(object='Gef'),
          definition=function(object){
            dat <- as.data.frameY(object, complete=TRUE)
            isShd=('Gef0d' %in% names(dat)) ##is there shadows?
            if (isShd) {
              shd <- with(dat, mean(1-Gefd/Gef0d))
              eff <- with(dat, mean(1-Gef0d/Gd))
            } else {
              shd <- 0
              eff <- with(dat, mean(1-Gefd/Gd))
            }
            result <- data.frame(id=c('Shadows', 'AoI'), values=c(shd, eff))
            result
          }
          )

setMethod('losses',
          signature=(object='ProdGCPV'),
          definition=function(object){
            DayOfMonth=c(31,28,31,30,31,30,31,31,30,31,30,31) ###OJO
            dat <- as.data.frameY(object, complete=TRUE)
            module0=object@module
            module0$CoefVT=0 ##No losses with temperature
            ## p0 <- prodGCPV(lat=object@lat, modeTrk=object@modeTrk,
            ##                modeRad='prev', prev=object,
            ##                module=module0, generator=object@generator,
            ##                inverter=object@inverter, effSys=object@effSys)
            ## p0Y <- as.data.frameY(p0)
            ## temp <- mean(1-dat$Edc/p0Y$Edc)
            Pg=object@generator$Pg
            Nm=1/sample2Hours(object@sample)
            datI <- as.zooI(object, complete=TRUE)
            if (object@type=='prom'){
              YfDC0=sum(monthlySum(datI$Vmpp*datI$Impp)/Pg*DayOfMonth)
              YfAC0=sum(monthlySum(datI$Pdc*datI$EffI)/Pg*DayOfMonth)
            } else {
              YfDC0 <- yearlySum(datI$Vmpp*datI$Impp)/Pg
              YfAC0 <- yearlySum(datI$Pdc*datI$EffI)/Pg
            }
            gen <- mean(1-YfDC0/dat$Gefd)
            YfDC <- dat$Edc/Pg*1000
            DC=mean(1-YfDC/YfDC0)
            inv=mean(1-YfAC0/YfDC)
            AC=mean(1-dat$Yf/YfAC0)
            result0 <- losses(as(object, 'Gef'))
            result1 <- data.frame(id=c('Generator', 'DC', 'Inverter', 'AC'),
                                  values=c(gen, DC, inv, AC))
            result <- rbind(result0, result1)
            result
          }
          )

###compareLosses

## compareLosses,ProdGCPV: no visible binding for global variable ‘name’
if(getRversion() >= "2.15.1") globalVariables(c('name'))

setGeneric('compareLosses', signature='...', function(...){standardGeneric('compareLosses')})

setMethod('compareLosses', 'ProdGCPV',
          definition=function(...){
            dots <- list(...)
            nms0 <- substitute(list(...))
            if (!is.null(names(nms0))){ ##estamos dentro de do.call
              nms <- names(nms0[-1])
            } else {
              nms <- as.character(nms0[-1])
            }
            foo <- function(object, label){
              yY <- losses(object)
              yY <- cbind(yY, name=label)
              yY
            }
            cdata <- mapply(FUN=foo, dots, nms, SIMPLIFY=FALSE)
            z <- do.call(rbind, cdata)
            z$id <- ordered(z$id, levels=c('Shadows', 'AoI', 'Generator', 'DC', 'Inverter', 'AC'))
            p <- dotplot(id~values*100, groups=name, data=z,
                         par.settings=solaR.theme, type='b',
                         auto.key=list(corner=c(0.95,0.2), cex=0.7), xlab='Losses (%)')
            print(p)
            return(z)
          }
          )
