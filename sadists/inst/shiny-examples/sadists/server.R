# Created: 2015.05.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

library(shiny)
library(ggplot2)
library(reshape2)
library(sadists)

text.size <<- 12   # sigh

# convert a string like "234, 12, 99" to an array of numerics c(234, 12, 99)
chars_to_num <- function(astr) {
	unlist(lapply(strsplit(astr,'\\s+|,'),as.numeric))
}

server <- function(input, output) {

	get_dpqr <- reactive({
		dpqr <- switch(input$distro,
						dnbeta = {
							argslist = list(df1=input$dnbeta_df1,df2=input$dnbeta_df2,
															ncp1=input$dnbeta_ncp1,ncp2=input$dnbeta_ncp2)
							list(d=function(x) do.call(ddnbeta,c(argslist,list(x=x))),
									 p=function(q) do.call(pdnbeta,c(argslist,list(q=q))),
									 q=function(p) do.call(qdnbeta,c(argslist,list(p=p))),
									 r=function(n) do.call(rdnbeta,c(argslist,list(n=n))))
						},
						dneta = {
							argslist = list(df=input$dneta_df,
															ncp1=input$dneta_ncp1,ncp2=input$dneta_ncp2)
							list(d=function(x) do.call(ddneta,c(argslist,list(x=x))),
									 p=function(q) do.call(pdneta,c(argslist,list(q=q))),
									 q=function(p) do.call(qdneta,c(argslist,list(p=p))),
									 r=function(n) do.call(rdneta,c(argslist,list(n=n))))
						},
						dnf = {
							argslist = list(df1=input$dnf_df1,df2=input$dnf_df2,
															ncp1=input$dnf_ncp1,ncp2=input$dnf_ncp2)
							list(d=function(x) do.call(ddnf,c(argslist,list(x=x))),
									 p=function(q) do.call(pdnf,c(argslist,list(q=q))),
									 q=function(p) do.call(qdnf,c(argslist,list(p=p))),
									 r=function(n) do.call(rdnf,c(argslist,list(n=n))))
						},
						dnt = {
							argslist = list(df=input$dnt_df,
															ncp1=input$dnt_ncp1,ncp2=input$dnt_ncp2)
							list(d=function(x) do.call(ddnt,c(argslist,list(x=x))),
									 p=function(q) do.call(pdnt,c(argslist,list(q=q))),
									 q=function(p) do.call(qdnt,c(argslist,list(p=p))),
									 r=function(n) do.call(rdnt,c(argslist,list(n=n))))
						},
						kprime = {
							argslist = list(a=input$kprime_a,b=input$kprime_b,
															v1=input$kprime_v1,v2=input$kprime_v2)
							list(d=function(x) do.call(dkprime,c(argslist,list(x=x))),
									 p=function(q) do.call(pkprime,c(argslist,list(q=q))),
									 q=function(p) do.call(qkprime,c(argslist,list(p=p))),
									 r=function(n) do.call(rkprime,c(argslist,list(n=n))))
						},
						lambdap = {
							argslist = list(df=input$lambdap_df,t=input$lambdap_t)
							list(d=function(x) do.call(dlambdap,c(argslist,list(x=x))),
									 p=function(q) do.call(plambdap,c(argslist,list(q=q))),
									 q=function(p) do.call(qlambdap,c(argslist,list(p=p))),
									 r=function(n) do.call(rlambdap,c(argslist,list(n=n))))
						},
						prodchisqpow = {
							argslist = list(df=chars_to_num(input$prodchisqpow_df),
															ncp=pmax(0,chars_to_num(input$prodchisqpow_ncp)),
															pow=chars_to_num(input$prodchisqpow_pow))
							list(d=function(x) do.call(dprodchisqpow,c(argslist,list(x=x))),
									 p=function(q) do.call(pprodchisqpow,c(argslist,list(q=q))),
									 q=function(p) do.call(qprodchisqpow,c(argslist,list(p=p))),
									 r=function(n) do.call(rprodchisqpow,c(argslist,list(n=n))))
						},
						proddnf = {
							argslist = list(df1=chars_to_num(input$proddnf_df1),
															df2=chars_to_num(input$proddnf_df2),
															ncp1=pmax(0,chars_to_num(input$proddnf_ncp1)),
															ncp2=pmax(0,chars_to_num(input$proddnf_ncp2)))
							list(d=function(x) do.call(dproddnf,c(argslist,list(x=x))),
									 p=function(q) do.call(pproddnf,c(argslist,list(q=q))),
									 q=function(p) do.call(qproddnf,c(argslist,list(p=p))),
									 r=function(n) do.call(rproddnf,c(argslist,list(n=n))))
						},
						sumchisqpow = {
							argslist = list(wts=chars_to_num(input$sumchisqpow_wts),
															df=pmax(3,chars_to_num(input$sumchisqpow_df)),
															ncp=pmax(0,chars_to_num(input$sumchisqpow_ncp)),
															pow=chars_to_num(input$sumchisqpow_pow))
							list(d=function(x) do.call(dsumchisqpow,c(argslist,list(x=x))),
									 p=function(q) do.call(psumchisqpow,c(argslist,list(q=q))),
									 q=function(p) do.call(qsumchisqpow,c(argslist,list(p=p))),
									 r=function(n) do.call(rsumchisqpow,c(argslist,list(n=n))))
						},
						sumlogchisq = {
							argslist = list(wts=chars_to_num(input$sumlogchisq_wts),
															df=pmax(3,chars_to_num(input$sumlogchisq_df)),
															ncp=pmax(0,chars_to_num(input$sumlogchisq_ncp)))
							list(d=function(x) do.call(dsumlogchisq,c(argslist,list(x=x))),
									 p=function(q) do.call(psumlogchisq,c(argslist,list(q=q))),
									 q=function(p) do.call(qsumlogchisq,c(argslist,list(p=p))),
									 r=function(n) do.call(rsumlogchisq,c(argslist,list(n=n))))
						},
						upsilon = {
							argslist = list(df=pmax(3,chars_to_num(input$upsilon_df)),
															t=chars_to_num(input$upsilon_t))
							list(d=function(x) do.call(dupsilon,c(argslist,list(x=x))),
									 p=function(q) do.call(pupsilon,c(argslist,list(q=q))),
									 q=function(p) do.call(qupsilon,c(argslist,list(p=p))),
									 r=function(n) do.call(rupsilon,c(argslist,list(n=n))))
						}
						)
		dpqr
	})

  sims <- reactive({
		dpqr <- get_dpqr()
    set.seed(input$randseed)
		rv <- sort(dpqr$r(input$nsamples))
		data <- data.frame(draws=rv,pvals=dpqr$p(rv))
		data
  })

  # dd plot the results.
  output$ddplot <- renderPlot({
		dpqr <- get_dpqr()
    data <- sims()

		do.log <- (min(data$draws) > 0) && ((max(data$draws) / min(data$draws)) > 20)

		# http://stackoverflow.com/a/5688125/164611
		p1 <- qplot(data$draws, geom = 'blank') +   
			geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
			stat_function(fun = dpqr$d, aes(colour = 'Theoretical')) +                       
			geom_histogram(aes(y = ..density..), alpha = 0.3) +                        
			scale_colour_manual(name = 'Density', values = c('red', 'blue')) +
			theme(text=element_text(size=text.size)) + 
			labs(x=paste0('draws from ',input$distro,' distribution'),
					 title="Density (tests dfunc)")
		#if (do.log)
			#p1 <- p1 + scale_x_log10()
		return(p1)
  })

  # qq plot the results.
  output$qqplot <- renderPlot({
		dpqr <- get_dpqr()
    data <- sims()

		do.log <- (min(data$draws) > 0) && ((max(data$draws) / min(data$draws)) > 20)

		# Q-Q plot
		p2 <- ggplot(data, aes(sample = draws)) + stat_qq(distribution=function(p) { dpqr$q(p) }) +
			geom_abline(slope=1,intercept=0,colour='red') + 
			theme(text=element_text(size=text.size)) + 
			labs(title="Q-Q plot (tests qfunc)")
		if (do.log)
			p2 <- p2 + scale_x_log10() + scale_y_log10()
		return(p2)
	})

  # pp plot the results.
  output$ppplot <- renderPlot({
    data <- sims()

		# empirical CDF of the p-values; should be uniform
		p3 <- ggplot(data, aes(sample = pvals)) + stat_qq(distribution=qunif) + 
			geom_abline(slope=1,intercept=0,colour='red') + 
			theme(text=element_text(size=text.size)) + 
			labs(title="P-P plot (tests pfunc)")
		return(p3)
	})
}
shinyServer(server)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
