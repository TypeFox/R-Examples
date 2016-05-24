# Created: 2015.05.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

library(shiny)

ui <- shinyUI(fluidPage(
  titlePanel("Simulations"),
  sidebarLayout(
    sidebarPanel(
      h3("parameters"),
      selectInput("distro", "Distribution:", 
									choices=c("dnbeta","dneta","dnf","dnt","kprime","lambdap",
														"prodchisqpow","proddnf","sumchisqpow","sumlogchisq","upsilon"),
									selected="upsilon",
									multiple=FALSE),
			conditionalPanel(
				condition = "input.distro == 'dnbeta'",
				numericInput("dnbeta_df1", "df1:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnbeta_df2", "df2:", min=1, max=Inf, value=100, step=0.1),
				numericInput("dnbeta_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dnbeta_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dneta'",
				numericInput("dneta_df", "df:", min=-Inf, max=Inf, value=50, step=0.1),
				numericInput("dneta_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dneta_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dnf'",
				numericInput("dnf_df1", "df1:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnf_df2", "df2:", min=1, max=Inf, value=100, step=0.1),
				numericInput("dnf_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dnf_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dnt'",
				numericInput("dnt_df", "df:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnt_ncp1", "ncp1:", min=-Inf, max=Inf, value=1, step=0.001),
				numericInput("dnt_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'kprime'",
				numericInput("kprime_a", "a:", min=0, max=Inf, value=5, step=0.1),
				numericInput("kprime_b", "b:", min=1, max=Inf, value=1, step=0.1),
				numericInput("kprime_v1", "v1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("kprime_v2", "v2:", min=5, max=Inf, value=20, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'lambdap'",
				numericInput("lambdap_df", "df:", min=3, max=Inf, value=10, step=0.1),
				numericInput("lambdap_t", "t:", min=-Inf, max=Inf, value=0, step=0.1)
			),
			conditionalPanel(
				condition = "input.distro == 'prodchisqpow'",
				helpText('Enter parameters separated by commas.',
								'Parameters are recycled against each other'),
				textInput("prodchisqpow_df", "df:", value="10,3,5"),
				textInput("prodchisqpow_ncp", "ncp:", value="0,2,1"),
				textInput("prodchisqpow_pow", "pow:", value="0.5,1,0.25")
			),
			conditionalPanel(
				condition = "input.distro == 'proddnf'",
				helpText('Enter parameters separated by commas.',
								'Parameters are recycled against each other'),
				textInput("proddnf_df1", "df1:", value="10,30,50,500"),
				textInput("proddnf_df2", "df2:", value="100,50,150,200"),
				textInput("proddnf_ncp1", "ncp1:", value="0,2,0,5"),
				textInput("proddnf_ncp2", "ncp2:", value="1,1,2,30")
			),
			conditionalPanel(
				condition = "input.distro == 'sumchisqpow'",
				helpText('Enter parameters separated by commas.',
								'Parameters are recycled against each other'),
				textInput("sumchisqpow_wts", "wts:", value="10,-20,1,-3"),
				textInput("sumchisqpow_df", "df:", value="100,50,150,200"),
				textInput("sumchisqpow_ncp", "ncp:", value="5,5,2,1"),
				textInput("sumchisqpow_pow", "ncp2:", value="1")
			),
			conditionalPanel(
				condition = "input.distro == 'sumlogchisq'",
				helpText('Enter parameters separated by commas.',
								'Parameters are recycled against each other'),
				textInput("sumlogchisq_wts", "wts:", value="1,-1,1,-1,1,-1"),
				textInput("sumlogchisq_df", "df:", value="100,200"),
				textInput("sumlogchisq_ncp", "ncp:", value="4,4,3,3,2,2")
			),
			conditionalPanel(
				condition = "input.distro == 'upsilon'",
				helpText('Enter parameters separated by commas.',
								'Parameters are recycled against each other'),
				textInput("upsilon_df", "df:", value="100,100,200,200,300,300"),
				textInput("upsilon_t", "t:", value="4,-4,10,-5,2,-3")
			),
      hr(),
      numericInput("nsamples", "Number of draws:", min = 200, max = 50000, value = 10000, step=100),
      numericInput("randseed", "Rand seed:", min = 1, max = .Machine$integer.max, value = 2015, step=1)
    ,width=3),
    mainPanel(
			tabsetPanel(
				tabPanel("d-d",plotOutput("ddplot")),
				tabPanel("q-q",plotOutput("qqplot")),
				tabPanel("p-p",plotOutput("ppplot"))
			)
		))
,title="Testing sadist distributions"))

shinyUI(ui)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
