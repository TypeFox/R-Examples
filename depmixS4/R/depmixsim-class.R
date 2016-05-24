# Classes for simulated mix and depmix models

setClass("mix.sim",
	contains="mix",
	representation(
		states="matrix"
	)
)

setClass("depmix.sim",
	contains="depmix",
	representation(
		states="matrix"
	)
)

setAs("mix.fitted","mix.sim",def=function(from) {
  as(as(from,"mix"),"mix.sim")
})

setAs("depmix.fitted","depmix.sim",def=function(from) {
  as(as(from,"mix"),"mix.sim")
})