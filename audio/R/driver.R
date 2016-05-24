load.audio.driver <- function(path) .Call("audio_load_driver", as.character(path), PACKAGE="audio")

audio.drivers <- function() .Call("audio_drivers_list", PACKAGE="audio")

set.audio.driver <- function(name) .Call("audio_use_driver", name, PACKAGE="audio")

current.audio.driver <- function() .Call("audio_current_driver", PACKAGE="audio")
