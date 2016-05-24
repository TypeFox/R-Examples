#
#	Rparallel.R
#Fri Jun 15 12:29:14 CEST 2012
#source('Rparallel.back.R');
library('tools');

#
#	<p> Lapply state reference classes
#
LapplyStateClass = setRefClass('LapplyState',
	fields = list( sequence = 'numeric', depth = 'numeric',
	runMode = 'logical', probeMode = 'logical', max_depth = 'numeric'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		sequence <<- 0;
		depth <<- 0;
		runMode <<- F;
		probeMode <<- F;
		.self$initFields(...);
		.self
	},
	# depth is defined by the number of times Lapply is recursively called
	depthInc = function() { depth <<- depth + 1; },
	depthDec = function() { depth <<- depth - 1; },
	# sequence is defined by the number of Lapplys that were started no matter how deeply nested
	sequenceInc = function() { sequence <<- sequence + 1; },
	sequenceDec = function() { sequence <<- sequence - 1; },
	isEqualTo = function(s) { depth == s$depth && sequence == s$sequence }
	#
	#	</p> methods
	#
	)
);
LapplyStateClass$accessors(names(LapplyStateClass$fields()));

LapplyProbeStateClass = setRefClass('LapplyProbeState',
	fields = list( elements = 'list' ),
	contains = 'LapplyState',
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		# initialize super class
		callSuper(...);
		# defaults
		elements <<- list();
		# overwrites
		.self$setProbeMode(T);
		.self
	},
	pushElements = function(es) {
		elements[[depth]] <<- if (length(elements) < depth) es else c(elements[[depth]], es);
		NULL
	},
	elementsCount = function(atDepth) {
		count = sum(if (atDepth > length(elements)) elements[[length(elements)]] else elements[[atDepth]]);
		count
	}
	#
	#	</p> methods
	#
	)
);
LapplyProbeStateClass$accessors(names(LapplyProbeStateClass$fields()));

LapplyRunStateClass = setRefClass('LapplyRunState',
	fields = list( chunkSize = 'numeric' ),
	contains = 'LapplyState',
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		# initialize super class
		callSuper(...);
		# overwrites
		.self$setRunMode(T);
		.self
	}
	#
	#	</p> methods
	#
	)
);
LapplyRunStateClass$accessors(names(LapplyRunStateClass$fields()));

#
#	<p> freezer classes
#

# copy functions code adapted from restorepoint R package
object.copy = function(obj) {
	# Dealing with missing values
	if (is.name(obj)) return(obj);
	obj_class = class(obj);

	copy =
		if ('environment' %in% obj_class) environment.copy(obj) else
		if (is.list(obj) && !(is.data.frame(obj))) list.copy(obj) else
		obj;
	return(copy)
}
list.copy = function(l)lapply(l, object.copy);
environment.copy = function(envir__)as.environment(eapply(envir__, object.copy));

# Freezer class stores individual calls for list elements iterated overwrites
# Also the structure of Lapply calls is stored to be able to re-associate results
#	with Lapply calls
# The isolation of individual calls allows for re-shuffeling of bundling calls for final
#	execution
LapplyFreezerClass = setRefClass('LapplyFreezer',
	fields = list( slots = 'list', calls = 'list', results = 'list' ),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		slots <<- list();
		calls <<- list();
		.self$initFields(...);
		.self
	},
	# depth is defined by the number of times Lapply is recursively called
	clear = function() {
		slots <<- list();
		calls <<- list();
		gc();
	},
	push = function(sequence, f, l, args, envir__ = parent.frame()) {
		# store by seqeunce id from LapplyState object
		Log(sprintf('Freezing %d invocations @ seq %d.', length(l), sequence), 5);
		slots[[as.character(sequence)]] <<- list(
			# definition of function called
			f = f,
			# number of list elements iterated
			N = length(l),
			# start index of result list in sequential order of calls
			start = sum(list.key(slots, 'N')) + 1
		);
		calls <<- c(calls, lapply(l, function(e) {
			callWithFunctionArgs(f, c(list(e), args), envir__ = envir__)
		}));
		NULL
	},
	N = function()length(calls),
	call = function(i)calls[[i]],

	# <p> results
	# for efficiency use nested structure
	pushResults = function(r){
		results[[length(results) + 1]] <<- r;
	},
	# in case results were stored in chunks (as lists) this method flattens the structure
	unlistResults = function() {
		results <<- unlist(results, recursive = F);
		NULL
	},
	finalizeResults = function() { NULL },
	# only to be called after finalizeResults
	resultsForSequence = function(s) {
		slot = slots[[as.character(s)]];
		results[slot$start : (slot$start + slot$N - 1)];
	}

	#
	#	</p> methods
	#
	)
);
LapplyFreezerClass$accessors(names(LapplyFreezerClass$fields()));

LapplyPersistentFreezerClass = setRefClass('LapplyPersistentFreezer',
	contains = 'LapplyFreezer',
	fields = list(),
	methods = list(
	finalizeResults = function() {
		callSuper();
	},
	resultsForSequence = function(s) {
		slot = slots[[as.character(s)]];
		seqCalls = slot$start : (slot$start + slot$N - 1);
		r = lapply(results[[length(results)]], function(r) {
			i = intersect(r$from:r$to, seqCalls);
			r1 = if (length(i) > 0) {
				r0 = frozenCallResults(r$file);
				r0[i - r$from + 1];
			} else NULL;
			r1
		});
		r = unlist.n(r, 1);
		r
	}

	)
);

# The sentinel stack records entry points into parallelization for the different sequential rampUps
#	occuring during parallelization the ramp down of a given Lapply is equivalent to the rampUp of
#	the ensueing parallelization
# As probing occurs for successively deeper levels, the first sequence number for a given level is stored
#	write-once and will represent the first lapply-loop to parallelize
#	As the depth of the rampDown is unclear, sequenceStop will be updated to represent the most current
#	sequence number of ongoing Lapply loops
# Recovering will than happen between sequence and sequenceStop at the recorded depth
LapplyExecutionStateClass = setRefClass('LapplyExecutionState',
	fields = list(
		# execution state
		sequenceNos = 'list', rampUp = 'numeric',
		# results
		freezerClass = 'character', freezers = 'list'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(freezerClass = 'LapplyFreezer', ...) {
		# initialize super class
		callSuper(freezerClass = freezerClass, ...);
		# defaults
		sequenceNos <<- list();
		rampUp <<- 1;
		# overwrites
		.self
	},
	# add a new sentinel in the parallelization stack
	addSentinel = function() {
		sequenceNos[[length(sequenceNos) + 1]] <<- list(depth = -1, sequence = -1, sequenceStop = -1);
	},
	# update only last element in stack (previous sentinels are fixed)
	# only record first sequence no for given rampUp and depth
	pushSequenceForRampUp = function(sequenceNo, depth) {
		N = length(sequenceNos);
		# if we probe for bigger depthes we re-record the starting sequence
		#	as parallelization will happen at that deeper level
		if (sequenceNos[[N]]$depth < depth) {
			sequenceNos[[N]] <<- list(depth = depth, sequence = sequenceNo);
		}
		# a new sequence number at the current maximal depth is recored as a potential stop of
		#	the current parallelization
		if (sequenceNos[[N]]$depth <= depth) {
#			sequenceNos[[rampUp + 1]] <<- merge.lists(sequenceNos[[rampUp + 1]],
#				list(sequenceStop = sequenceNo));
			sequenceNos[[N]]$sequenceStop <<- sequenceNo;
			Log(sprintf('new sentinel stop: %d', sequenceNo), 6);
		}
		NULL
	},
	# the currentSentinel is the one to skip, which comes from the previous cursor position
	currentSentinel = function() {
		sequenceNos[[rampUp]]
#		if (rampUp == 1 || rampUp - 1 > length(sequenceNos))
#			list(depth = -1, sequence = -1, sequenceStop = -1) else
#			sequenceNos[[rampUp - 1]]
	},
	incCursor = function() {
		rampUp <<- rampUp + 1;
		.self$currentSentinel()
	},
	resetCursor = function() {
		rampUp <<- 1;
	},

	#
	#	functional methods
	#

	rampUpForeFront = function()length(sequenceNos),
	# detect range where results need to be recovered (thawed from the freezer)
	# the latest state (stack position N) is nascent and ignored
	#	there is currently probing or parallelization going on
	checkAgainstState = function(state) {
		#N = length(sequenceNos);
		N = .self$rampUpForeFront();
		sentinel = .self$currentSentinel();
		r = (N > rampUp &&
			state$sequence >= sentinel$sequence &&
			state$sequence <= sentinel$sequenceStop &&
			state$depth == sentinel$depth);
		#Log(sprintf('sentinel check %d [rampup %d]', r,  Lapply_executionState__$rampUp));
		r
	},
	skipToRampDown = function(state) {
		N = length(sequenceNos);
		sentinel = .self$currentSentinel();
#if (rampUp > 0 && state$depth > 1) browser();
		r = (N > rampUp && state$sequence <= sentinel$sequenceStop);
		#Log(sprintf('sentinel check %d [rampup %d]', r,  Lapply_executionState__$rampUp));
		r
	},
	# <N> after probing length of sequenceNos is increased by one
	#	when reaching this point we want to parallelize
	isLastRampUp = function() {
		rampUp == length(sequenceNos)
	},
	# <N> tb called after the state has been processed
	#	if the last sequence was processed the cursor is advanded to proceed to the next rampUp
	adjustCursor = function(state) {
		sentinel = .self$currentSentinel();
		if (.self$checkAgainstState(state) && state$sequence == sentinel$sequenceStop)
			.self$incCursor();
	},
	#
	#	freezer methods
	#
	currentFreezer = function() {
		if (rampUp > length(freezers)) {
			freezers[[rampUp]] <<- getRefClass(freezerClass)$new();
		}
		freezers[[rampUp]]
	}

	#
	#	</p> methods
	#
	)
);
LapplyExecutionStateClass$accessors(names(LapplyExecutionStateClass$fields()));


#
#	<p> core parallize functions
#

if (!exists('parallelize_env')) parallelize_env <- new.env();

# force_rerun instructs backends to ignore state-retaining files and re-run all computations
parallelize_initialize = Lapply_initialize = function(Lapply_config = Lapply_config_default,
	stateClass = 'LapplyState', backend = 'local', freezerClass = 'LapplyFreezer', ...,
	force_rerun = FALSE, sourceFiles = NULL, parallel_count = NULL) {
	# <p> check for turning off
	if (backend == 'off') {
		parallelize_setEnable(F);
		return(NULL);
	} else parallelize_setEnable(T);
	# <p> misc setup
	Log.setLevel(firstDef(Lapply_config$logLevel, Log.level(), 4));
	parallelize_setEnable(T);
	# <p> config
	sourceFiles = c(Lapply_config$sourceFiles, Lapply_config$backends[[backend]]$sourceFiles, sourceFiles);
	backendConfig = merge.lists(
		Lapply_backendConfig_default,
		Lapply_config$backends[[backend]],
		list(force_rerun = force_rerun),
		list(sourceFiles = sourceFiles)
	);
	Lapply_config = merge.lists(
		Lapply_config_default,
		Lapply_config,
		list(backend = backend, backendConfig = backendConfig, parallel_count = parallel_count)
	);
	Lapply_setConfig(Lapply_config);
	# <p> backend
	if (exists('Lapply_backend__',  envir = parallelize_env)) rm('Lapply_backend__', envir = parallelize_env);
	freezerClass = firstDef(backendConfig$freezerClass, freezerClass);
	# <p> iteration states
	state = getRefClass(stateClass)$new(...);
	assign('Lapply__', state, envir = parallelize_env);
	assign('Lapply_executionState__', LapplyExecutionStateClass$new(freezerClass = freezerClass),
		envir = parallelize_env);
	NULL
}
parallelize_initializeBackendWithCall = function(call_, Lapply_config) with(Lapply_config, {
	# heuristic to get original function name if not supplied
	#functionName = firstDef(call_$name, deparse(sys.call(-6)[[2]][[1]]));
	functionName = firstDef(call_$name, deparse(sys.call(-6)[[2]]));
	#signature = md5sumString(sprintf('%s%s', tag, deparse(.f)));
	signature = md5sumString(sprintf('%s%s', functionName, backend));
	Log(sprintf('parallelize signature %s', signature), 5);

	backendClass = sprintf('ParallelizeBackend%s', uc.first(firstDef(backendConfig$backend, backend)));
	#backendConfig = rget(sprintf('%s_config__', backendClass), default = list());
	assign('Lapply_backend__', new(backendClass, config = backendConfig, signature = signature),
		envir = parallelize_env);
})

Lapply_initializeState = function(stateClass = 'LapplyState', ...) {
	state = getRefClass(stateClass)$new(...);
	assign('Lapply__', state, envir = parallelize_env);
}
Lapply_initialze_probing = function() {
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply_executionState__$resetCursor();
}

Lapply_setConfig = function(config) {
	assign('Lapply_globalConfig__', config, envir = parallelize_env);
}
Lapply_getConfig = function() {
	get('Lapply_globalConfig__', envir = parallelize_env);
	#Lapply_globalConfig__
}

#
#	</p> Lapply state reference classes
#


#
#	<p> S3 classes
#

Lapply_error = function() {
	Lapply__ = get('Lapply__', envir = parallelize_env);
	e = structure(list(state = Lapply__$copy()), class =  c('Lapply_error', 'simpleError'));
	e
}
Lapply_error.as.character = function(e) {
	msg = structure(sprintf('Lapply stopped at sequence %d, depth %d', e$state$sequence, e$state$depth),
		error = e);
	msg
}

#
#	</p> S3 classes
#

#
#	<p> error handling
#

Throw = function(msg, state) {
	assign('Global_error_state__', state, envir = parallelize_env);
	stop(msg);
}

Catch = function(result, errorClass, handler) {
	doCallHandler = class(result) == 'try-error' &&
		any(class(get('Global_error_state__', envir = parallelize_env)) == errorClass);
	r = if (doCallHandler) {
		errorState = get('Global_error_state__', envir = parallelize_env);
		handler(errorState);
	} else NULL;
	r = list(result = r, didCall = doCallHandler);
	r
}

Try = function(expr, catch = list(), silent = T, setClass = F) {
	r = try(expr, silent = silent);
	didCall = F;
	if (exists('Global_error_state__', envir = parallelize_env)) {
		for (i in 1:length(catch)) {
			errorClass = names(catch)[i];
			r0 = Catch(r, errorClass, catch[[i]]);
			if (r0$didCall) {
				r = r0$result;
				if (setClass) {
					if (is.null(r)) r = integer(0);
					class(r) = errorClass;
				}
			}
			didCall = didCall || r0$didCall;
		}
		remove('Global_error_state__', envir = parallelize_env);
	}
	if (!didCall && class(r) == 'try-error') stop(r[1]);
	r
}

#
#	</p> error handling
#

Lapply_config_default = list(
	max_depth = 2, parallel_count = 32, parallel_stack = 10,
	provideChunkArgument = F, offline = F, stateDir = '.',
	wait_interval = 30
);
Lapply_backendConfig_default = list(doNotReschedule = F, doSaveResult = F);

Lapply_do = function(l, .f, ..., Lapply_config, Lapply_chunk = 1, envir__) {
	f_ = function(e, ...)do.call(.f, list(e, ...), envir = envir__);
	r = if (Lapply_config$provideChunkArgument)
		lapply(l, f_, Lapply_chunk = Lapply_chunk, ...) else
		lapply(l, f_, ...);
	r
}

Lapply_probeDepth = function(l, .f, ..., Lapply_config, Lapply_chunk = 1, envir__) {
	# <p> probe deeper levels if Lapply_depth not yet reached
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	Lapply__$pushElements(length(l));
	Log(sprintf('Lapply: Adding %d elements @depth %d.', length(l), Lapply__$depth), 5);
#if (Lapply_executionState__$rampUp == 2 && Lapply__$depth == 2) browser();
	r = if (Lapply__$max_depth > Lapply__$depth) {
		probeWrapper = function(e, ...) {
			Try(
				.f(e, ...), catch = list(Lapply_error = function(e) {
					Log(sprintf('caught escape from depth %d', e$depth));
					# <N> the escape left the Lapply state stale
					Lapply__$depthDec();
					e
				})
			)
		};
		Lapply_do(l, probeWrapper, ...,
			Lapply_config = Lapply_config, Lapply_chunk = Lapply_chunk, envir__ = envir__);
	} else list(Lapply_error());
	# pushSequenceForRampUp records depending of its current state
	#	see implementation thereof
	Lapply_executionState__$pushSequenceForRampUp(Lapply__$sequence, Lapply__$depth);
	# generate escape from ramp-down
	#Throw('lapply probe escape', Lapply_error());
	if (any(as.vector(sapply(r, class)) == 'Lapply_error')) Throw('lapply probe escape', Lapply_error());
	# reached when no parallelization happened
	r
}

# from determines the starting point of probing
Lapply_probe = function(call_, Lapply_config) with(Lapply_config,  {
 	depths = 1:max_depth;
	Log(sprintf("Lapply_probe: depths to probe: %s", join(depths)), 5);
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	# probing will determine a new sentinel
	Lapply_executionState__$addSentinel();
	r = NULL;
	for (i in depths) {
		# reset state: first rampUp, cursor to beginning
		Lapply_executionState__$resetCursor();
		Lapply_initializeState('LapplyProbeState', max_depth = i);
		Lapply__ = get('Lapply__', envir = parallelize_env);
		# probing function
		Log(sprintf('Lapply: probing depth %d.', i), 5);
		r = Try(Do.call(call_$fct, call_$args, envir = call_$envir, envirArgs = call_$envirArgs),
			catch = list(Lapply_error = function(e) {
				Log('final catch', 6);
				e
			}));
		# no parallelization found, real result already computed
		if (all(class(r) != 'Lapply_error')) break;
		# compute registered parallelization 
		count = Lapply__$elementsCount(i);
		Log(sprintf('Lapply: registered %d parallel jobs @depth %d.', count, i), 5);
		# <p> determine stopping condition
		rampUp = Lapply_executionState__$getRampUp();
		# specific counts per rampUp
		this_parallel_count =
			if (rampUp > length(parallel_count)) parallel_count[1] else parallel_count[rampUp];
		if (count >= this_parallel_count) break;
		# break if no parallelism was found
	}
	# in case of parallelization Lapply_errors were thrown. these are re-thrown at higher levels such
	#	that return values are nested lists of Lapply_errors. We detect this case by probing for a return
	#	list and checking classes of members
	if (any(as.vector(sapply(r, class)) == 'Lapply_error')) {
		r = Lapply_error();
	}
	r
})

Lapply_parallelize = function(l, .f, ..., Lapply_config, envir__) {
	Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	r = if (Lapply__$depth == Lapply__$max_depth) {
		lapply_dispatch(Lapply_backend__, l, .f, ..., envir__ = envir__);
		# <i><N> raise exception only if asynchroneous
		Throw('lapply run escape', Lapply_error());
	} else {
		Log(sprintf('entering parallelization at depth %d', Lapply__$depth), 6);
		Lapply_do(l, function(e)Try(.f(e, ...),
			catch = list(Lapply_error = function(e){
				Log('Lapply_do catch', 6);
				Lapply__$depthDec();	# <A> balance depth
				e
			})),
		..., Lapply_config = Lapply_config, envir__ = envir__);
	}
	r
}

# excute code for the given rampUp
#	Lapply_depth: depth at which to parallelize
Lapply_run = function(call_, Lapply_depth, Lapply_config) {
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply_executionState__$resetCursor();
	# reset state cursor
	Lapply_initializeState('LapplyRunState', max_depth = Lapply_depth);

	r = Try(Do.call(call_$fct, call_$args, envir = call_$envir, envirArgs = call_$envirArgs),
		catch = list(Lapply_error = function(e)Log('final run catch', 5)));
	Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
	lapply_dispatchFinalize(Lapply_backend__);
	r
}

Lapply_recoverState = function(sequence) {
	Lapply__ = get('Lapply__', envir = parallelize_env);
	Log(sprintf('Recovering state for sequence %d, depth %d.', Lapply__$sequence, Lapply__$depth), 5);
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();
	r = freezer$resultsForSequence(sequence);
	r
}

.lapply = Lapply = Lapply_backup = function(l, .f, ...,
	Lapply_config = Lapply_getConfig(),
	#Lapply_local = rget('Lapply_local', envir = parallelize_env, default = T),
	Lapply_local = Lapply_config$local,
	Lapply_chunk = 1) {

	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	# <p> Lapply__ state
	Lapply__$depthInc();
	envir__ = parent.frame();
	#	starting new Lapply
	Lapply__$sequenceInc();
	# <p> handle cases
	Log(sprintf("Sequence %d.", Lapply__$sequence), 6);
	r = if (Lapply_executionState__$checkAgainstState(Lapply__)) {
		r = Lapply_recoverState(Lapply__$sequence);
		Lapply_executionState__$adjustCursor(Lapply__);
		r
	#	<p> Lapply within range of sentinel but not at max_depth
	} else if (Lapply_executionState__$skipToRampDown(Lapply__)) {
		Log(sprintf("Skipping sequence %d.", Lapply__$sequence), 4);
		# <p> either do regular lapply or skip parallelisation to required rampUp
		Lapply_do(l, .f, ..., Lapply_config = Lapply_config, Lapply_chunk = Lapply_chunk, envir__ = envir__);
	#	<p> probe for degree of parallelism
	} else if (Lapply__$probeMode) {
		Lapply_probeDepth(l, .f, ...,
			Lapply_config = Lapply_config, Lapply_chunk = Lapply_chunk, envir__ = envir__);
	#	<p> parallelization
	} else if (Lapply__$runMode) {
		Lapply_parallelize(l, .f, ..., Lapply_config = Lapply_config, envir__ = envir__);
	#	<p> local mode
	} else {
		Lapply_do(l, .f, ..., Lapply_config = Lapply_config, Lapply_chunk = Lapply_chunk, envir__ = envir__);
	};
	# <p> Lapply__ state
	Lapply__$depthDec();

	r
}

Sapply = sapply;
Sapply_backup = function(X, FUN, ...) {
	r = Lapply(X, FUN, ...);
	r0 = sapply(r, function(e)e);
}

Apply_margin_error = 'wrong MARGIN argument supplied to Apply';
Apply = apply;
Apply_backup = function(X, MARGIN, FUN, ...) {
	r = if (length(MARGIN) == 1) {
		extractor = if (MARGIN == 1) function(X, i)X[i, ] else
			if (MARGIN == 2) function(X, i)X[, i] else stop(Apply_margin_error);
		r0 = Lapply(1:dim(X)[MARGIN], function(i, ..., Apply_object__, Apply_FUN__, Apply_extractor__) {
			Apply_FUN__(Apply_extractor__(Apply_object__, i), ...)
		} , ..., Apply_object__ = X, Apply_FUN__ = FUN, Apply_extractor__ = extractor);
		r = sapply(r0, function(e)e);
		r
	} else if (length(MARGIN) == 2 && all(MARGIN == 1:2)) {
		extractor = function(X, tuple)X[tuple[1], tuple[2]];
		els = apply(merge(data.frame(row = 1:dim(X)[1]), data.frame(col = 1:dim(X)[2])), 1, as.list);
		r0 = Lapply(els, function(i, ..., Apply_object__, Apply_FUN__, Apply_extractor__) {
			Apply_FUN__(Apply_extractor__(Apply_object__, unlist(i)), ...)
		} , ..., Apply_object__ = X, Apply_FUN__ = FUN, Apply_extractor__ = extractor);
		r = sapply(r0, function(e)e);
		if (is.vector(r)) r = matrix(r, ncol = dim(X)[1]);
		r
	} else {
		stop(Apply_margin_error);
	}
	r
}

parallelizeStep = function(call_, Lapply_config) {
	# probe parallelism
	r = Lapply_probe(call_, Lapply_config = Lapply_config);
	# no parallelization was possible
	if (all(class(r) != 'Lapply_error')) return(r);
	# run computation for this rampUp sequence
	Lapply__ = get('Lapply__', envir = parallelize_env);
	Lapply_run(call_, Lapply_depth = Lapply__$max_depth, Lapply_config = Lapply_config);
	Lapply_error();
}

# tag: allow to uniquify the signature for multiple calls to the same function
parallelizeOfflineStep = function(call_, Lapply_config) with(Lapply_config, {
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
	Log(sprintf('parallelize rampUp %d', Lapply_executionState__$rampUpForeFront()), 5);
	#statePath = sprintf('%s/.parallelize%s.RData', stateDir, signature);
	if (Lapply_executionState__$rampUpForeFront() == 0) {
		initScheduling(Lapply_backend__, call_);
	} else {
		restoreParallelizationState(Lapply_backend__);
	}
	#r = parallelizeStep(.f, ..., Lapply_config = Lapply_config);
	r = performParallelizationStep(Lapply_backend__, call_, Lapply_config = Lapply_config);
	saveParallelizationState(Lapply_backend__);
	if (any(class(r) == 'Lapply_error')) {
		if (!backendConfig$doNotReschedule)
			scheduleNextParallelization(Lapply_backend__, call_);
	} else {
		r = finalizeParallelization(Lapply_backend__, r);
	}
	r
})

parallelize_dummy = function(.f, ..., Lapply_config = NULL).f(...);
parallelize_call_dummy = function(.call, Lapply_config = NULL)eval(.call);

parallelize_internal = function(call_, Lapply_local = rget('Lapply_local', default = F),
	parallelize_wait = T) {
	Lapply_config = merge.lists(Lapply_getConfig(), list(local = Lapply_local));
	r = if (Lapply_local) {
		# <!><i> setEnable(F), re-enable afterwards
		do.call(call_$fct, call_$args, envir = call_$envir);
	} else {
		Lapply_setConfig(Lapply_config);
		parallelize_initializeBackendWithCall(call_, Lapply_config = Lapply_config);
		Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
		r = parallelize_backend(Lapply_backend__, call_);
		# this is a delegating backend (only supported in offline mode)
		if (parallelize_wait && Lapply_backend__@offline) {
			while (pollParallelization(Lapply_backend__)$continue) Sys.sleep(Lapply_config$wait_interval);
			r = getResult(Lapply_backend__);
		}
		r
	}
	r
}

parallelize = parallelize_backup = function(.f, ..., Lapply_local = rget('Lapply_local', default = FALSE),
	parallelize_wait = TRUE) {
	call_ = list(fct = .f, args = list(...), envir = parent.frame(), name = as.character(sys.call()[[2]]));
	parallelize_internal(call_, Lapply_local = Lapply_local, parallelize_wait = parallelize_wait);
}


parallelize_call = parallelize_call_backup = function(.call, ..., parallelize_wait = TRUE) {
	call_  = encapsulateCall(sys.call()[[2]], ..., envir__ = parent.frame());
	parallelize_internal(call_, ..., parallelize_wait = parallelize_wait);
}

#
#	<p> utility functions
#

tempcodefile = function(fcts) {
	fctNames = as.character(as.list(sys.call()[[2]])[-1]);
	# create source file with code
	code = join(sapply(fctNames,
		function(name) {
			def = join(deparse(get(name)), sep = "\n");
			code = sprintf('%s = %s', name, def);
			code
		}), sep = "\n");
	codeFile = tempfile();
	# windows specific code
	codeFile = gsub('([\\])', '/', codeFile);
	writeFile(codeFile, code);
	codeFile
}

