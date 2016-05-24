#
#	Rparallel.back.R
#Sun Jul 15 10:48:17 UTC 2012

#
#	<p> general documentation
#

# online vs offline mode: online means that rampUps are computed in one go, whereas offline backends compute one rampUp for a single invocation
# delegating backends are backends that forward execution to another offline backend
#	example: OGSremote -> OGS

#
#	<p> generic interface
#

setGeneric("isSynchroneous", function(self) standardGeneric("isSynchroneous"));
setGeneric("lapply_dispatch", function(self, l, f, ...) standardGeneric("lapply_dispatch"));
setGeneric("lapply_dispatchFinalize", function(self) standardGeneric("lapply_dispatchFinalize"));
setGeneric("lapply_results", function(self, r) standardGeneric("lapply_results"));
# parallelize function as customized by the backend
setGeneric('parallelize_backend', function(self, call_) standardGeneric('parallelize_backend'));
#	scheduling
setGeneric('initScheduling',
	function(self, call_) standardGeneric('initScheduling'));
setGeneric('performParallelizationStep',
	function(self, call_, Lapply_config) standardGeneric('performParallelizationStep'));
setGeneric('finalizeParallelization',
	function(self, r) standardGeneric('finalizeParallelization'));

setGeneric('saveParallelizationState',
	function(self) standardGeneric('saveParallelizationState'));
setGeneric('restoreParallelizationState',
	function(self) standardGeneric('restoreParallelizationState'));
setGeneric('scheduleNextParallelization',
	function(self, call_) standardGeneric('scheduleNextParallelization'));
setGeneric('pollParallelization',
	function(self, options) standardGeneric('pollParallelization'));
setGeneric('getResult',
	function(self) standardGeneric('getResult'));

#
#	<p> default class
#

setClass('ParallelizeBackend',
	representation = list(
		config = 'list', offline = 'logical', signature = 'character'
	),
	prototype = list(config = list(), offline = F, signature = '')
);
setMethod('initialize', 'ParallelizeBackend', function(.Object, config = list(), signature = '') {
	.Object@config = config;
	.Object@signature = signature;
	if (!is.null(config$offline)) .Object@offline = config$offline;
	.Object
});

#
#	<p> default class implementation
#

setMethod('isSynchroneous', 'ParallelizeBackend', function(self) { return(T); });
# use envir__ to evaluate ...
setMethod('lapply_dispatch', 'ParallelizeBackend', function(self, l, f, ..., envir__ = parent.frame()) { 
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();
	args = eval(list(...), envir = envir__);
	freezer$push(Lapply__$sequence, f, l, args, envir__ = environment(f));
	NULL
});
setMethod('lapply_dispatchFinalize', 'ParallelizeBackend', function(self) { 
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();
	parallelize_setEnable(F);
	r = lapply(1:freezer$N(), function(i) {
		call = freezer$call(i);
		call = callEvalArgs(call);
		r = Do.call(call$f, call$args, envir = call$envir);
	});
	freezer$finalizeResults();
	parallelize_setEnable(T);
	r
});
setMethod('lapply_results', 'ParallelizeBackend', function(self, r) { 
	stop('ParallelizeBackend: result retrieval only supported for asynchroneous backends.');
});
setMethod('parallelize_backend', 'ParallelizeBackend', function(self, call_) {
	with(Lapply_getConfig(), if (self@offline) {
		parallelizeOfflineStep(call_, Lapply_config = Lapply_getConfig());
	} else {
		Lapply_initialze_probing();
		for (i in 1:parallel_stack) {
			r = performParallelizationStep(self, call_, Lapply_config = Lapply_getConfig());
			if (all(class(r) != 'Lapply_error')) break;
		}
		r
	});
});
setMethod('performParallelizationStep', 'ParallelizeBackend', function(self, call_, Lapply_config) {
	parallelizeStep(call_, Lapply_config = Lapply_config);
});
setMethod('finalizeParallelization', 'ParallelizeBackend', function(self, r)r);
setMethod('pollParallelization', 'ParallelizeBackend',
	function(self, options = list())list(continue = F, message = '')
);

#
#		<p> parallelization state
#

# <A> running in '.' will not create sub-directory
#	used by remoting computations and already changing to remote stateDir
parallelizationStatePath = function(self, tag = '', ..., ext = '.RData') {
	tagStr = sprintf(tag, ...);
	path = if (self@config$stateDir == '.')
		sprintf('./%s%s', tagStr, ext) else
		sprintf('%s/parallelization_%s/%s%s', self@config$stateDir, self@signature, tagStr, ext);
	Log(sprintf('parallelization path: %s', path), 7);
	path
}
parallelizationStateObjects = c(
	'Lapply_globalConfig__', 'Lapply__', 'Lapply_executionState__', 'Lapply_backend__'
);
saveParallelizationStatePath = function(self, path = NULL) {
	if (is.null(path)) path = parallelizationStatePath(self, 'state');
	Log(sprintf('Saving state to %s', path), 5);
	parallelizationStateObjects = names(as.list(parallelize_env));
	Save(parallelizationStateObjects, file = path, symbolsAsVectors = T, envir = parallelize_env);
}
restoreParallelizationStatePath = function(self, path = NULL) {
	if (is.null(path)) path = parallelizationStatePath(self, 'state');
	Load(file = path, envir = parallelize_env);
}

setMethod('initScheduling', 'ParallelizeBackend', function(self, call_) {
	stateDir = parallelizationStatePath(self, '', ext = '');
	Log(sprintf('State dir: %s', stateDir), 5);
	Dir.create(stateDir, recursive = T);
	saveParallelizationStatePath(self);
});
setMethod('saveParallelizationState', 'ParallelizeBackend', function(self) {
	saveParallelizationStatePath(self);
});
setMethod('restoreParallelizationState', 'ParallelizeBackend', function(self) {
	restoreParallelizationStatePath(self);
});
setMethod('scheduleNextParallelization', 'ParallelizeBackend', function(self, call_) {
	NULL
});
setMethod('getResult', 'ParallelizeBackend', function(self) {
	if (self@config$doSaveResult)
		r = get(Load(file = parallelizationStatePath(self, 'result'))[1]) else
		stop(sprintf('result was not saved for signature %s', self@signature));
});

#
#	<p> local execution
#

setClass('ParallelizeBackendLocal',
	contains = 'ParallelizeBackend',
	representation = list(),
	prototype = list()
);
setMethod('initialize', 'ParallelizeBackendLocal', function(.Object, config, ...) {
	.Object = callNextMethod(.Object, config = config, ...);
	Dir.create(config$path, recursive = T);
	.Object
});
setMethod('lapply_dispatchFinalize', 'ParallelizeBackendLocal', function(self) { 
	Log(sprintf('Dispatching locally, tmp: %s', self@config$path), 5);
	parallelize_setEnable(F);
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();
	r = lapply(1:freezer$N(), function(i) {
		call = freezer$call(i);
		call = callEvalArgs(call);
		r = Do.call(call$f, call$args, envir = call$envir);
		freezer$pushResults(r);
		r
	});
	freezer$finalizeResults();
	save(r, file = sprintf('%s/sequence-%d.RData', self@config$path, Lapply__$sequence));
	parallelize_setEnable(T);
	NULL
});

#
#	<p> SNOW execution
#

setClass('ParallelizeBackendSnow',
	contains = 'ParallelizeBackend',
	representation = list(),
	prototype = list()
);
setMethod('initialize', 'ParallelizeBackendSnow', function(.Object, config, ...) {
	.Object = callNextMethod(.Object, config = config, ...);
	args = List_(config[c('sourceFiles', 'localNodes', 'splitN', 'libraries')], rm.null = T);
	args$libraries = c(args$libraries, 'parallelize.dynamic');
	do.call('specifyCluster', args);
	.Object
});
setMethod('lapply_dispatchFinalize', 'ParallelizeBackendSnow', function(self) { 
	Log(sprintf('Dispatching locally, tmp: %s', self@config$path), 5);
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();
	r = clapply(freezer$getCalls(), function(call) {
		parallelize_setEnable(F);
		call = callEvalArgs(call);
		Do.call(call$f, call$args, envir = call$envir)
	});
	freezer$pushResults(r);
	freezer$unlistResults();
	freezer$finalizeResults();
	NULL
});

#
#	<p> OGS execution
#

#
#	ParallelizeBackendOGS S4 class
#

.ParallelizeBackendOGSstateClass = setRefClass('ParallelizeBackendOGSstate',
	fields = list( steps = 'list', chunks = 'list', logPath = 'character' ),
	methods = list(
	initialize = function(...) {
		steps <<- list();
		chunks <<- list();
		logPath <<- '';
		.self
	},
	log = function() { if (logPath != '') save(.self, file = logPath); },
	pushStep = function(jid) {
		steps[[length(steps) + 1]] <<- jid;
		.self$log();
	},
	pushChunks = function(jids) {
		chunks[[length(chunks) + 1]] <<- jids;
		.self$log();
	},
	chunksJids = function() { if (!length(chunks)) c() else chunks[[length(chunks)]]; },
	setLogPath = function(path) {
		logPath <<-path;
		if (file.exists(logPath)) file.remove(logPath);
	}
	)
);
.ParallelizeBackendOGSstateClass$accessors(names(.ParallelizeBackendOGSstateClass$fields()));

setClass('ParallelizeBackendOGS',
	contains = 'ParallelizeBackend',
	representation = list(jids = 'ParallelizeBackendOGSstate'),
	prototype = list(jids = .ParallelizeBackendOGSstateClass$new())
);
.ParallelizeBackendOGSDefaultConfig = list(
	qsubOptions = '--queue all.q'
);
setMethod('initialize', 'ParallelizeBackendOGS', function(.Object, config, ...) {
	# <p> super-class
	config = merge.lists(.ParallelizeBackendOGSDefaultConfig, config);
	.Object = callNextMethod(.Object, config = config, ...);
	# <p> OGS initialization
	Log('initializing OGS', 6);
	.Object@offline = T;
	# <p> RNG
	RNGkind("L'Ecuyer-CMRG");
	set.seed(as.integer(Sys.time()));

	# <p> jid state
	.Object@jids$setLogPath(parallelizationStatePath(.Object, 'jids'));
	.Object
});

setMethod('initScheduling', 'ParallelizeBackendOGS', function(self, call_) {
	callNextMethod(self);
	# <p> dir initialization
	dir = parallelizationStatePath(self, tag = '', ext = '');
	Dir.create(dir, recursive = T);
	# <p> initialize files
	sentinelPath = parallelizationStatePath(self, 'sentinel');
	if (file.exists(sentinelPath)) file.remove(sentinelPath);
});

.parallelizationStepOGS = function(call_, pathHandover) {
	# <!> potential race condition with scheduleNextParallelization
	r0 = get(Load(file = pathHandover, Load_sleep = 5)[1]);
	Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
	Lapply_backend__@jids$pushStep(r0$jid);
	parallelize_setEnable(T);	# default is off
	parallelizeOfflineStep(call_, Lapply_config = Lapply_getConfig());
}

.freezeCallOGS = function(self, ..f, ...,
	freeze_file = tempfile(), freeze_control = list(), waitForJids = c(),
	patterns = 'qsub', cwd = NULL, ssh_host = 'localhost', ssh_source_file = NULL,
	qsubPath = parallelizationStatePath(self, 'qsub', ext = ''), envir = NULL) {

	path = freezeCall(freeze_f = ..f, ...,
		freeze_file = freeze_file, freeze_save_output = T, freeze_control = freeze_control,
		freeze_envir = NULL, freeze_objects = 'parallelize_env');
	wrap = frozenCallWrap(path, freeze_control);
	qsubOptions = sprintf('%s --outputDir %s %s',
		self@config$qsubOptions,
		qs(qsubPath),
		if (!length(waitForJids)) '' else sprintf('--waitForJids %s', paste(waitForJids, collapse = ','))
	);
	r = System(wrap, 5, patterns = patterns, qsubOptions = qsubOptions, cwd = cwd,
		ssh_host = ssh_host, ssh_source_file = ssh_source_file);
	r
}

# we use the freeze/thaw mechanism and a handover such that restoring the state would
#	destroy handover changes, the saving still occurs for tracking purposes
setMethod('restoreParallelizationState', 'ParallelizeBackendOGS', function(self) {
	NULL
});

setMethod('scheduleNextParallelization', 'ParallelizeBackendOGS', function(self, call_) {
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	c = Lapply_getConfig();
	freeze_control = list(
		sourceFiles = self@config$sourceFiles,
		libraries = self@config$libraries,
		objects = parallelizationStateObjects,
		logLevel = Log.level()
	);
	path = parallelizationStatePath(self, 'rampUp:%03d', Lapply_executionState__$rampUp);
	pathHandover = parallelizationStatePath(self, 'rampUp:%03d_handover', Lapply_executionState__$rampUp);
	# <i> gather information from previous step
	#qacct -j 257
	# new path for each rampUp due to potential race condition
	r0 = .freezeCallOGS(self, ..f = .parallelizationStepOGS, call_,
		# .parallelizationStepOGS
		pathHandover = pathHandover,
		# freeze
		freeze_file = path, freeze_control = freeze_control, waitForJids = self@jids$chunksJids())
	save(r0, file = pathHandover);
	r0
});

setMethod('lapply_dispatchFinalize', 'ParallelizeBackendOGS', function(self) { 
	Log(sprintf('Dispatching locally, tmp: %s', self@config$path), 5);
	Lapply_executionState__ = get('Lapply_executionState__', envir = parallelize_env);
	Lapply__ = get('Lapply__', envir = parallelize_env);
	freezer = Lapply_executionState__$currentFreezer();

	# <p> setup
	c = Lapply_getConfig();
	freeze_control = list(
		sourceFiles = self@config$sourceFiles,
		libraries = self@config$libraries,
		objects = parallelizationStateObjects,
		logLevel = Log.level()
	);
	calls = freezer$getCalls();
	# <p> split up calls
	N = c$parallel_count;	# No of splits
	idcs = splitListIndcs(length(calls), N);

	f0 = function(calls) {
		parallelize_setEnable(F);
		lapply(calls, function(call)try(do.call(call$f, call$args, envir = call$envir)));
	}
	r = lapply(1:dim(idcs)[1], function(.i) {
		path = parallelizationStatePath(self, 'sequence:%03d_chunk:%05d', Lapply__$sequence, .i);
		mycalls = calls[idcs[.i, 1]:idcs[.i, 2]];	# force evaluation
		mycalls = lapply(mycalls, callEvalArgs);
		r = .freezeCallOGS(self, f0, calls = mycalls, freeze_file = path, freeze_control = freeze_control);
		r = c(r, list(file = path, from = idcs[.i, 1], to = idcs[.i, 2]));
		r
	});
	self@jids$pushChunks(list.kp(r, 'jid', do.unlist = T));
	freezer$pushResults(r);
	#freezer$unlistResults();
	freezer$finalizeResults();
	NULL
});

setMethod('finalizeParallelization', 'ParallelizeBackendOGS', function(self, r) {
	if (self@config$doSaveResult)
		save(r, file = parallelizationStatePath(self, 'result'));
	sentinel = list(signature = self@signature);
	save(sentinel, file = parallelizationStatePath(self, 'sentinel'));
	r
});

.progressStat = function(jidsTasks, i, jidsRunning) {
	jidsTask = if (length(jidsTasks) < i) NULL else jidsTasks[[i]];
	jidsPending = intersect(jidsTask, jidsRunning);
	N = length(jidsTask);
	Npending = length(jidsPending);
	r = list(N = N, Npending = Npending, Ncomplete = N - Npending, complete = 1 - Npending / N);
	r
}
.stdProgressFormat = list(title = '%-30s', N = '%4d', progress = '%25s', Perc = '%3.0f%%');
progressString = function(stat, title = 'Task', format = .stdProgressFormat, NanString = '----') {
	format = merge.lists(.stdProgressFormat, format);
	L = nchar(sprintf(format$progress, '-'));	# length progress bar
	progressBar = if (is.nan(stat$complete)) sprintf('%-*s', L, 'count pending') else
		paste(c(rep('#', round(stat$complete * L, 0)),
			rep('.', round((1 - stat$complete) * L, 0))), collapse = '');
	values = list(title = title, Perc = floor(100 * stat$complete), N = stat$N, progress = progressBar);
	r = unlist(nlapply(format, function(n) {
		if (is.nan(values[[n]])) NanString else sprintf(format[[n]], values[[n]])
	}));
	r = paste(r, collapse = ' ');
	r
		
}

.pollJids = function(...) {
	qstat = System("qstat -u \\* -xml | xml sel -t -m '//JB_job_number' -v 'text()' -o ' '", 5,
		..., return.output = T);
	jids = fetchRegexpr('(\\d+)', qstat$output, captures = T);
	jids
}

.pollMessageRaw = function(jids, qstat_jids) {
	N = max(length(jids$steps), length(jids$chunks));
	msg = as.vector(sapply(1:N, function(i) {
		psc = .progressStat(jids$chunks, i, qstat_jids);
		pss = .progressStat(jids$steps, i, qstat_jids);
		c(
			progressString(psc, title = sprintf('  Parallelization %d', i)),
			progressString(pss, title = sprintf('Rampdown %d', i))
		)
	}));
	msg
}

.pollMessage = function(msg, continue) {
	header = paste(rep('-', 79), collapse = '');
	conclusion = if (continue) 'Further scheduling pending' else 'Computation complete';
	#messageRaw = paste(msg, collapse = "\n");
	#message = paste(c(header, messageRaw, header, conclusion, '', ''), collapse = "\n");
	message = c(header, msg, header, conclusion);
	message
}

setMethod('pollParallelization', 'ParallelizeBackendOGS', function(self, options = list()) {
	continue = !file.exists(parallelizationStatePath(self, 'sentinel'));
	# <p> fetch jids
	qstat_jids = .pollJids();
	# <p> restore state locally
	Load(file = parallelizationStatePath(self, 'state'));
	# <p> raw message
	Lapply_backend__ = get('Lapply_backend__', envir = parallelize_env);
	message = .pollMessageRaw(Lapply_backend__@jids, qstat_jids);
	# <p> refine
	message = .pollMessage(message, continue);
	#			=~ m{(\d+)}sog)
	r = list(continue = continue, message = message);
	r
});

#
#	ParallelizeBackendOGSremote S4 class
#

.ParallelizeBackendOGSremoteDefaultConfig = list(
	remote = 'localhost:parallelize_projects'
);
setClass('ParallelizeBackendOGSremote',
	contains = 'ParallelizeBackend',
	representation = list(jids = 'ParallelizeBackendOGSstate'),
	prototype = list(jids = .ParallelizeBackendOGSstateClass$new())
);
setMethod('initialize', 'ParallelizeBackendOGSremote', function(.Object, config, ...) {
	# <p> super-class
	config = merge.lists(.ParallelizeBackendOGSDefaultConfig, config);
	.Object = callNextMethod(.Object, config = config, ...);
	# <p> OGS initialization
	Log('initializing OGS for remote execution', 6);
	.Object@offline = T;
	# restart on other host
	.Object
});

.remoteConfigForOGSremote = function(stateDir = '.') {
	Lapply_remote_config = Lapply_getConfig();
	backendConfig = merge.lists(
		Lapply_remote_config$backendConfig,
		list(backend = 'OGS', stateDir = stateDir, logLevel = Log.level())
	);
	Lapply_remote_config$backends[[Lapply_remote_config$backend]] = 
		Lapply_remote_config$backendConfig = backendConfig;
	Lapply_remote_config
}
.OGSremoteFile = function(self, tag = '', ext = '.RData') {
	Lapply_remote_config = .remoteConfigForOGSremote(stateDir = self@config$remote);
	remoteDummy = new('ParallelizeBackendOGS', config =
		Lapply_remote_config$backendConfig, signature = self@signature);
	remoteDir = parallelizationStatePath(remoteDummy, tag = tag, ext = ext);
	remoteDir
}
.OGSremoteWorkingDir = function(self).OGSremoteFile(self, tag = '', ext = '')


setMethod('initScheduling', 'ParallelizeBackendOGSremote', function(self, call_) {
	callNextMethod(self);
	r = with(self@config, {
	# <p> check starting sentinel
	sentinelPath = parallelizationStatePath(self, 'OGSremote_sentinel');
	if (file.exists(sentinelPath) && !self@config$force_rerun) {
		Log(sprintf('Signature %s already scheduled.', self@signature), 5);
		return(NULL);
	}
# 	# prevent further parallelize calls from re-initializing
# 	c = Lapply_getConfig();
# 	c$backendConfig$force_rerun = F;
# 	Lapply_setConfig(c);

	# <p> establish start sentinel
	sentinel = list(signature = self@signature);
	save(sentinel, file = sentinelPath);

	# <p> setup remote environment
	#remoteDir = sprintf('%s/%s', remote, self@signature);
	remoteDir = .OGSremoteWorkingDir(self);
	sp = splitPath(remoteDir, ssh = T);
	Log(sprintf('setting up parallelization step in dir %s', remoteDir), 5);
	ignore.shell = Log.level() < 5;
	Dir.create(remoteDir, recursive = T, ignore.shell = ignore.shell);
	# either copy source files or explicitely spcified copyFiles (important for dirs);
	copyFiles = if (length(self@config$copyFiles) > 0) copyFiles else sourceFiles;
	File.copy(copyFiles, remoteDir, ignore.shell = ignore.shell, recursive = T);
	# clear jids
	File.remove(.OGSremoteFile(self, 'jids'));

	# <p> create remote wrappers
	parallelize_remote = function(call_, Lapply_config) {
		parallelize_initialize(Lapply_config = Lapply_config, backend = Lapply_config$backend);
		r = parallelize_internal(call_, parallelize_wait = F);
	};
	# <p> start rampup on remote host
	freeze_control = list(
		sourceFiles = self@config$sourceFiles,
		libraries = self@config$libraries,
		logLevel = Log.level(),
		freeze_relative = T
	);
	remoteConfig = .remoteConfigForOGSremote(stateDir = '.');
	r = .freezeCallOGS(self, parallelize_remote,
		# parallelize_remote
		callEvalArgs(call_), Lapply_config = remoteConfig,
		# freeze
		freeze_control = freeze_control,
		freeze_file = sprintf('%s/rampUp:000.RData', remoteDir),
		# System
		patterns = c('cwd', 'qsub', 'ssh'),
		cwd = sp$path, ssh_host = sp$userhost, qsubPath = sprintf('%s/qsub', sp$path),
		ssh_source_file = self@config$ssh_source_file);
	});
	self@jids$pushStep(r$jid);
	r
});

# instead of doing something here, we poll the remote backend
setMethod('performParallelizationStep', 'ParallelizeBackendOGSremote',
	function(self, call_, Lapply_config) {
	# prevent from completing computation, result has to be gathered by polling
	Lapply_error();

	if (0) {
	stop('ParallelizeBackendOGSremote backend is a delegating backend and does not perform parallelization itself. Use the following to monitor this backend in a loop:
		r = NULL;
		p = pollParallelization(self);
		if (!p$continue) r = getResult(Lapply_backend__);
		r
	');
	}
});

.catVectorAsLine = function(message, width = options('width')$width) {
	messagePadded = sapply(message, function(line) sprintf('%s%*s', line, width - nchar(line) - 1, ' '));
	messageInALine = paste(messagePadded, collapse = '');
	cat(messageInALine);
	flush.console();
	cat("\r");
}

.catVector = function(message, width = options('width')$width, clear = T, padLines = 40) {
	if (clear) cat(paste(rep("\n", 100), collapse = ''));
	cat(paste(c(message, ''), collapse = "\n"));
	if (padLines > 0) cat(paste(rep("\n", padLines), collapse = ''));
}

setMethod('pollParallelization', 'ParallelizeBackendOGSremote', function(self,
	options = list(printProgress = T)) {
	# <p> overwrite backend configuration
	remote_config = .remoteConfigForOGSremote();
	jidFile = .OGSremoteFile(self, 'jids');
	jids = get(Load(file = jidFile, Load_sleep = 30, Load_retries = 60)[[1]]);
	qstat_jids = .pollJids(patterns = 'ssh',
		ssh_host = splitPath(jidFile, ssh = T)$userhost, ssh_source_file = self@config$ssh_source_file);
	print(jids); print(qstat_jids);
	message = .pollMessageRaw(jids, qstat_jids);
	# <p> check for completion
	continue = !File.exists(.OGSremoteFile(self, 'sentinel'));
	# <p> add rampup
	message = c(
		progressString(.progressStat(self@jids$steps, 1, qstat_jids), title = 'Rampup 1'),
		message
	);
	# <p> refine
	message = .pollMessage(message, continue);
	.catVector(message);
	r = list(message = message, continue = continue);
	r
});

setMethod('lapply_dispatchFinalize', 'ParallelizeBackendOGSremote',
	function(self) { NULL });

setMethod('getResult', 'ParallelizeBackendOGSremote', function(self) {
	r = get(Load(file = .OGSremoteFile(self, 'result'))[1]);
	r
});

