library(rjson);

Tuple = setRefClass("Tuple",
                    fields = list(
                      #The tuple's id - this is a string to support languages lacking 64-bit precision
                      id="character",
                      #The id of the component that created this tuple
                      #XXX integer?
                      comp="character",
                      #The id of the stream this tuple was emitted to
                      #XXX integer?
                      stream="character",
                      #The id of the task that created this tuple
                      #XXX integer?
                      task="character",
                      #All input values in this tuple
                      input="list",
                      #All output values for this tuple
                      output="vector",
                      #All tuples used to create field @tuple.out
                      anchors="vector"
                      
                      #parse="function"
                    ));

Tuple$methods(
  parse = function(json=character) {
    t = fromJSON(json);
    .self$id=as.character(t$id);
    .self$comp=as.character(t$comp);
    .self$stream=as.character(t$stream);
    .self$task=as.character(t$task);
    .self$input=t$tuple;
    .self$output=vector(mode="character");
    .self$anchors=vector(mode="character");
    .self;
  }
);

Storm = setRefClass("Storm",
                    fields = list(
                      tuple = "Tuple",
                      setup = "list",
                      lambda="function"
                    )
);

#Storm$methods(
#  lambda = function(s=Storm) {
#    s$log(c("skipping tuple id='",s$tuple$id,"'"));
#  }
#);
#Storm$lambda = function(s) {
#  s$log(c("skipping tuple id='",s$tuple$id,"'"));
#};

Storm$methods(
  initialize = function() {
    #.self$lambda = Storm$lambda;
    .self$lambda = function(s) {
      s$log(c("skipping tuple id='",s$tuple$id,"'"));
    };
    .self;
  }
);

Storm$methods(
  run = function() {
    buf = "";

    x.stdin = file("stdin");
    open(x.stdin);

    cat(paste('{"pid": ',Sys.getpid(),'}',"\nend\n",sep=""));
    cat('{"command": "emit", "anchors": [], "tuple": ["bolt initializing"]}\nend\n');

    while (TRUE) {
      rl = as.character(readLines(con=x.stdin,n=1,warn=FALSE));
      if (length(rl) == 0) {
        close(x.stdin);
        break;
      }
      if (rl == "end" || rl == "END") {
        #print(buf);
        .self$process_json(buf);
        buf = "";
      }
      else {
        buf = paste(buf,rl,"\n",sep="");
      }
    }
  }
);

Storm$methods(
  process_json = function(json=character) {
    t = NULL;
    if (is.list(json)) {
      t = fromJSON(unlist(json));
    }
    else {
      t = fromJSON(json);
    }
    # we assume that it is a tuple if the "tuple" field is not empty.
    if (is.list(t)) {
    if (!is.null(t$tuple)) {
      .self$process_tuple(json);
    }
    # we assume that it is a setup if the "pidDir" field is not empty.
    else if (!is.null(t$pidDir)) {
      .self$setup = t;
      .self$process_setup(t);
    }
    # otherwise, not sure what it is.  log an error and stop.
    else {
#TODO      .self$log(paste("unrecognized JSON:\n'",json,"'\n",sep=""));
      #stop("unrecognized JSON: ",json);
    }
    }
  }
);

Storm$methods(
  process_tuple = function(json=character) {
    #print(Tuple);
    #tuple = getRefClass("Tuple")$new();
    tt = Tuple$new();
    tt$parse(json);
    .self$tuple = tt;
    .self$lambda(.self);      
    #    print(.self$input.tuple);
  }
);
Storm$methods(
  process_setup = function(json=character) {
    file.create(paste(json$pidDir,"/",Sys.getpid(),sep=""))
    #    print(.self$input.setup);
  }
);

Storm$methods(
  log = function(msg=character) {
    cat(c(
      '{',
      '"command": "log",',
      '"msg": "',msg,'"',
      '}\n',
      'end\n'
    ),sep="");
  }
);
Storm$methods(
  ack = function(tuple=Tuple) {
    cat(c(
      '{\n',
      '\t"command": "ack",\n',
      '\t"id": "',tuple$id,'"\n',
      '}\n',
      'end\n'
    ),sep="");
  }
);
Storm$methods(
  fail = function(tuple=Tuple) {
    cat(c(
      '{',
      '"command": "fail",',
      '"id": "',tuple$id,'"',
      '}\n',
      'end\n'
    ),sep="") 
  }
);
Storm$methods(
  emit = function(tuple=Tuple) {
    t.prefix="[";
    t.suffix="]";
    #if (length(tuple$output) == 1) {
    #  t.prefix="[";
    #  t.suffix="]";
    #}
    a.prefix="[";
    a.suffix="]";
    #if (length(tuple$anchors) == 1) {
    #  a.prefix="[";
    #  a.suffix="]";
    #}
    cat(c(
      '{',
      '"command": "emit", ',
      
      #TODO enable reliable tuples, see
      #https://github.com/nathanmarz/storm/wiki/Guaranteeing-message-processing
      #'"id": "',as.character(runif(1) * 2**64 - 1),'", ',
      
      #TODO enable multi-anchoring, see:
      #https://github.com/nathanmarz/storm/wiki/Guaranteeing-message-processing      
      '"anchors": ',a.prefix,toJSON(tuple$id),a.suffix,', ',

      #TODO enable id of stream tuple is emitted to, see:
      #https://github.com/nathanmarz/storm/wiki/Multilang-protocol
      #'"stream": "',tuple$stream,'", ',

      #TODO enable "emit direct", see:
      #https://github.com/nathanmarz/storm/wiki/Multilang-protocol
      #'"task": "',tuple$task,'", ',

      #'"tuple": ',t.prefix,toJSON(tuple$output),t.suffix,
      '"tuple": ',toJSON(tuple$output),
      '}\n',
      'end\n'
    ),sep="");
  }
);
