#' @import R6
library(R6)
Deserializer <- R6Class(
  'Deserializer',
  public = list(
    bytes = NULL,
    offset = NULL,
    object = NULL,
    initialize = function(bytes, offset){
      if (!is.raw(bytes)) stop("bytes must be of type raw")
      if (is.null(offset)){
        self$offset = 1
      } else {
        if (!is.numeric(offset)) stop("offset must be of type numeric")
        self$offset = offset 
      }
      self$bytes = bytes
      version = self$readObject()
      if (version != TSON_SPEC_VERSION) stop(paste0("TSON spec version mismatch, found ", version, " expected " , TSON_SPEC_VERSION))
      self$object = self$readObject()
    },
    readType = function() {
      type = self$bytes[self$offset]
      self$offset = self$offset + 1
      return (type)
    },
    readObject = function(){
      type = self$readType()
      if (type == NULL_TYPE){
        return(NULL)
      } else if (type == STRING_TYPE){
        return(self$readString())
      } else if (type == INTEGER_TYPE){
        return(self$readInteger())
      } else if (type == DOUBLE_TYPE){
        return(self$readDouble())
      } else if (type == BOOL_TYPE){
        return(self$readBool())
      } else if (type == LIST_TYPE){
        return(self$readList())
      } else if (type == MAP_TYPE){
        return(self$readMap())
      } else if (type == LIST_UINT8_TYPE){
        return(self$readUint8List())
      } else if (type == LIST_UINT16_TYPE){
        return(self$readUint16List())
      } else if (type == LIST_UINT32_TYPE){
        return(self$readUint32List())
      } else if (type == LIST_INT8_TYPE){
        return(self$readInt8List())
      } else if (type == LIST_INT16_TYPE){
        return(self$readInt16List())
      } else if (type == LIST_INT32_TYPE){
        return(self$readInt32List())
      } else if (type == LIST_INT64_TYPE){
        return(self$readInt64List())
      } else if (type == LIST_FLOAT32_TYPE){
        return(self$readFloat32List())
      } else if (type == LIST_FLOAT64_TYPE){
        return(self$readFloat64List())
      } else if (type == LIST_STRING_TYPE){
        return(self$readStringList())
      } 
    },
    readString = function(){
      bytes = self$bytes[self$offset:length(self$bytes)]
      object = readBin(bytes, character(), n=1)
      self$offset = self$offset + nchar(object) + 1
      return (object)
    },
    readInteger = function(){
      bytes = self$bytes[self$offset:length(self$bytes)]
      object = readBin(bytes, integer(), n=1)
      self$offset = self$offset + 4
      return (object)
    },
    readDouble = function(){
      bytes = self$bytes[self$offset:length(self$bytes)]
      object = readBin(bytes, double(), n=1)
      self$offset = self$offset + 8
      return (object)
    },
    readBool = function(){
      bytes = self$bytes[self$offset:length(self$bytes)]
      object = readBin(bytes, integer(), n=1, size=1)
      self$offset = self$offset + 1
      return (object == 1)
    },
    readLength = function(){
      len = self$readInteger()
      if (len < 0) stop("Length must be greater or equals to 0")
      return (len)
    },
    readList = function(){
      len = self$readLength()
      list = NULL
      if (len > 0){
        list = lapply(seq(1,len), function(i) self$readObject())
      } else {
        list = list()
      }
      return (list)
    },
    readMap = function(){
      len = self$readLength()
      list = NULL
      if (len > 0){
        keys = list()
        list = list()
        for (i in 1:len){
          keys[[i]] = self$readObject()
          value = self$readObject()
          if (is.null(value)){
            list[i] = list(NULL)
          } else {
            list[[i]] = value
          }
        }
        names(list) = keys
      } else {
        list = list()
      }
      return (list)
    },
    readTypedList = function(what, size, signed){
      len = self$readLength()
      bytes = self$bytes[self$offset:length(self$bytes)]
      object = readBin(bytes, what, n=len, size=size, signed = signed, endian =  "little")
      self$offset = self$offset + (len * size)
      return (object)
    },
    readUint8List = function() self$readTypedList(integer(), 1 , FALSE),
    readUint16List = function() self$readTypedList(integer(), 2, FALSE),
    readUint32List = function() self$readTypedList(integer(), 4, TRUE),
    readInt8List = function() self$readTypedList(integer(), 1 , TRUE),
    readInt16List = function() self$readTypedList(integer(), 2, TRUE),
    readInt32List = function() self$readTypedList(integer(), 4, TRUE),
    readFloat32List = function() self$readTypedList(double(), 4, TRUE),
    readFloat64List = function() self$readTypedList(double(), 8, TRUE),
    readStringList = function(){
      lengthInBytes = self$readLength()
      bytes = self$bytes[self$offset:(self$offset+lengthInBytes-1)]
      n = 0
      for (i in bytes){
        if (i == 0) n = n + 1
      }
      object = readBin(bytes, character(), n=n)
      self$offset = self$offset + lengthInBytes
      return (object)
    }
  )
)