GDataStreamByteOrder<-c("big-endian" = 0,
	"little-endian" = 1,
	"host-endian" = 2)
storage.mode(GDataStreamByteOrder) <- 'integer'
class(GDataStreamByteOrder) <- 'enums' 

GDataStreamNewlineType<-c("lf" = 0,
	"cr" = 1,
	"cr-lf" = 2,
	"any" = 3)
storage.mode(GDataStreamNewlineType) <- 'integer'
class(GDataStreamNewlineType) <- 'enums' 

GFileAttributeType<-c("invalid" = 0,
	"string" = 1,
	"byte-string" = 2,
	"boolean" = 3,
	"uint32" = 4,
	"int32" = 5,
	"uint64" = 6,
	"int64" = 7,
	"object" = 8)
storage.mode(GFileAttributeType) <- 'integer'
class(GFileAttributeType) <- 'enums' 

GFileAttributeStatus<-c("unset" = 0,
	"set" = 1,
	"error-setting" = 2)
storage.mode(GFileAttributeStatus) <- 'integer'
class(GFileAttributeStatus) <- 'enums' 

GFileType<-c("unknown" = 0,
	"regular" = 1,
	"directory" = 2,
	"symbolic-link" = 3,
	"special" = 4,
	"shortcut" = 5,
	"mountable" = 6)
storage.mode(GFileType) <- 'integer'
class(GFileType) <- 'enums' 

GFileMonitorEvent<-c("changed" = 0,
	"changes-done-hint" = 1,
	"deleted" = 2,
	"created" = 3,
	"attribute-changed" = 4,
	"pre-unmount" = 5,
	"unmounted" = 6)
storage.mode(GFileMonitorEvent) <- 'integer'
class(GFileMonitorEvent) <- 'enums' 

GIOErrorEnum<-c("failed" = 0,
	"not-found" = 1,
	"exists" = 2,
	"is-directory" = 3,
	"not-directory" = 4,
	"not-empty" = 5,
	"not-regular-file" = 6,
	"not-symbolic-link" = 7,
	"not-mountable-file" = 8,
	"filename-too-long" = 9,
	"invalid-filename" = 10,
	"too-many-links" = 11,
	"no-space" = 12,
	"invalid-argument" = 13,
	"permission-denied" = 14,
	"not-supported" = 15,
	"not-mounted" = 16,
	"already-mounted" = 17,
	"closed" = 18,
	"cancelled" = 19,
	"pending" = 20,
	"read-only" = 21,
	"cant-create-backup" = 22,
	"wrong-etag" = 23,
	"timed-out" = 24,
	"would-recurse" = 25,
	"busy" = 26,
	"would-block" = 27,
	"host-not-found" = 28,
	"would-merge" = 29,
	"failed-handled" = 30)
storage.mode(GIOErrorEnum) <- 'integer'
class(GIOErrorEnum) <- 'enums' 

GPasswordSave<-c("never" = 0,
	"for-session" = 1,
	"permanently" = 2)
storage.mode(GPasswordSave) <- 'integer'
class(GPasswordSave) <- 'enums' 

GMountOperationResult<-c("handled" = 0,
	"aborted" = 1,
	"unhandled" = 2)
storage.mode(GMountOperationResult) <- 'integer'
class(GMountOperationResult) <- 'enums' 

GSeekType<-c("cur" = 0,
	"set" = 1,
	"end" = 2)
storage.mode(GSeekType) <- 'integer'
class(GSeekType) <- 'enums' 

GEmblemOrigin<-c("unknown" = 0,
	"device" = 1,
	"livemetadata" = 2,
	"tag" = 3)
storage.mode(GEmblemOrigin) <- 'integer'
class(GEmblemOrigin) <- 'enums' 

GDriveStartFlags<-c("none" = 0)
storage.mode(GDriveStartFlags) <- 'integer'
class(GDriveStartFlags) <- 'enums' 

GDriveStartStopType<-c("unknown" = 0,
	"shutdown" = 1,
	"network" = 2,
	"multidisk" = 3,
	"password" = 4)
storage.mode(GDriveStartStopType) <- 'integer'
class(GDriveStartStopType) <- 'enums' 

GSocketFamily<-c("invalid" = 0,
	"unix" = 1,
	"ipv4" = 2,
	"ipv6" = 3)
storage.mode(GSocketFamily) <- 'integer'
class(GSocketFamily) <- 'enums' 

GSocketType<-c("invalid" = 0,
	"stream" = 1,
	"datagram" = 2,
	"seqpacket" = 3)
storage.mode(GSocketType) <- 'integer'
class(GSocketType) <- 'enums' 

GSocketProtocol<-c("unknown" = -1,
	"default" = 0,
	"tcp" = 6,
	"udp" = 17,
	"sctp" = 132)
storage.mode(GSocketProtocol) <- 'integer'
class(GSocketProtocol) <- 'enums' 

GIOCondition<-c("in" = 0,
	"out" = 1,
	"pri" = 2,
	"err" = 3,
	"hup" = 4,
	"nval" = 5)
storage.mode(GIOCondition) <- 'integer'
class(GIOCondition) <- 'enums' 

GAppInfoCreateFlags<-c("one" = 1,
	"eeds-terminal" = 2)
storage.mode(GAppInfoCreateFlags) <- 'numeric'
class(GAppInfoCreateFlags) <- 'flags' 

GFileAttributeInfoFlags<-c("none" = 1,
	"copy-with-file" = 2,
	"copy-when-moved" = 4)
storage.mode(GFileAttributeInfoFlags) <- 'numeric'
class(GFileAttributeInfoFlags) <- 'flags' 

GFileQueryInfoFlags<-c("ne" = 1,
	"follow-symlinks" = 2)
storage.mode(GFileQueryInfoFlags) <- 'numeric'
class(GFileQueryInfoFlags) <- 'flags' 

GFileCreateFlags<-c("none" = 1,
	"private" = 2)
storage.mode(GFileCreateFlags) <- 'numeric'
class(GFileCreateFlags) <- 'flags' 

GMountMountFlags<-c("none" = 1)
storage.mode(GMountMountFlags) <- 'numeric'
class(GMountMountFlags) <- 'flags' 

GMountUnmountFlags<-c("none" = 1,
	"force" = 2)
storage.mode(GMountUnmountFlags) <- 'numeric'
class(GMountUnmountFlags) <- 'flags' 

GFileCopyFlags<-c("none" = 1,
	"overwrite" = 2,
	"backup" = 4,
	"nofollow-symlinks" = 8,
	"all-metadata" = 16,
	"no-fallback-for-move" = 32)
storage.mode(GFileCopyFlags) <- 'numeric'
class(GFileCopyFlags) <- 'flags' 

GFileMonitorFlags<-c("none" = 1,
	"watch-mounts" = 2)
storage.mode(GFileMonitorFlags) <- 'numeric'
class(GFileMonitorFlags) <- 'flags' 

GAskPasswordFlags<-c("need-password" = 1,
	"need-username" = 2,
	"need-domain" = 4,
	"saving-supported" = 8,
	"anonymous-supported" = 16)
storage.mode(GAskPasswordFlags) <- 'numeric'
class(GAskPasswordFlags) <- 'flags' 

GOutputStreamSpliceFlags<-c("none" = 1,
	"close-source" = 2,
	"close-target" = 4)
storage.mode(GOutputStreamSpliceFlags) <- 'numeric'
class(GOutputStreamSpliceFlags) <- 'flags' 

