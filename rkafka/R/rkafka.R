#Licensed to the Apache Software Foundation (ASF) under one
#or more contributor license agreements.  See the NOTICE file
#distributed with this work for additional information
#regarding copyright ownership.  The ASF licenses this file
#to you under the Apache License, Version 2.0 (the
#                                              "License"); you may not use this file except in compliance
#with the License.  You may obtain a copy of the License at
#
#http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing,
#software distributed under the License is distributed on an
#"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#KIND, either express or implied.  See the License for the
#specific language governing permissions and limitations
#under the License.
# Author:Shruti Gupta
###############################################################################
.onLoad <- function(libname, pkgname) {

	#specifying package location
  #.jaddClassPath("inst/java/rkafka-1.0-jar-with-dependencies.jar")
	.jpackage(pkgname, lib.loc = libname)

}

#function to create Producer

#Definition of the arguments
#   /**
#     * @param metadataBrokerList:String
#   *            !!Mandatory list of brokers used for bootstrapping knowledge
#   *            about the rest of the cluster format: host1:port1,host2:port2
#   *            ... default:localhost:9092
#   * 
#     * @param producerType:String
#   *            !!Mandatory specifies whether the messages are sent
#   *            asynchronously (async) or synchronously (sync) default:sync
#   * 
#     * @param compressionCodec:String
#   *            !!Mandatory specify the compression codec for all data
#   *            generated: none , gzip, snappy. default:none
#   * 
#     * @param serializerClass:String
#   *            !!Mandatory message encoder
#   *            default:kafka.serializer.StringEncoder
#   * 
#     * @param partitionerClass:String
#   *            --Optional name of the partitioner class for partitioning
#   *            events; default partition spreads data randomly
#               default:NULL
#   * 
#     * @param compressedTopics:String
#   *            --Optional allow topic level compression
#   * #               default:NULL
#     * @param queueBufferingMaxTime:String
#   *            --Optional(for Async Producer only) maximum time, in
#   *            milliseconds, for buffering data on the producer queue
#   * #               default:NULL
#     * @param queueBufferingMaxMessages:String
#   *            --Optional(for Async Producer only) the maximum size of the
#   *            blocking queue for buffering on the producer
#   *#               default:NULL
#     * @param queueEnqueueTimeoutTime:String
#   *            --Optional(for Async Producer only) 0: events will be enqueued
#   *            immediately or dropped if the queue is full -ve: enqueue will
#   *            block indefinitely if the queue is full +ve: enqueue will
#   *            block up to this many milliseconds if the queue is full
#   *#               default:NULL
#     * @param batchNumMessages:String
#   *            --Optional(for Async Producer only) the number of messages
#   *            batched at the producer
#   * #               default:NULL
#     * @return returns a Properties Object containing properties for the
#   *         Producer, to be passed to MuProducer class
#   */
rkafka.createProducer = function(metadataBrokerList,producerType="sync",compressionCodec="none",
		serializerClass="kafka.serializer.StringEncoder",partitionerClass="NULL",compressedTopics="NULL",
		queueBufferingMaxTime="NULL",
		queueBufferingMaxMessages="NULL",
		queueEnqueueTimeoutTime="NULL",
		batchNumMessages="NULL")
{
	producerType=as.character(producerType);
	compressionCodec=as.character(compressionCodec);
	serializerClass=as.character(serializerClass);
	partitionerClass=as.character(partitionerClass);
	queueBufferingMaxTime=as.character(queueBufferingMaxTime);
	queueBufferingMaxMessages=as.character(queueBufferingMaxMessages);
	queueEnqueueTimeoutTime=as.character(queueEnqueueTimeoutTime);
	batchNumMessages=as.character(batchNumMessages);
	
	# Create an object of the producer properties class
	producerProperties <- .jnew("com/musigma/producer/ProducerProperties")
	
	#set Properties from passed arguments and receive Properties Object
	producerPropertiesObj <- .jcall(producerProperties,"Ljava/util/Properties;","setProducerProperties",metadataBrokerList,producerType,compressionCodec,serializerClass,partitionerClass,compressedTopics,queueBufferingMaxTime,queueBufferingMaxMessages,queueEnqueueTimeoutTime,batchNumMessages)
	
	#Creating a producer
	producer<- .jnew("com/musigma/producer/MuProducer",producerPropertiesObj)
	return(producer)
}

#Function for producer to send message
#Definition of the arguments
#   /**
#     * @param producer:producer(Java object)
#   *            !!Mandatory: Producer through which messages are to be sent
#   * 
#   * @param topicName:String
#   *            !!Mandatory: Topic to which messages are to be sent. If topicName doesn't exist, new topic is created
#   * 
#     * @param ip:String
#   *            !!Mandatory: ip on which producer is running
#   * 
#     * @param message:String
#   *            !!Mandatory: message to be sent
#   */
rkafka.send <-function(producer, topicName, ip, message)
{
	topicName <- as.character(topicName)
	ip <- as.character(ip)
	message <- as.character(message)
	
	.jcall(producer,"V","sendMessage", topicName, ip, message)
	print("INFO:Remember to close the producer after done sending messages");
}

#function to shut down the producer
#Definition of the arguments
#   /**
#     * @param producer:producer(Java object)
#   *            !!Mandatory: Producer which is to be terminated
rkafka.closeProducer <-function(producer)
{
	.jcall(producer,"V","close")
}

#function to create consumer

#/**
#		* Definition of arguments
#* 
#@param zookeeperConnect
#*            !!Mandatory:Zookeeper connection string comma separated
#*            host:port pairs, each corresponding to a zk server. e.g.
#*            "127.0.0.1:3000,127.0.0.1:3001,127.0.0.1:3002"
#*			  default:"127.0.0.1:2181"
#* @param groupId
#*            !!Mandatory:consumer group id default:test-consumer-group
#* @param zookeeperConnectionTimeoutMs
#*            !!Mandatory:timeout in ms for connecting to zookeeper
#*            default:100000
#* @param consumerTimeoutMs
#*            !!Mandatory:Throw a timeout exception to the consumer if no
#*            message is available for consumption after the specified
#*            interval default:10000
#*@param topicName: 
#*            !!Mandatory: Name of the topic from where to read messages
#* @param autoCommitEnable
#*            --Optional:default:true If true, periodically commit to
#*            ZooKeeper the offset of messages already fetched by the
#*            consumer. This committed offset will be used when the process
#*            fails as the position from which the new consumer will begin.
#* @param autoCommitIntervalMs
#*            --Optional:default:60*1000 The frequency in ms that the
#*            consumer offsets are committed to zookeeper.
#* @param autoOffsetReset
#*            --Optional:default:largest 
#             * smallest : automatically reset the offset to the smallest offset 
#             largest : automatically: reset the offset to the largest offset 
#             anything else: throw exception to the consumer
#*/

rkafka.createConsumer<- function(zookeeperConnect,topicName,groupId="test-consumer-group",zookeeperConnectionTimeoutMs="100000",consumerTimeoutMs="10000",autoCommitEnable="NULL",autoCommitInterval="NULL",autoOffsetReset="NULL"){
	
	zookeeperConnect=as.character(zookeeperConnect)
	topicName=as.character(topicName)
	groupId=as.character(groupId)
	zookeeperConnectionTimeoutMs=as.character(zookeeperConnectionTimeoutMs)
	consumerTimeoutMs=as.character(consumerTimeoutMs) 
	autoCommitEnable=as.character(autoCommitEnable)
	autoCommitInterval=as.character(autoCommitInterval)
	autoOffsetReset=as.character(autoOffsetReset)

  
  Consumer<-.jnew("com/musigma/consumer/MuConsumer")
  print(Consumer)
  .jcall(Consumer,"V","CreateConsumer",zookeeperConnect,groupId,zookeeperConnectionTimeoutMs,consumerTimeoutMs,autoCommitEnable,autoCommitInterval,autoOffsetReset)
 success= .jcall(Consumer,"Z","startConsumer",topicName)
 
 if(success=="false"){
   stop("Consumer didn't start propertly")
   
 }

  return(Consumer)
}

#/**
#		* Reads messages from the topic passed as parameter while creating the consumer object one at a time.Waits for a time
#* specified by consumer timeout property to check for new messages, else returns ""
#		  @param ConsumerObj: Consumer through which messages are to be read(Java Object)	
#*		  @return String: next available message
#*/

rkafka.read<-function(ConsumerObj)
{
	message=.jcall(ConsumerObj,"Ljava/lang/String;","tail")
	print("INFO: Remember to close the consumer after done reading messages")
  return(message)
	
}
#/**
#  	* Reads messages from the topic passed as parameter while creating the consumer object in batch.Waits for a time
#* specified by consumer timeout property and returns all the messages it read.
#		  @param ConsumerObj: Consumer through which messages are to be read(Java Object)	
#*		  @return String[]: Array of Strings
#*/
rkafka.readPoll<-function(ConsumerObj){
  messages=.jcall(ConsumerObj,"[Ljava/lang/String;","poll")
  print("INFO: Remember to close the consumer after done reading messages")
  return(messages)
}

#function to shut down the consumer
rkafka.closeConsumer<-function(ConsumerObj){
	.jcall(ConsumerObj,"V","close")
}
#function to start Simple Consumer
# /**
#   * Creates a KAFKA simple consumer
# * @param kafkaServerURL : URL of the KAFKA server
# * @param kafkaServerPort :Port number of the KAFKA server
# * @param connectionTimeOut : Connection Timeout in ms
# * @param kafkaProducerBufferSize: Buffer size
# * @param clientId : ID Of the client
# */
rkafka.createSimpleConsumer<-function(kafkaServerURL, kafkaServerPort,connectionTimeOut, kafkaProducerBufferSize, clientId){
  MuSimpleConsumer<-.jnew("com/musigma/consumer/MuSimpleConsumer")
  
  kafkaServerURL<-as.character(kafkaServerURL)
  kafkaServerPort<-as.character(kafkaServerPort)
  connectionTimeOut<-as.character(connectionTimeOut)
  kafkaProducerBufferSize<-as.character(kafkaProducerBufferSize)
  clientId<-as.character(clientId)
  .jcall(MuSimpleConsumer,"V","CreateSimpleConsumer",kafkaServerURL,kafkaServerPort,connectionTimeOut,kafkaProducerBufferSize,clientId)
  
  return(MuSimpleConsumer)
}

# /**Function for consumer to read messages from topic
# * @param topicName: Name of the topic from where to read messages
# * @param partition: Partition Number
# * @param Offset :Offset Number
# * @param msgReadSize : Size of the message to be read**/
rkafka.receiveFromSimpleConsumer=function(SimpleConsumerObj,topicName,partition,Offset,msgReadSize){
  
  topicName<-as.character(topicName)
  partition<-as.character(partition)
  Offset<-as.character(Offset)
  msgReadSize<-as.character(msgReadSize)
  .jcall(SimpleConsumerObj,"V","receive",topicName,partition,Offset,msgReadSize)

}
#Function to return messages read by Simple Consumer one at a time
rkafka.readFromSimpleConsumer=function(SimpleConsumerObj){
  msg=.jcall(SimpleConsumerObj,"Ljava/lang/String;","read")
  return(msg)
}
# /**
#   * Close the simple consumer
# */
rkafka.closeSimpleConsumer=function(SimpleConsumer){
  .jcall(SimpleConsumer,"V","close")
}


