/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.solr.store.blob.process;

import java.lang.invoke.MethodHandles;

import org.apache.solr.core.CoreContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Utils related to blob background process e.g. init/shutdown of blob's background processes
 */
public class BlobProcessUtil {
  
  private static final Logger log = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private volatile CorePullerFeeder runningFeeder = null;
  
  /*
   * Start the Blob store async core pull machinery
   */
  public BlobProcessUtil(CoreContainer coreContainer) {
    this(coreContainer, new CorePullerFeeder(coreContainer));
  }
  
  /*
   * Start the Blob store async core pull machinery
   */
  public BlobProcessUtil(CoreContainer coreContainer, CorePullerFeeder cpf) {    
    runningFeeder = initializeCorePullerFeeder(cpf);
  }
  
  public CorePullerFeeder getCorePullerFeeder() {
    return runningFeeder;
  }
  
  /**
   * Shutdown background blob puller process
   */
  public void shutdown() {
    shutdownCorePullerFeeder();
  }
  
  /**
   * Initializes the CorePullerFeeder and starts running thread
   * @param cpf CorePullerFeeder
   * @return CorePullerFeeder 
   */
  private CorePullerFeeder initializeCorePullerFeeder(CorePullerFeeder cpf) {
    Thread t = new Thread(cpf);
    t.setName("blobPullerFeeder-" + t.getName());
    t.start();
    
    log.info("CorePullerFeeder initialized : " + t.getName());
    
    return cpf;
  }

  /**
   * Closes down the CorePullerFeeder thread
   */
  private void shutdownCorePullerFeeder() {
    final CoreSyncFeeder rf = runningFeeder;
    if (rf != null) {
      log.info("Shutting down CorePullerFeeder");
      rf.close();
    }
    runningFeeder = null;
  }
}
