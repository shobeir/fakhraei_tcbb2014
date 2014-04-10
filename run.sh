#!/bin/bash
NOW=$(date +"%F")
NOWT=$(date +"%H-%M-%S")

echo --- Output will be saved "in" PSLOutput.$NOW.$NOWT.txt ---

mvn compile

mvn dependency:build-classpath -Dmdep.outputFile=classpath.out

java -Xms4G -Xmx8G -cp ./target/classes:`cat classpath.out` edu.umd.cs.psl.fakhraei_tcbb2014.cross_validation_triads &> PSLOutput.$NOW.$NOWT.txt
