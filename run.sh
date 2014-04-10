#!/bin/bash
# Use this script to run the experiments.

NOW=$(date +"%F")
NOWT=$(date +"%H-%M-%S")
echo --- Output will be saved "in" PSLOutput.$NOW.$NOWT.txt ---

mvn compile

mvn dependency:build-classpath -Dmdep.outputFile=classpath.out

# Change the Java virtual machine memory based on your machine
# Cross validation experiemtns with triads
java -Xms4G -Xmx8G -cp ./target/classes:`cat classpath.out` edu.umd.cs.psl.fakhraei_tcbb2014.cross_validation_triads &> PSLOutput.$NOW.$NOWT.txt

# Cross validation experiemtns with triads and tetrads
#java -Xms4G -Xmx8G -cp ./target/classes:`cat classpath.out` edu.umd.cs.psl.fakhraei_tcbb2014.cross_validation_triads_and_tetrads &> PSLOutput.$NOW.$NOWT.txt

# New interaction predictions experiemtns with triads
#java -Xms4G -Xmx8G -cp ./target/classes:`cat classpath.out` edu.umd.cs.psl.fakhraei_tcbb2014.new_interactions_triads &> PSLOutput.$NOW.$NOWT.txt

# New interaction predictions experiemtns with triads and tetrads
#java -Xms4G -Xmx8G -cp ./target/classes:`cat classpath.out` edu.umd.cs.psl.fakhraei_tcbb2014.new_interactions_triads_and_tetrads &> PSLOutput.$NOW.$NOWT.txt

