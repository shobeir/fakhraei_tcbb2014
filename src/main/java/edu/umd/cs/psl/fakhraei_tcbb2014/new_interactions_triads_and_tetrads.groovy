/* Implemented by Shobeir Fakhraei for experimental evaluation in the following paper:
 *
 *     @article{fakharei2014network,
 *     title={Network-Based Drug-Target Interaction Prediction with Probabilistic Soft Logic},
 *     author={Fakhraei, Shobeir and Huang, Bert and Raschid, Louiqa and Getoor, Lise},
 *     journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
 *     year={2014},
 *     }
 *
 * This file contains a PSL model for drug-target interaction predictions.
 * Please refer to http://psl.umiacs.umd.edu/ for help on how to install and use PSL.
 * 
 * 
 * The file contains the experiments with new interations for traids with k=5 and tetrads with k=15.
 * 
 */

package edu.umd.cs.psl.fakhraei_tcbb2014;

import edu.umd.cs.psl.application.inference.MPEInference
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.MaxLikelihoodMPE
import edu.umd.cs.psl.config.*
import edu.umd.cs.psl.core.*
import edu.umd.cs.psl.core.inference.*
import edu.umd.cs.psl.database.DataStore
import edu.umd.cs.psl.database.Database
import edu.umd.cs.psl.database.DatabasePopulator
import edu.umd.cs.psl.database.Partition
import edu.umd.cs.psl.database.rdbms.RDBMSDataStore
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver.Type
import edu.umd.cs.psl.evaluation.result.*
import edu.umd.cs.psl.evaluation.statistics.RankingScore
import edu.umd.cs.psl.evaluation.statistics.SimpleRankingComparator
import edu.umd.cs.psl.groovy.*
import edu.umd.cs.psl.model.argument.ArgumentType
import edu.umd.cs.psl.model.argument.GroundTerm
import edu.umd.cs.psl.model.argument.Variable
import edu.umd.cs.psl.model.atom.QueryAtom
import edu.umd.cs.psl.model.parameters.PositiveWeight
import edu.umd.cs.psl.model.predicate.Predicate
import edu.umd.cs.psl.ui.loading.*
import edu.umd.cs.psl.util.database.*

// Setting the config file parameters
ConfigManager cm = ConfigManager.getManager();
ConfigBundle dtBundle = cm.getBundle("fakhraei_tcbb2014");

// Settings the experiments parameters
today = new Date();
double initialWeight = 5;
int folds = 10;
int numDrugs = 315;
int numTargets = 250;
boolean sq = true;
boolean doWeightLearning = true; //Set to false to avoid weight learning

// Setting the data path
// Change these to change the blocking parameter K
def base_dir = 'data'+java.io.File.separator;
def triad_similarity_dir = base_dir+'Similarities_K5_0.5'+java.io.File.separator;
def tetrad_similarity_dir = base_dir+'Similarities_K15_0.5'+java.io.File.separator;
def interactions_dir = base_dir+'NewInteractions'+java.io.File.separator;

System.out.println "Triad Similarity Folder: "+triad_similarity_dir;
System.out.println "Tetrad Similarity Folder: "+tetrad_similarity_dir;
System.out.println "Interactions Folds Folder: "+interactions_dir+" \n";

// Creating the variables to save the results of each fold
double[] AUC_Folds = new double[folds];
double[] AUPR_P_Folds = new double[folds];
double[] AUPR_N_Folds = new double[folds];


// Setting up the database
String dbpath = "./testdb"+today.getDate()+""+today.getHours()+""+today.getMinutes()+""+today.getSeconds();
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, true), dtBundle);

// Creating the PSL Model
// ======================

PSLModel m = new PSLModel(this,data)

// Defining Predicates
// Adding predicates for target-target similarities to be used in triads
m.add predicate: "targetSimilarity_dist",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "targetSimilarity_GO",		types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "targetSimilarity_seq",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
// Adding predicates for target-target similarities to be used in tetrads
// This is only necessary if one wants to use different similarites of tetrad (e.g. triads with k=15 and tetrads with k=5)
m.add predicate: "tetrad_targetSimilarity_dist",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_targetSimilarity_GO",		types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_targetSimilarity_seq",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// Adding predicates for drug-drug similarities to be used in triads
m.add predicate: "drugSimilarity_ATCHier",			types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "drugSimilarity_chemical",			types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "drugSimilarity_ligandJaccard",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "drugSimilarity_newCMapJaccard", 	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "drugSimilarity_SideEffect", 		types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
// Adding predicates for drug-drug similarities to be used in tetrads
// This is only necessary if one wants to use different similarites of tetrad (e.g. triads with k=15 and tetrads with k=5)
m.add predicate: "tetrad_drugSimilarity_ATCHier",			types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_drugSimilarity_chemical",			types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_drugSimilarity_ligandJaccard",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_drugSimilarity_newCMapJaccard", 	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate: "tetrad_drugSimilarity_SideEffect", 		types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

// Adding predicate for interactions
m.add predicate: "interacts_DT",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
// Adding predicate to hide the held out interactions in weight learning
m.add predicate: "ignored_interacts_DT",	types: [ArgumentType.UniqueID, ArgumentType.UniqueID] // For cross validation

// Defining Rules
// Setting up variables for rules to reset the weights in every fold only necessary when one needs to use the same rule in different folds
def triad_1, triad_2, triad_3, triad_4, triad_5, triad_6, triad_7, triad_8, prior;
def tetrad_1, tetrad_2, tetrad_3, tetrad_4, tetrad_5, tetrad_6, tetrad_7, tetrad_8, tetrad_9, tetrad_10, tetrad_11, tetrad_12, tetrad_13, tetrad_14, tetrad_15;

// Triad Rules
// Triads for drugs
triad_1 = m.add rule : (~ignored_interacts_DT(Drug1,Target) & ~ignored_interacts_DT(Drug2,Target) &  drugSimilarity_ATCHier(Drug1,Drug2) 		& interacts_DT(Drug1,Target) & (Drug1-Drug2) ) >> interacts_DT(Drug2,Target),  weight : initialWeight, squared: sq
triad_2 = m.add rule : (~ignored_interacts_DT(Drug1,Target) & ~ignored_interacts_DT(Drug2,Target) &  drugSimilarity_chemical(Drug1,Drug2) 		& interacts_DT(Drug1,Target) & (Drug1-Drug2) ) >> interacts_DT(Drug2,Target),  weight : initialWeight, squared: sq
triad_3 = m.add rule : (~ignored_interacts_DT(Drug1,Target) & ~ignored_interacts_DT(Drug2,Target) &  drugSimilarity_ligandJaccard(Drug1,Drug2) 	& interacts_DT(Drug1,Target) & (Drug1-Drug2) ) >> interacts_DT(Drug2,Target),  weight : initialWeight, squared: sq
triad_4 = m.add rule : (~ignored_interacts_DT(Drug1,Target) & ~ignored_interacts_DT(Drug2,Target) &  drugSimilarity_newCMapJaccard(Drug1,Drug2) 	& interacts_DT(Drug1,Target) & (Drug1-Drug2) ) >> interacts_DT(Drug2,Target),  weight : initialWeight, squared: sq
triad_5 = m.add rule : (~ignored_interacts_DT(Drug1,Target) & ~ignored_interacts_DT(Drug2,Target) &  drugSimilarity_SideEffect(Drug1,Drug2) 		& interacts_DT(Drug1,Target) & (Drug1-Drug2) ) >> interacts_DT(Drug2,Target),  weight : initialWeight, squared: sq

// Triads for targets
triad_6 = m.add rule : (~ignored_interacts_DT(Drug,Target1) & ~ignored_interacts_DT(Drug,Target2) &  targetSimilarity_dist(Target1,Target2) 		& interacts_DT(Drug,Target1) & (Target1-Target2) ) >> interacts_DT(Drug,Target2),  weight : initialWeight, squared: sq
triad_7 = m.add rule : (~ignored_interacts_DT(Drug,Target1) & ~ignored_interacts_DT(Drug,Target2) &  targetSimilarity_GO(Target1,Target2) 		& interacts_DT(Drug,Target1) & (Target1-Target2) ) >> interacts_DT(Drug,Target2),  weight : initialWeight, squared: sq
triad_8 = m.add rule : (~ignored_interacts_DT(Drug,Target1) & ~ignored_interacts_DT(Drug,Target2) &  targetSimilarity_seq(Target1,Target2) 		& interacts_DT(Drug,Target1) & (Target1-Target2) ) >> interacts_DT(Drug,Target2),  weight : initialWeight, squared: sq

// Tetrads Rules
tetrad_1 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) 	& ~ignored_interacts_DT(Drug2,Target2) & tetrad_drugSimilarity_ATCHier(Drug1,Drug2)	& tetrad_targetSimilarity_dist(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_2 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) & tetrad_drugSimilarity_ATCHier(Drug1,Drug2) 	& tetrad_targetSimilarity_GO(Target1,Target2)		& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_3 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) & tetrad_drugSimilarity_ATCHier(Drug1,Drug2) 	& tetrad_targetSimilarity_seq(Target1,Target2)		& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq

tetrad_4 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_chemical(Drug1,Drug2)	& tetrad_targetSimilarity_dist(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_5 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_chemical(Drug1,Drug2)	& tetrad_targetSimilarity_GO(Target1,Target2)		& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_6 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_chemical(Drug1,Drug2)	& tetrad_targetSimilarity_seq(Target1,Target2)		& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq

tetrad_7 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_ligandJaccard(Drug1,Drug2) & tetrad_targetSimilarity_dist(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_8 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_ligandJaccard(Drug1,Drug2) & tetrad_targetSimilarity_GO(Target1,Target2)		& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_9 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_ligandJaccard(Drug1,Drug2) & tetrad_targetSimilarity_seq(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq

tetrad_10 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_newCMapJaccard(Drug1,Drug2) & tetrad_targetSimilarity_dist(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_11 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_newCMapJaccard(Drug1,Drug2) & tetrad_targetSimilarity_GO(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_12 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_newCMapJaccard(Drug1,Drug2) & tetrad_targetSimilarity_seq(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq

tetrad_13 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_SideEffect(Drug1,Drug2) & tetrad_targetSimilarity_dist(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_14 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_SideEffect(Drug1,Drug2) & tetrad_targetSimilarity_GO(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq
tetrad_15 = m.add rule : ( ~interacts_DT(Drug1,Target2) & ~interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target2) & ~ignored_interacts_DT(Drug2,Target1) & ~ignored_interacts_DT(Drug1,Target1) & ~ignored_interacts_DT(Drug2,Target2) &  tetrad_drugSimilarity_SideEffect(Drug1,Drug2) & tetrad_targetSimilarity_seq(Target1,Target2)	& interacts_DT(Drug1,Target1) & (Drug1-Drug2) & (Target1-Target2) ) >> interacts_DT(Drug2,Target2),  weight : initialWeight, squared: sq

// Prior
prior = m.add rule : ~ignored_interacts_DT(Drug,Target) >> ~interacts_DT(Drug,Target),  weight : initialWeight, squared: sq


// Printing the model
System.out.println m;


// Loading the data
// ================

// Creating the partition to read the data
Partition readSimilarities =  new Partition(1);
Partition write = new Partition(2);

// Reading triad target similarities from file
for (Predicate p : [targetSimilarity_dist,targetSimilarity_GO,targetSimilarity_seq])
{
	System.out.println "Reading " + p.getName();
	insert = data.getInserter(p, readSimilarities)
	InserterUtils.loadDelimitedDataTruth(insert, triad_similarity_dir+"PSL_"+p.getName()+".txt")
}

// Reading tetrad target similarities from file
for (Predicate p : [tetrad_targetSimilarity_dist,tetrad_targetSimilarity_GO,tetrad_targetSimilarity_seq])
{
	System.out.println "READING " + p.getName();
	insert = data.getInserter(p, readSimilarities)
	InserterUtils.loadDelimitedDataTruth(insert, tetrad_similarity_dir+"PSL_"+p.getName().substring(7,p.getName().length())+".txt")
}


// Reading triad drug similarities from file
for (Predicate p : [drugSimilarity_ATCHier, drugSimilarity_chemical, drugSimilarity_ligandJaccard, drugSimilarity_newCMapJaccard, drugSimilarity_SideEffect])
{
	System.out.println "Reading " + p.getName();
	insert = data.getInserter(p, readSimilarities)
	InserterUtils.loadDelimitedDataTruth(insert, triad_similarity_dir+"PSL_"+p.getName()+".txt")
}

// Reading tetrad drug similarities from file
for (Predicate p : [tetrad_drugSimilarity_ATCHier, tetrad_drugSimilarity_chemical, tetrad_drugSimilarity_ligandJaccard, tetrad_drugSimilarity_newCMapJaccard, tetrad_drugSimilarity_SideEffect])
{
	System.out.println "READING " + p.getName();
	insert = data.getInserter(p, readSimilarities)
	InserterUtils.loadDelimitedDataTruth(insert, tetrad_similarity_dir+"PSL_"+p.getName().substring(7,p.getName().length())+".txt")
}

// Setting which predicates are closed
Set <Predicate>closedPredicates = [targetSimilarity_dist, targetSimilarity_GO, targetSimilarity_seq, drugSimilarity_ATCHier, drugSimilarity_chemical, drugSimilarity_ligandJaccard, drugSimilarity_newCMapJaccard, drugSimilarity_SideEffect, ignored_interacts_DT];
Set <Predicate>closedPredicatesLinks=[interacts_DT];
Set <Predicate>closedPredicatesAll=[interacts_DT];
closedPredicatesAll.addAll(closedPredicates);

//Making the sets to use in database population
Variable Drug = new Variable("Drug");
Set<GroundTerm> drugGroundings = new HashSet<GroundTerm>();
for (int i = 1; i <= numDrugs; i++)
{
	drugGroundings.add(data.getUniqueID(i));
}

Variable Target = new Variable("Target");
Set<GroundTerm> targetGroundings = new HashSet<GroundTerm>();
for (int i = 1; i <= numTargets; i++)
{
	targetGroundings.add(data.getUniqueID(i));
}

Map<Variable, Set<GroundTerm>> substitutions = new HashMap<Variable, Set<GroundTerm>>();
substitutions.put(Drug, drugGroundings);
substitutions.put(Target, targetGroundings);


// Cross Validation
// ================

for (int k=0;k<folds;k++)
{
	System.out.println("\n-------------------");
	TimeNow = new Date();
	System.out.println("Fold "+(k+1).toString()+" Start: "+TimeNow);

	// Setting up the partitions for final model
	Partition writeCVLabels =  new Partition(666+(k*10)); // Labels for the Cross-validation hold-outs
	Partition writeCVTest =  new Partition(667+(k*10)); // Partition to write the predictions in
	Partition readCVTrain =  new Partition(668+(k*10)); // Observed training data for the training (i.e., all data minus hold-outs)

	//Setting up the the partitions for weight learning
	Partition writeWLLabels =  new Partition(669+(k*10)); // Labels for the weight Learning
	Partition readWLTrain =  new Partition(671+(k*10)); // Training data for Weight Learning
	Partition writeWLTest =  new Partition(672+(k*10)); // Partition to write prediction in for Weight Learning

	// Setting up the inserters
	insertWLTrain = data.getInserter(interacts_DT, readWLTrain);
	insertWLLabels = data.getInserter(interacts_DT, writeWLLabels);
	insertCVIgnored = data.getInserter(ignored_interacts_DT, readWLTrain);

	insertCVTrain = data.getInserter(interacts_DT, readCVTrain);
	insertCVLabels = data.getInserter(interacts_DT, writeCVLabels);

	// Reading the interactions and setting the data for current fold

	CVHoldoutFold = k+1; // current hold out fold
	WLHoldoutFold = k+2; // current hold out fold for weight learning
	if (k==9) WLHoldoutFold = 1;

	System.out.print "\nReading INTERACTS_DT files for fold "+(k+1).toString()+" ";

	// Reading all the other folds as training data
	for (int j=1;j<=folds;j++)
	{
		System.out.print ".";
		if ((j!=CVHoldoutFold) && (j!=WLHoldoutFold))
		{
			InserterUtils.loadDelimitedDataTruth(insertCVTrain, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+j+'.txt');
			InserterUtils.loadDelimitedDataTruth(insertWLTrain, interactions_dir+'PSL_NewInteractions_Positive_Fold_'+j+'.txt');
			InserterUtils.loadDelimitedDataTruth(insertWLTrain, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+j+'.txt');
		}
	}

	System.out.print ".";
	// Adding all the positive interactions as training data for the final model
	InserterUtils.loadDelimitedDataTruth(insertCVTrain, interactions_dir+'INTERACTS_DT_T.txt');
	
	// Adding the weight learning held-out to the final training data
	InserterUtils.loadDelimitedDataTruth(insertCVTrain, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+WLHoldoutFold+'.txt');

	// Reading the weight learning held-out as the the labels
	// The positive (known) interaction are only held out for weight learning not for the final model
	// "PSL_NewInteractions_Positive_Fold_*" contains the observed positive interactions ("INTERACTS_DT_T.txt") divided into 10 folds (samples).
	InserterUtils.loadDelimitedDataTruth(insertWLLabels, interactions_dir+'PSL_NewInteractions_Positive_Fold_'+WLHoldoutFold+'.txt');
	InserterUtils.loadDelimitedDataTruth(insertWLLabels, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+WLHoldoutFold+'.txt');

	System.out.print ".";
	// Reading the held out as labels for weight learning and also into the ignored_interacts_DT.
	// Weight learning will not be able to see these labels because of they will be ignored and the body of their rules will be 0.
	InserterUtils.loadDelimitedDataTruth(insertWLLabels, interactions_dir+'PSL_NewInteractions.txt');
	InserterUtils.loadDelimitedDataTruth(insertWLLabels, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+CVHoldoutFold+'.txt');
	InserterUtils.loadDelimitedDataTruth(insertCVIgnored, interactions_dir+'PSL_NewInteractions.txt');
	InserterUtils.loadDelimitedDataTruth(insertCVIgnored, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+CVHoldoutFold+'.txt');

	// Reading the held out as labels for the final model
	// All the new interaction and a sample of Negatives are held for the final model
	InserterUtils.loadDelimitedDataTruth(insertCVLabels, interactions_dir+'PSL_NewInteractions.txt');
	InserterUtils.loadDelimitedDataTruth(insertCVLabels, interactions_dir+'PSL_NewInteractions_Negetive_Fold_'+CVHoldoutFold+'.txt');

	// Populating the dataset; the unobserved interactions should be given to the model so it knows to what links to provide predictions for.
	Database dbPopTemp1 = data.getDatabase(writeWLTest, closedPredicates, readWLTrain);
	DatabasePopulator dbPop0 = new DatabasePopulator(dbPopTemp1);
	dbPop0.populate(new QueryAtom(interacts_DT, Drug, Target), substitutions);
	dbPopTemp1.close();

	Database dbPopTemp2 = data.getDatabase(writeCVTest, closedPredicates ,readCVTrain);
	DatabasePopulator dbPop2 = new DatabasePopulator(dbPopTemp2);
	dbPop2.populate(new QueryAtom(interacts_DT, Drug, Target), substitutions);
	dbPopTemp2.close();

	System.out.println "\n";

	// Weight Learning
	// ===============

	if (doWeightLearning)
	{
		System.out.println "Weight Learning fold: "+(k+1).toString();

		//Reseting all the weights in each fold
		triad_1.setWeight(new PositiveWeight(initialWeight));
		triad_2.setWeight(new PositiveWeight(initialWeight));
		triad_3.setWeight(new PositiveWeight(initialWeight));
		triad_4.setWeight(new PositiveWeight(initialWeight));
		triad_5.setWeight(new PositiveWeight(initialWeight));
		triad_6.setWeight(new PositiveWeight(initialWeight));
		triad_7.setWeight(new PositiveWeight(initialWeight));
		triad_8.setWeight(new PositiveWeight(initialWeight));

		tetrad_1.setWeight(new PositiveWeight(initialWeight));
		tetrad_2.setWeight(new PositiveWeight(initialWeight));
		tetrad_3.setWeight(new PositiveWeight(initialWeight));
		tetrad_4.setWeight(new PositiveWeight(initialWeight));
		tetrad_5.setWeight(new PositiveWeight(initialWeight));
		tetrad_6.setWeight(new PositiveWeight(initialWeight));
		tetrad_7.setWeight(new PositiveWeight(initialWeight));
		tetrad_8.setWeight(new PositiveWeight(initialWeight));
		tetrad_9.setWeight(new PositiveWeight(initialWeight));
		tetrad_10.setWeight(new PositiveWeight(initialWeight));
		tetrad_11.setWeight(new PositiveWeight(initialWeight));
		tetrad_12.setWeight(new PositiveWeight(initialWeight));
		tetrad_13.setWeight(new PositiveWeight(initialWeight));
		tetrad_14.setWeight(new PositiveWeight(initialWeight));
		tetrad_15.setWeight(new PositiveWeight(initialWeight));

		prior.setWeight(new PositiveWeight(initialWeight));

		// the actual weight learning happens here
		Database dbWLTrain = data.getDatabase(writeWLTest, closedPredicates, readWLTrain, readSimilarities);
		Database dbWLLabels = data.getDatabase(writeWLLabels, closedPredicatesAll);

		MaxLikelihoodMPE wLearn = new MaxLikelihoodMPE(m,dbWLTrain,dbWLLabels,dtBundle);
		wLearn.learn();

		dbWLTrain.close();
		dbWLLabels.close();

		// Printing the new weights
		System.out.println m;
	};

	// Inferring
	// =========

	System.out.println "Inferring fold: "+(k+1).toString();

	Database dbCVTrain = data.getDatabase(writeCVTest, closedPredicates ,readCVTrain, readSimilarities);

	MPEInference mpe = new MPEInference(m, dbCVTrain, dtBundle);
	FullInferenceResult result = mpe.mpeInference();
	mpe.close();
	mpe.finalize();

	dbCVTrain.close();

	timeNow = new Date();
	System.out.println("\nFold "+(k+1).toString()+" End: "+timeNow);
	System.out.println("-------------------\n");

	//Begin Evaluate
	//==============

	System.out.println "Evaluating fold: "+(k+1).toString();

	def LabelsDB = data.getDatabase(writeCVLabels, closedPredicatesAll)
	Database PredictionsDB = data.getDatabase(new Partition(6000+k), writeCVTest)
	def comparator = new SimpleRankingComparator(PredictionsDB)
	comparator.setBaseline(LabelsDB)

	// Choosing what metrics to report
	def metrics = [RankingScore.AUPRC, RankingScore.NegAUPRC,  RankingScore.AreaROC]
	double [] score = new double[metrics.size()]

	try {
		for (int i = 0; i < metrics.size(); i++) {
			comparator.setRankingScore(metrics.get(i))
			score[i] = comparator.compare(interacts_DT)
		}
		//Storing the performance values of the current fold
		AUPR_P_Folds[k]=score[0];
		AUPR_N_Folds[k]=score[1];
		AUC_Folds[k]=score[2];

		System.out.println("\nArea under positive-class PR curve: " + score[0])
		System.out.println("Area under negetive-class PR curve: " + score[1])
		System.out.println("Area under ROC curve: " + score[2])
	}
	catch (ArrayIndexOutOfBoundsException e) {
		System.out.println("No evaluation data! Terminating!");
	}

	PredictionsDB.close();
	LabelsDB.close();

	// Cleaning up
	data.deletePartition(writeCVLabels);
	data.deletePartition(writeCVTest);
	data.deletePartition(readCVTrain);
	data.deletePartition(writeWLLabels);
	data.deletePartition(readWLTrain);
	data.deletePartition(writeWLTest);

}
// End of cross validation

System.out.println("\n========================");

// Generating Performance Mean and STDev

double FinalAUC=0,FinalAUPR_N=0,FinalAUPR_P=0;
double STDAUC=0, STDAUPR_N=0, STDAUPR_P=0;

for (int i=0; i<folds; i++)
{
	FinalAUC+=AUC_Folds[i];
	FinalAUPR_N+=AUPR_N_Folds[i];
	FinalAUPR_P+=AUPR_P_Folds[i];
}

FinalAUC=FinalAUC/folds;
FinalAUPR_N=FinalAUPR_N/folds;
FinalAUPR_P=FinalAUPR_P/folds;

for (int i=0; i<folds; i++)
{
	STDAUC=STDAUC + Math.pow(FinalAUC-AUC_Folds[i],2);
	STDAUPR_N=STDAUPR_N + Math.pow(FinalAUPR_N-AUPR_N_Folds[i],2);
	STDAUPR_P=STDAUPR_P + Math.pow(FinalAUPR_P-AUPR_P_Folds[i],2);
}

STDAUC=Math.sqrt(STDAUC/folds);
STDAUPR_N=Math.sqrt(STDAUPR_N/folds);
STDAUPR_P=Math.sqrt(STDAUPR_P/folds);

System.out.println("\nFinal Area under positive-class PR curve: " + FinalAUPR_P +" +/- "+ STDAUPR_P)
System.out.println("Final Area under negative-class PR curve: " + FinalAUPR_N +" +/- "+ STDAUPR_N)
System.out.println("Final Area under ROC curve: " + FinalAUC +" +/- "+ STDAUC)
System.out.println(" ");
System.out.println "========================";
System.out.println "Done! ";
