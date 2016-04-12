/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:        DataMethod
//- Description:  Class implementation
//- Owner:        Mike Eldred

#include "DataMethod.hpp"
#include "dakota_data_io.hpp"
#include "pecos_global_defs.hpp"
#ifdef HAVE_OPTPP
#include "globals.h"
#endif


namespace Dakota {

DataMethodRep::DataMethodRep():
  methodName(DEFAULT_METHOD), subMethod(SUBMETHOD_DEFAULT),
  methodOutput(NORMAL_OUTPUT), maxIterations(-1), maxFunctionEvaluations(1000),
  speculativeFlag(false), methodUseDerivsFlag(false),
  convergenceTolerance(1.e-4), constraintTolerance(0.), methodScaling(false),
  numFinalSolutions(0),
  // Meta-iterators
  iteratorServers(0), procsPerIterator(0), // 0 defaults to detect user spec
  iteratorScheduling(DEFAULT_SCHEDULING), hybridLSProb(0.1),
  //hybridProgThresh(0.5),
  concurrentRandomJobs(0),
  // Local surrogate-based opt/NLS
  softConvLimit(0), // dummy value -> method-specific default
  surrBasedLocalLayerBypass(false),      surrBasedLocalTRInitSize(0.4),
  surrBasedLocalTRMinSize(1.0e-6),       surrBasedLocalTRContractTrigger(0.25),
  surrBasedLocalTRExpandTrigger(0.75),   surrBasedLocalTRContract(0.25),
  surrBasedLocalTRExpand(2.0), surrBasedLocalSubProbObj(ORIGINAL_PRIMARY),
  surrBasedLocalSubProbCon(ORIGINAL_CONSTRAINTS),
  surrBasedLocalMeritFn(AUGMENTED_LAGRANGIAN_MERIT),
  surrBasedLocalAcceptLogic(FILTER),     surrBasedLocalConstrRelax(NO_RELAX),
  // Global surrogate-based opt/NLS
  surrBasedGlobalReplacePts(false),
  // Branch and bound
  //branchBndNumSamplesRoot(0), branchBndNumSamplesNode(0),
  // DL_SOLVER
  dlLib(NULL), //dlDetails(""),
  // NPSOL
  verifyLevel(-1), functionPrecision(1.e-10), lineSearchTolerance(0.9),
  // NL2SOL: Real values of -1. ==> use NL2SOL default
  absConvTol(-1.), xConvTol(-1.), singConvTol(-1.), singRadius(-1.),
  falseConvTol(-1.), initTRRadius(-1.), covarianceType(0), regressDiag(false),
  // OPT++
  // searchMethod default is null since "trust_region" is preferred for 
  // unconstrained opt., whereas "line_search" is preferred for bc opt.
  gradientTolerance(0.0001), maxStep(1.e+3), 
#if HAVE_OPTPP
  meritFn(OPTPP::ArgaezTapia),
#else
  meritFn(NULL),
#endif
  stepLenToBoundary(-1.), centeringParam(-1.), // dummy defaults (see SNLLBase)
  searchSchemeSize(32),
  // APPSPACK
  initStepLength(1.0), contractStepLength(0.5), threshStepLength(0.01),
  meritFunction("merit2_squared"), constrPenalty(1.0), smoothFactor(0.0),
  // COLINY
  constantPenalty(false), globalBalanceParam(-1.),
  localBalanceParam(-1.), maxBoxSize(-1.), minBoxSize(-1.),
  //boxDivision("major_dimension"), // leave empty string as default
  showMiscOptions(false), mutationAdaptive(true),
  // These attributes must replicate the Coliny defaults due to Coliny 
  // member fn. structure:
  mutationRate(1.0),
  mutationScale(0.1), contractFactor(0.5),
  // These attributes replicate COLINY defaults due to convenience:
  totalPatternSize(0), // simplifies maxEvalConcurrency calculation
  // These attributes define DAKOTA defaults which may differ from Coliny
  // defaults:
  solnTarget(-DBL_MAX), // COLINY default of 1.e-5 can cause premature term.
  // These attributes use dummy defaults which are used to trigger conditional
  // option processing (since we don't want to get out of synch with SGOPT's 
  // defaults which may change).  The dummy defaults should be values which 
  // are _not_ reasonable user inputs.
  mutationMinScale(-1.), initDelta(-1.), threshDelta(-1.),
  newSolnsGenerated(-9999), numberRetained(-9999),
  expansionFlag(true), // default = on, no_expansion spec turns off
  expandAfterSuccess(0), contractAfterFail(0), mutationRange(-9999),
  randomizeOrderFlag(false), //betaSolverName(""),
  // JEGA
  numCrossPoints(2), numParents(2), numOffspring(2), //convergenceType(""),
  fitnessLimit(6.0), shrinkagePercent(0.9), percentChange(0.1),
  numGenerations(15), nichingType("null_niching"), numDesigns(100),
  postProcessorType("null_postprocessor"), logFile("JEGAGlobal.log"),
  printPopFlag(false),
  // JEGA/COLINY
  constraintPenalty(-1.), crossoverRate(-1.), //crossoverType(""),
  initializationType("unique_random"),
  //mutationType(""), replacementType(""), fitnessType(""),
  populationSize(50), //flatFile(),
  // NOMAD
  historyFile("mads_history"), displayFormat("bbe obj"), 
  vns(0.0), neighborOrder(1), showAllEval(false),
  // NCSU 
  volBoxSize(-1.),
  // DDACE
  numSymbols(0),mainEffectsFlag(false),
  // FSUDace
  numTrials(10000), latinizeFlag(false), volQualityFlag(false),
  fixedSequenceFlag(false), //default is variable sampling patterns
  //initializationType("grid"), trialType("random"),
  // COLINY, JEGA, NonD, & DACE
  randomSeed(0),
  // NonD & DACE
  numSamples(0), fixedSeedFlag(false), previousSamples(0), vbdFlag(false),
  vbdDropTolerance(-1.),backfillFlag(false), pcaFlag(false),
  percentVarianceExplained(0.95),
  // NonD
  vbdOrder(0), covarianceControl(DEFAULT_COVARIANCE), rngName("mt19937"),
  refinementType(Pecos::NO_REFINEMENT), refinementControl(Pecos::NO_CONTROL),
  nestingOverride(Pecos::NO_NESTING_OVERRIDE),
  growthOverride(Pecos::NO_GROWTH_OVERRIDE), expansionType(EXTENDED_U),
  piecewiseBasis(false), expansionBasisType(Pecos::DEFAULT_BASIS),
  cubIntOrder(USHRT_MAX), collocationRatio(0.), collocRatioTermsOrder(1.),
  regressionType(Pecos::DEFAULT_REGRESSION), lsRegressionType(DEFAULT_LS),
  regressionL2Penalty(0.), crossValidation(false), //adaptedBasisInitLevel(0),
  adaptedBasisAdvancements(3), normalizedCoeffs(false), tensorGridFlag(false),
  //expansionSampleType("lhs"),
  sampleType(SUBMETHOD_DEFAULT), reliabilitySearchType(MV),
  integrationRefine(NO_INT_REFINE), refineSamples(0),
  distributionType(CUMULATIVE), responseLevelTarget(PROBABILITIES),
  responseLevelTargetReduce(COMPONENT), emulatorSamples(0), emulatorOrder(0),
  emulatorType(NO_EMULATOR), mcmcType("dram"), standardizedSpace(false),
  adaptPosteriorRefine(false), logitTransform(false),
  preSolveMethod(SUBMETHOD_DEFAULT), proposalCovUpdates(0),
  fitnessMetricType("predicted_variance"), batchSelectionType("naive"),
  batchSize(0), calibrateErrorMode(CALIBRATE_NONE), numChains(3), numCR(3),
  crossoverChainPairs(3), grThreshold(1.2), jumpStep(5), lipschitzType("local"),
  // Parameter Study
  numSteps(0), pstudyFileFormat(TABULAR_ANNOTATED), pstudyFileActive(false), 
  // Verification
  refinementRate(2.),
  // Point import/export files
  importBuildFormat(TABULAR_ANNOTATED),  importBuildActive(false),
  importApproxFormat(TABULAR_ANNOTATED), importApproxActive(false),
  exportApproxFormat(TABULAR_ANNOTATED), exportMCMCFormat(TABULAR_ANNOTATED),
  referenceCount(1)
{ }


void DataMethodRep::write(MPIPackBuffer& s) const
{
  s << idMethod << modelPointer << methodOutput << maxIterations
    << maxFunctionEvaluations << speculativeFlag << methodUseDerivsFlag
    << convergenceTolerance << constraintTolerance << methodScaling
    << numFinalSolutions << linearIneqConstraintCoeffs << linearIneqLowerBnds
    << linearIneqUpperBnds << linearIneqScaleTypes << linearIneqScales
    << linearEqConstraintCoeffs << linearEqTargets << linearEqScaleTypes
    << linearEqScales << methodName << subMethod << subMethodName
    << subModelPointer << subMethodPointer;

  // Meta-iterators
  s << iteratorServers << procsPerIterator << iteratorScheduling
    << hybridMethodNames << hybridModelPointers << hybridMethodPointers
  //<< hybridProgThresh
    << hybridGlobalMethodName << hybridGlobalModelPointer
    << hybridGlobalMethodPointer << hybridLocalMethodName
    << hybridLocalModelPointer << hybridLocalMethodPointer << hybridLSProb
  //<< branchBndNumSamplesRoot << branchBndNumSamplesNode
    << concurrentRandomJobs << concurrentParameterSets;

  // Surrogate-based
  s << softConvLimit << surrBasedLocalLayerBypass
    << surrBasedLocalTRInitSize << surrBasedLocalTRMinSize
    << surrBasedLocalTRContractTrigger << surrBasedLocalTRExpandTrigger
    << surrBasedLocalTRContract << surrBasedLocalTRExpand
    << surrBasedLocalSubProbObj << surrBasedLocalSubProbCon
    << surrBasedLocalMeritFn << surrBasedLocalAcceptLogic
    << surrBasedLocalConstrRelax << surrBasedGlobalReplacePts;

  // DL_SOLVER
  s << dlDetails;

  // NPSOL
  s << verifyLevel << functionPrecision << lineSearchTolerance;

  // NL2SOL
  s << absConvTol << xConvTol << singConvTol << singRadius << falseConvTol
    << initTRRadius << covarianceType << regressDiag;

  // OPT++
  s << searchMethod << gradientTolerance << maxStep << meritFn
    << stepLenToBoundary << centeringParam << searchSchemeSize;

  // APPSPACK
  s << initStepLength << contractStepLength << threshStepLength << meritFunction
    << constrPenalty << smoothFactor;

  // COLINY
  s << constraintPenalty << constantPenalty << globalBalanceParam
    << localBalanceParam << maxBoxSize << minBoxSize << boxDivision
    << mutationAdaptive << showMiscOptions << miscOptions << solnTarget
    << crossoverRate << mutationRate << mutationScale << mutationMinScale
    << initDelta << threshDelta << contractFactor << newSolnsGenerated
    << numberRetained << expansionFlag << expandAfterSuccess
    << contractAfterFail << mutationRange << totalPatternSize
    << randomizeOrderFlag << selectionPressure << replacementType
    << crossoverType << mutationType << exploratoryMoves << patternBasis;

  // COLINY + APPSPACK
  s << evalSynchronize;

  // JEGA
  s << numCrossPoints << numParents << numOffspring << fitnessType
    << convergenceType << percentChange << numGenerations << fitnessLimit
    << shrinkagePercent << nichingType << nicheVector << numDesigns
    << postProcessorType << distanceVector;

  // JEGA/COLINY
  s << initializationType << flatFile << logFile << populationSize
    << printPopFlag;

  // NCSU 
  s << volBoxSize;

  // DDACE
  s << numSymbols << mainEffectsFlag;

  // FSUDace 
  s << latinizeFlag << volQualityFlag << sequenceStart << sequenceLeap
    << primeBase << numTrials << trialType;

  // COLINY, NonD, DACE, & JEGA
  s << randomSeed;

  // MADS
  s << historyFile << displayFormat << vns << neighborOrder << showAllEval;

  // NonD & DACE
  s << numSamples << fixedSeedFlag << fixedSequenceFlag << previousSamples
    << vbdFlag << vbdDropTolerance << backfillFlag << pcaFlag
    << percentVarianceExplained;

  // NonD
  s << vbdOrder << covarianceControl << rngName << refinementType
    << refinementControl << nestingOverride << growthOverride << expansionType
    << piecewiseBasis << expansionBasisType << expansionOrder
    << expansionSamples << expansionSampleType << quadratureOrder
    << sparseGridLevel << anisoDimPref << cubIntOrder << collocationPoints
    << collocationRatio << collocRatioTermsOrder << regressionType
    << lsRegressionType << regressionNoiseTol << regressionL2Penalty
    << crossValidation //<< adaptedBasisInitLevel
    << adaptedBasisAdvancements << normalizedCoeffs << pointReuse
    << tensorGridFlag << tensorGridOrder << importExpansionFile
    << exportExpansionFile << sampleType << reliabilitySearchType
    << reliabilityIntegration << integrationRefine << refineSamples
    << distributionType << responseLevelTarget << responseLevelTargetReduce
    << responseLevels << probabilityLevels << reliabilityLevels
    << genReliabilityLevels << emulatorSamples << emulatorOrder << emulatorType
    << mcmcType << standardizedSpace << adaptPosteriorRefine << logitTransform
    << preSolveMethod << proposalCovType << proposalCovUpdates
    << proposalCovInputType << proposalCovData << proposalCovFile
    << fitnessMetricType << batchSelectionType << batchSize
    << calibrateErrorMode << hyperPriorAlphas << hyperPriorBetas
    << numChains << numCR << crossoverChainPairs
    << grThreshold << jumpStep << lipschitzType << dataDistType 
    << dataDistCovInputType << dataDistMeans << dataDistCovariance
    << dataDistFile << posteriorDensityExportFilename
    << posteriorSamplesExportFilename << posteriorSamplesImportFilename
    << generatePosteriorSamples << evaluatePosteriorDensity;
    ;

  // Parameter Study
  s << finalPoint << stepVector << numSteps << stepsPerVariable << listOfPoints
    << pstudyFilename << pstudyFileFormat << pstudyFileActive
    << varPartitions;

  // Verification
  s << refinementRate;
 
  // Point import/export files
  s << importBuildPtsFile  << importBuildFormat  << importBuildActive
    << importApproxPtsFile << importApproxFormat << importApproxActive
    << exportApproxPtsFile << exportApproxFormat << exportMCMCPtsFile
    << exportMCMCFormat;
}


void DataMethodRep::read(MPIUnpackBuffer& s)
{
  s >> idMethod >> modelPointer >> methodOutput >> maxIterations
    >> maxFunctionEvaluations >> speculativeFlag >> methodUseDerivsFlag
    >> convergenceTolerance >> constraintTolerance >> methodScaling
    >> numFinalSolutions >> linearIneqConstraintCoeffs >> linearIneqLowerBnds
    >> linearIneqUpperBnds >> linearIneqScaleTypes >> linearIneqScales
    >> linearEqConstraintCoeffs >> linearEqTargets >> linearEqScaleTypes
    >> linearEqScales >> methodName >> subMethod >> subMethodName
    >> subModelPointer >> subMethodPointer;

  // Meta-iterators
  s >> iteratorServers >> procsPerIterator >> iteratorScheduling
    >> hybridMethodNames >> hybridModelPointers >> hybridMethodPointers
  //>> hybridProgThresh
    >> hybridGlobalMethodName >> hybridGlobalModelPointer
    >> hybridGlobalMethodPointer >> hybridLocalMethodName
    >> hybridLocalModelPointer >> hybridLocalMethodPointer >> hybridLSProb
  //>> branchBndNumSamplesRoot >> branchBndNumSamplesNode
    >> concurrentRandomJobs >> concurrentParameterSets;

  // Surrogate-based
  s >> softConvLimit >> surrBasedLocalLayerBypass
    >> surrBasedLocalTRInitSize >> surrBasedLocalTRMinSize
    >> surrBasedLocalTRContractTrigger >> surrBasedLocalTRExpandTrigger
    >> surrBasedLocalTRContract >> surrBasedLocalTRExpand
    >> surrBasedLocalSubProbObj >> surrBasedLocalSubProbCon
    >> surrBasedLocalMeritFn >> surrBasedLocalAcceptLogic
    >> surrBasedLocalConstrRelax >> surrBasedGlobalReplacePts;
  //>> branchBndNumSamplesRoot >> branchBndNumSamplesNode

  // DL_SOLVER
  s >> dlDetails;

  // NPSOL
  s >> verifyLevel >> functionPrecision >> lineSearchTolerance;

  // NL2SOL
  s >> absConvTol >> xConvTol >> singConvTol >> singRadius >> falseConvTol
    >> initTRRadius >> covarianceType >> regressDiag;

  // OPT++
  s >> searchMethod >> gradientTolerance >> maxStep >> meritFn
    >> stepLenToBoundary >> centeringParam >> searchSchemeSize;

  // APPSPACK
  s >> initStepLength >> contractStepLength >> threshStepLength >> meritFunction
    >> constrPenalty >> smoothFactor;

  // COLINY
  s >> constraintPenalty >> constantPenalty >> globalBalanceParam
    >> localBalanceParam >> maxBoxSize >> minBoxSize >> boxDivision
    >> mutationAdaptive >> showMiscOptions >> miscOptions >> solnTarget
    >> crossoverRate >> mutationRate >> mutationScale >> mutationMinScale
    >> initDelta >> threshDelta >> contractFactor >> newSolnsGenerated
    >> numberRetained >> expansionFlag >> expandAfterSuccess
    >> contractAfterFail >> mutationRange >> totalPatternSize
    >> randomizeOrderFlag >> selectionPressure >> replacementType
    >> crossoverType >> mutationType >> exploratoryMoves >> patternBasis;

  // COLINY + APPSPACK
  s >> evalSynchronize;

  // JEGA
  s >> numCrossPoints >> numParents >> numOffspring >> fitnessType
    >> convergenceType >> percentChange >> numGenerations >> fitnessLimit
    >> shrinkagePercent >> nichingType >> nicheVector >> numDesigns
    >> postProcessorType >> distanceVector;

  // JEGA/COLINY
  s >> initializationType >> flatFile >> logFile >> populationSize
    >> printPopFlag;

  // NCSU 
  s >> volBoxSize;

  // DDACE
  s >> numSymbols >> mainEffectsFlag;

  // FSUDace 
  s >> latinizeFlag >> volQualityFlag >> sequenceStart >> sequenceLeap
    >> primeBase >> numTrials >> trialType;

  // COLINY, NonD, DACE, & JEGA
  s >> randomSeed;

  // MADS
  s >> historyFile >> displayFormat >> vns >> neighborOrder >> showAllEval;

  // NonD & DACE
  s >> numSamples >> fixedSeedFlag >> fixedSequenceFlag >> previousSamples
    >> vbdFlag >> vbdDropTolerance >> backfillFlag >> pcaFlag 
    >> percentVarianceExplained;

  // NonD
  s >> vbdOrder >> covarianceControl >> rngName >> refinementType
    >> refinementControl >> nestingOverride >> growthOverride >> expansionType
    >> piecewiseBasis >> expansionBasisType >> expansionOrder
    >> expansionSamples >> expansionSampleType >> quadratureOrder
    >> sparseGridLevel >> anisoDimPref >> cubIntOrder >> collocationPoints
    >> collocationRatio >> collocRatioTermsOrder >> regressionType
    >> lsRegressionType >> regressionNoiseTol >> regressionL2Penalty
    >> crossValidation //>> adaptedBasisInitLevel
    >> adaptedBasisAdvancements >> normalizedCoeffs >> pointReuse
    >> tensorGridFlag >> tensorGridOrder >> importExpansionFile
    >> exportExpansionFile >> sampleType >> reliabilitySearchType
    >> reliabilityIntegration >> integrationRefine >> refineSamples
    >> distributionType >> responseLevelTarget >> responseLevelTargetReduce
    >> responseLevels >> probabilityLevels >> reliabilityLevels
    >> genReliabilityLevels >> emulatorSamples >> emulatorOrder >> emulatorType
    >> mcmcType >> standardizedSpace >> adaptPosteriorRefine >> logitTransform
    >> preSolveMethod >> proposalCovType >> proposalCovUpdates
    >> proposalCovInputType >> proposalCovData >> proposalCovFile
    >> fitnessMetricType >> batchSelectionType >> batchSize
    >> calibrateErrorMode  >> hyperPriorAlphas >> hyperPriorBetas
    >> numChains >> numCR >> crossoverChainPairs
    >> grThreshold >> jumpStep >> lipschitzType >> dataDistType 
    >> dataDistCovInputType >> dataDistMeans >> dataDistCovariance
    >> dataDistFile >> posteriorDensityExportFilename
    >> posteriorSamplesExportFilename >> posteriorSamplesImportFilename
    >> generatePosteriorSamples >> evaluatePosteriorDensity;

  // Parameter Study
  s >> finalPoint >> stepVector >> numSteps >> stepsPerVariable >> listOfPoints
    >> pstudyFilename >> pstudyFileFormat >> pstudyFileActive
    >> varPartitions;

  // Verification
  s >> refinementRate;

  // Point import/export files
  s >> importBuildPtsFile  >> importBuildFormat  >> importBuildActive
    >> importApproxPtsFile >> importApproxFormat >> importApproxActive
    >> exportApproxPtsFile >> exportApproxFormat >> exportMCMCPtsFile
    >> exportMCMCFormat;
}


void DataMethodRep::write(std::ostream& s) const
{
  s << idMethod << modelPointer << methodOutput << maxIterations
    << maxFunctionEvaluations << speculativeFlag << methodUseDerivsFlag
    << convergenceTolerance << constraintTolerance << methodScaling
    << numFinalSolutions << linearIneqConstraintCoeffs << linearIneqLowerBnds
    << linearIneqUpperBnds << linearIneqScaleTypes << linearIneqScales
    << linearEqConstraintCoeffs << linearEqTargets << linearEqScaleTypes
    << linearEqScales << methodName << subMethod << subMethodName
    << subModelPointer << subMethodPointer;

  // Meta-iterators
  s << iteratorServers << procsPerIterator << iteratorScheduling
    << hybridMethodNames << hybridModelPointers << hybridMethodPointers
  //<< hybridProgThresh
    << hybridGlobalMethodName << hybridGlobalModelPointer
    << hybridGlobalMethodPointer << hybridLocalMethodName
    << hybridLocalModelPointer << hybridLocalMethodPointer << hybridLSProb
  //<< branchBndNumSamplesRoot << branchBndNumSamplesNode
    << concurrentRandomJobs << concurrentParameterSets;

  // Surrogate-based
  s << softConvLimit << surrBasedLocalLayerBypass
    << surrBasedLocalTRInitSize << surrBasedLocalTRMinSize
    << surrBasedLocalTRContractTrigger << surrBasedLocalTRExpandTrigger
    << surrBasedLocalTRContract << surrBasedLocalTRExpand
    << surrBasedLocalSubProbObj << surrBasedLocalSubProbCon
    << surrBasedLocalMeritFn << surrBasedLocalAcceptLogic
    << surrBasedLocalConstrRelax << surrBasedGlobalReplacePts;
  //<< branchBndNumSamplesRoot << branchBndNumSamplesNode

  // DL_SOLVER
  s << dlDetails;

  // NPSOL
  s << verifyLevel << functionPrecision << lineSearchTolerance;

  // NL2SOL
  s << absConvTol << xConvTol << singConvTol << singRadius << falseConvTol
    << initTRRadius << covarianceType << regressDiag;

  // OPT++
  s << searchMethod << gradientTolerance << maxStep << meritFn
    << stepLenToBoundary << centeringParam << searchSchemeSize;

  // APPSPACK
  s << initStepLength << contractStepLength << threshStepLength << meritFunction
    << constrPenalty << smoothFactor;

  // COLINY
  s << constraintPenalty << constantPenalty << globalBalanceParam
    << localBalanceParam << maxBoxSize << minBoxSize << boxDivision
    << mutationAdaptive << showMiscOptions << miscOptions << solnTarget
    << crossoverRate << mutationRate << mutationScale << mutationMinScale
    << initDelta << threshDelta << contractFactor << newSolnsGenerated
    << numberRetained << expansionFlag << expandAfterSuccess
    << contractAfterFail << mutationRange << totalPatternSize
    << randomizeOrderFlag << selectionPressure << replacementType
    << crossoverType << mutationType << exploratoryMoves << patternBasis;

  // COLINY + APPSPACK
  s << evalSynchronize;

  // JEGA
  s << numCrossPoints << numParents << numOffspring << fitnessType
    << convergenceType << percentChange << numGenerations << fitnessLimit
    << shrinkagePercent << nichingType << nicheVector << numDesigns
    << postProcessorType << distanceVector;

  // JEGA/COLINY
  s << initializationType << flatFile << logFile << populationSize
    << printPopFlag;

  // NCSU 
  s << volBoxSize;

  // DDACE
  s << numSymbols << mainEffectsFlag;

  // FSUDace 
  s << latinizeFlag << volQualityFlag << sequenceStart << sequenceLeap
    << primeBase << numTrials << trialType;

  // COLINY, NonD, DACE, & JEGA
  s << randomSeed;

  // MADS
  s << historyFile << displayFormat << vns << neighborOrder << showAllEval;

  // NonD & DACE
  s << numSamples << fixedSeedFlag << fixedSequenceFlag << previousSamples
    << vbdFlag << vbdDropTolerance << backfillFlag << pcaFlag 
    << percentVarianceExplained;

  // NonD
  s << vbdOrder << covarianceControl << rngName << refinementType
    << refinementControl << nestingOverride << growthOverride << expansionType
    << piecewiseBasis << expansionBasisType << expansionOrder
    << expansionSamples << expansionSampleType << quadratureOrder
    << sparseGridLevel << anisoDimPref << cubIntOrder << collocationPoints
    << collocationRatio << collocRatioTermsOrder << regressionType
    << lsRegressionType << regressionNoiseTol << regressionL2Penalty
    << crossValidation //<< adaptedBasisInitLevel
    << adaptedBasisAdvancements << normalizedCoeffs << pointReuse
    << tensorGridFlag << tensorGridOrder << importExpansionFile
    << exportExpansionFile << sampleType << reliabilitySearchType
    << reliabilityIntegration << integrationRefine << refineSamples
    << distributionType << responseLevelTarget << responseLevelTargetReduce
    << responseLevels << probabilityLevels << reliabilityLevels
    << genReliabilityLevels << emulatorSamples << emulatorOrder << emulatorType
    << mcmcType << standardizedSpace << adaptPosteriorRefine << logitTransform
    << preSolveMethod << proposalCovType << proposalCovUpdates
    << proposalCovInputType << proposalCovData << proposalCovFile
    << fitnessMetricType << batchSelectionType << batchSize
    << calibrateErrorMode << hyperPriorAlphas << hyperPriorBetas
    << numChains << numCR << crossoverChainPairs
    << grThreshold << jumpStep << lipschitzType << dataDistType 
    << dataDistCovInputType << dataDistMeans << dataDistCovariance
    << dataDistFile << posteriorDensityExportFilename
    << posteriorSamplesExportFilename << posteriorSamplesImportFilename
    << generatePosteriorSamples << evaluatePosteriorDensity;

  // Parameter Study
  s << finalPoint << stepVector << numSteps << stepsPerVariable << listOfPoints
    << pstudyFilename << pstudyFileFormat << pstudyFileActive
    << varPartitions;

  // Verification
  s << refinementRate;

  // Point import/export files
  s << importBuildPtsFile  << importBuildFormat  << importBuildActive
    << importApproxPtsFile << importApproxFormat << importApproxActive
    << exportApproxPtsFile << exportApproxFormat << exportMCMCPtsFile
    << exportMCMCFormat;
}


DataMethod::DataMethod(): dataMethodRep(new DataMethodRep())
{
#ifdef REFCOUNT_DEBUG
  Cout << "DataMethod::DataMethod(), dataMethodRep referenceCount = "
       << dataMethodRep->referenceCount << std::endl;
#endif
}


DataMethod::DataMethod(const DataMethod& data_method)
{
  // Increment new (no old to decrement)
  dataMethodRep = data_method.dataMethodRep;
  if (dataMethodRep) // Check for an assignment of NULL
    ++dataMethodRep->referenceCount;

#ifdef REFCOUNT_DEBUG
  Cout << "DataMethod::DataMethod(DataMethod&)" << std::endl;
  if (dataMethodRep)
    Cout << "dataMethodRep referenceCount = " << dataMethodRep->referenceCount
	 << std::endl;
#endif
}


DataMethod& DataMethod::operator=(const DataMethod& data_method)
{
  if (dataMethodRep != data_method.dataMethodRep) { // normal case: old != new
    // Decrement old
    if (dataMethodRep) // Check for NULL
      if ( --dataMethodRep->referenceCount == 0 ) 
	delete dataMethodRep;
    // Assign and increment new
    dataMethodRep = data_method.dataMethodRep;
    if (dataMethodRep) // Check for NULL
      ++dataMethodRep->referenceCount;
  }
  // else if assigning same rep, then do nothing since referenceCount
  // should already be correct

#ifdef REFCOUNT_DEBUG
  Cout << "DataMethod::operator=(DataMethod&)" << std::endl;
  if (dataMethodRep)
    Cout << "dataMethodRep referenceCount = " << dataMethodRep->referenceCount
	 << std::endl;
#endif

  return *this;
}


DataMethod::~DataMethod()
{
  if (dataMethodRep) { // Check for NULL
    --dataMethodRep->referenceCount; // decrement
#ifdef REFCOUNT_DEBUG
    Cout << "dataMethodRep referenceCount decremented to "
         << dataMethodRep->referenceCount << std::endl;
#endif
    if (dataMethodRep->referenceCount == 0) {
#ifdef REFCOUNT_DEBUG
      Cout << "deleting dataMethodRep" << std::endl;
#endif
      delete dataMethodRep;
    }
  }
}

} // namespace Dakota
