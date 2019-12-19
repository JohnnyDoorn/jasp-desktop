#ifndef ENGINEDEFINITIONS_H
#define ENGINEDEFINITIONS_H

#include "enumutilities.h"

DECLARE_ENUM(engineState,			initializing, idle, analysis, filter, rCode, computeColumn, moduleRequest, pauseRequested, paused, resuming, stopRequested, stopped, logCfg, settings);
DECLARE_ENUM(performType,			init, run, abort, saveImg, editImg, rewriteImgs);
DECLARE_ENUM(analysisResultStatus,	validationError, fatalError, imageSaved, imageEdited, imagesRewritten, complete, inited, running, changed, waiting);
DECLARE_ENUM(moduleStatus,			initializing, installNeeded, installModPkgNeeded, loadingNeeded, unloadingNeeded, readyForUse, error);
DECLARE_ENUM(engineAnalysisStatus,	empty, toInit, initing, inited, toRun, running, changed, complete, error, exception, aborted, stopped, saveImg, editImg, rewriteImgs, synchingData);

#endif // ENGINEDEFINITIONS_H
