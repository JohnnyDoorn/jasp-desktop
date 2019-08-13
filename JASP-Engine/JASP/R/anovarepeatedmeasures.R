#
# Copyright (C) 2013-2018 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
AnovaRepeatedMeasures <- function(jaspResults, dataset = NULL, options) {
  
  numericVariables <- c(unlist(options$repeatedMeasuresCells),unlist(options$covariates))
  numericVariables <- numericVariables[numericVariables != ""]
  factorVariables <- c(unlist(options$betweenSubjectFactors))
  factorVariables <- factorVariables[factorVariables != ""]
  
  if (is.null(dataset)) 
    dataset <- .readDataSetToEnd(columns.as.numeric=numericVariables, columns.as.factor=factorVariables, 
                                 exclude.na.listwise=c(numericVariables, factorVariables))
  
  longData <- .BANOVAreadRManovaData(dataset, options)

  ready <- length(options$repeatedMeasuresCells) > 1 &&  length(options$withinModelTerms) > 0
  
  ## Error checking
  # .rmAnovaCheckErrors(dataset, options, ready)
  .BANOVAerrorhandling(longData, options, "RM-ANOVA")
  
  ## Create container 
  rmAnovaContainer <- .getRMAnovaContainer(jaspResults)

  ## Compute All Anova Results
  .rmAnovaComputeResultsContainer(rmAnovaContainer, longData, options, ready)
  
  ## Create Within Subjects Effects Table
  .rmAnovaWithinSubjectsTable(rmAnovaContainer, dataset, options, ready)  
  
  ## Create Between Subjects Effects Table
  .rmAnovaBetweenSubjectsTable(rmAnovaContainer, dataset, options, ready)
  
  ## Create Reference Grid
  .referenceGrid(rmAnovaContainer, options, ready)
  
  ## Create Sphericity Table
  .rmAnovaAssumptionsContainer(rmAnovaContainer, dataset, options, ready)
  
  # Descriptives
  options[["credibleInterval"]] <- 0.95
  # plotDepends <- c("plotHorizontalAxis", "plotSeparateLines", "plotSeparatePlots", "plotErrorBars", "usePooledStandErrorCI")
  .BANOVAdescriptives(rmAnovaContainer, longData, options, list(noVariables=FALSE), "RM-ANOVA")
  
  ## Create Post Hoc Tables
  .resultsPostHoc(rmAnovaContainer, dataset, options, ready)

  .rmAnovaContrastTable(rmAnovaContainer, dataset, longData, options, ready)
  return()
  
  .resultsContrasts(rmAnovaContainer, dataset, options, ready)
  
  return()
  
  
  
  ## Create Simple Effects Table
  result <- .rmAnovaSimpleEffects(dataset, options, perform, fullModel, results[["withinSubjectsEffects"]], 
                                  results[["betweenSubjectsEffects"]], status, singular, stateSimpleEffects)
  
  
  # Create Friedman Table
  result <- .rmAnovaFriedman(dataset, fullModel, options, perform, status, singular, stateFriedman)
  
  # Create Conover Table
  result <- .rmAnovaConoverTable(dataset, options, perform, anovaModel$fullModel, status, stateConover, singular)
  
  

  ## Create Marginal Means Tables
  result <- .rmAnovaMarginalMeansTable(dataset, options, perform, status, fullModel)
  
  ## Create Marginal Means via Bootstrapping Tables
  result <- .rmAnovaMarginalMeansBootstrappingTable(dataset, options, perform, status, fullModel)
  
  
  return()
}

.getRMAnovaContainer <- function(jaspResults) {
  
  if (!is.null(jaspResults[["rmAnovaContainer"]])) {
    
    anovaContainer <- jaspResults[["rmAnovaContainer"]]
    
  } else {
    
    anovaContainer <- createJaspContainer(title = "Repeated Measures ANOVA")
    # we set the dependencies on the container, this means that all items inside the container automatically have these dependencies
    anovaContainer$dependOn(c("withinModelTerms", "betweenModelTerms", "repeatedMeasuresCells", 
                              "repeatedMeasuresFactors", "covariates", "sumOfSquares"))
    jaspResults[["rmAnovaContainer"]] <- anovaContainer
  }
  
  return(anovaContainer)
  
}

.rmAnovaCheckErrors <- function(dataset, options, ready) {
  if (!ready) return()
  
  modelTerms <- unlist(options$withinModelTerms, recursive = FALSE)
  betweenModelTerms <- options$betweenModelTerms
  
  for(betweenTerm in rev(betweenModelTerms)) {
    .hasErrors(
      dataset = dataset, 
      type = c('observations', 'variance', 'infinity'),
      all.target = c(options$repeatedMeasuresCells, options$covariates), 
      all.grouping = betweenTerm,
      observations.amount = "< 2", 
      exitAnalysisIfErrors = TRUE)
  }
  
  for(betweenTerm in rev(betweenModelTerms)) {
    .hasErrors(
      dataset = dataset, 
      type = c('infinity', 'factorLevels'),
      all.target = betweenTerm, 
      factorLevels.amount  = "< 2",
      exitAnalysisIfErrors = TRUE)
  }
  
}

.rmAnovaCheck <- function(dataset, options, perform) {
  
  error <- FALSE
  errorMessage <- NULL
  ready <- "" %in% options$repeatedMeasuresCells == FALSE && length(options$withinModelTerms) > 0
  
  if (ready && perform == "run") {
    
    components <- unique(unlist(options$betweenSubjectFactors))
    independentsWithLessThanTwoLevels <- c()
    
    for (component in components) {
      
      nLevels <- length(levels(dataset[[ .v(component) ]]))
      
      if (nLevels < 2)
        independentsWithLessThanTwoLevels <- c(independentsWithLessThanTwoLevels, component)
    }
    
    if (length(independentsWithLessThanTwoLevels) > 0) {
      
      error <- TRUE
      if(length(independentsWithLessThanTwoLevels) == 1) {
        errorMessage <- paste("Factor: <em>", independentsWithLessThanTwoLevels, "</em>, contains fewer than two levels.", sep="")
      } else {
        errorMessage <- paste("Factors: <em>", paste(independentsWithLessThanTwoLevels, collapse=",", sep=""), "</em>, contain fewer than two levels.", sep="")
      }
    }
    
    repeatedMeasuresData <- list()
    for(i in options$repeatedMeasuresCells) {
      repeatedMeasuresData[[i]] <- dataset[[.v(i)]]
    }
    infiniteRM <- unlist(lapply(repeatedMeasuresData,function(x)sum(is.infinite(x)) > 0))
    
    if (!is.null(infiniteRM) && sum(infiniteRM) > 0) {
      
      error <- TRUE
      if(sum(infiniteRM) == 1) {
        errorMessage <- paste("The repeated measure: <em>", options$repeatedMeasuresCells[infiniteRM], "</em>, contains infinite values.", sep="")
      } else {
        errorMessage <- paste("The repeated measures: <em>", paste(options$repeatedMeasuresCells[infiniteRM], collapse=", "), "</em>, contain infinite values.", sep="")
      }
    }
    
    covariatesData <- list()
    for(i in options$covariates) {
      covariatesData[[i]] <- dataset[[.v(i)]]
    }
    infiniteCov <- unlist(lapply(covariatesData,function(x)sum(is.infinite(x)) > 0))
    
    if (!is.null(infiniteCov) && sum(infiniteCov) > 0) {
      
      error <- TRUE
      if(sum(infiniteCov) == 1) {
        errorMessage <- paste("The covariate: <em>", options$covariates[infiniteCov], "</em>, contains infinite values.", sep="")
      } else {
        errorMessage <- paste("The covariates: <em>", paste(options$covariates[infiniteCov], collapse=", "), "</em>, contain infinite values.", sep="")
      }
    }
    
    allNames <- unlist(lapply(options[['repeatedMeasuresFactors']], function(x) x$name)) # Factornames 
    for(factorName in allNames){
      if (any(factorName %in% options$betweenSubjectFactors )) {
        error <- TRUE
        errorMessage <- paste("Please choose a name for the RM factors that differs from those for the 
		                          between subjects factors.", sep="")
      }
    } 
    
  }
  
  list(ready=ready, error=error, errorMessage=errorMessage)
}

.rmModelFormula <- function(options) {
  
  termsRM.base64 <- c()
  termsRM.normal <- c()
  
  for (term in options$withinModelTerms) {
    
    components <- unlist(term$components)
    termRM.base64 <- paste(.v(components), collapse=":", sep="")
    termRM.normal <- paste(components, collapse=" \u273B ", sep="")
    
    termsRM.base64 <- c(termsRM.base64, termRM.base64)
    termsRM.normal <- c(termsRM.normal, termRM.normal)
  }
  
  termsBS.base64 <- c()
  termsBS.normal <- c()
  
  for (term in options$betweenModelTerms) {
    
    components <- unlist(term$components)
    termBS.base64 <- paste(.v(components), collapse=":", sep="")
    termBS.normal <- paste(components, collapse=" \u273B ", sep="")
    
    termsBS.base64 <- c(termsBS.base64, termBS.base64)
    termsBS.normal <- c(termsBS.normal, termBS.normal)
  }
  
  terms.base64 <- list()
  terms.normal <- list()
  terms.base64[[1]] <- termsBS.base64
  terms.normal[[1]] <- termsBS.normal
  
  for (i in 1:length(termsRM.base64)) {
    if (is.null(termsBS.base64)) {
      terms.base64[[i+1]] <- termsRM.base64[i]
      terms.normal[[i+1]] <- termsRM.normal[i]
    } else if (!is.null(termsRM.base64)){
      terms.base64[[i+1]] <- c(termsRM.base64[i], paste(termsRM.base64[i], termsBS.base64, sep = ":"))
      terms.normal[[i+1]] <- c(termsRM.normal[i], paste(termsRM.normal[i], termsBS.normal, sep = " \u273B "))
    }
  }
  
  main <- paste("(",paste(unlist(terms.base64), collapse=" + "),")", sep="")
  termsBS <- paste("(",paste(termsBS.base64, collapse=" + "),")", sep="")
  errorRM <- paste("Error(",paste("Xc3ViamVjdA/(", termsRM.base64, ")",sep="", collapse=" + "),")",sep="")
  
  if (is.null(termsBS.base64) && is.null(termsRM.base64)) {
    model.def <- XZGVwZW5kZW50 ~ 1
  } else if (is.null(termsBS.base64)) {
    model.def <- paste("XZGVwZW5kZW50", "~", paste(main, errorRM, sep=" + "))
  } else if (is.null(termsRM.base64)) {
    model.def <- paste("XZGVwZW5kZW50", "~", main)
  } else {
    model.def <- paste("XZGVwZW5kZW50", "~", paste(main, errorRM, termsBS, sep=" + "))
  }
  
  list(model.def = model.def, terms.normal = terms.normal, terms.base64 = terms.base64, termsRM.normal = termsRM.normal, termsRM.base64 = termsRM.base64)
}

.rmAnovaComputeResultsContainer <- function(rmAnovaContainer, dataset, options, ready) {
  if (!ready) return()
  
  # Take results from state if possible
  if (!is.null(rmAnovaContainer[["anovaResult"]]))
    return()
  
  rmAnovaResult <- .rmAnovaComputeResults(dataset, options)
  
  if (rmAnovaResult[["tryResult"]] == "try-error") {
    rmAnovaContainer$setError("Some parameters are not estimable, most likely due to empty cells of the design.")
    return()
  }
  
  # Save model to state
  rmAnovaContainer[["anovaResult"]] <- createJaspState(object = rmAnovaResult)
}

.rmAnovaComputeResults <- function(dataset, options, status) {
  
  modelDef <- .rmModelFormula(options)
  model.formula <- as.formula(modelDef$model.def)

  variables <- unlist(c(options$betweenSubjectFactors, lapply(options$repeatedMeasuresFactors, function(x) x$name)))
  
  for (i in variables)
    dataset[[.v(i)]] <- .v(dataset[[.v(i)]])
  
  options(contrasts=c("contr.sum","contr.poly"))
  
  # set these options once for all afex::aov_car calls,
  # this ensures for instance that afex::aov_car always returns objects of class afex_aov.
  afex::afex_options(
    check_contrasts = TRUE, correction_aov = "GG", 
    emmeans_model = "univariate", es_aov = "ges", factorize = TRUE, 
    lmer_function = "lmerTest", method_mixed = "KR", return_aov = "afex_aov", 
    set_data_arg = FALSE, sig_symbols = c(" +", " *", " **", " ***"), type = 3
  )
  
  # Computations:
  if (options$sumOfSquares == "type1") {
    
    tryResult <- try({
      result <- stats::aov(model.formula, data=dataset)
      summaryResultOne <- summary(result, expand.split = FALSE)
    
      result <- afex::aov_car(model.formula, data=dataset, type= 3, factorize = FALSE)
      summaryResult <- summary(result)

      # Reformat the results to make it consistent with types 2 and 3
      model <- as.data.frame(unclass(summaryResult$univariate.tests))

      for (mySub in unlist(summaryResultOne, recursive = FALSE)) {
        for(term in trimws(rownames(mySub)[-nrow(mySub)])) {
          model[term, "Sum Sq"] <- mySub[term, "Sum Sq"]
          model[term, "num Df"] <- mySub[term, "Df"]
          model[term, "F value"] <- mySub[term, "F value"]
          model[term, "Pr(>F)"] <- mySub[term, "Pr(>F)"]
          model[term, "Error SS"] <-  mySub["Residuals", "Sum Sq"]
          model[term, "den Df"] <- mySub["Residuals", "Df"]
        }
      }
    })
    
  } else if (options$sumOfSquares == "type2") {
    
    tryResult <- try({
      result <- afex::aov_car(model.formula, data=dataset, type= 2, factorize = FALSE)
      summaryResult <- summary(result)
      model <- as.data.frame(unclass(summaryResult$univariate.tests))
    })
    
  } else {
    
    tryResult <- try({
      result <- afex::aov_car(model.formula, data=dataset, type= 3, factorize = FALSE)
      summaryResult <- summary(result)
      model <- as.data.frame(unclass(summaryResult$univariate.tests))
    })
    
  }

  if (class(tryResult) == "try-error") {
    return(list(tryResult == "try-error"))
  }

  # Now we reformat the results table some more to make it flow with jaspResults later
  interceptRow <- model["(Intercept", ]
  model <- model[-1,]
  rownames(model) <- trimws(rownames(model))
  model[["isWithinTerm"]] <- model[[".isNewGroup"]] <- numeric(nrow(model))

  sortedModel <- model
  cases <- unlist(sapply(modelDef$terms.base64, function(x) x[[1]]))
  residualResults <- sortedModel[cases, ]

  nextNewGroup <- 0
  counter <- 1
  for (modelTerm in modelDef$terms.base64) {
    
    if (!is.null(modelTerm)) {
      isWithin <- any(modelTerm %in% modelDef$termsRM.base64)
      indices <- .mapAnovaTermsToTerms(modelTerm, rownames(model))
      nextNewGroup <- c(TRUE, rep(FALSE, length(indices) - 1))
      sortedModel[indices, ] <- model[indices, ]
      sortedModel[indices, c(".isNewGroup", "isWithinTerm")] <- c(nextNewGroup, rep(isWithin, length(indices)))
  
      residualResults[modelTerm[[1]], ] <- c(model[indices[1],  c("Error SS", "den Df")], rep(NA, 4), 0, isWithin) 
    }

  }

  sortedModel[["case"]] <- unlist(modelDef$terms.normal)
  rownames(sortedModel) <- unlist(modelDef$terms.base64)
  rownames(residualResults) <- cases
  sortedModel[["Mean Sq"]] <- sortedModel[["Sum Sq"]] / sortedModel[["num Df"]]
  residualResults[["Mean Sq"]] <- residualResults[["Sum Sq"]] / residualResults[["num Df"]]
  residualResults[["case"]] <- "Residuals"

  # Now we calculate effect sizes
  SSr <- sortedModel[["Error SS"]]
  MSr <- SSr/sortedModel[["den Df"]]
  
  sortedModel[["eta"]] <- sortedModel[["Sum Sq"]] / (sum(sortedModel[["Sum Sq"]]) + sum(residualResults[["Sum Sq"]]))
  sortedModel[["etaPart"]] <- sortedModel[["Sum Sq"]] / (sortedModel[["Sum Sq"]] + SSr)
  sortedModel[["genEta"]]<- result[["anova_table"]][["ges"]]
  
  n <- interceptRow[["den Df"]] + 1
  MSb <- interceptRow[["Error SS"]] / (n-1)
  MSm <- sortedModel[["Mean Sq"]]
  df <- sortedModel[["num Df"]]
  
  omega <- (df / (n * (df + 1)) * (MSm - MSr)) / (MSr + ((MSb - MSr) / (df + 1)) + 
                                                                     (df / (n * (df + 1))) * (MSm - MSr))
  sortedModel[["omega"]] <- sapply(omega, max, 0)
  # Now we include the results from the corrections
  corrections <- summaryResult$pval.adjustments
  withinAnovaTable <- ggTable <- hfTable <- subset(sortedModel, isWithinTerm == 1)
  withinIndices <- .mapAnovaTermsToTerms(rownames(sortedModel), rownames(corrections))
  
  withinAnovaTable[["correction"]] <- "None"
  
  ggTable[["num Df"]] <- withinAnovaTable[["num Df"]] * corrections[withinIndices, "GG eps"]
  ggTable[["Mean Sq"]] <-  withinAnovaTable[["Sum Sq"]] / ggTable[["num Df"]]
  ggTable[["den Df"]] <- withinAnovaTable[["den Df"]] * corrections[withinIndices, "GG eps"]
  ggTable[["Pr(>F)"]] <-  pf(withinAnovaTable[["F value"]], ggTable[["num Df"]], 
                             ggTable[["den Df"]], lower.tail = FALSE)
  ggTable[["correction"]] <-  "Greenhouse-Geisser"
  ggTable[[".isNewGroup"]] <- FALSE 
  
  hfTable[["Mean Sq"]] <- withinAnovaTable[["Sum Sq"]] / hfTable[["num Df"]]
  hfTable[["num Df"]] <-  withinAnovaTable[["num Df"]] * corrections[, "HF eps"]
  hfTable[["den Df"]] <- withinAnovaTable[["den Df"]] * corrections[, "HF eps"]
  hfTable[["Pr(>F)"]] <-  pf(withinAnovaTable[["F value"]], hfTable[["num Df"]], 
                             hfTable[["den Df"]], lower.tail = FALSE)
  hfTable[["correction"]] <-  "Huyn-Feldt"
  hfTable[[".isNewGroup"]] <- FALSE 
 
  wResidualResults <- wResidualResultsGG <- wResidualResultsHF <- subset(residualResults, isWithinTerm == 1)
  wResidualResults[["correction"]] <- "None"
  
  residualIndices <- .mapAnovaTermsToTerms(rownames(wResidualResults), rownames(corrections))
  wResidualResultsGG[["num Df"]] <- wResidualResults[["num Df"]] * corrections[residualIndices, "GG eps"]
  wResidualResultsGG[["Mean Sq"]] <- wResidualResults[["Sum Sq"]] / wResidualResultsGG[["num Df"]]
  wResidualResultsGG[["correction"]] <-  "Greenhouse-Geisser"
  
  wResidualResultsHF[["num Df"]] <- wResidualResults[["num Df"]] * corrections[residualIndices, "HF eps"]
  wResidualResultsHF[["Mean Sq"]] <- wResidualResults[["Sum Sq"]] / wResidualResultsHF[["num Df"]]
  wResidualResultsHF[["correction"]] <-  "Huyn-Feldt"
  
  withinAnovaTable <- cbind(withinAnovaTable, corrections[withinIndices, ], 
                            summaryResult$sphericity.tests[withinIndices, ])
  # Makes lists with results and insert Intercept row back in anova table
  withinAnovaTableCollection <- list("None" = withinAnovaTable, "Huyn-Feldt" = hfTable, "Greenhouse-Geisser" = ggTable)
  
  for (tableIndex in .indices(withinAnovaTableCollection)) {
    withinAnovaTableCollection[[tableIndex]][["VovkSellkeMPR"]] <- .VovkSellkeMPR(withinAnovaTableCollection[[tableIndex]][["Pr(>F)"]])
  }
  sortedModel[["VovkSellkeMPR"]] <- .VovkSellkeMPR(sortedModel[["Pr(>F)"]])
  
  
  wResidualResultsList <- list("None" = wResidualResults, "Huyn-Feldt" = wResidualResultsHF, 
                               "Greenhouse-Geisser" = wResidualResultsGG)
  # sortedModel["Intercept", ] <-  c(interceptRow, 0, 0, "Residuals", rep(NA, ncol(sortedModel)-3))

  return(list(anovaResult = sortedModel, 
              residualTable = wResidualResultsList, 
              withinAnovaTable = withinAnovaTableCollection,
              assumptionResult = cbind(summaryResult$sphericity.tests, summaryResult$pval.adjustments), 
              fullModel = result, 
              tryResult = "tryResult"))
}

.mapAnovaTermsToTerms <- function(oneTerms, twoTerms) {
  nTerms <- length(oneTerms)
  indices <- numeric(nTerms)
  counter <- 1
  
  for (i in 1:nTerms) {
    splitFirst <- strsplit(oneTerms[[i]],":")[[1]]
    for (j in 1:length(twoTerms)) {
      matchedTerms <- match(splitFirst, strsplit(twoTerms[[j]],":")[[1]])
      if ((length(strsplit(twoTerms[[j]],":")[[1]]) == length(splitFirst)) && !any(is.na(matchedTerms))) {
        indices[counter] <- j
        counter <- counter + 1
      }
    }
  }
  
  return(indices)
}

.identicalTerms <- function(x,y) {
  
  equalLength <- length(x) == length(y)
  equalTerms <- all(x %in% y)
  
  if(equalLength && equalTerms) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

.rmAnovaBetweenSubjectsTable <- function(rmAnovaContainer, dataset, options, ready) {
  if(!is.null(rmAnovaContainer[["betweenTable"]]) || length(options$betweenSubjectFactors) == 0)
    return()
  
  betweenTable <- createJaspTable(title = "Between Subjects Effects")
  
  betweenTable$addColumnInfo(title = "Cases", name = "case", type = "string" )
  betweenTable$addColumnInfo(title = "Sum of Squares", name = "Sum Sq", type = "number")
  betweenTable$addColumnInfo(title = "df", name = "num Df", type = "integer")
  betweenTable$addColumnInfo(title = "Mean Square", name = "Mean Sq", type = "number")
  betweenTable$addColumnInfo(title = "F", name = "F value", type = "number")
  betweenTable$addColumnInfo(title = "p", name = "Pr(>F)", type = "number")
  
  if (options$VovkSellkeMPR) {
    betweenTable$addColumnInfo(title = "VS-MPR\u002A", name = "VovkSellkeMPR", type = "number")
    betweenTable$addFootnote(message = .messages("footnote", "VovkSellkeMPR"), symbol = "\u002A")
  }
  
  if (options$effectSizeEstimates) {
    
    if (options$effectSizeEtaSquared) 
      betweenTable$addColumnInfo(title = "\u03B7\u00B2", name = "eta", type = "number")
    
    if (options$effectSizePartialEtaSquared) 
      betweenTable$addColumnInfo(title = "\u03B7\u00B2\u209A", name = "etaPart", type = "number")
    
    if (options$effectSizeGenEtaSquared) 
      betweenTable$addColumnInfo(title = "\u03B7\u00B2<sub>G</sub>", name = "genEta", type = "number")
    
    if (options$effectSizeOmegaSquared) 
      betweenTable$addColumnInfo(title = "\u03C9\u00B2", name = "omega", type = "number")
    
  }
  
  betweenTable$showSpecifiedColumnsOnly <- TRUE
  
  # set the type footnote already
  typeFootnote <- switch(options$sumOfSquares,
                         type1 = "Type I Sum of Squares",
                         type2 = "Type II Sum of Squares",
                         type3 = "Type III Sum of Squares")
  betweenTable$addFootnote(message = typeFootnote, symbol = "<em>Note.</em>")
  
  rmAnovaContainer[["betweenTable"]] <- betweenTable
  
  if (!ready) {
    return()
  }
  
  result <- rmAnovaContainer[["anovaResult"]]$object$anovaResult
  result <- result[result$isWithinTerm == 0, ]
  result["Residuals", ] <- NA
  result["Residuals", "case"] <- "Residuals"
  result["Residuals", "num Df"] <- result[["den Df"]][1]
  result["Residuals", "Sum Sq"] <- result[["Error SS"]][1]
  result["Residuals", "Mean Sq"] <- result[["Error SS"]][1] / result[["den Df"]][1]
  
  betweenTable$setData(result)
  return()
  
  
  #     if (options$effectSizeEstimates) {
  #       
  #       if (sum(gsub(" ", "", row.names(resultTable), fixed = TRUE) == "Residuals") > 0) {
  #         
  #         SSt <- sum(resultTable[,"Sum Sq"])
  #         SSr <- resultTable["Residuals","Sum Sq"]
  #         MSr <- SSr/resultTable["Residuals","Df"]
  # 
  #         row[["eta"]] <- SS / SSt
  #         row[["partialEta"]] <- SS / (SS + SSr)
  #         row[["genEta"]] <- fullModel[["anova_table"]][modelTermsCase, "ges"]
  #         
  #         omega <- (SS - (df * MSr)) / (SSt + MSr) 
  # 
  #         if (omega < 0) {
  #           row[["omega"]] <- 0
  #         } else {
  #           row[["omega"]] <- omega
  #         }
  #         
  # 
  # # } else {	# if sum of squares is equal to type 2 or 3
  #   
  #          # if (result[index,"Error SS"] > 0) {
  #         # Is between subjects factor
  #         
  #         SSr <- result[index,"Error SS"]
  #         SSt <- sum(result[indices,"Sum Sq"]) + SSr
  #         MSr <- SSr/result[index,"den Df"]
  # 
  #         row[["eta"]] <- SS / SSt
  #         row[["partialEta"]] <- SS / (SS + SSr)
  #         row[["genEta"]] <- fullModel[["anova_table"]][modelTermsCase, "ges"]
  #         
  #         omega <- (SS - (df * MSr)) / (SSt + MSr)
  # 
  #         if (omega < 0) {
  #           row[["omega"]] <- 0
  #         } else {
  #           row[["omega"]] <- omega
  #         }
  # 
  return()
  
}

.rmAnovaWithinSubjectsTable <- function(rmAnovaContainer, dataset, options, ready) {
  if (!is.null(rmAnovaContainer[["withinAnovaTable"]]))
    return()
  
  anovaTable <- createJaspTable(title = "Within Subjects Effects", position = 1)
  rmAnovaContainer[["withinAnovaTable"]] <- anovaTable
  anovaTable$showSpecifiedColumnsOnly <- TRUE
  
  if (options$sphericityCorrections) {
    corrections <- c("None", "Greenhouse-Geisser", "Huyn-Feldt")[c(options$sphericityNone, 
                                                                   options$sphericityGreenhouseGeisser,
                                                                   options$sphericityHuynhFeldt)]
  } else 
    corrections <- "None"
  
  anovaTable$addColumnInfo(title = "Cases", name = "case", type = "string", combine = TRUE)
  
  dfType <- "integer" # Make df an integer unless corrections are applied
  if ((length(corrections) > 1 || any(!"None" %in% corrections)) && options$sphericityCorrections) {
    anovaTable$addColumnInfo(title = "Sphericity Correction", name = "correction", type = "string", combine = TRUE)
    dfType <- "number"
  }
  
  anovaTable$addColumnInfo(title = "Sum of Squares", name = "Sum Sq", type = "number")
  anovaTable$addColumnInfo(title = "df", name = "num Df", type = dfType)
  anovaTable$addColumnInfo(title = "Mean Square", name = "Mean Sq", type = "number")
  anovaTable$addColumnInfo(title = "F", name = "F value", type = "number")
  anovaTable$addColumnInfo(title = "p", name = "Pr(>F)", type = "number")
  
  if (options$VovkSellkeMPR) {
    anovaTable$addColumnInfo(title = "VS-MPR\u002A", name = "VovkSellkeMPR", type = "number")
    anovaTable$addFootnote(message = .messages("footnote", "VovkSellkeMPR"), symbol = "\u002A")
  }
  
  if (options$effectSizeEstimates) {
    
    if (options$effectSizeEtaSquared) {
      anovaTable$addColumnInfo(title = "\u03B7\u00B2", name = "eta", type = "number")
    }
    
    if (options$effectSizePartialEtaSquared) {
      anovaTable$addColumnInfo(title = "\u03B7\u00B2\u209A", name = "etaPart", type = "number")
    }
    
    if(options$effectSizeGenEtaSquared) {
      anovaTable$addColumnInfo(name="genEta", type="number", title="\u03B7\u00B2<sub>G</sub>")
    }
    
    if (options$effectSizeOmegaSquared) {
      anovaTable$addColumnInfo(title = "\u03C9\u00B2", name = "omega", type = "number")
    }
    
  }
  
  typeFootnote <- switch(options$sumOfSquares,
                         type1 = "Type I Sum of Squares",
                         type2 = "Type II Sum of Squares",
                         type3 = "Type III Sum of Squares")
  anovaTable$addFootnote(message = typeFootnote, symbol = "<em>Note.</em>")
  
  if (!ready)
    return()
  
  withinResults <- rmAnovaContainer[["anovaResult"]]$object$withinAnovaTable
  residualResults <- rmAnovaContainer[["anovaResult"]]$object$residualTable
  modelTerms <- .rmModelFormula(options)$termsRM.base64
  allCases <- rownames(withinResults[[1]])
  addResidualAfter <- allCases[.mapAnovaTermsToTerms(modelTerms, allCases) + 1]
  
  for (case in allCases) {
    
    for (cor in corrections)
      anovaTable$addRows(withinResults[[cor]][case, ])
    
    if (case %in% modelTerms) {
      currentCase <- case
    } else if (case %in% addResidualAfter) {
      for (cor in corrections)
        anovaTable$addRows(residualResults[[cor]][currentCase, ])
    }
  }
  
  
  return()
}

.rmAnovaAssumptionsContainer <- function(rmAnovaContainer, dataset, options, ready) {
  if (!is.null(rmAnovaContainer[["assumptionsContainer"]]))
    return()
  
  assumptionsContainer <- createJaspContainer(title = "Assumption Checks",
                                              dependencies = c("homogeneityTests", "sphericityTests"))
  
  rmAnovaContainer[["assumptionsContainer"]] <- assumptionsContainer

  if (options$homogeneityTests == TRUE)  
    .rmAnovaLevenesTable(rmAnovaContainer, dataset, options, ready)
  
  if (options$sphericityTests == TRUE) 
    .rmAnovaSphericityTable(rmAnovaContainer, dataset, options, ready)
  
  return()
}

.rmAnovaLevenesTable <- function(rmAnovaContainer, dataset, options, ready) {
  if (!is.null(rmAnovaContainer[["rmAnovaLevenesTable"]]) || options$homogeneityTests == FALSE)
    return()

  rmAnovaLevenesTable <- createJaspTable("Test for Equality of Variances (Levene's)")
  rmAnovaContainer[["assumptionsContainer"]][["rmAnovaLevenesTable"]] <- rmAnovaLevenesTable
  
  rmAnovaLevenesTable$addColumnInfo(name="case", type="string", title="")
  rmAnovaLevenesTable$addColumnInfo(name="F", type="number")
  rmAnovaLevenesTable$addColumnInfo(name="df1", type="number")
  rmAnovaLevenesTable$addColumnInfo(name="df2", type="number")
  rmAnovaLevenesTable$addColumnInfo(name="p", type="number")
  
  if (options$VovkSellkeMPR) {
    rmAnovaLevenesTable$addColumnInfo(title = "VS-MPR\u002A", name = "VovkSellkeMPR", type = "number")
    rmAnovaLevenesTable$addFootnote(message = .messages("footnote", "VovkSellkeMPR"), symbol = "\u002A")
  }
  rmAnovaLevenesTable$showSpecifiedColumnsOnly <- TRUE

  if (!ready)
    return()
  
  rmAnovaLevenesTable$setExpectedSize(length(options$repeatedMeasuresCells))
  if (length(options$betweenModelTerms) == 0) {
    rmAnovaLevenesTable$setError("Cannot perform homogeneity tests because there are no between subjects factors specified.")
    return()
  }
  
  for (i in .indices(options$repeatedMeasuresCells)) {
    
    interaction <- paste(.v(options$betweenSubjectFactors), collapse=":", sep="")
    if (length(options$covariates) > 0 ) {
      
      covterms <- paste(.v(options$covariates), collapse="+", sep="")
      combterms <- paste(c(interaction,covterms), collapse="+", sep="")
      levene.def <- paste(.v(options$repeatedMeasuresCells[i]), "~", combterms)
      
    } else {
      
      levene.def <- paste(.v(options$repeatedMeasuresCells[i]), "~", interaction)
      
    }
    
    levene.formula <- as.formula(levene.def)
    
    dummyAov <- aov(levene.formula, data = dataset, qr = T)
    resids <- abs(dummyAov$residuals)
    levene.def <- paste("resids", "~", interaction)
    levene.formula <- as.formula(levene.def)
    
    r <- summary(aov(levene.formula, dataset))
    error <- base::tryCatch(summary(aov(levene.formula, dataset)),error=function(e) e, warning=function(w) w)
    
    row <- data.frame(case = options$repeatedMeasuresCells[i],
                      F = r[[1]]$`F value`[1], 
                      df1 = r[[1]]$Df[1],
                      df2 = r[[1]]$Df[2],
                      p=r[[1]]$`Pr(>F)`[1],
                      ".isNewGroup" = i == 1,
                      VovkSellkeMPR = .VovkSellkeMPR(r[[1]]$`Pr(>F)`[1]))
    
    rmAnovaLevenesTable$addRows(row)
    
  }
  
  return()
}

.rmAnovaSphericityTable <- function(rmAnovaContainer, dataset, options, ready) {
  if (!is.null(rmAnovaContainer[["sphericityTable"]]) || !options$sphericityTests)
    return()

  sphericityTable <- createJaspTable("Test of Sphericity")
  
  sphericityTable$addColumnInfo(name="case", type="string", title="")
  sphericityTable$addColumnInfo(name="Test statistic", type="number", title="Mauchly's W")
  sphericityTable$addColumnInfo(name="approxChi", type="number", title="Approx. \u03A7\u00B2")
  sphericityTable$addColumnInfo(name="dfSphericity", type="integer")
  sphericityTable$addColumnInfo(name="p-value", type="number")
  sphericityTable$addColumnInfo(name="GG eps", type="number", title="Greenhouse-Geisser \u03B5")
  sphericityTable$addColumnInfo(name="HF eps", type="number", title="Huynh-Feldt \u03B5")
  sphericityTable$addColumnInfo(name="LB", type="number", title="Lower Bound \u03B5")

  sphericityTable$showSpecifiedColumnsOnly <- TRUE
  
  rmAnovaContainer[["assumptionsContainer"]][["sphericityTable"]] <- sphericityTable
  
  if(!ready)
    return()
  
  .approxChi <- function(df, n, W){
    d <- 1 - (2*df^2 + df + 2) / (6*df*(n-1))
    -(n-1)*d*log(W)
  }
  
  assumptionResult <- rmAnovaContainer[["anovaResult"]]$object$assumptionResult
  anovaResult <- rmAnovaContainer[["anovaResult"]]$object$withinAnovaTable$None
  
  df <- anovaResult[["num Df"]]
  anovaResult[["dfSphericity"]] <- df * (df + 1) / 2 - 1
  n  <- anovaResult[["den Df"]] / df + 1
  
  anovaResult[["approxChi"]] <- .approxChi(df, n, anovaResult[["Test statistic"]])
  anovaResult[["LB"]] <- 1 / df
  anovaResult[["HF eps"]] <- sapply(anovaResult[["HF eps"]], min, 1)
  anovaResult[[".isNewGroup"]] <- FALSE
  
  sphericityTable$setData(anovaResult)
  
  return()
}

.referenceGrid <- function (rmAnovaContainer, options, ready) {
  if (!is.null(rmAnovaContainer[["referenceGrid"]]) || !ready)
    return()
  
  fullModel <- rmAnovaContainer[["anovaResult"]]$object$fullModel 

  referenceGridList <- list()
  variables <- unlist(c(lapply(options$betweenModelTerms, 
                               function(x) {
                                 if (length(x$components) == 1) {
                                   return (.v(x$components))
                                 } else {
                                   return(NULL)
                                 }
                               }), lapply(options$withinModelTerms,
                                          function(x) {
                                            if (length(x$components) == 1) {
                                              return (.v(x$components))
                                            } else {
                                              return(NULL)
                                            }
                                          })
  ))
  
  postHocVariables <- unlist(options$postHocTestsVariables, recursive = FALSE)
  variablesPost <- unname(sapply(postHocVariables, function(x) paste(.v(x), collapse = ":")))
  
  variables <- union(variables, variablesPost)
  
  for (var in variables) {
    formula <- as.formula(paste("~", var))
    referenceGrid <- emmeans::emmeans(fullModel, formula)
    
    referenceGridList[[var]] <- referenceGrid
    
  }
  
  rmAnovaContainer[["referenceGrid"]] <- createJaspState(object = referenceGridList)
  
  return()
}

.resultsPostHoc <- function (rmAnovaContainer, dataset, options, ready) {
  if(!is.null(rmAnovaContainer[["postHocStandardContainer"]]) || length(options$postHocTestsVariables) ==0)
    return()
  
  postHocContainer <- createJaspContainer(title = "Post Hoc Tests")
  postHocContainer$dependOn(c("postHocTestsVariables", "postHocTestEffectSize", "postHocTestsBonferroni", 
                                      "postHocTestsHolm", "postHocTestsScheffe", "postHocTestsTukey",
                                      "postHocFlagSignificant", "confidenceIntervalsPostHoc", 
                                      "confidenceIntervalIntervalPostHoc", "postHocTestPooledError"))
  
  rmAnovaContainer[["postHocStandardContainer"]] <- postHocContainer
  
  postHocVariables <- unlist(options$postHocTestsVariables, recursive = FALSE)
  postHocVariablesListV <- unname(lapply(postHocVariables, .v))
  variables <- unname(sapply(postHocVariables, function(x) paste(.v(x), collapse = ":")))
  
  options$postHocTestsSidak <- FALSE
  for (postHocVarIndex in 1:length(postHocVariables)) {
    
    thisVarName <- paste(postHocVariables[[postHocVarIndex]], collapse = " \u273B ")
    interactionTerm <- length(postHocVariables[[postHocVarIndex]]) > 1
    
    postHocContainer[[variables[postHocVarIndex]]] <- .createPostHocStandardTable(thisVarName, interactionTerm, options)
  }
  
  if (!ready)
    return()
  
  referenceGrid <- rmAnovaContainer[["referenceGrid"]]$object
  fullModel <- rmAnovaContainer[["anovaResult"]]$object$fullModel
  

  postHocData <- fullModel$data$wide
  factorNamesV <- colnames(postHocData) # Names to use to refer to variables in data
  # Because there are multiple names for each variable in JASP, one of the things the following code does is make sure to get the correct naming
  # and refer to the correct actual variable. The different names are the actual name of the variable, the name the user gives in jasp for the lvel and factor, 
  # and also the name that JASP gives to it, which is a concatenation of "Level#_Level#', where the first refers to the factor and second to the level. 
 
  rmFactorIndex <- 1
  allNames <- unlist(lapply(options$repeatedMeasuresFactors, function(x) x$name)) # Factornames 

  for (var in variables) {
    
    # Results using the Bonferroni method
    resultPostHoc <- summary(pairs(referenceGrid[[var]], adjust="bonferroni"), 
                          infer = TRUE, level = options$confidenceIntervalIntervalPostHoc)
    
    resultPostHoc[["bonferroni"]] <- resultPostHoc[["p.value"]] 
    
    resultPostHoc[["tukey"]] <-  summary(pairs(referenceGrid[[var]], adjust="tukey"))[["p.value"]]
    
    resultPostHoc[["scheffe"]] <-  summary(pairs(referenceGrid[[var]], adjust="scheffe"))[["p.value"]]
    
    resultPostHoc[["holm"]] <-  summary(pairs(referenceGrid[[var]], adjust="holm"))[["p.value"]]
    
    comparisons <- strsplit(as.character(resultPostHoc$contrast), " - ")
    
    orderOfTerms <- unlist(lapply(options$repeatedMeasuresFactors, function(x) x$name))
    indexofOrderFactors <- match(allNames,orderOfTerms)

    if (any(var == .v(allNames))) {     ## If the variable is a repeated measures factor

      if (!options$postHocTestPooledError) {

        levelsOfThisFactor <- unlist(lapply(options$repeatedMeasuresFactors[rmFactorIndex], function(x) x$levels)) # Levels within Factor
        numberOfLevels <- length(unique(levelsOfThisFactor))
        splitNames <- unlist(lapply(strsplit(factorNamesV,  split = "_"), function(x) x[indexofOrderFactors[rmFactorIndex]]))
        listVarNamesToLevel <- list()  # create a list of vectors of variable names, used to group the dataset for the post-hoc t-tests
        
        for(i in 1:numberOfLevels){
          
          listVarNamesToLevel[[i]] <- factorNamesV[(splitNames %in% .v(levelsOfThisFactor[i]))]  
        
        }
        
        countr <- 1
        for (k in 1:numberOfLevels) {  ### Loop over all the levels within factor and do pairwise t.tests on them
          
          for (i in .seqx(k+1, numberOfLevels)) {
            
            bonfAdjustCIlevel <- 1-((1-options$confidenceIntervalIntervalPostHoc)/choose(numberOfLevels, 2))
            tResult <- t.test(rowMeans(postHocData[listVarNamesToLevel[[k]]]),rowMeans(postHocData[listVarNamesToLevel[[i]]]),
                              paired= TRUE, var.equal = FALSE, conf.level = bonfAdjustCIlevel)
            resultPostHoc[["estimate"]][countr] <- tResult$estimate
            resultPostHoc[["t.ratio"]][countr] <- tResult$statistic
            resultPostHoc[["SE"]][countr] <- tResult$estimate / tResult$statistic
            resultPostHoc[["p.value"]][countr] <- tResult$p.value
            resultPostHoc[["lower.CL"]][countr] <- tResult$conf.int[1]
            resultPostHoc[["upper.CL"]][countr] <- tResult$conf.int[2]
            countr <- countr + 1
            
          }
          
        }
        resultPostHoc[["bonferroni"]] <- p.adjust(resultPostHoc[["p.value"]], method = "bonferroni")
        resultPostHoc[["holm"]] <- p.adjust(resultPostHoc[["p.value"]], method = "holm")
      }
      resultPostHoc[["scheffe"]] <- "."
      resultPostHoc[["tukey"]] <-  "."
      
      if (options$postHocTestsScheffe || options$postHocTestsTukey) {
        cors <- paste(c("Tukey", "Scheffe")[c(options$postHocTestsTukey, options$postHocTestsScheffe)], collapse = " and ")
        
        postHocContainer[[var]]$addFootnote(paste(cors, "corrected p-values are not appropriate for repeated", 
                                                   "measures post-hoc tests (Maxwell, 1980; Field, 2012)."))
      }
      rmFactorIndex <- rmFactorIndex + 1
    }

    
    resultPostHoc[['cohenD']] <- resultPostHoc[['t.ratio']] / sqrt(nrow(dataset))

    resultPostHoc[["contrast_A"]] <- lapply(comparisons, function(x) paste(.unv(strsplit(x[[1]], ",")[[1]]), 
                                                                           collapse = ", "))
    resultPostHoc[["contrast_B"]] <- lapply(comparisons, function(x) paste(.unv(strsplit(x[[2]], ",")[[1]]), 
                                                                           collapse = ", "))

    if (options$confidenceIntervalsPostHoc & nrow(resultPostHoc) > 1)
      postHocContainer[[var]]$addFootnote(
        message = gsub(x = attr(resultPostHoc, "mesg")[3], "Conf-level", "Confidence interval"),
        symbol = "<i>Note.</i>")

    avFootnote <- attr(resultPostHoc, "mesg")[1]
    avTerms <- .unv(strsplit(gsub(avFootnote, pattern = "Results are averaged over the levels of: ", replacement = ""), 
                             ", ")[[1]])
    
    postHocContainer[[var]]$addFootnote(
      message = paste0("Results are averaged over the levels of: ", paste(avTerms, collapse = ", ")),
      symbol = "<i>Note.</i>")

    for (pCorrection in c("bonferroni", "scheffe", "tukey", "holm")) {
      if (options$postHocFlagSignificant) {
        for (i in 3:1) {
          signifComparisons <- rownames(resultPostHoc)[resultPostHoc[[pCorrection]] < c(0.05, 0.01, 0.001)[i]]
          if (length(signifComparisons) > 0) {
            colNames <- rep(pCorrection, length(signifComparisons))
            postHocContainer[[var]]$addFootnote(message = "p < .05, ** p < .01, *** p < .001", 
                                                colNames = colNames, rowNames = signifComparisons,
                                                symbol = paste0(rep("*", i), collapse = ""))
          }
        }
      }
    }
    
    postHocContainer[[var]]$setData(resultPostHoc)
  }
  
  return()
}

.createPostHocStandardTable <- function(myTitle, interactionTerm, options, makeBootstrapTable = FALSE) {
  
  preTitle <- if (!makeBootstrapTable) "Post Hoc Comparisons - " else "Bootstrapped Post Hoc Comparisons - "
  postHocTable <- createJaspTable(title = paste0(preTitle, myTitle))
  
  # postHocTable$addColumnInfo(name="contrast", title="Comparison", type="string")
  postHocTable$addColumnInfo(name="contrast_A", title=" ", type="string", combine = TRUE)
  postHocTable$addColumnInfo(name="contrast_B", title=" ", type="string")
  
  postHocTable$addColumnInfo(name="estimate", title="Mean Difference", type="number")
  
  if (options$confidenceIntervalsPostHoc || makeBootstrapTable) {
    
    if (makeBootstrapTable) {
      thisOverTitle <- paste0(options$confidenceIntervalIntervalPostHoc * 100, "% bca\u002A CI")
    } else {
      thisOverTitle <- paste0(options$confidenceIntervalIntervalPostHoc * 100, "% CI for Mean Difference")
    }
    
    postHocTable$addColumnInfo(name="lower.CL", type = "number", title = "Lower",
                               overtitle = thisOverTitle)
    postHocTable$addColumnInfo(name="upper.CL", type = "number", title = "Upper",
                               overtitle = thisOverTitle)
  } 
  
  postHocTable$addColumnInfo(name="SE", type="number")
  
  if (makeBootstrapTable) {
    
    postHocTable$addColumnInfo(name="bias", type="number")
    
  } 
  
  postHocTable$addColumnInfo(name="t.ratio", title="t", type="number")
  
  if (options$postHocTestEffectSize & !interactionTerm) {
    postHocTable$addColumnInfo(name="cohenD", title="Cohen's d", type="number")
    postHocTable$addFootnote("Cohen's d does not correct for multiple comparisons.")
  }
  
  if (options$postHocTestsTukey)
    postHocTable$addColumnInfo(name="tukey", title="p<sub>tukey</sub>", type="number")
  
  if (options$postHocTestsScheffe)
    postHocTable$addColumnInfo(name="scheffe", title="p<sub>scheffe</sub>", type="number")
  
  if (options$postHocTestsBonferroni)
    postHocTable$addColumnInfo(name="bonferroni", title="p<sub>bonf</sub>", type="number")
  
  if (options$postHocTestsHolm)
    postHocTable$addColumnInfo(name="holm",title="p<sub>holm</sub>", type="number")
  
  if(options$postHocTestsSidak)
    postHocTable$addColumnInfo(name="sidak",title="p<sub>sidak</sub>", type="number")
  
  
  postHocTable$showSpecifiedColumnsOnly <- TRUE
  
  return(postHocTable)
}

.resultsContrasts <- function(rmAnovaContainer, dataset, options, ready) {
  if(!ready)
    return()
  
  referenceGrid <- rmAnovaContainer[["referenceGrid"]]$object

  resultsContrasts <- list()
  datasetLong <- .shortToLong(dataset, options$repeatedMeasuresFactors, options$repeatedMeasuresCells, options$betweenSubjectFactors)
  
  contrastTypes <- c("none", "deviation", "simple", "Helmert", "repeated", "difference", "polynomial")
  
  for (contrast in options$contrasts) {
    
    if (! .v(contrast$variable) %in% names(referenceGrid)) {
      next
    }
    
    resultsContrasts[[contrast$variable]] <- list()
    
    column <- datasetLong[[.v(contrast$variable)]]
    
    for(contrastType in contrastTypes) {
      
      contrastMatrix <- .rmAnovaCreateContrast(column, contrastType)
      contrCoef <- lapply(as.data.frame(contrastMatrix), as.vector)
      names(contrCoef) <- .v(.anovaContrastCases(column, contrastType))
      
      if (contrastType == "none") {
        r <- NULL
      } else {
        r <- emmeans::contrast(referenceGrid[[.v(contrast$variable)]], contrCoef)
      }
      resultsContrasts[[contrast$variable]][[contrastType]] <- summary(r)
    }
  }
  
  return(resultsContrasts)
}

.rmAnovaCreateContrast <- function (column, contrast.type) {
  
  levels <- levels(column)
  n.levels <- length(levels)
  
  contr <- NULL
  
  if (contrast.type == "none") {
    
    options(contrasts = c("contr.sum","contr.poly"))
    contr <- NULL
    
  } else if (contrast.type == "deviation") {
    
    contr <- matrix(0,nrow = n.levels, ncol = n.levels - 1)
    for (i in 2:n.levels) {
      contr[,(i-1)] <-  -1 / n.levels
      contr[i,(i-1)] <- (n.levels - 1) / n.levels
    }
    
  } else if (contrast.type == "simple") {
    
    contr <- matrix(0,nrow = n.levels, ncol = n.levels - 1)
    for (i in 1:n.levels-1) {
      contr[c(1,i+1),i]<- c(1,-1) * -1
    }
    
  } else if (contrast.type == "Helmert") {
    
    contr <- contr.helmert(levels) 
    contr <- apply(contr, 2, function(x){ x/max(abs(x))})
    contr <- matrix(rev(contr), ncol = ncol(contr), nrow = nrow(contr))
    
  } else if (contrast.type == "repeated") {
    
    contr <- matrix(0,nrow = n.levels, ncol = n.levels - 1)
    
    for (i in 1:(n.levels-1)) {
      contr[i,i] <- 1
      contr[i+1,i] <- -1
    }
    
  } else if (contrast.type == "difference") {
    
    contr <- contr.helmert(levels) 
    contr <- apply(contr, 2, function(x){ x/max(abs(x))})
    
  } else if (contrast.type == "polynomial") {
    
    contr <- contr.poly(levels)
  }
  
  if (! is.null(contr)) {
    dimnames(contr) <- list(NULL, 1:dim(contr)[2])
  }
  
  contr
}

.rmAnovaContrastTable <- function(rmAnovaContainer, dataset, longDataset, options, ready) {
  if (!is.null(rmAnovaContainer[["contrastContainer"]]) || all(grepl("none", options$contrasts)))
    return()
  
  contrastContainer <- createJaspContainer(title = "Contrast Tables")
  contrastContainer$dependOn(c("contrasts", "contrastAssumeEqualVariance", "confidenceIntervalIntervalContrast", 
                               "confidenceIntervalsContrast"))
  
  createContrastTable <- function(myTitle, options) {
    
    contrastTable <- createJaspTable(title = myTitle)
    contrastTable$addColumnInfo(name = "Comparison", type = "string")
    contrastTable$addColumnInfo(name = "estimate", "Estimate", type = "number")
    
    if (options$confidenceIntervalsContrast) {
      
      thisOverTitle <- paste0(options$confidenceIntervalIntervalContrast * 100, "% CI for Mean Difference")
      contrastTable$addColumnInfo(name="lower.CL", type = "number", title = "Lower",
                                  overtitle = thisOverTitle, )
      contrastTable$addColumnInfo(name="upper.CL", type = "number", title = "Upper",
                                  overtitle = thisOverTitle)
      
    } 
    
    contrastTable$addColumnInfo(name = "SE", type = "number")
    
    dfType <- if (options$contrastAssumeEqualVariance) "integer" else "number"
    contrastTable$addColumnInfo(name = "df", type = dfType)
    contrastTable$addColumnInfo(name = "t.ratio", title = "t", type = "number")
    contrastTable$addColumnInfo(name = "p.value", title = "p", type = "number")
    
    contrastTable$showSpecifiedColumnsOnly <- TRUE
    
    return(contrastTable)
  }
  
  
  
  for (contrast in options$contrasts) {
    
    if (contrast$contrast != "none") {
      
      contrastType <- unlist(strsplit(contrast$contrast, ""))
      contrastType[1] <- toupper(contrastType[1])
      contrastType <- paste0(contrastType, collapse = "")
      
      myTitle <- paste0(contrastType, " Contrast", " - ",  contrast$variable)
      contrastContainer[[paste0(contrast$contrast, "Contrast_",  contrast$variable)]] <- createContrastTable(myTitle, options)
    }
    
  }
  
  rmAnovaContainer[["contrastContainer"]] <- contrastContainer
  
  if (!ready) 
    return()  
  
  referenceGrid <- rmAnovaContainer[["referenceGrid"]]$object
  
  for (contrast in options$contrasts) {
    
    if (contrast$contrast != "none") {
      column <- longDataset[[.v(contrast$variable)]]
      contrastMatrix <- .rmAnovaCreateContrast(column, contrast)
      contrCoef <- lapply(as.data.frame(contrastMatrix), as.vector)
      names(contrCoef) <- .v(.anovaContrastCases(column, contrast$contrast))
      
      r <- emmeans::contrast(referenceGrid[[.v(contrast$variable)]], contrCoef)

      contrastContainer[[paste0(contrast$contrast, "Contrast_",  contrast$variable)]]$addFootnote(
        message = paste0("Results are averaged over the levels of: ", paste(.unv(r@misc$avgd.over), collapse = ", ")),
        symbol = "<i>Note.</i>")
      
      r <- cbind(r, confint(r, level = options$confidenceIntervalIntervalContrast)[,5:6])
      r[["Comparison"]] <- .unv(r[["contrast"]])
      
      # New feature - verify To do!!!!
      if (options$contrastAssumeEqualVariance == FALSE) {
        newDF <- do.call(data.frame, tapply(longDataset[["XZGVwZW5kZW50"]], longDataset[[.v(contrast$variable)]], cbind))
        
        allTestResults <- list()
        
        for (coefIndex in 1:length(contrCoef)) {
          allTestResults[[coefIndex]] <- t.test(as.matrix(newDF) %*% contrCoef[[coefIndex]])
        }
        
        r[["t.ratio"]] <- sapply(allTestResults, function(x) x[["statistic"]])
        r[["df"]] <- sapply(allTestResults, function(x) x[["parameter"]])
        r[["SE"]] <- sapply(allTestResults, function(x) x[["estimate"]] /  x[["statistic"]])
        r[["p.value"]] <- sapply(allTestResults, function(x) x[["p.value"]])
      }
      
      contrastContainer[[paste0(contrast$contrast, "Contrast_",  contrast$variable)]]$setData(r)
    }
  }
 
}

.rmAnovaMarginalMeansTable <- function(dataset, options, perform, status, fullModel = NULL) {
  
  if (is.null(options$marginalMeansTerms))
    return (list(result=NULL, status=status))
  
  
  terms <- options$marginalMeansTerms
  # the following adds automatically interactions of repeated measures terms with between subject terms
  # (a workaround the qml/ui stuff)
  repeatedMeasuresFactors <- sapply(options$repeatedMeasuresFactors, function(fac) fac$name)
  repeatedMeasuresTerms <- sapply(terms, function(term) any(term$components %in% repeatedMeasuresFactors))
  if(any(!repeatedMeasuresTerms)){
    termsInteract <- list()
    for(rmterm in which(repeatedMeasuresTerms)){
      for(bsterm in which(!repeatedMeasuresTerms)){
        termsInteract <- c(termsInteract, list(list(components = c(terms[[bsterm]]$components, terms[[rmterm]]$components))))
      }
    }
    terms <- c(terms, termsInteract)
  }
  terms.base64 <- c()
  terms.normal <- c()
  
  for (term in terms) {
    
    components <- unlist(term)
    term.base64 <- paste(.v(components), collapse=":", sep="")
    term.normal <- paste(components, collapse=" \u273B ", sep="")
    
    terms.base64 <- c(terms.base64, term.base64)
    terms.normal <- c(terms.normal, term.normal)
  }
  
  marginalMeans <- list()
  
  for (i in .indices(terms.base64)) {
    
    result <- list()
    
    result[["title"]] <- paste("Marginal Means - ",terms.normal[i], sep="")
    result[["name"]] <- paste("marginalMeans_",gsub("\u273B","*",gsub(" ", "", terms.normal[i], fixed=TRUE), fixed=TRUE), sep="")
    
    fields <- list()
    
    for(j in .indices(unlist(terms[[i]])))
      fields[[j]] <- list(name=unlist(terms[[i]])[[j]], type="string", combine=TRUE)
    
    fields[[length(fields) + 1]] <- list(name="Marginal Mean", type="number", format="sf:4;dp:3")
    fields[[length(fields) + 1]] <- list(name="SE", type="number", format="sf:4;dp:3")
    fields[[length(fields) + 1]] <- list(name="Lower", type="number", format="sf:4;dp:3", overTitle="95% CI")
    fields[[length(fields) + 1]] <- list(name="Upper", type="number", format="sf:4;dp:3", overTitle="95% CI")
    
    footnotes <- .newFootnotes()
    
    if(options$marginalMeansCompareMainEffects) {
      fields[[length(fields) + 1]] <- list(name="t", type="number", format="sf:4;dp:3")
      fields[[length(fields) + 1]] <- list(name="p", type="number", format="dp:3;p:.001")
      
      if(options$marginalMeansCIAdjustment == "bonferroni") {
        .addFootnote(footnotes, text = "Bonferroni CI adjustment", symbol = "<em>Note.</em>")
      } else if(options$marginalMeansCIAdjustment == "sidak") {
        .addFootnote(footnotes, text = "Sidak CI adjustment", symbol = "<em>Note.</em>")
      }
    }
    
    result[["schema"]] <- list(fields=fields)
    
    termsTemp <- as.vector(terms[[i]])
    
    lvls <- list()
    
    for (variable in unlist(termsTemp)) {
      
      whichRMFactor <- unlist(lapply(options[['repeatedMeasuresFactors']], 
                                     FUN = function(x){x$name == variable}))
      if (any(whichRMFactor)) {
        lvls[[.v(variable)]] <- .v(options$repeatedMeasuresFactors[[which(whichRMFactor == TRUE)]]$levels)
      } else {
        whichBSFactor <- variable %in% options[['betweenSubjectFactors']]
        lvls[[.v(variable)]] <- .v(levels(dataset[[.v(options$betweenSubjectFactors[[which(whichBSFactor == TRUE)]])]]))
      }
      
    }
    
    cases <- rev(expand.grid(rev(lvls)))
    cases <- as.data.frame(apply(cases,2,as.character))
    
    nRows <- dim(cases)[1]
    nCol <- dim(cases)[2]
    
    if (perform == "run" && status$ready && status$error == FALSE)  {
      
      formula <- as.formula(paste("~", terms.base64[i]))
      
      if(options$marginalMeansCIAdjustment == "bonferroni") {
        adjMethod <- "bonferroni"
      } else if(options$marginalMeansCIAdjustment == "sidak") {
        adjMethod <- "sidak"
      } else {
        adjMethod <- "none"
      }
      
      r <- summary(emmeans::lsmeans(fullModel, formula), adjust = adjMethod, infer = c(TRUE,TRUE))
      
      rows <- list()
      
      for(k in 1:nRows) {
        
        row <- list()
        
        for(j in 1:nCol) {
          row[[ .unv(colnames(cases)[j]) ]] <- .unv(cases[k,j])
        }
        
        if(nCol > 1) {
          index <- apply(r[,1:nCol], 1, function(x) all(x==cases[k,]))
        } else {
          index <- k
        }
        
        row[["Marginal Mean"]] <- .clean(r$lsmean[index])
        row[["SE"]] <- .clean(r$SE[index])
        row[["Lower"]] <- .clean(r$lower.CL[index])
        row[["Upper"]] <- .clean(r$upper.CL[index])
        
        if(options$marginalMeansCompareMainEffects) {
          row[["t"]] <- .clean(r$t.ratio[index])
          row[["p"]] <- .clean(r$p.value[index])
        }
        
        if(cases[k,nCol] == lvls[[ nCol ]][1]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[k]] <- row
        
      }
      
      result[["data"]] <- rows
      result[["status"]] <- "complete"
      
    } else {
      
      rows <- list()
      
      for(k in 1:nRows) {
        
        row <- list()
        
        for(j in 1:nCol)
          row[[ .unv(colnames(cases)[j]) ]] <- .unv(cases[k,j])
        
        row[["Marginal Mean"]] <- "."
        row[["SE"]] <- "."
        row[["Lower"]] <- "."
        row[["Upper"]] <- "."
        
        if(options$marginalMeansCompareMainEffects) {
          row[["t"]] <- "."
          row[["p"]] <- "."
        }
        
        if(cases[k,nCol] == lvls[[ nCol ]][1]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[k]] <- row
        
      }
      
      result[["data"]] <- rows
    }
    
    result[["footnotes"]] <- as.list(footnotes)
    
    if (status$error)
      result[["error"]] <- list(error="badData")
    
    marginalMeans[[i]] <- result
    
  }
  
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    stateMarginalMeans <- marginalMeans
    
  } else {
    
    stateMarginalMeans <- NULL
    
  }
  list(result=marginalMeans, status=status, stateMarginalMeans=stateMarginalMeans)
}

.rmAnovaMarginalMeansBootstrappingTable <- function(dataset, options, perform, status, fullModel = NULL) {
  
  if (is.null(options$marginalMeansTerms))
    return (list(result=NULL, status=status))
  
  terms <- options$marginalMeansTerms
  # the following adds automatically interactions of repeated measures terms with between subject terms
  # (a workaround the qml/ui stuff)
  repeatedMeasuresFactors <- sapply(options$repeatedMeasuresFactors, function(fac) fac$name)
  repeatedMeasuresTerms <- sapply(terms, function(term) any(term$components %in% repeatedMeasuresFactors))
  if(any(!repeatedMeasuresTerms)){
    termsInteract <- list()
    for(rmterm in which(repeatedMeasuresTerms)){
      for(bsterm in which(!repeatedMeasuresTerms)){
        termsInteract <- c(termsInteract, list(list(components = c(terms[[bsterm]]$components, terms[[rmterm]]$components))))
      }
    }
    terms <- c(terms, termsInteract)
  }
  terms.base64 <- c()
  terms.normal <- c()
  
  for (term in terms) {
    
    components <- unlist(term)
    term.base64 <- paste(.v(components), collapse=":", sep="")
    term.normal <- paste(components, collapse=" \u273B ", sep="")
    
    terms.base64 <- c(terms.base64, term.base64)
    terms.normal <- c(terms.normal, term.normal)
  }
  
  marginalMeans <- list()
  
  if(length(terms.base64) > 0 && perform == "run" && status$ready && status$error == FALSE){
    ticks <- options[['marginalMeansBootstrappingReplicates']] * length(terms.base64)
    progress <- .newProgressbar(ticks = ticks, callback = callback, response = TRUE)
  }
  
  for (i in .indices(terms.base64)) {
    
    result <- list()
    
    result[["title"]] <- paste("Bootstrapped Marginal Means - ",terms.normal[i], sep="")
    result[["name"]] <- paste("marginalMeans_",gsub("\u273B","*",gsub(" ", "", terms.normal[i], fixed=TRUE), fixed=TRUE), sep="")
    
    fields <- list()
    
    for(j in .indices(unlist(terms[[i]])))
      fields[[j]] <- list(name=unlist(terms[[i]])[[j]], type="string", combine=TRUE)
    
    fields[[length(fields) + 1]] <- list(name="Marginal Mean", type="number", format="sf:4;dp:3")
    fields[[length(fields) + 1]] <- list(name="Bias", type="number", format="sf:4;dp:3")
    fields[[length(fields) + 1]] <- list(name="SE", type="number", format="sf:4;dp:3")
    fields[[length(fields) + 1]] <- list(name="Lower", type="number", format="sf:4;dp:3", overTitle="95% bca\u002A CI")
    fields[[length(fields) + 1]] <- list(name="Upper", type="number", format="sf:4;dp:3", overTitle="95% bca\u002A CI")
    
    footnotes <- .newFootnotes()
    
    .addFootnote(footnotes, symbol = "<em>Note.</em>",
                 text = paste0("Bootstrapping based on ", options[['marginalMeansBootstrappingReplicates']], " replicates."))
    .addFootnote(footnotes, symbol = "<em>Note.</em>",
                 text = "Marginal Means estimate is based on the median of the bootstrap distribution.")
    .addFootnote(footnotes, symbol = "\u002A",
                 text = "Bias corrected accelerated.")
    
    result[["schema"]] <- list(fields=fields)
    
    termsTemp <- as.vector(terms[[i]])
    
    lvls <- list()
    
    for (variable in unlist(termsTemp)) {
      
      whichRMFactor <- unlist(lapply(options[['repeatedMeasuresFactors']], 
                                     FUN = function(x){x$name == variable}))
      if (any(whichRMFactor)) {
        lvls[[.v(variable)]] <- .v(options$repeatedMeasuresFactors[[which(whichRMFactor == TRUE)]]$levels)
      } else {
        whichBSFactor <- variable %in% options[['betweenSubjectFactors']]
        lvls[[.v(variable)]] <- .v(levels(dataset[[.v(options$betweenSubjectFactors[[which(whichBSFactor == TRUE)]])]]))
      }
      
    }
    
    cases <- rev(expand.grid(rev(lvls)))
    cases <- as.data.frame(apply(cases,2,as.character))
    
    nRows <- dim(cases)[1]
    nCol <- dim(cases)[2]
    
    if (perform == "run" && status$ready && status$error == FALSE)  {
      
      formula <- as.formula(paste("~", terms.base64[i]))
      
      .bootstrapMarginalMeans <- function(data, indices, options){
        pr <- progress()
        response <- .optionsDiffCheckBootstrapRMAnovaMarginalMeans(pr, options)
        
        if(response$status == "changed" || response$status == "aborted")
          stop("Bootstrapping options have changed")
        
        resamples <- data[indices, , drop=FALSE]
        
        anovaModelBoots <- .rmAnovaComputeResults(resamples, options, TRUE) # refit model
        
        modelBoots <- anovaModelBoots$fullModel
        singularBoots <- anovaModelBoots$singular
        r <- suppressMessages( # to remove clutter
          summary(emmeans::lsmeans(modelBoots, formula), infer = c(FALSE,FALSE))
        )
        
        if(length(r$lsmean) == nRows){ # ensure that the bootstrap has all levels
          return(r$lsmean)
        } else {
          return(rep(NA, nRows))
        }
      }
      
      bootstrapMarginalMeans <- try(boot::boot(data = dataset, statistic = .bootstrapMarginalMeans, 
                                               R = options[["marginalMeansBootstrappingReplicates"]],
                                               options = options), silent = TRUE)
      if(inherits(bootstrapMarginalMeans, "try-error") && 
         identical(attr(bootstrapMarginalMeans, "condition")$message, "Bootstrapping options have changed"))
        return("Bootstrapping options have changed")
      
      bootstrapMarginalMeans.summary <- summary(bootstrapMarginalMeans)
      ci.fails <- FALSE
      bootstrapMarginalMeans.ci <- t(sapply(1:nrow(bootstrapMarginalMeans.summary), function(index){
        res <- try(boot::boot.ci(boot.out = bootstrapMarginalMeans, conf = 0.95, type = "bca",
                                 index = index)[['bca']][1,4:5])
        if(!inherits(res, "try-error")){
          return(res)
        } else if(identical(attr(res, "condition")$message, "estimated adjustment 'a' is NA")){
          ci.fails <<- TRUE
          return(c(NA, NA))
        } else{
          return(res)
        }
      }))
      
      if(ci.fails){
        .addFootnote(footnotes,
                     symbol = "<i>Note.</i>", 
                     text = "Some confidence intervals could not be computed. Possibly too few bootstrap replicates.")
      }
      
      bootstrapMarginalMeans.summary[,"lower.CL"] <- bootstrapMarginalMeans.ci[,1]
      bootstrapMarginalMeans.summary[,"upper.CL"] <- bootstrapMarginalMeans.ci[,2]
      
      # the next chunk of code ensures that the rows in bootstrap
      # table are in the same order as the rows in object cases
      getModelCases <- summary(emmeans::lsmeans(fullModel, formula), infer = c(FALSE,FALSE))
      getModelCases <- getModelCases[,names(cases), drop = FALSE]
      names(getModelCases) <- .unv(names(getModelCases))
      r <- as.data.frame(bootstrapMarginalMeans.summary)
      r <- cbind(getModelCases, r)
      
      rows <- list()
      
      for(k in 1:nRows) {
        
        row <- list()
        
        for(j in 1:nCol) {
          row[[ .unv(colnames(cases)[j]) ]] <- .unv(cases[k,j])
        }
        
        if(nCol > 1) {
          index <- apply(r[,1:nCol], 1, function(x) all(x==cases[k,]))
        } else {
          index <- k
        }
        
        row[["Marginal Mean"]] <- .clean(r$bootMed[index])
        row[["Bias"]] <- .clean(r$bootBias[index])
        row[["SE"]] <- .clean(r$bootSE[index])
        row[["Lower"]] <- .clean(r$lower.CL[index])
        row[["Upper"]] <- .clean(r$upper.CL[index])
        
        
        if(cases[k,nCol] == lvls[[ nCol ]][1]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[k]] <- row
        
      }
      
      result[["data"]] <- rows
      result[["status"]] <- "complete"
      
    } else {
      
      rows <- list()
      
      for(k in 1:nRows) {
        
        row <- list()
        
        for(j in 1:nCol)
          row[[ .unv(colnames(cases)[j]) ]] <- .unv(cases[k,j])
        
        row[["Marginal Mean"]] <- "."
        row[["Bias"]] <- "."
        row[["SE"]] <- "."
        row[["Lower"]] <- "."
        row[["Upper"]] <- "."
        
        if(cases[k,nCol] == lvls[[ nCol ]][1]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[k]] <- row
        
      }
      
      result[["data"]] <- rows
    }
    
    result[["footnotes"]] <- as.list(footnotes)
    
    if (status$error)
      result[["error"]] <- list(error="badData")
    
    marginalMeans[[i]] <- result
    
  }
  
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    stateMarginalMeans <- marginalMeans
    
  } else {
    
    stateMarginalMeans <- NULL
    
  }
  list(result=marginalMeans, status=status, stateMarginalMeansBoots=stateMarginalMeans)
}

.optionsDiffCheckBootstrapRMAnovaMarginalMeans <- function(response, options) {
  if(response$status == "changed"){
    change <- .diff(options, response$options)
    
    if(change$repeatedMeasuresCells || change$covariates || change$betweenSubjectFactors ||
       change$withinModelTerms || change$marginalMeansTerms ||
       change$marginalMeansBootstrapping || change$marginalMeansBootstrappingReplicates)
      return(response)
    
    response$status <- "ok"
  }
  
  return(response)
}

.rmAnovaSimpleEffects <- function(dataset, options, perform, fullModel, fullAnovaTableWithin, 
                                  fullAnovaTableBetween, status, singular, stateSimpleEffects) {
  
  if (identical(options$simpleFactor, "") | identical(options$moderatorFactorOne, "")) {
    return (list(result=NULL, status=status))
  }
  
  
  terms <- c(options$moderatorFactorOne,options$moderatorFactorTwo)
  terms.base64 <- c()
  terms.normal <- c()
  simpleFactor.base64 <- .v(options[['simpleFactor']])
  simpleFactor <- options[['simpleFactor']]
  moderatorFactorOne <- options[['moderatorFactorOne']]
  moderatorFactorTwo <- options[['moderatorFactorTwo']]
  
  for (term in terms) {
    
    components <- unlist(term)
    term.base64 <- paste(.v(components), collapse=":", sep="")
    term.normal <- paste(components, collapse=" \u273B ", sep="")
    
    terms.base64 <- c(terms.base64, term.base64)
    terms.normal <- c(terms.normal, term.normal)
  }
  
  simpleEffectsTable <- list()
  simpleEffectsTable[["title"]] <- paste("Simple Main Effects - ", simpleFactor, sep = "")
  
  fields <- list(
    list(name="ModOne", type="string", combine = TRUE, title = paste0("Level of ", terms.normal[1])),
    list(name="ModTwo", type="string", combine = TRUE, title = paste0("Level of ", terms.normal[2])),
    list(name="SumSquares", type="number", format="sf:4;dp:3", title = "Sum of Squares"),
    list(name="df", type="integer", title = "df"),
    list(name="MeanSquare", type="number", format="sf:4;dp:3", title = "Mean Square"),
    list(name="F", type="number", format="sf:4;dp:3", title = "F"),
    list(name="p", type="number", format="dp:3;p:.001", title = "p"))
  
  # If there is no second moderator, so remove from table:
  if (identical(options$moderatorFactorTwo, ""))
    fields <- fields[-2]
  
  footnotes <- .newFootnotes()
  
  simpleEffectsTable[["schema"]] <- list(fields=fields)
  
  if (perform == "run" && status$ready && status$error == FALSE && !isTryError(fullModel))  {
    
    isMixedAnova <-   length(options[['betweenSubjectFactors']]) > 0
    isSimpleFactorWithin <- !simpleFactor %in% unlist(options[['betweenModelTerms']] )
    isModeratorOneWithin <- !moderatorFactorOne %in% unlist(options[['betweenModelTerms']] )
    isModeratorTwoWithin <- !moderatorFactorTwo %in% unlist(options[['betweenModelTerms']] )
    errorIndexWithin <- length(fullAnovaTableWithin$data)
    errorIndexBetween <-  length(fullAnovaTableBetween$data)
    ssWithin <- fullAnovaTableWithin$data[[errorIndexWithin]]$SS
    ssBetween <- fullAnovaTableBetween$data[[errorIndexBetween]]$SS
    dfWithin <- fullAnovaTableWithin$data[[errorIndexWithin]]$df
    dfBetween <- fullAnovaTableBetween$data[[errorIndexBetween]]$df
    fullAnovaTable <- fullAnovaTableWithin
    tableCounter <- 1
    if(isMixedAnova & isSimpleFactorWithin){
      fullAnovaMS <- ssWithin / dfWithin
      fullAnovaDf <- dfWithin
    } else if(isMixedAnova & !isSimpleFactorWithin){
      fullAnovaMS <- ssBetween / dfBetween
      fullAnovaDf <- dfBetween
    } else {
      fullAnovaMS <- fullAnovaTable$data[[length(fullAnovaTable$data)]]$MS 
      fullAnovaDf <- fullAnovaTable$data[[length(fullAnovaTable$data)]]$df
    }
    
    
    simpleEffectRows <- list()
    rows <- list()
    
    termsBothModerators <- as.vector(c(moderatorFactorOne, moderatorFactorTwo))
    if(identical(moderatorFactorTwo, "")) {
      termsBothModerators <- termsBothModerators[1]
    }
    
    lvls <- list()
    factors <- list()
    
    for (variable in termsBothModerators) {
      if (variable %in% unlist(options[['withinModelTerms']])) {
        whichFactor <- unlist(lapply(options[['repeatedMeasuresFactors']], 
                                     FUN = function(x){x$name == variable}))
        lvls[[variable]] <- options$repeatedMeasuresFactors[[which(whichFactor == TRUE)]]$levels 
      } else if (variable %in% unlist(options[['betweenModelTerms']])) {
        thisFactor <- dataset[[ .v(variable) ]]
        factors[[length(factors)+1]] <- thisFactor
        lvls[[variable]] <- levels(thisFactor)
      }
    }
    #
    allNames <- unlist(lapply(options[['repeatedMeasuresFactors']], function(x) x$name)) # Factornames 
    
    # make separate covariate dataframe to avoid mismatching dataframe names
    covDataset <- dataset[.v(options[['covariates']])]
    wideDataset <- fullModel$data$wide[,(-1)]
    
    covariatesInModel <- ifelse(length(options[['covariates']]) > 0, TRUE, FALSE)
    if (covariatesInModel) {
      dataset <- dataset[, names(dataset) != (.v(options[['covariates']]))]
      wideDataset <- wideDataset[ -match(covDataset, wideDataset)]
    }
    # Following steps to make sure column ordering is the same for the datasets
    
    factorNamesV <- colnames(wideDataset)
    withinFactorNamesV <- factorNamesV[!(.unv(factorNamesV) %in% options[['betweenSubjectFactors']])] 
    
    withinOrder <- match(wideDataset, dataset)
    betweenOrder <- match(names(wideDataset), names(dataset))
    betweenOrder[is.na(betweenOrder)] <- withinOrder[!is.na(withinOrder)]
    fullOrder <- betweenOrder
    
    dataset <- dataset[fullOrder]
    # orderOfTerms <- unlist(options[['withinModelTerms']][[length(options[['withinModelTerms']])]]$components)
    orderOfTerms <- unlist(lapply(options$repeatedMeasuresFactors, function(x) x$name))
    
    indexofOrderFactors <- match(allNames,orderOfTerms)
    
    for (level in lvls[[1]]) {
      # For each level of the first moderator factor, take a subset of the dataset, and adjust the options object
      # Suboptions is the same as options, except that the first moderator factor has been removed as a predictor 
      # (because each subset only has one level of that factor). The same procedure is applied to the second moderator, if specified.
      
      subOptions <- options
      # subOptions[['covariates']] <- list()
      # Prepare the options object to handle the contional dataset
      # Based on whether the first moderator variable is within or between
      if (isModeratorOneWithin) {
        rmFactorIndex <- which(lapply(options[['repeatedMeasuresFactors']], 
                                      FUN = function(x){x$name == moderatorFactorOne}) == TRUE)
        splitNames <- unlist(lapply(strsplit(factorNamesV,  split = "_"), 
                                    FUN = function(x) x[indexofOrderFactors[rmFactorIndex]]))
        splitWithinNames <- unlist(lapply(strsplit(withinFactorNamesV,  split = "_"), 
                                          FUN = function(x) x[indexofOrderFactors[rmFactorIndex]]))
        subDataset <- dataset[, ((splitNames %in% .v(level)) | (names(dataset)) %in% .v(options[['betweenSubjectFactors']]))]
        subCovDataset <- covDataset
        subFactorNamesV <- factorNamesV[((splitNames %in% .v(level)) | (names(dataset)) %in% .v(options[['betweenSubjectFactors']]))]
        whichFactorsBesidesModerator <- !unlist(lapply((options[['withinModelTerms']]), FUN = function(x){grepl(moderatorFactorOne, x)}))
        subOptions[['withinModelTerms']] <- options[['withinModelTerms']][whichFactorsBesidesModerator]
        subOptions[['repeatedMeasuresFactors']] <- options[['repeatedMeasuresFactors']][-rmFactorIndex]
        subOptions[['repeatedMeasuresCells']]<- options$repeatedMeasuresCells[(splitWithinNames %in% .v(level))]
      } else {
        subDataset <- subset(dataset, dataset[terms.base64[1]] == level)
        if (covariatesInModel) {
          subCovDataset  <- subset(covDataset, dataset[terms.base64[1]] == level)
        } else {
          subCovDataset <- covDataset[dataset[terms.base64[1]] == level,] 
        }
        subFactorNamesV <- factorNamesV
        whichTermsBesidesModerator <- !unlist(lapply((options[['betweenModelTerms']]), FUN = function(x){grepl(moderatorFactorOne, x)}))
        whichFactorsBesidesModerator <- !unlist(lapply((options[['betweenSubjectFactors']]), FUN = function(x){grepl(moderatorFactorOne, x)}))
        
        subOptions[['betweenModelTerms']] <- options[['betweenModelTerms']][whichTermsBesidesModerator]
        subOptions[['betweenSubjectFactors']] <- options[['betweenSubjectFactors']][whichFactorsBesidesModerator]
      }
      areSimpleFactorCellsDropped <- ifelse(isSimpleFactorWithin, FALSE, (nrow(unique(subDataset[simpleFactor.base64])) <  
                                                                            nrow(unique(dataset[simpleFactor.base64]))))
      model <- NULL
      singular <- FALSE
      if (identical(moderatorFactorTwo, "")) {
        # This means there is only one moderator variable (i.e., one factor on which to condition)
        newGroup <- ifelse( level == lvls[[1]][1], TRUE, FALSE )
        if (nrow(subDataset) < 3 || areSimpleFactorCellsDropped ) {
          row <- list("ModOne"=level, "SumSquares"=".", "df"=".", "MeanSquare"=".", "F"=".", "p"=".", ".isNewGroup" = newGroup)
          .addFootnote(footnotes, text = paste0("Not enough observations in cell ", level, " of ", subOptions$moderatorFactorOne), 
                       symbol = "<em>Note.</em>")
        } else {
          if (length(subOptions[['repeatedMeasuresFactors']]) > 0) {
            # There are still multiple levels of RM factors, so proceed with conditional RM ANOVA
            anovaModel <- .rmAnovaComputeResults(cbind(subDataset, subCovDataset), subOptions, status = status)
            modelSummary <- anovaModel$model
            if (subOptions$sumOfSquares != "type1") {
              modOneIndex <- which(row.names(modelSummary) == .v(simpleFactor))
              df <- modelSummary[modOneIndex,'num Df']
              SS <- modelSummary[modOneIndex,'Sum Sq']  
              if (!options$poolErrorTermSimpleEffects) {
                fullAnovaMS <- modelSummary[modOneIndex,'Error SS'] / modelSummary[modOneIndex,'den Df']
                fullAnovaDf <- modelSummary[modOneIndex,'den Df']
              }
            } else {
              modelSummary <- modelSummary[[-1]][[1]]
              modOneIndex <- which(row.names(modelSummary) == .v(simpleFactor))
              df <- modelSummary[modOneIndex,'Df']
              SS <- modelSummary[modOneIndex,'Sum Sq']   
              if (!options$poolErrorTermSimpleEffects) {
                fullAnovaMS <- modelSummary['Residuals','Mean Sq']
                fullAnovaDf <- modelSummary['Residuals','Df']
              }
            }
            
          } else {
            # There is only one level of RM factor left, so proceed with conditional simple ANOVA
            subOptionsSimpleAnova <- subOptions
            subOptionsSimpleAnova['fixedFactors']  <- list(subOptions[['betweenSubjectFactors']])
            subOptionsSimpleAnova['modelTerms'] <- list(subOptions[['betweenModelTerms']])
            subOptionsSimpleAnova['dependent'] <-  options$repeatedMeasuresCells[splitWithinNames %in% .v(level)]
            
            anovaModel <- .anovaModel(cbind(subDataset, subCovDataset), options = subOptionsSimpleAnova)
            model <- anovaModel$model
            modelSummary <- summary(model)[[1]]
            modOneIndex <- which(subOptionsSimpleAnova[['fixedFactors']] == simpleFactor)
            df <- modelSummary$Df[modOneIndex]
            SS <- modelSummary$`Sum Sq`[modOneIndex]
            if (!options$poolErrorTermSimpleEffects) {
              fullAnovaMS <- modelSummary$`Mean Sq`[length(modelSummary$`Mean Sq`)]
              fullAnovaDf <- modelSummary$Df[length(modelSummary$Df)]
            }
          }
          MS <- SS / df
          F <- MS / fullAnovaMS
          p <- pf(F, df, fullAnovaDf, lower.tail = FALSE)
          row <- list("ModOne"=level, "SumSquares"=SS, "df"=df, "MeanSquare"=MS, "F"=F, "p"=p, ".isNewGroup" = newGroup)
        }
        simpleEffectRows[[length(simpleEffectRows) + 1]] <- row
      } else {
        # This means there are two moderator variables (i.e., two factors on which to condition)
        for (levelTwo in lvls[[2]]) {
          newGroup <- ifelse( levelTwo == lvls[[2]][1], TRUE, FALSE )
          # Prepare the options object to handle the contional dataset
          # Based on whether the second moderator variable is within or between
          if (isModeratorTwoWithin) {
            rmFactorIndex <- which(lapply(options[['repeatedMeasuresFactors']], 
                                          FUN = function(x){x$name == moderatorFactorTwo}) == TRUE)
            splitNames <- unlist(lapply(strsplit(subFactorNamesV,  split = "_"), 
                                        FUN = function(x) x[indexofOrderFactors[rmFactorIndex]]))
            subWithinFactorNamesV <- subFactorNamesV[ !(names(subDataset) %in% .v(options[['betweenSubjectFactors']]))]
            splitWithinNames <- unlist(lapply(strsplit(subWithinFactorNamesV,  split = "_"), 
                                              FUN = function(x) x[indexofOrderFactors[rmFactorIndex]]))
            subSubDataset <- subDataset[, (splitNames %in% .v(levelTwo)) | (names(subDataset) %in% .v(subOptions[['betweenSubjectFactors']]))]
            subSubCovDataset <- subCovDataset
            subSubOptions <- subOptions
            whichFactorsBesidesModerator <- !unlist(lapply((subOptions[['withinModelTerms']]), FUN = function(x){grepl(moderatorFactorTwo, x)}))
            subSubOptions[['withinModelTerms']] <- subOptions[['withinModelTerms']][whichFactorsBesidesModerator]
            whichFactorToRemove <- which(lapply(subOptions[['repeatedMeasuresFactors']], 
                                                FUN = function(x){x$name == moderatorFactorTwo}) == TRUE)
            subSubOptions[['repeatedMeasuresFactors']] <- subOptions[['repeatedMeasuresFactors']][-whichFactorToRemove]
            subSubOptions[['repeatedMeasuresCells']]<- subOptions$repeatedMeasuresCells[(splitWithinNames %in% .v(levelTwo))]
            
          } else {
            subSubDataset <- subset(subDataset, subDataset[terms.base64[2]] == levelTwo)
            if (covariatesInModel) {
              subSubCovDataset  <- subset(subCovDataset, subDataset[terms.base64[2]] == levelTwo)
            } else {
              subSubCovDataset <- subCovDataset[subDataset[terms.base64[2]] == levelTwo,] 
            }
            subSubOptions <- subOptions
            whichFactorsBesidesModerator <- !unlist(lapply((subOptions[['betweenModelTerms']]), FUN = function(x){grepl(moderatorFactorTwo, x)}))
            subSubOptions[['betweenModelTerms']] <- subOptions[['betweenModelTerms']][whichFactorsBesidesModerator]
          }
          
          areSimpleFactorCellsDropped <- ifelse(isSimpleFactorWithin, FALSE, (nrow(unique(subSubDataset[simpleFactor.base64])) <  
                                                                                nrow(unique(dataset[simpleFactor.base64]))))
          
          if (nrow(subSubDataset) < 3 || areSimpleFactorCellsDropped) {
            row <- list("ModOne"=level, "ModTwo" = levelTwo, "SumSquares"=".", "df"=".", "MeanSquare"=".", "F"=".", "p"=".", ".isNewGroup" = newGroup)
            .addFootnote(footnotes, text = paste0("Not enough observations in cell ",level, " of ", moderatorFactorOne, " and ",
                                                  levelTwo," of ", moderatorFactorTwo), symbol = "<em>Note.</em>")
          } else {
            if (length(subSubOptions[['repeatedMeasuresFactors']]) > 0) {
              # There are still multiple levels of RM factors, so proceed with conditional RM ANOVA
              anovaModel <- .rmAnovaComputeResults(cbind(subSubDataset, subSubCovDataset), subSubOptions, status = status)
              modelSummary <- anovaModel$model
              if (subSubOptions$sumOfSquares != "type1") {
                modTwoIndex <- which(row.names(modelSummary) == .v(simpleFactor))
                df <- modelSummary[modTwoIndex,'num Df']
                SS <- modelSummary[modTwoIndex,'Sum Sq']   
                if (!options$poolErrorTermSimpleEffects) {
                  fullAnovaMS <- modelSummary[modTwoIndex,'Error SS'] / modelSummary[modTwoIndex,'den Df']
                  fullAnovaDf <- modelSummary[modTwoIndex,'den Df']
                }
              } else {
                modelSummary <- modelSummary[[-1]][[1]]
                modTwoIndex <- which(row.names(modelSummary) == .v(simpleFactor))
                df <- modelSummary[modTwoIndex,'Df']
                SS <- modelSummary[modTwoIndex,'Sum Sq']   
                if (!options$poolErrorTermSimpleEffects) {
                  fullAnovaMS <- modelSummary['Residuals','Mean Sq']
                  fullAnovaDf <- modelSummary['Residuals','Df']
                }
              }
              
            } else {
              # There is only one level of RM factor left, so proceed with conditional simple ANOVA
              subSubOptionsSimpleAnova <- subSubOptions
              subSubOptionsSimpleAnova['fixedFactors']  <- subSubOptions[['betweenSubjectFactors']]
              subSubOptionsSimpleAnova['modelTerms'] <- list(subSubOptions[['betweenModelTerms']])
              subSubOptionsSimpleAnova['dependent'] <-  subSubOptions[['repeatedMeasuresCells']]
              anovaModel <- .anovaModel(cbind(subSubDataset, subSubCovDataset), options = subSubOptionsSimpleAnova)
              model <- anovaModel$model
              modelSummary <- summary(model)[[1]]
              modTwoIndex <- which(subSubOptionsSimpleAnova[['fixedFactors']] == simpleFactor)
              df <- modelSummary$Df[modTwoIndex]
              SS <- modelSummary$`Sum Sq`[modTwoIndex]
              if (!options$poolErrorTermSimpleEffects) {
                fullAnovaMS <- modelSummary$`Mean Sq`[length(modelSummary$`Mean Sq`)]
                fullAnovaDf <- modelSummary$Df[length(modelSummary$Df)]
              }
            }
            
            MS <- SS / df
            F <- MS / fullAnovaMS
            p <- pf(F, df, fullAnovaDf, lower.tail = FALSE)
            row <- list("ModOne"=level, "ModTwo" = levelTwo, "SumSquares"=SS, "df"=df, "MeanSquare"=MS, "F"=F, "p"=p, ".isNewGroup" = newGroup)
          }
          simpleEffectRows[[length(simpleEffectRows) + 1]] <- row
        }
      }
      
      
    }
    
    simpleEffectsTable[["data"]] <- simpleEffectRows
    
    if (options$sumOfSquares == "type1") {
      
      .addFootnote(footnotes, text = "Type I Sum of Squares", symbol = "<em>Note.</em>")
      
    } else if (options$sumOfSquares == "type2") {
      
      .addFootnote(footnotes, text = "Type II Sum of Squares", symbol = "<em>Note.</em>")
      
    } else if (options$sumOfSquares == "type3") {
      
      .addFootnote(footnotes, text = "Type III Sum of Squares", symbol = "<em>Note.</em>")
      
    }
    
  } else {
    if(options$sumOfSquares == "type1" ) {
      .addFootnote(footnotes, text = "Simple effects not yet available for type 1 SS.", 
                   symbol = "<em>Note.</em>")  }
    simpleEffectsTable[["data"]]  <- list(list("ModOne"=terms.normal, "SumSquares"=".", "df"=".", "MeanSquare"=".", "F"=".", "p"=".", ".isNewGroup" = TRUE))
  }
  
  simpleEffectsTable[["footnotes"]] <- as.list(footnotes)
  simpleEffectsTable[["status"]] <- "complete"
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    stateSimpleEffects <- simpleEffectsTable
    
  } else {
    
    stateSimpleEffects <- NULL
    
  }
  simpleEffectsTable[["citation"]] <- list(
    "Howell, D. C. (2002). Statistical Methods for Psychology (8th. ed.). Pacific Grove, CA: Duxberry. "
  )
  
  list(result=simpleEffectsTable, status=status, stateSimpleEffects=stateSimpleEffects)
}

.rmAnovaFriedman <- function(dataset, fullModel, options, perform, status, singular, stateFriedman) {
  
  if (length(options$friedmanWithinFactor) == 0)
    return (list(result=NULL, status=status))
  
  withinTerms <- options$friedmanWithinFactor
  betweenTerm <- options$friedmanBetweenFactor
  
  withinTerms.base64 <- .v(withinTerms)
  betweenTerms.base64 <- .v(betweenTerm)
  
  result <- list()
  
  if( any(!(withinTerms %in% unlist(options$withinModelTerms))) | 
      (betweenTerm %in% unlist(options$withinModelTerms)) ) {
    status$error <- TRUE
    status$errorMessage <- "Please specify appropriate terms for the Friedman/Durbin test."
    result[["error"]] <- list(errorType="badData", errorMessage=status$errorMessage)
  }
  
  result[["title"]] <- paste("Friedman Test")
  result[["name"]] <- paste("friedmanTable")
  
  fields <- list()
  fields[[length(fields) + 1]] <- list(name="Factor", type="string")
  fields[[length(fields) + 1]] <- list(name="Chi-Squared", type="number", format="sf:4;dp:3")
  fields[[length(fields) + 1]] <- list(name="df", type="integer")
  fields[[length(fields) + 1]] <- list(name="p", type="number", format="dp:3;p:.001")
  fields[[length(fields) + 1]] <- list(name="Kendall's W", type="number", format="sf:4;dp:3")
  # fields[[length(fields) + 1]] <- list(name="Kendall's W corr.", type="number", format="sf:4;dp:3")
  
  footnotes <- .newFootnotes()
  
  rows <- list()
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    longData <- fullModel$data$long
    
    if (identical(betweenTerm, "")) {
      betweenTerms.base64 <- 'subject'
    }
    
    for (i in 1:length(withinTerms)) {
      
      groups <- as.factor(longData[, withinTerms.base64[i]])
      blocks <- as.factor(longData[, betweenTerms.base64])
      y <- longData[, 'dependent']
      
      useDurbin <- any(table(groups, blocks) != 1)
      
      t <- nlevels(groups)
      b <- nlevels(blocks)
      r <- unique(table(groups))
      k <- unique(table(blocks))
      
      
      if (length(r) == 1 & length(k) == 1) {
        rankPerBlock <- unlist(tapply(y, blocks, rank))
        rankPerGroup <- unlist(tapply(y, groups, rank))    
        
        rankJ <- tapply(rankPerBlock, groups, sum)    
        
        sumRanks <- sum(rankPerBlock^2)
        cVal <- (b * k * (k + 1)^2) / 4
        dVal <- sum(rankJ^2) - r * cVal
        testStatOne <- (t - 1) / (sumRanks - cVal) * dVal
        testStatTwo <- (testStatOne / (t - 1)) / ((b * k - b - testStatOne) / (b * k - b - t + 1))
        
        ## Code from PMCMRplus
        dfChi <- t - 1
        dfOneF <- k - 1
        dfTwoF <- b * k - b - t + 1 
        pValOne <- pchisq(testStatOne, dfChi, lower.tail = FALSE)
        pValTwo <- pf(testStatTwo, dfOneF, dfTwoF, lower.tail = FALSE)
        
        # Kendall W
        rankMatrixRM <- matrix(rankPerGroup, ncol = t)
        rowSumsMatrix <- rowSums(rankMatrixRM)
        nTies <- unlist(apply(rankMatrixRM, 2, function(x) {
          tab <- table(x)
          tab[tab > 1] }))
        nTies <- sum(nTies^3 - nTies)
        kendallW <- (sum(rowSumsMatrix^2) - sum(rowSumsMatrix)^2 / b) / (t^2 * (b^3 - b) / 12)
        kendallWcor <-(sum(rowSumsMatrix^2) - sum(rowSumsMatrix)^2 / b) / ((t^2 * (b^3 - b) - t * nTies) / 12)
        
        row <- list()
        
        row[["Factor"]] <- withinTerms[i]
        row[["Chi-Squared"]] <- .clean(testStatOne)
        row[["df"]] <- .clean(dfChi)
        row[["p"]] <- .clean(pValOne)
        row[["Kendall's W"]] <- .clean(kendallWcor)
        # row[["Kendall's W corr."]] <- .clean(kendallWcor)
        
        
        if (useDurbin) {
          result[["title"]] <- "Durbin Test"
          
          row[["F"]] <- .clean(testStatTwo)
          row[["df num"]] <- .clean(dfOneF)
          row[["df den"]] <- .clean(dfTwoF)
          row[["pF"]] <-.clean(pValTwo)
          
          if (i == 1) {
            fields[[length(fields) + 1]] <- list(name="F", type="number", format="sf:4;dp:3")
            fields[[length(fields) + 1]] <- list(name="df num", type="integer")
            fields[[length(fields) + 1]] <- list(name="df den", type="integer")
            fields[[length(fields) + 1]] <- list(name="pF", title="p<sub>F</sub>",type="number", format="dp:3;p:.001")
          }
          
        } 
        
        rows[[i]] <- row
        
      } else {
        status$error <- TRUE
        status$errorMessage <- "Specified ANOVA design is not balanced."
        result[["error"]] <- list(errorType="badData", errorMessage=status$errorMessage)
        
        row <- list()
        row[["Factor"]] <- "."
        row[["Statistic"]] <- "."
        row[["df"]] <- "."
        row[["p"]] <- "."
        
        rows[[i]] <- row
      }
    }
  } else {
    
    row <- list()
    row[["Factor"]] <- "."
    row[["Statistic"]] <- "."
    row[["df"]] <- "."
    row[["p"]] <- "."
    
    rows[[1]] <- row
  }
  
  
  result[["data"]] <- rows
  result[["status"]] <- "complete"
  
  result[["schema"]] <- list(fields=fields)
  result[["footnotes"]] <- as.list(footnotes)
  
  
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    stateFriedman <- result
    
  } else {
    
    stateFriedman <- NULL
    
  }
  
  list(result=result, status=status, stateFriedman=stateFriedman)
}

.rmAnovaConoverTable <- function(dataset, options, perform, fullModel, status, stateConover, singular) {
  
  if (options$conoverTest == FALSE | identical(options$friedmanWithinFactor, ""))
    return (list(result=NULL, status=status))
  
  groupingVariables <- unlist(options$friedmanWithinFactor)
  blockingVar <- ifelse( identical(options$friedmanBetweenFactor, ""), "subject", .v(options$friedmanBetweenFactor))
  
  conoverTables <- list()
  
  for (groupingVar in groupingVariables) {
    
    conoverTable <- list()
    
    conoverTable[["title"]] <- paste("Conover's Post Hoc Comparisons - ", groupingVar, sep="")
    conoverTable[["name"]] <- paste("conoverTest_", groupingVar, sep="")
    
    fields <- list(
      list(name="(I)",title="", type="string", combine=TRUE),
      list(name="(J)",title="", type="string"),
      list(name="t",  title="T-Stat", type="number", format="sf:4;dp:3"),
      list(name="df", type="integer"),
      list(name="wA", title="W<sub>i</sub>", type="number", format="sf:4;dp:3"),
      list(name="wB", title="W<sub>j</sub>", type="number", format="sf:4;dp:3"),
      list(name="pval", title="p", type="number", format="dp:3;p:.001"),
      list(name="bonferroni", title="p<sub>bonf</sub>", type="number", format="dp:3;p:.001"),
      list(name="holm",title="p<sub>holm</sub>", type="number", format="dp:3;p:.001")
    )
    
    conoverTable[["schema"]] <- list(fields=fields)
    
    rows <- list()
    
    if (perform == "run" && status$ready && status$error == FALSE) {
      
      longData <- fullModel$data$long
      
      groups <- as.factor(longData[, .v(groupingVar)])
      blocks <- as.factor(longData[, blockingVar])
      y <- longData[, 'dependent']
      
      groupNames <- .unv(levels(groups))
      ## Code from PMCMRplus
      t <- nlevels(groups)
      b <- nlevels(blocks)
      r <- unique(table(groups))
      k <- unique(table(blocks)) 
      rij <- unlist(tapply(y, blocks, rank))
      Rj <- unname(tapply(rij, groups, sum))
      
      df <- b * k - b - t + 1
      
      S2 <- 1 / ( 1 * t -1 ) * (sum(rij^2) - t * b * ( t + 1)^2 / 4)
      T2 <- 1 / S2 * (sum(Rj) - b * ((t + 1) / 2)^2)
      A <- S2 * (2 * b * (t - 1)) / ( b * t - t - b + 1)
      B <- 1 - T2 / (b * (t- 1))
      denom <- sqrt(A) * sqrt(B)
      
      for (i in 1:t) {
        
        for (j in .seqx(i+1, t)) {
          
          row <- list("(I)"=groupNames[[i]], "(J)"=groupNames[[j]])
          
          diff <-  abs(Rj[i] - Rj[j]) 
          tval <- diff / denom
          pval <- 2 * pt(q = abs(tval), df = df, lower.tail=FALSE)
          
          row[["t"]] <- .clean(tval)
          row[["wA"]]  <- .clean(Rj[i])
          row[["wB"]] <- .clean(Rj[j])
          row[["pval"]] <- .clean(pval)
          row[["bonferroni"]] <- .clean(pval)
          row[["holm"]] <- .clean(pval)
          row[["df"]] <- .clean(df)
          
          conoverTable[["status"]] <- "complete"
          rows[[length(rows)+1]] <- row
          
        }
        
        if (length(rows) == 0)  {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
      }
      
      allP <- unlist(lapply(rows, function(x) x$p))
      allBonf <- p.adjust(allP, method = "bonferroni")
      allHolm <- p.adjust(allP, method = "holm")
      
      for (p in 1:length(rows)) {
        rows[[p]][['bonferroni']] <- .clean(allBonf[p])
        rows[[p]][['holm']] <- .clean(allHolm[p])
      }
      
    } else {
      row <- list("(I)"= ".", "(J)"= ".")
      row[["t"]] <- "."
      row[["wA"]]  <- "."
      row[["wB"]] <- "."
      row[["pval"]] <- "."
      row[["bonferroni"]] <- "."
      row[["holm"]] <- "."
      row[["df"]] <- "."
      
      rows[[length(rows)+1]] <- row
    }
    
    conoverTable[["data"]] <- rows
    
    conoverTables[[length(conoverTables)+1]] <- conoverTable
  }
  
  if (perform == "run" && status$ready && status$error == FALSE)  {
    
    stateConover <- conoverTables
    
  } else {
    
    stateConover <- NULL
    
  }
  
  list(result=conoverTables, status=status, stateConover=stateConover)
}

.rmAnovaDescriptivesTable <- function(dataset, options, perform, status, stateDescriptivesTable) {
  
  if (options$descriptives == FALSE)
    return(list(result=NULL, status=status))
  
  rmFactors <- c()
  rmLevels <- list()
  
  for (i in .indices(options$repeatedMeasuresFactors)) {
    
    rmFactors[i] <- options$repeatedMeasuresFactors[[i]]$name
    rmLevels[[i]] <- options$repeatedMeasuresFactors[[i]]$levels
    
  }
  
  bsFactors <- c()
  bsLevels <- list()
  
  for (i in .indices(options$betweenSubjectFactors)) {
    
    bsFactors[i] <- options$betweenSubjectFactors[i]
    bsLevels[[i]] <- levels(dataset[[ .v(options$betweenSubjectFactors[i]) ]])
    
  }
  
  factors <- c(rmFactors, bsFactors)
  lvls <- c(rmLevels, bsLevels)
  
  descriptives.table <- list()
  
  descriptives.table[["title"]] <- "Descriptives"
  
  fields <- list()
  
  for (variable in factors) {
    
    name <- paste(".", variable, sep="")  # in case variable is "Mean", "SD" or "N"
    fields[[length(fields)+1]] <- list(name=name, type="string", title=variable, combine=TRUE)
    
  }
  
  fields[[length(fields)+1]] <- list(name="Mean", type="number", format="sf:4;dp:3")
  fields[[length(fields)+1]] <- list(name="SD", type="number", format="sf:4;dp:3")
  fields[[length(fields)+1]] <- list(name="N", type="integer")
  
  descriptives.table[["schema"]] <- list(fields=fields)
  
  cases <- rev(expand.grid(rev(lvls)))
  
  namez <- unlist(factors)
  column.names <- paste(".", namez, sep="")
  
  if (length(factors) > 0) {
    
    rows <- list()
    
    if (perform == "run" && status$ready && status$error == FALSE) {
      
      dataset <- .shortToLong(dataset, options$repeatedMeasuresFactors, options$repeatedMeasuresCells, options$betweenSubjectFactors)
      
      for (i in 1:dim(cases)[1]) {
        
        row <- list()
        
        for (j in 1:dim(cases)[2])
          row[[ column.names[[j]] ]] <- as.character(cases[i, j])
        
        sub  <- eval(parse(text=paste("dataset$", .v(namez), " == \"", row, "\"", sep="", collapse=" & ")))
        
        data <- base::subset(dataset, sub, select="dependent")[[1]]
        
        N <- base::length(data)
        
        row[["N"]] <- N
        
        if (N == 0) {
          
          row[["Mean"]] <- ""
          row[["SD"]]   <- ""
          
        } else if (N == 1) {
          
          row[["Mean"]] <- data
          row[["SD"]]   <- ""
          
        } else {
          
          row[["Mean"]] <- base::mean(data)
          row[["SD"]]   <- stats::sd(data)
        }
        
        if(cases[i,dim(cases)[2]] == lvls[[ dim(cases)[2] ]][[1]]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[i]] <- row
      }
      
    } else {
      
      for (i in 1:dim(cases)[1]) {
        
        row <- list()
        
        for (j in 1:dim(cases)[2])
          row[[ column.names[[j]] ]] <- as.character(cases[i, j])
        
        if(cases[i,dim(cases)[2]] == lvls[[ dim(cases)[2] ]][[1]]) {
          row[[".isNewGroup"]] <- TRUE
        } else {
          row[[".isNewGroup"]] <- FALSE
        }
        
        rows[[i]] <- row
      }
    }
    
    descriptives.table[["data"]] <- rows
    
    if (perform == "run" && status$ready && status$error == FALSE)
      descriptives.table[["status"]] <- "complete"
  }
  
  if (status$error)
    descriptives.table[["error"]] <- list(error="badData")
  
  if (perform == "run" && status$ready && status$error == FALSE) {
    
    stateDescriptivesTable <- descriptives.table
    
  } else {
    
    stateDescriptivesTable <- NULL
    
  }
  
  list(result=descriptives.table, status=status, stateDescriptivesTable=stateDescriptivesTable)
}

.rmAnovaDescriptivesPlot <- function(dataset, options, perform, status, stateDescriptivesPlot) {
  
  descriptivesPlotList <- list()
  
  if (perform == "run" && status$ready && !status$error && options$plotHorizontalAxis != "") {
    
    dataset <- .shortToLong(dataset, options$repeatedMeasuresFactors, options$repeatedMeasuresCells, options$betweenSubjectFactors)
    
    groupVars <- c(options$plotHorizontalAxis, options$plotSeparateLines, options$plotSeparatePlots)
    groupVars <- groupVars[groupVars != ""]
    groupVarsV <- .v(groupVars)
    
    betweenSubjectFactors <- groupVars[groupVars %in% options$betweenSubjectFactors]
    repeatedMeasuresFactors <- groupVars[groupVars %in% sapply(options$repeatedMeasuresFactors,function(x)x$name)]
    
    usePooledSE <- ifelse(is.null(options$usePooledStandErrorCI), FALSE, options$usePooledStandErrorCI)
    
    if (length(repeatedMeasuresFactors) == 0) {
      
      summaryStat <- .summarySE(as.data.frame(dataset), measurevar = "dependent", groupvars = .v(betweenSubjectFactors),
                                conf.interval = options$confidenceIntervalInterval, na.rm = TRUE, .drop = FALSE, errorBarType = options$errorBarType, 
                                usePooledSE=usePooledSE)
      
    } else {
      
      summaryStat <- .summarySEwithin(as.data.frame(dataset), measurevar="dependent", betweenvars=.v(betweenSubjectFactors), withinvars=.v(repeatedMeasuresFactors),
                                      idvar="subject", conf.interval=options$confidenceIntervalInterval, na.rm=TRUE, .drop=FALSE, errorBarType=options$errorBarType, 
                                      usePooledSE=usePooledSE)
      
    }
    
    if ( options$plotHorizontalAxis != "" ) {
      colnames(summaryStat)[which(colnames(summaryStat) == .v(options$plotHorizontalAxis))] <- "plotHorizontalAxis"
    }
    
    if ( options$plotSeparateLines != "" ) {
      colnames(summaryStat)[which(colnames(summaryStat) == .v(options$plotSeparateLines))] <- "plotSeparateLines"
    }
    
    if ( options$plotSeparatePlots != "" ) {
      colnames(summaryStat)[which(colnames(summaryStat) == .v(options$plotSeparatePlots))] <- "plotSeparatePlots"
    }
    
    base_breaks_x <- function(x){
      b <- unique(as.numeric(x))
      d <- data.frame(y=-Inf, yend=-Inf, x=min(b), xend=max(b))
      list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE, size = 1))
    }
    
    base_breaks_y <- function(x, plotErrorBars){
      if (plotErrorBars) {
        ci.pos <- c(x[,"dependent"], x[,"ciLower"],x[,"ciUpper"])
        b <- pretty(ci.pos)
        d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
        list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE, size = 1),
             ggplot2::scale_y_continuous(breaks=c(min(b),max(b))))
      } else {
        b <- pretty(x[,"dependent"])
        d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
        list(ggplot2::geom_segment(data=d, ggplot2::aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE, size = 1),
             ggplot2::scale_y_continuous(breaks=c(min(b),max(b))))
      }
    }
    
    if (options$plotSeparatePlots != "") {
      subsetPlots <- levels(summaryStat[,"plotSeparatePlots"])
      nPlots <- length(subsetPlots)
    } else {
      nPlots <- 1
    }
    
    for (i in 1:nPlots) {
      
      descriptivesPlot <- list()
      
      if (options$plotSeparateLines != "") {
        
        descriptivesPlot[["width"]] <- options$plotWidthDescriptivesPlotLegend
        descriptivesPlot[["height"]] <- options$plotHeightDescriptivesPlotLegend
        descriptivesPlot[["custom"]] <- list(width="plotWidthDescriptivesPlotLegend", height="plotHeightDescriptivesPlotLegend")
        
      } else {
        
        descriptivesPlot[["width"]] <- options$plotWidthDescriptivesPlotNoLegend
        descriptivesPlot[["height"]] <- options$plotHeightDescriptivesPlotNoLegend
        descriptivesPlot[["custom"]] <- list(width="plotWidthDescriptivesPlotNoLegend", height="plotHeightDescriptivesPlotNoLegend")
        
      }
      
      if (options$plotSeparatePlots != "") {
        summaryStatSubset <- subset(summaryStat,summaryStat[,"plotSeparatePlots"] == subsetPlots[i])
      } else {
        summaryStatSubset <- summaryStat
      }
      
      if(options$plotSeparateLines == "") {
        
        p <- ggplot2::ggplot(summaryStatSubset, ggplot2::aes(x=plotHorizontalAxis,
                                                             y=dependent,
                                                             group=1))
        
      } else {
        
        p <- ggplot2::ggplot(summaryStatSubset, ggplot2::aes(x=plotHorizontalAxis,
                                                             y=dependent,
                                                             group=plotSeparateLines,
                                                             shape=plotSeparateLines,
                                                             fill=plotSeparateLines))
        
      }
      
      if (options$plotErrorBars) {
        
        pd <- ggplot2::position_dodge(.2)
        p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=ciLower,
                                                     ymax=ciUpper),
                                        colour="black", width=.2, position=pd)
        
      } else {
        
        pd <- ggplot2::position_dodge(0)
        
      }
      
      p <- p + ggplot2::geom_line(position=pd, size = .7) +
        ggplot2::geom_point(position=pd, size=4) +
        ggplot2::scale_fill_manual(values = c(rep(c("white","black"),5),rep("grey",100))) +
        ggplot2::scale_shape_manual(values = c(rep(c(21:25),each=2),21:25,7:14,33:112)) +
        ggplot2::scale_color_manual(values = rep("black",200)) +
        ggplot2::ylab(options$labelYAxis) +
        ggplot2::xlab(options$plotHorizontalAxis) +
        ggplot2::labs(shape=options$plotSeparateLines, fill=options$plotSeparateLines) +
        base_breaks_y(summaryStat, options$plotErrorBars) +
        base_breaks_x(summaryStatSubset[,"plotHorizontalAxis"])
      
      nrowsInLegend <- min(10, nlevels(as.factor(summaryStatSubset[["plotSeparateLines"]])))
      guide <- ggplot2::guide_legend(nrow = nrowsInLegend)
      p <- p + ggplot2::guides(fill = guide, shape = guide, color = guide)
      
      p <- JASPgraphs::themeJasp(p, legend.position="right")
      
      if (nPlots > 1) {
        descriptivesPlot[["title"]] <- paste(options$plotSeparatePlots,": ",subsetPlots[i], sep = "")
      } else {
        descriptivesPlot[["title"]] <- "Descriptives Plot"
      }
      
      if (options$plotSeparateLines != "") {
        
        # image <- .beginSaveImage(options$plotWidthDescriptivesPlotLegend, options$plotHeightDescriptivesPlotLegend)
        content <- .writeImage(width = options$plotWidthDescriptivesPlotLegend, 
                               height = options$plotHeightDescriptivesPlotLegend, 
                               plot = p, obj = TRUE)
        
      } else {
        
        # image <- .beginSaveImage(options$plotWidthDescriptivesPlotNoLegend, options$plotHeightDescriptivesPlotNoLegend)
        content <- .writeImage(width = options$plotWidthDescriptivesPlotNoLegend, 
                               height = options$plotHeightDescriptivesPlotNoLegend, 
                               plot = p, obj = TRUE)
        
      }
      
      # content <- .endSaveImage(image)
      
      descriptivesPlot[["convertible"]] <- TRUE
      descriptivesPlot[["obj"]] <- content[["obj"]]
      descriptivesPlot[["data"]] <- content[["png"]]
      
      # descriptivesPlot[["data"]] <- content
      descriptivesPlot[["status"]] <- "complete"
      
      descriptivesPlotList[[i]] <- descriptivesPlot
      
    }
    
    stateDescriptivesPlot <- descriptivesPlotList
    
  } else if (options$plotHorizontalAxis != "") {
    
    if (options$plotSeparatePlots != "") {
      
      repeatedMeasuresNames <- sapply(options$repeatedMeasuresFactors, function(x) x$name)
      repeatedMeasuresLevels <- lapply(options$repeatedMeasuresFactors, function(x) x$levels)
      
      if (sum(options$plotSeparatePlots == repeatedMeasuresNames) > 0) {
        
        index <- which(options$plotSeparatePlots == repeatedMeasuresNames)
        nPlots <- length(unlist(repeatedMeasuresLevels[[index]]))
        
      } else {
        
        nPlots <- length(levels(dataset[[ .v(options$plotSeparatePlots) ]]))
        
      }
      
    } else {
      
      nPlots <- 1
      
    }
    
    for (i in 1:nPlots) {
      
      descriptivesPlot <- list()
      
      if (nPlots == 1) {
        descriptivesPlot[["title"]] <- "Descriptives Plot"
      } else {
        descriptivesPlot[["title"]] <- ""
      }
      
      if (options$plotSeparateLines != "") {
        
        descriptivesPlot[["width"]] <- options$plotWidthDescriptivesPlotLegend
        descriptivesPlot[["height"]] <- options$plotHeightDescriptivesPlotLegend
        descriptivesPlot[["custom"]] <- list(width="plotWidthDescriptivesPlotLegend", height="plotHeightDescriptivesPlotLegend")
        
      } else {
        
        descriptivesPlot[["width"]] <- options$plotWidthDescriptivesPlotNoLegend
        descriptivesPlot[["height"]] <- options$plotHeightDescriptivesPlotNoLegend
        descriptivesPlot[["custom"]] <- list(width="plotWidthDescriptivesPlotNoLegend", height="plotHeightDescriptivesPlotNoLegend")
        
      }
      
      descriptivesPlot[["data"]] <- ""
      
      if (status$error)
        descriptivesPlot[["error"]] <- list(errorType="badData")
      
      descriptivesPlotList[[i]] <- descriptivesPlot
    }
    
    stateDescriptivesPlot <- NULL
    
  }
  
  list(result=descriptivesPlotList, status=status, stateDescriptivesPlot=stateDescriptivesPlot)
}

.summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE, 
                       errorBarType="confidenceInterval", usePooledSE=FALSE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  # First aggregate over unused RM factors, if desired:
  if (usePooledSE & measurevar == "dependent") {
    data <- plyr::ddply(data, c("subject", groupvars), plyr::summarise, dependent = mean(dependent))
    names(data)[which(names(data) == "dependent")] <- measurevar
  } else if (usePooledSE & measurevar == "dependent_norm") {
    data <- plyr::ddply(data, c("subject", groupvars), plyr::summarise, dependent = mean(dependent_norm))
    names(data)[which(names(data) == "dependent")] <- measurevar
  }
  
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  
  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  if (errorBarType == "confidenceInterval") {
    
    datac$ciLower <- datac[,measurevar] - datac[,"ci"]
    datac$ciUpper <- datac[,measurevar] + datac[,"ci"]
    
  } else {
    
    datac$ciLower <- datac[,measurevar] - datac[,"se"]
    datac$ciUpper <- datac[,measurevar] + datac[,"se"]
    
  }
  
  return(datac)
}

.normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL, na.rm=FALSE, .drop=TRUE) {
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- plyr::ddply(data, c(idvar, betweenvars), .drop=.drop,
                               .fun = function(xx, col, na.rm) {
                                 c(subjMean = mean(xx[,col], na.rm=na.rm))
                               },
                               measurevar,
                               na.rm
  )
  
  
  
  # Put the subject means with original data
  data <- base::merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

.summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL, idvar=NULL, na.rm=FALSE, 
                             conf.interval=.95, .drop=TRUE, errorBarType="confidenceInterval", usePooledSE=FALSE) {
  
  # Get the means from the un-normed data
  datac <- .summarySE(data, measurevar, groupvars=c(betweenvars, withinvars), na.rm=na.rm, 
                      conf.interval=conf.interval, .drop=.drop, errorBarType=errorBarType, usePooledSE=usePooledSE)
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  datac$ciLower <- NULL
  datac$ciUpper <- NULL
  
  # Norm each subject's data
  ndata <- .normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- .summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars), na.rm=na.rm, conf.interval=conf.interval, .drop=.drop, errorBarType=errorBarType,
                       usePooledSE=usePooledSE)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  # Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels, FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  if (errorBarType == "confidenceInterval") {
    
    ndatac$ciLower <- datac[,measurevar] - ndatac[,"ci"]
    ndatac$ciUpper <- datac[,measurevar] + ndatac[,"ci"]
    
  } else {
    
    ndatac$ciLower <- datac[,measurevar] - ndatac[,"se"]
    ndatac$ciUpper <- datac[,measurevar] + ndatac[,"se"]
    
  }
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

