#require caret package
library(caret)

#vote on the best number of features

#' Iterated Recursive Feature Elimination
#'
#' For each possible size of optimal feature subset (OFS) performs multiple interations of RFE algorithm, using random forest as a classifier and accuracy as a performance metrics.
#' Each iteration 'votes' on the features it has selected.
#' For every N considered as possible OFS size, N features with top votes number are selected
#'
#'
#' @param trainSet Dataframe used for Random Forrest cassifier training
#' @param testSet Dataframe used for Random Forrest testing - accuracy calculation. Not used if useCV == T. Must have the same column names as trainSet
#' @param initFeatures A vector containing a subset of trainSet column names to choose optmial features from
#' @param classLab A name of trainSet column containing class variable
#' @param checkNFeatures Range of possible OFS sizes to consider (the algorithm will consider 1:checkNFeatures sizes)
#' @param votingIterations A number of times RFE algorithm will be iterated, each time with different, random hyperparameters
#' @param useCV Whether to use a cross validation (T) or Test set (F) for accuracy evaluation
#' @param nfolds Number of folds to use for cross validation (Used only if useCV == T)
#' @param initRandomState Initial random state to use for model hyperparameters selection
#'
#' @return Returns a list containing three elements:
#' 'accuracyPerNFeatures' - a dataframe containing average model accuracy for each of 1:checkNFeatures OFS sizes
#' 'votesPerN'- a Dataframe containing 'voting results' for each OFS size
#' 'topFeaturesPerN' - a list containing OFS for each (1:checkNFeatures) OFS size
#' @export
#'
#' @examples
iteratedRFE <- function(trainSet, testSet, initFeatures, classLab, checkNFeatures = 25, votingIterations = 100000, useCV = F, nfolds = 10, initRandomState = 42 ) {

  set.seed(initRandomState)

  #prepare output data structures
  resAcc <- c(rep(0, checkNFeatures))
  resVotes <- data.frame(matrix(0, nrow = length(initFeatures), ncol = checkNFeatures), row.names = initFeatures)
  for(i in 1:checkNFeatures) colnames(resVotes)[i] <- toString(i)
  resTop <- list()


  for (i in 1:votingIterations) {

    if(useCV == F) {
      params <- rfeControl(functions = rfFuncs, saveDetails = T)
      iter <- rfeIter(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), testX = testSet[, initFeatures], testY = as.factor(testSet[, classLab]), sizes = 1:checkNFeatures,
                 metric = "Accuracy", rfeControl = params)

      for(j in 1:checkNFeatures) {
        tmp <- iter$pred[iter$pred$Variables == j, ]

        acc <- length(which(tmp$pred == tmp$obs)) / nrow(tmp) #calculate and add accuracy
        resAcc[j] <- resAcc[j] + acc

        selected <- iter$finalVariables[[j+1]]$var
        numb <- iter$finalVariables[[j+1 ]]$Variables[1]

        resVotes[selected, numb] <- resVotes[selected, numb] + 1
      }


    }
    else {

      seeds <- vector(mode = "list", length = nfolds + 1) # add random seeds for cross validation
      for(i in 1:nfolds) seeds[[i]] <- sample.int(1000000000, checkNFeatures + 1)
      seeds[nfolds + 1] <- sample.int(1000000000, 1)

      params <- rfeControl(functions = rfFuncs, number = nfolds, saveDetails = T)
      iter <- rfe(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), sizes = 1:checkNFeatures, rfeControl = params)

      for(j in 1:checkNFeatures) {
        tmp <- iter$variables[iter$variables$Variables == j, ]
        for(k in tmp$var) resVotes[k, j] <- resVotes[k, j] + 1 # increase a voting score for each fold

        resAcc[j] <- resAcc[j] + iter$results[iter$results$Variables == j, "Accuracy"]

      }
    }
  }

  resAcc <- resAcc / votingIterations #make average accuracy

  for(i in 1:ncol(resVotes)) resTop[[i]] <- rownames(resVotes[order(-resVotes[, i])[1:i], ])

  returning <- list(data.frame(resAcc), resVotes, resTop)
  names(returning) <- c("accuracyPerNFeatures", "votesPerN", "topFeaturesPerN")
  return(returning)

}
