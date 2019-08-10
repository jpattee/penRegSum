#' @title Function to obtain quasi-correlation estimates given 
#' 
#' @details Function requires univariate effect size estimates, standard errors, and sample sizes, candidate polygenic risk score estimates, and a reference panel.
#' 
#' @param testBetas Univariate effect sizes for the out-of-sample data
#' @param testSes Univariate standard errors for the effect sizes of the out-of-sample data
#' @param N Sample sizes corresponding to the univariate effect size estimates. Can be a constant or a vector.
#' @param penalizedBetas The polygenic risk scores to be tested on the out-of-sample data. May be a matrix.
#' @param refPanel PLINK stem for binary file to be used as the reference panel
#' @param allele1 Vector of effect alleles for the penalizedBetas. Corresponds to the fifth column of a PLINK .bim file.
#' @param allele2 Vector of reference alleles for the penalizedBetas. Corresponds to the sixth column of a PLINK .bim file.
#' @param standardized Set to true if the coefficient estimates for penalizedBetas are standardized. Note that elastSum and tlpSum output standardized estimates. 
#' @param trainBetas Vector of univariate effect sizes for the data used to estimate the polygenic risk score. Used to estimate variance of training phenotype. Unneccessary if standardized = FALSE.
#' @param trainBetas Vector of univariate standard errors for the data used to estimate the polygenic risk score. Used to estimate variance of training phenotype. Unneccessary if standardized = FALSE.
#' @param trainN Sample size of the training data. Can be a constant or a vector.
#' @export

quasicors <- function(testBetas, testSes, N, penalizedBetas, refPanel, allele1 = NULL, allele2 = NULL, standardized = TRUE, trainBetas = NULL, trainSes = NULL, trainN = NULL){
  if(length(N) == 1) N = rep(N, length(betas))
  
  penalizedBetas = as.matrix(penalizedBetas)
  
  bim = read.table(paste0(refPanel,".bim"), stringsAsFactors = FALSE)

  if(!is.null(allele1)){
    flippedInd = which(bim$V5==allele2 & bim$V5==allele1)
    matchInd = which(bim$V5==allele1 & bim$V5==allele2)
    testBetas[flippedInd] = testBetas[flippedInd] * -1
    penalizedBetas[flippedInd,] = penalizedBetas[flippedInd,] * -1
  }

  p = nrow(bim)

  fam = read.table(paste0(refPanel,".fam"))
  nFam = nrow(fam)

  genoMat=genotypeMatrix(paste0(refPanel,".bed"),N=nFam,P=p,integer(0),integer(0),integer(0),integer(0),1)

  matList = vector("list", length(unique(bim$V1)))
  matListInd = 1
  for(i in unique(bim$V1)){
    tempInd = which(bim$V1 == i)
    matList[[matListInd]] = cov(genoMat[,tempInd])
    matListInd = matListInd+1
  }
  V = bdiag(matList)

  sds = normalize(genoMat)
  rm(genoMat)


  #estimate yty if not standardized
  if(standardized){
    tempEst = rep(0, length(trainSes))
    for(i in 1:length(trainSes)){
      tempEst[i] = trainN[i] * sds[i]^2 * trainSes[i]^2 + sds[i]^2*trainBetas[i]^2
    }
    ytyEst = median(tempEst,na.rm = TRUE)
  }
  
  #estimate yty for the test data
  tempEstTest = rep(0, length(testSes))
  for(i in 1:length(testSes)){
    tempEstTest = N[i] * sds[i]^2 * testSes[i]^2 + sds[i]^2*testBetas[i]^2
  }
  ytyEstTest = median(tempEstTest,na.rm = TRUE)


  quasicorrelations = rep(0, ncol(penalizedBetas))
  num = rep(0, ncol(penalizedBetas))
  denom = rep(0, ncol(penalizedBetas))
  for(i in 1:length(quasicorrelations)){
    penalizedBetasTemp = penalizedBetas[,i]
    if(standardized) penalizedBetasTemp = penalizedBetas[,i]*sqrt(ytyEst)/sds
    num[i] = sum((sds^2)*testBetas*penalizedBetasTemp)
    denom[i] = sqrt(t(penalizedBetasTemp)%*%V%*%penalizedBetasTemp)
    quasicorrelations[i] = num[i]/(denom[i]*sqrt(ytyEstTest))
  }
  toReturn = structure(list(quasicorrelations = quasicorrelations, bxy = num, bxxb = denom^2))
  return(toReturn)
  #' @return a list with the following
  #' \item{quasicorrelations} Quasi-correlation values for each of the input polygenic risk scores.
  #' \item{bxy} The estimate of the dot product of the univariate out-of-sample effect size estimates and the candidate polygenic risk score: b'x'y.
  #' \item{bxxb} The estimate of the variance of the predicted phenotype for the out of sample data: b'x'xb.
}