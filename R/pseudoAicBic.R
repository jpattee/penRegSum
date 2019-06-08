#' @title Function to obtain pseudo-AIC and pseudo-BIC estimates given summary statistics and a reference panel
#' 
#' @details Function requires univariate effect size estimates, standard errors, and sample sizes, candidate polygenic risk score estimates, and a reference panel.
#' 
#' @param penalizedBetas Polygenic risk score effect size estimates. These are the estimates for which we want to approximate AIC and BIC. May be a matrix.
#' @param betas Univariate effect size estimates.
#' @param ses Standard errors of the univariate effect size estimates.
#' @param N Sample sizes corresponding to the univariate effect size estimates. Can be a constant.
#' @param refPanel Stem of PLINK binary file
#' @param sigSqReg The matrix regularization parameter for the estimation of the residual variance. Default is .2.
#' @param sseReg The matrix regularization parameter for the estimation of the SSE. Default is .1.
#' @param sigSqInd Index of SNPs to be used for estimation of residual variance. If null, default is to pick the n/5 SNPs with largest univariate effect sizes, where n is the sample size of the reference panel.
#' @export

pseudoAicBic <- function(penalizedBetas, betas, ses, N, refPanel, sigSqReg = .2, sseReg = .1, sigSqInd = NULL){

  tuneBim = read.table("tuning.bim",stringsAsFactors = FALSE)
  p = nrow(tuneBim)
  N = nrow(trainFam)

  genoMat=genotypeMatrix("testing-1.bed",N=1560,P=p,integer(0),integer(0),integer(0),integer(0),1)

  flippedInd=which(trainBim$V5==tuneBim$V6 & trainBim$V5==tuneBim$V6)

  n = nrow(genoMat)

  betas[flippedInd,] = betas[flippedInd,]

  covMat = cov(genoMat)

  sds=normalize(genoMat)
  sds = sds[,1]
  rm(genoMat)

  #FOR NOW
  SSEvec = NULL
  qVec = NULL
  aicVec = NULL
  bicVec = NULL
  bxxb = NULL
  bxy = NULL

  tempEst = rep(0, nrow(sumStats))
  for(i in 1:nrow(sumStats)){
    tempEst[i] = N *sds[i]^2 * ses[i,l]^2 + sds[i]^2*betas[i,l]^2
  }
  ytyEst = median(tempEst)

  xtyEst = sds^2 * betas

  if(is.null(sigSqInd)){
    frac = floor(n / 5)
    quant = 1 - frac / length(betas)
    if(quant < 0) quant = 0
    cutoff = quantile(abs(betas), quant)
    sigSqInd = which(abs(betas) > cutoff)
  }
  
  #calculate sigma^2
  xtyTemp = xtyEst[sigSqInd]
  xtxTemp = covMat[sigSqInd,sigSqInd]
  diag(xtxTemp) = diag(xtxTemp) + sigSqReg
  xtxInvTemp = solve(xtxTemp)
  qTemp = length(sigSqInd)
  sigSqTilde = (ytyEst - t(xtyTemp)%*%xtxInvTemp%*%xtyTemp)*median(N) / (median(N) - qTemp)

  for(k in 1:ncol(penalizedBetas)){
    penalizedBetasTemp = penalizedBetas[,k]

    penalizedBetasTemp = penalizedBetasTemp *sqrt(ytyEst) /sds


    qElast = sum(penalizedBetasTemp!=0)
    qVec = c(qVec,qElast)
  
    xtxWeightDiag = covMat
    diag(xtxWeightDiag) = diag(xtxWeightDiag) + sseReg
    bxxbTemp = t(penalizedBetasTemp)%*%covMat%*%penalizedBetasTemp
    bxxbWeight = t(penalizedBetasTemp)%*%xtxWeightDiag%*%penalizedBetasTemp
    bxxb = c(bxxb, bxxbTemp[1,1])
  
    bxyTemp = t(penalizedBetasTemp)%*%xtyEst
    bxy = c(bxy,bxyTemp[1,1])

    SSEest = ytyEst - (2 * bxyTemp) + bxxbWeight
    SSEvec = c(SSEvec,SSEest[1,1])

    logLik = -SSEest*median(N)/ (2 *  sigSqTilde)

    aicTemp = (2 * qElast) - (2 * logLik)
    bicTemp = (log(median(N)) * qElast) - (2 * logLik)

    aicVec = c(aicVec, aicTemp[1,1])
    bicVec = c(bicVec, bicTemp[1,1])
  }

  toReturn = structure(list(aic = aicVec, bic=bicVec, SSE = SSEvec, q = qVec, bxxb = bxxb, bxy = bxy, sigSqTilde = sigSqTilde))
  return(toReturn)
  #' @return a list with the following
  #' \item{aic} Estimates of AIC for the candidate polygenic risk scores.
  #' \item{bic} Estimates of BIC for the candidate polygenic risk scores.
  #' \item{SSE} Estimates of SSE for the candidate polygenic risk scores.
  #' \item{q} Number of active parameters for the candidate polygenic risk scores.
  #' \item{bxxb} Estimates of the b'x'xb quantity for the candidate polygenic risk scores
  #' \item{bxy} Estimates of b'x'y for the candidate polygenic risk scores.
  #' \item{sigSqTilde} Estimate of the residual variance.
}