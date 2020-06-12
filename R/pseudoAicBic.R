#' @title Function to obtain pseudo-AIC and pseudo-BIC estimates given summary statistics and a reference panel
#' 
#' @details Function requires univariate effect size estimates, standard errors, and sample sizes, candidate polygenic risk score estimates, and a reference panel.
#' 
#' @param penalizedBetas Polygenic risk score effect size estimates. These are the estimates for which we want to approximate AIC and BIC. May be a matrix.
#' @param betas Univariate effect size estimates.
#' @param ses Standard errors of the univariate effect size estimates.
#' @param N Sample sizes corresponding to the univariate effect size estimates. Can be a constant or a vector.
#' @param refPanel Stem of PLINK binary file
#' @param sigSqReg The matrix regularization parameter for the estimation of the residual variance. Default is .2.
#' @param sseReg The matrix regularization parameter for the estimation of the SSE. Default is .2.
#' @param sigSqInd Index of SNPs to be used for estimation of residual variance. If null, default is to pick the n/5 SNPs with largest univariate effect sizes, where n is the sample size of the reference panel.
#' @param allele1 Vector of effect alleles for the betas and penalizedBetas. Corresponds to the fifth column of a PLINK .bim file.
#' @param allele2 Vector of reference alleles for the betas and penalizedBetas. Corresponds to the sixth column of a PLINK .bim file.
#' @param standardized Set to true if the coefficient estimates for penalizedBetas are standardized. Note that elastSum and tlpSum output standardized estimates.
#' @param extract Vector of the indices of SNPs to keep. If null, will keep all SNPs.
#' @param ldBlocks Location of file specifying independent LD Blocks to be used. File should be in the BED file format. If null, estimation is done by chromosome.
#' @export

pseudoAicBic <- function(penalizedBetas, betas, ses, N, refPanel, sigSqReg = .2, sseReg = .2, sigSqInd = NULL, allele1 = NULL, allele2 = NULL, standardized = TRUE, extract = NULL, ldBlocks = NULL){

  if(length(N) == 1) N = rep(N, length(betas))
  penalizedBetas = as.matrix(penalizedBetas)
  
  bim = read.table(paste0(refPanel,".bim"),stringsAsFactors = FALSE)
  P = nrow(bim)
  fam = read.table(paste0(refPanel,".fam"))
  n = nrow(fam)

  genoMat=genotypeMatrix(paste0(refPanel,".bed"),N=n,P=P,integer(0),integer(0),integer(0),integer(0),1)

  flippedInd = NULL
  
  if(!is.null(allele1)){
    flippedInd = which(bim$V5==allele2 & bim$V5==allele1)
    matchInd = which(bim$V5==allele1 & bim$V5==allele2)
    betas[flippedInd] = betas[flippedInd] * -1
  }
  
  if(is.null(extract)) extract=c(1:P)
  
  if(length(extract)<P) P = length(extract)
  
  genoMat=genoMat[,extract]
  bim = bim[extract,]
  
  if(!is.null(ldBlocks)){
    ldDat = read.table(ldBlocks, header = TRUE)
    ldDat[,1] = as.character(sub("^chr","",ldDat[,1], ignore.case = TRUE))
  }
  
  #if the LD blocks file is null, do estimation by chromosome
  if(is.null(ldBlocks)){
    vec1 = c(1:22)
    vec2 = rep(0, 22)
    vec3 = rep(Inf, 22)
    ldDat = cbind.data.frame(vec1,vec2,vec3)
  }
  
  matList = vector("list")
  matListInd = 1
  indexList = vector("list")
  
  for(i in 1:nrow(ldDat)){
    curChr = ldDat[i,1]
    curMin = ldDat[i,2]
    curMax = ldDat[i,3]
    curSnps = bim$V2[bim$V4 > curMin & bim$V4 <= curMax & bim$V1 == curChr]
    tempInd = which(bim$V2%in%curSnps)
    if(length(tempInd) > 0){
      tempMat = cov(genoMat[,tempInd, drop = FALSE])
      diag(tempMat) = diag(tempMat) + sseReg
      matList[[matListInd]] = tempMat
      indexList[[matListInd]] = tempInd
      matListInd = matListInd+1
    }
  }

  sds=normalize(genoMat)
  sds = sds[,1]
  rm(genoMat)

  SSEvec = NULL
  qVec = NULL
  aicVec = NULL
  bicVec = NULL
  bxxb = NULL
  bxy = NULL

  tempEst = rep(0, length(ses))
  for(i in 1:length(ses)){
    tempEst[i] = N[i] *sds[i]^2 * ses[i]^2 + sds[i]^2*betas[i]^2
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
  
  xtxMatList = vector("list")
  matListInd = 1
  for(g in 1:length(matList)){
    tempInd = indexList[[g]]
    matchInd = which(tempInd%in%sigSqInd)
    if(length(matchInd) > 0){
      tempMat = matList[[g]]
      diag(tempMat) = diag(tempMat) - sseReg + sigSqReg
      xtxMatList[[matListInd]] = tempMat[matchInd,matchInd]
      matListInd = matListInd + 1
    }
  }
  xtxTemp = bdiag(xtxMatList)
  xtxInvTemp = solve(xtxTemp)
  qTemp = length(sigSqInd)
  sigSqTilde = (ytyEst - t(xtyTemp)%*%xtxInvTemp%*%xtyTemp)*median(N) / (median(N) - qTemp)

  for(k in 1:ncol(penalizedBetas)){
    
    penalizedBetasTemp = penalizedBetas[,k]
    
    if(length(flippedInd) > 0) penalizedBetasTemp[flippedInd] = penalizedBetasTemp[flippedInd]*-1

    if(standardized){
      penalizedBetasTemp = penalizedBetasTemp *sqrt(ytyEst) /sds
      penalizedBetasTemp[sds == 0] = 0
    } 

    qElast = sum(penalizedBetasTemp!=0)
    qVec = c(qVec,qElast)
  
    bxxbSum = 0
    for(g in 1:length(matList)){
      tempInd = indexList[[g]]
      bxxbInc = t(penalizedBetasTemp[tempInd])%*%matList[[g]]%*%penalizedBetasTemp[tempInd]
      bxxbSum = bxxbSum + bxxbInc
    }

    bxxbWeight = bxxbSum
    bxxb = c(bxxb, bxxbWeight[1,1])
  
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
  #' \item{aic}Estimates of AIC for the candidate polygenic risk scores.
  #' \item{bic}Estimates of BIC for the candidate polygenic risk scores.
  #' \item{SSE}Estimates of SSE for the candidate polygenic risk scores.
  #' \item{q}Number of active parameters for the candidate polygenic risk scores.
  #' \item{bxxb}Estimates of the variance of the predicted phenotype resulting from the polygenic risk score applied to the reference panel: b'x'xb.
  #' \item{bxy}Estimates of the dot product of the polygenic risk score and the univariate effect size estimates: b'x'y.
  #' \item{sigSqTilde} Estimate of the residual variance.
}