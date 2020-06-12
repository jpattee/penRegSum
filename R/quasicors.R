#' @title Function to obtain quasi-correlation estimates given 
#' 
#' @details Function requires univariate effect size estimates, standard errors, and sample sizes, candidate polygenic risk score estimates, and a reference panel.
#' 
#' @param testBetas Univariate effect sizes for the out-of-sample data
#' @param testSes Univariate standard errors for the effect sizes of the out-of-sample data
#' @param N Sample sizes corresponding to the testBetas. Can be a constant or a vector.
#' @param penalizedBetas The polygenic risk scores to be tested on the out-of-sample data. May be a matrix.
#' @param refPanel PLINK stem for binary file to be used as the reference panel
#' @param allele1 Vector of effect alleles for the penalizedBetas. Corresponds to the fifth column of a PLINK .bim file.
#' @param allele2 Vector of reference alleles for the penalizedBetas. Corresponds to the sixth column of a PLINK .bim file.
#' @param standardized Set to true if the coefficient estimates for penalizedBetas are standardized. Note that elastSum and tlpSum output standardized estimates. 
#' @param trainVar Phenotypic variance for data used to estimate the penalizedBetas. Only necessary if standardized = TRUE
#' @param allele1.test Vector of effect alleles for the testBetas. Corresponds to the fifth column of a PLINK .bim file.
#' @param allele2.test Vector of reference alleles for the testBetas. Corresponds to the sixth column of a PLINK .bim file.
#' @param ldBlocks Location of file specifying independent LD Blocks to be used. File should be in the BED file format. If null, estimation is done by chromosome.
#' @export

quasicors <- function(testBetas, testSes, N, penalizedBetas, refPanel, allele1 = NULL, allele2 = NULL, standardized = TRUE, trainVar = NULL, allele1.test = NULL, allele2.test = NULL){
  if(length(N) == 1) N = rep(N, length(betas))
  
  penalizedBetas = as.matrix(penalizedBetas)
  
  bim = read.table(paste0(refPanel,".bim"), stringsAsFactors = FALSE)

  if(!is.null(allele1)){
    flippedInd = which(bim$V5==allele2 & bim$V6==allele1)
    penalizedBetas[flippedInd,] = penalizedBetas[flippedInd,] * -1
  }
  
  if(!is.null(allele1.test)){
    flippedIndTest = which(bim$V5==allele2.test & bim$V6==allele1.test)
    testBetas[flippedIndTest] = testBetas[flippedIndTest] * -1
  }

  p = nrow(bim)

  fam = read.table(paste0(refPanel,".fam"))
  nFam = nrow(fam)

  genoMat=genotypeMatrix(paste0(refPanel,".bed"),N=nFam,P=p,integer(0),integer(0),integer(0),integer(0),1)

  #If LD blocks are specified, format the file and proceed
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
  
  for(i in 1:nrow(ldDat)){
    curChr = substring(ldDat[i,1],4)
    curMin = ldDat[i,2]
    curMax = ldDat[i,3]
    curSnps = bim$V2[bim$V4 > curMin & bim$V4 <= curMax & bim$V1 == curChr]
    tempInd = which(bim$V2%in%curSnps)
    if(length(tempInd) > 0){
      matList[[matListInd]] = cov(genoMat[,tempInd, drop = FALSE])
      matListInd = matListInd+1
    }
  }
  V = bdiag(matList)

  sds = normalize(genoMat)
  rm(genoMat)

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
    if(standardized) penalizedBetasTemp = penalizedBetas[,i]*sqrt(trainVar)/sds
    num[i] = sum((sds^2)*testBetas*penalizedBetasTemp)
    denom[i] = sqrt(t(penalizedBetasTemp)%*%V%*%penalizedBetasTemp)
    quasicorrelations[i] = num[i]/(denom[i]*sqrt(ytyEstTest))
  }
  toReturn = structure(list(quasicorrelations = quasicorrelations, bxy = num, bxxb = denom^2))
  return(toReturn)
}