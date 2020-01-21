#' @title Function to obtain elastic net effect size estimates given summary statistics and a reference panel
#' 
#' @details Function finds the minimum of the elastic net objective function 
#' 
#' @param cors A vector of correlations, representing the univariate SNP-phenotype associations from the summary statistic data.
#' @param bfile The stem of a PLINK binary file to be used as a reference panel.
#' @param lambdas Vector of values of lambda.
#' @param alphas Vector of values of alpha.
#' @param s Vector of values for s.
#' @param thr Convergence threshold.
#' @param init Initial values for \eqn{\beta}.
#' @param maxIter Maximum number of iterations.
#' @param extract Vector of the indices of SNPs to keep. If null, will keep all SNPs.
#' @param ldBlocks Location of file specifying independent LD Blocks to be used. File should be in the BED file format. If null, estimation is done by chromosome.
#' @export

elastSum <- function(cors, bfile, lambdas, alphas, s=0.5, thr=1e-4, init=NULL, maxIter=1000, extract=NULL, ldBlocks = NULL){
  
  fam=read.table(paste0(bfile,".fam"))
  N=nrow(fam)
  
  bim=read.table(paste0(bfile,".bim"))
  P=nrow(bim)
  
  genoMat=genotypeMatrix(paste0(bfile,".bed"),N=N,P=P,integer(0),integer(0),integer(0),integer(0),1)
  
  sds=normalize2(genoMat)
  genoMat[is.na(genoMat)]=0
  
  lassoBetasFull=NULL
  yhatFull=NULL
  convergedFull=NULL
  
  if(is.null(extract)) extract=c(1:P)
  
  if(sum(extract!=FALSE)!=P) P=P-sum(extract==FALSE)
  if(length(extract)<P) P = length(extract)
  
  if(is.null(init)) init = rep(0,P)
  
  genoMat=genoMat[,extract]
  
  sMult=rep(0, length(s))
  for(i in 1:length(s)){
    if(i==1) sMult[i] = 1-s[i]
    else sMult[i] = (1-s[i])/(1-s[i-1])
  }
  
  
  temp = rep(lambdas, each = length(alphas))
  alphas = rep(alphas, length(lambdas))
  lambdas = temp

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
  
  permList = NULL
  
  for(b in 1:nrow(ldDat)){
    curChr = substring(ldDat[b,1],4)
    curMin = ldDat[b,2]
    curMax = ldDat[b,3]
    curSnps = tuneBim$V2[tuneBim$V4 > curMin & tuneBim$V4 <= curMax & tuneBim$V1 == curChr]
    curInd = which(tuneBim$V2%in%curSnps)
    permList = c(permList,curInd)
    if(length(curInd) > 0){
      corTemp = cors[curInd]
      lassoBetasBlock = NULL
      yhatBlock = NULL
      converged = NULL
      zeroSD = which(sds[curInd] == 0)
      for(l in 1:length(s)){
        genoMatTemp=genoMat[,curInd, drop=FALSE]*sqrt(1 - s[l])
        diagVec=rep(1-s[l],length(curInd))
        diagVec[zeroSD]=0
        denom=diagVec+s[l]
        for(i in 1:length(lambdas)){
          maxDiff=Inf
          iterator=0
          lambda=lambdas[i]
          alpha=alphas[i]
          lassoBetas=init[curInd]
          yhat=rep(0, nrow(genoMatTemp))
          while(maxDiff>thr & iterator<maxIter){
            iterator=iterator+1
            lassoBetasOld=lassoBetas
            for(j in 1:length(lassoBetas)){
              xj=lassoBetas[j]
              lassoBetas[j]=0
              tempEst=diagVec[j]*xj + corTemp[j] - sum(genoMatTemp[,j]*yhat)
              if(zeroSD[j]==1) tempEst=0
              if(abs(tempEst)>(lambda*alpha)){
                if(tempEst>0){newEst=(tempEst-(lambda*alpha))/denom[j]}
                if(tempEst<0){newEst=(tempEst+(lambda*alpha))/denom[j]}
              }
              if(abs(tempEst)<(lambda*alpha)){newEst=0}
              lassoBetas[j]=newEst/(1+lambda*(1-alpha))
              if(denom[j]==0){lassoBetas[j]=0}
              if(lassoBetas[j]!=xj){
                del=lassoBetas[j]-xj
                yhat=yhat+(del*genoMatTemp[,j])
              }
            }
            maxDiff=max(abs(lassoBetasOld-lassoBetas))
          }
          converged=rbind(converged,c(lambda,alpha,s[l],iterator!=maxIter))
          lassoBetasBlock=cbind(lassoBetasBlock,lassoBetas)
          yhatBlock=cbind(yhatFull,cbind(yhatBlock,yhat))
        }
      }
      lassoBetasFull = rbind(lassoBetasFull,lassoBetasBlock)
      if(!is.null(yhatFull)) yhatFull = yhatFull+yhatBlock
      if(is.null(yhatFull)) yhatFull = yhatBlock
      if(!is.null(convergedFull)) convergedFull[converged[,4]==0,4] = 0
      if(is.null(convergedFull)) convergedFull = converged
    }
  }
  
  ###Account for out-of-order SNPs in the BED file
  lassoBetasFull[permList,] = lassoBetasFull
  
  ###generating labeling for the output
  tempS = rep(s, each = length(lambdas))
  tempLambdas = rep(lambdas, length(s))
  tempAlphas = rep(alphas, length(s))
  
  toReturn=structure(list(lambdas=tempLambdas, alphas=tempAlphas, s = tempS, beta=lassoBetasFull, converged=converged, pred=yhatFull))
  return(toReturn)
  #' @return a list with the following
  #' \item{lambdas}Vector of lambda values corresponding to columns of output files 'beta' and 'pred'
  #' \item{alphas}Vector of alpha values corresponding to columns of output files 'beta' and 'pred'
  #' \item{s}Vector of s values corresponding to columns of output files 'beta' and 'pred'
  #' \item{beta}Matrix of estimated coefficients.
  #' \item{converged}Matrix of convergence indicators with following format: column 1 is lambda value, column 2 is alpha value, column 3 is s value, column 4 is 1 if there was convergence.
  #' \item{pred}Matrix of predicted phenotypes.
}