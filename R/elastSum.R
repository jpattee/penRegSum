#' @title Function to obtain elastic net effect size estimates given summary statistics and a reference panel
#' 
#' @details Function finds the minimum of the elastic net objective function 
#' 
#' @param cor A vector of correlations, representing the univariate SNP-phenotype associations from the summary statistic data.
#' @param bfile The stem of a PLINK binary file to be used as a reference panel.
#' @param lambdas Vector of values of lambda.
#' @param alphas Vector of values of alpha.
#' @param s Vector of values for s.
#' @param thr Convergence threshold.
#' @param init Initial values for \eqn{\beta}.
#' @param maxIter Maximum number of iterations.
#' @param extract Vector of the indices of SNPs to keep. If null, will keep all SNPs.
#' @export

elastSum <- function(cor, bfile, lambdas, alphas, s=0.5, thr=1e-4, init=NULL, maxIter=1000, extract=NULL){
  
  fam=read.table(paste0(bfile,".fam"))
  N=nrow(fam)
  
  bim=read.table(paste0(bfile,".bim"))
  P=nrow(bim)
  
  genoMat=genotypeMatrix(paste0(bfile,".bed"),N=N,P=P,integer(0),integer(0),integer(0),integer(0),1)
  
  sds=normalize2(genoMat)
  genoMat[is.na(genoMat)]=0
  
  lassoBetasFull=NULL
  yhatFull=NULL
  converged=NULL
  lambdas=lambdas
  alphas=alphas
  maxIter=maxIter
  s=s
  
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

  for(l in 1:length(s)){
    genoMat=genoMat*sqrt(sMult[l])
    diagVec=rep(1-s[l],ncol(genoMat))
    diagVec[sds==0]=0
    denom=diagVec+s[l]
    for(i in 1:length(lambdas)){
      maxDiff=Inf
      iterator=0
      lambda=lambdas[i]
      alpha=alphas[i]
      lassoBetas=init
      yhat=genoMat%*%lassoBetas
      while(maxDiff>thr & iterator<maxIter){
        lassoBetasOld=lassoBetas
        for(j in 1:length(lassoBetas)){
          xj=lassoBetas[j]
          lassoBetas[j]=0
          tempEst=diagVec[j]*xj + cor[j] - sum(genoMat[,j]*yhat)
          if(sds[j]==0) tempEst=0
          if(abs(tempEst)>(lambda*alpha)){
            if(tempEst>0){newEst=(tempEst-(lambda*alpha))/denom[j]}
            if(tempEst<0){newEst=(tempEst+(lambda*alpha))/denom[j]}
          }
          if(abs(tempEst)<(lambda*alpha)){newEst=0}
          lassoBetas[j]=newEst/(1+lambda*(1-alpha))
          if(denom[j]==0){lassoBetas[j]=0}
          if(lassoBetas[j]!=xj){
            del=lassoBetas[j]-xj
            yhat=yhat+(del*genoMat[,j])
          }
        }
        iterator=iterator+1
        maxDiff=max(abs(lassoBetasOld-lassoBetas))
      }
      converged=rbind(converged,c(lambda,alpha,s[l],iterator!=maxIter))
      lassoBetasFull=cbind(lassoBetasFull,c(lambda, alpha, s[l], lassoBetas))
      yhatFull=cbind(yhatFull,c(lambda,alpha,s[l],yhat))
    }
  }
  toReturn=structure(list(lambdas=lambdas,alphas=alphas,beta=lassoBetasFull[4:nrow(lassoBetasFull),],converged=converged,pred=yhatFull[4:nrow(yhatFull),]))
  return(toReturn)
  #' @return a list with the following
  #' \item{lambdas} Same as lambdas input.
  #' \item{alphas} Same as alphas input.
  #' \item{beta} Matrix of estimated coefficients.
  #' \item{converged} Matrix of convergence indicators with following format: column 1 is lambda value, column 2 is alpha value, column 3 is s value, column 4 is 1 if there was convergence.
  #' \item{pred} Matrix of predicted phenotypes.
}