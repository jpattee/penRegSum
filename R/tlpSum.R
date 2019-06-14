#' @title Function to obtain TLP effect size estimates given summary statistics and a reference panel
#' 
#' @details Function finds the minimum of the TLP objective function 
#' 
#' @param cor A vector of correlations, representing the univariate SNP-phenotype associations from the summary statistic data.
#' @param bfile The stem of a PLINK binary file to be used as a reference panel.
#' @param lambdas Vector of values of lambda.
#' @param taus Vector of values of tau.
#' @param s Vector of values for s.
#' @param thr Convergence threshold.
#' @param init Initial values for \eqn{\beta}.
#' @param maxIter Maximum number of iterations.
#' @param extract Vector of the indices of SNPs to keep. If null, will keep all SNPs.
#' @export

tlpSum <- function(cor, bfile, lambdas, taus, s=0.5, thr=1e-4, init=NULL, maxIter=1000, extract=NULL){

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
  taus=taus
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

  betas=cor

  for(k in 1:length(s)){
    genoMat=genoMat*sqrt(sMult[k])
    diagVec=rep(1-s[k],ncol(genoMat))
    diagVec[sds==0]=0
    denom=diagVec+s[k]
    for(i in 1:length(lambdas)){
      for(l in 1:length(taus)){
        maxDiff=Inf
        iterator=0
        lambda=lambdas[i]
        tau=taus[l]
        lassoBetas=init
        yhat=genoMat%*%lassoBetas
        while(maxDiff>1e-4 & iterator<maxIter){
          iterator=iterator+1
          lassoBetasOld=lassoBetas
          lambdaWts=abs(lassoBetas)<=tau
          for(j in c(1:length(lassoBetas))){
            lambdaTmp=lambda*lambdaWts[j]
            if(iterator==1) lambdaTmp=lambda
            xj=lassoBetas[j]
            lassoBetas[j]=0
            tempEst=diagVec[j]*xj + betas[j] - sum(genoMat[,j]*yhat)
            if(sds[j]==0) tempEst=0
            if(abs(tempEst)>lambdaTmp){
              tempSign=sign(tempEst)
              tempVal=abs(tempEst)-lambdaTmp
              temp2=tempVal*tempSign
              lassoBetas[j]=temp2
            }
            if(lassoBetas[j]!=xj){
              del=lassoBetas[j]-xj
              yhat=yhat+(del*genoMat[,j])
            }
          }
          maxDiff=max(abs(lassoBetasOld-lassoBetas))
        }
        converged=rbind(converged,c(lambda,taus[l],s[k],iterator!=maxIter))
        lassoBetasFull=cbind(lassoBetasFull,c(lambda,tau,s[k],lassoBetas))
        yhatFull=cbind(yhatFull,c(lambda,tau,s[k],yhat))
      }
    }
  }
  toReturn=structure(list(lambdas=lambdas,taus=taus,beta=lassoBetasFull[4:nrow(lassoBetasFull),],converged=converged,pred=yhatFull[4:nrow(yhatFull),]))
  return(toReturn)
  #' @return a list with the following
  #' \item{lambdas} Same as lambdas input.
  #' \item{taus} Same as taus input.
  #' \item{beta} Matrix of estimated coefficients.
  #' \item{converged} Matrix of convergence indicators with following format: column 1 is lambda value, column 2 is tau value, column 3 is s value, column 4 is 1 if there was convergence.
  #' \item{pred} Matrix of predicted phenotypes.
}

