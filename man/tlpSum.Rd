% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tlpSum.R
\name{tlpSum}
\alias{tlpSum}
\title{Function to obtain TLP effect size estimates given summary statistics and a reference panel}
\usage{
tlpSum(
  cors,
  bfile,
  lambdas,
  taus,
  s = 0.5,
  thr = 1e-04,
  init = NULL,
  maxIter = 1000,
  extract = NULL,
  ldBlocks = NULL,
  corBim = NULL
)
}
\arguments{
\item{cors}{A vector of correlations, representing the univariate SNP-phenotype associations from the summary statistic data.}

\item{bfile}{The stem of a PLINK binary file to be used as a reference panel.}

\item{lambdas}{Vector of values of lambda.}

\item{taus}{Vector of values of tau.}

\item{s}{Vector of values for s.}

\item{thr}{Convergence threshold.}

\item{init}{Initial values for \eqn{\beta}.}

\item{maxIter}{Maximum number of iterations.}

\item{extract}{Vector of the indices of SNPs to keep. If null, will keep all SNPs.}

\item{ldBlocks}{Location of file specifying independent LD Blocks to be used. File should be in the BED file format. If null, estimation is done by chromosome.}

\item{corBim}{Optional PLINK formatted '.bim' file that specifies SNP labels and locations for 'cors' vector. If non-null, will check for SNPs missing in the reference panel and handle ensuing estimation.}
}
\description{
Function to obtain TLP effect size estimates given summary statistics and a reference panel
}
\details{
Function finds the minimum of the TLP objective function
}
