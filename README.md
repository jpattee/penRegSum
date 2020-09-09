# penRegSum
'penRegSum' is an R package that implements penalized regression methods for genetic data using summary statistics and a reference panel. This R package can be thought of as a close companion to the LassoSum package: https://github.com/tshmak/lassosum. The LassoSum documentation may prove helpful in implementing our package, as the use cases and requirements are very similar.

penRegSum implements TLP and elastic net penalized regression models within the framework of summary statistics and reference panels. Note that the lassoSum model is a special case of the elastic net model. This package also implements the pseudo-AIC and pseudo-BIC methods for model selection in situations where validation data is not available. Also included is an implementation of the quasi-correlation method, which assesses the predictive accuracy of a polygenic risk score applied to out-of-sample data comprised of summary statistics.

The tlpSum and elastSum functions require the following inputs: a vector of SNP-wise correlations, a stem of a PLINK binary file to use as a reference panel, one or more values for the tuning parameter &lambda;, and one or more values for the tuning parameter &tau;. Optional parameters include one or more values of tuning parameter s (default is .5), a convergence threshold that determines when the coordinate descent algorithm terminates (defaults to 1e-4), initialization for the effect estimates (if a 'warm start' is desired), maximum number of iterations before the algorithm declares failure to converge (default is 1000), an index of SNPs to extract from the reference panel (defaults to all SNPs), a '.bed' file format denoting LD blocks to use in estimation (default is by chromosome), and a PLINK '.bim' file denoting the SNPs in the 'cor' vector, which is used to determine if any SNPs are missing from the reference panel. Output is a list with the following times (and a brief description):  
lambdas: list of &lambda; values corresponding to the columns of 'beta' and 'pred'. Not the same as the input.  
taus (or alphas) : list of &tau; (or &alpha;) values corresponding to the columns of 'beta' and 'pred'. Not the same as the input.  
s: list of s values corresponding to the columns of 'beta' and 'pred'. Not the same as the input.  
beta: matrix of effect size estimates for each of the SNPs. Each column corresponds to a unique combination of tuning parameters. Beta values are scaled as correlations - i.e. they correspond to a standardized y and x.  
converged: matrix of values representing whether the model converged for a given set of tuning parameters. Columns 1-3 display tuning parameter values, column for is an indicator that equals 1 if the model converged.  
pred: matrix of predicted phenotypes for the tuning data. Note that these are appropriately scaled, i.e. assuming minor allele coding for SNPs and rescaled (or 'unstandardized') betas.

The pseudoAicBic function requires the following inputs: a polygenic risk score model, univariate &beta; and corresponding standard error estimates for the training data, the size N of the training data (can be a vector of SNPs come from analyses of different sample size), and a stem of a PLINK binary file to use as a reference panel. Optional parameters include controls for the penalization of the estimated covariance matrix in estimation of SSE and the residual variance (default is to add .2 to the diagonal for both), an index of SNPs to use for residual variance estimation (defaults to most strongly associated n/5 SNPs, where n is the sample size of the reference panel), vectors of reference and alternative alleles for the polygenic risk score and training betas (if you want the program to handle flipped alleles for you), a true/false indicator that denotes whether the polygenic risk score is scaled as correlations (if 'TRUE', program will rescale estimates), a vector of SNP indices to extract from the reference panel, and an LD blocks file in '.bim' format.

The quasicors function requires the following inputs: effect size estimates and SEs from the 'testing' data, sample size of the testing data, a polygenic risk score model, and a reference panel (specified again as the stem of a PLINK binary file). Optional parameters include allele vectors specifying the reference and alternative alleles for the polygenic risk score, a standardization indicator, the variance of the phenotype for the training data (required if standardized == TRUE), and allele vectors specifying reference and alternative alleles for the testing data.

Note that much of the C++ code in this package is repurposed from the LassoSum R package. Generous credit to the authors of the LassoSum package is reiterated here. 
