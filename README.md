# penRegSum
'penRegSum' is an R package that implements penalized regression methods for genetic data using summary statistics and a reference panel. This R package can be thought of as a close companion to the LassoSum package: https://github.com/tshmak/lassosum. The LassoSum documentation may prove helpful in implementing our package, as the use cases and requirements are very similar.

penRegSum implements TLP and elastic net penalized regression models within the framework of summary statistics and reference panels. Note that the lassoSum model is a special case of the elastic net model. This package also implements the pseudo-AIC and pseudo-BIC methods for model selection in situations where validation data is not available. Also included is an implementation of the quasi-correlation method, which assesses the predictive accuracy of a polygenic risk score applied to out-of-sample data comprised of summary statistics.

Note that much of the C++ code in this package is repurposed from the LassoSum R package.
