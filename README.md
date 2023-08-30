# CPSR
Heterogeneity and high-dimensionality of time-series datasets bring computational burden for model fitting. This article introduces a regression model combining heterogeneity with homogeneity 
to accommodate the change points and time-invarying high-dimensional sparse parameters. A double penalization strategy leveraging the adaptive group fused lasso and adaptive lasso technique was 
developed to detect the Change Points and achieve Sparsity Recovery (CPSR). The proposed approach accelerates the computational procedure of estimating the number and positions of change points 
and identifying significant covariates for varaible selection. Theoretical analyses confirm that the estimated proportions of structural breaks converge towards the true values with high probability. 
Furthermore, this work also explores the asymptotic distribution of the post-lasso estimators of all regression parameters. The comprehensive Monte Carlo simulations showcase the efficacy of the 
proposed methodology in finite sample scenarios, and a real data example is provided for illustration. The application has been coded into an R package named ``CPSR'' to facilitate wide use.

\section{Installation}
