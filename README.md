# Bayesian-Run-estimates
The program reads historical daily weir-tower escapement data and estimateds missing passages, applying Hierarchical Bayesian apporach. 
In this model daily esscapement/run is modeled as log-noraml n(t) = a*(0.5*((ln(t) - ln(mu))/b)^2), where a, mu, b are hierachical parameters. 
The model will estimate missing daily passage as well as 95% Credible Interval. 
The code are written for WINBUS/OpenBUGS, JAG, and RSTAN.  


