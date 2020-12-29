# quiQI_simulations
Some simulations for the quiQI project

The point of the code is to simulate some signals the try to estimate its variance-covariance components and the parameters of the GLM. Two cases are considered:

- one-sample t-test with an age regressor, i.e. one group and one regressor;
- two-sample t-test without any regressor, i.e. two groups with unequal variance.

There are a few questions that need to be answered:

1. when calculating the WLS solution of the GLM, the residuals are not mean centered anymore. Is that an issue?
2. which components should be used for the modelization of the variance-covariance matrix? And do the rescaling and/or setting the minimal value of each component to zero matter?
3. what about the positivity constraints on covariance parameters in ReML?

