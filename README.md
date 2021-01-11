# quiQI_simulations
Some simulations for the quiQI project

The point of the code is to simulate some signals the try to estimate its variance-covariance components and the parameters of the GLM. Two cases are considered:

- one-sample t-test with an age regressor, i.e. one group and one regressor;
- two-sample t-test without any regressor, i.e. two groups with unequal variance.

### Questions

There are a few questions that need to be answered:

1. when calculating the WLS solution of the GLM, the residuals are not mean centred anymore. Is that an issue?
2. which components should be used for the modelling of the variance-covariance matrix? And do the rescaling and/or setting the minimal value of each component to zero matter?
3. what about the positivity constraints on covariance parameters in ReML?
4. what is the model linking the quality index (QI) and the individual signal?

### Answers?

Some answers to the above questions?

#### For #1. 
In SPM after the estimation of the GLM, with the whitening of data and design matrix, the residuals are not mean centred anymore, so it looks like this is not important.
#### For #3. 
The SPM code in `spm_est_non_sphericity`  relies on  `spm_reml`, i.e. NO positivity constraints.
#### For #2. 
There is some form of regularisation in `spm_reml`, so not sure about the effect of scaling some components or setting the smallest value to zero. This will really depend on how we think the QIndex should account for the variance in the residuals...

#### For #4.

Here in the simulation, the assumption is that the a random value, drawn from a 0-centred Normal distribution, is added to the true signal and the variance of this N-distribution is proportionate to the QI. Thus this should "just" mean noisier data for subjects with larger QI but no biasing effect.

:arrow_forward: Could there be another way of representing this? Is the normal distribution a good assumption?

---

## Further comments

Following some quick checks on the true signal/residuals:

1. When creating the noisy data, from some true (=noise free) signal, we add some random numbers drawn from a mean-centered normal distribution. This constitute the "true residuals" but this does NOT imply that the mean of the "true residuals" is exactly zero. Of course, if we generate a very large dataset, this mean will tend to zero... So we may be confusing the estimate of the residual mean with its true mean.
2. When we calculate the WLS of our problem, i.e. an OLS of the whitened data and design matrix,  we get some estimates of beta. These `$\hat{\beta}$` are for our original problem, so why bother with the whitened one? And for the residuals, the true one is not centered around zero anyway, so here we are just hoping to get a better estimate of these residuals.
3. Here I/we should probably have a look at the significance of the estimated parameters and check that it improves when we better handle/model the residuals variance.