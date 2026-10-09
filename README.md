# Effect of a feedback loop on clinical prediction model performance

This repository contains the code and results for the paper "Effect of a feedback loop on clinical prediction model performance" by Samantha Pacynko, Matthew Sperrin, and David A. Jenkins.

## Abstract

Clinical prediction models (CPMs) that guide treatment-initiation, and therefore affect patient outcomes, may create feedback loops when updated using data from patients whose care they influenced. Treatment can alter predictor-outcome relationships and outcome prevalence, potentially reducing model performance. This simulation study aims to quantify this effect using setting parameters based on cardiovascular disease prediction and intervention and assess if including a treatment variable in the model mitigates any performance degradation. In the simulation, models were developed to predict the probability of a binary adverse outcome, guide treatment allocation, and then updated using data containing patients in whom the model was applied. Scenarios varied treatment inclusion as a predictor in the model and the presence of unmeasured confounding. Performance was assessed with calibration-in-the-large, calibration slope, c-statistic, and correct risk assignment relative to treatment threshold. Without unmeasured confounding, excluding treatment as a predictor resulted in miscalibration, but including treatment solves this problem completely. With unmeasured confounding, performance deteriorated even when treatment was included. These results agree with previous findings and illustrate the likelihood that models that guide interventions may experience performance deterioration due to the difficulty of recording treatment and avoiding unmeasured confounders. Causal inference methods combined with risk prediction may help address this problem. Further research is needed on monitoring techniques to allow detection of this issue post-model implementation.

## Code and results

### Update methods

There are three folders, split by model update method:
- Full refit
- Recalibration
- Intercept

These folders each contain the code and associated results for a different model update method. For more details on these methods, please refer to the following paper:

Su T-L, Jaki T, Hickey GL, Buchan I, Sperrin M. A review of statistical updating methods for clinical prediction models. Statistical Methods in Medical Research. 2018;27(1):185-197. doi:10.1177/0962280215626466

In the section "2.1 Regression coefficients updating", the following methods correspond to our simulation in the above folders:
- a - intercept only
- b - recalibration
- e & f - full refit

With regards to setting up and running the simulation, there is no difference between these update methods.

### Setting parameters

#### Overview

The following parameters can be altered:
- prop.y - the prevalence of the outcome in an untreated population
- x1.mean.original - the mean value of X, a continuous predictor for the outcome
- x1.sd - the standard deviation of X, a continuous predictor for the outcome
- b1 - the coefficient for X
- theta.original - the risk threshold for treatment allocation
- risk - the average risk reduction observed when all individuals receive treatment
- n - the number of individuals per cycle
- c - the number of cycles
- r - the number of repetitions per scenario
- u.prob.original - the prevalence of individuals with the binary confounder U. If the unmeasured confounder is not to be included, this can be set to 0
- u.coef - the coefficient for binary confounder U
- u0.below.prob.original - the probability of someone receiving treatment if they do not have U and their predicted risk is below the threshold
- u0.above.prob.original - the probability of someone receiving treatment if they do not have U and their predicted risk is above the threshold
- u1.below.prob.original - the probability of someone receiving treatment if they do have U and their predicted risk is below the threshold
- u1.above.prob.original - the probability of someone receiving treatment if they do have U and their predicted risk is above the threshold


Data drift can be introduced to the simulation with the following parameters:
- x1.adjustment - the mean value of X will change by this value at the beginning of each cycle, before the next set of data is generated. This will naturally lead to a change in prevalence 
- theta.adjustment - the risk threshold for treatment will change by this value at the beginning of each cycle
- u.prob.adjustment - the prevalence of individuals with the binary confounder U will change by this value at the beginning of each cycle. This will naturally lead to a change in prevalence
- u.treatment.drift - controls if drift is included for the relationship between U and effect of treatment. Set to 1 if this drift should be included or 0 if it should be excluded. If it is excluded, the probability of individuals receiving treatment is equal to the the treatment probabilities previously set. If it is included, the probability of individuals receiving treatment will be equal to the u0 values for all individuals; the probability of individuals with U receiving treatment will then incrementally change until they reach the specified u1 values at the last cycle.

#### Defaults - CVD example

Here are the default values used for the above parameters:
- prop.y - 0.07
- b1 - 1
- x1.mean.original - 0
- x1.adjustment - 0 if no drift; 0.1 if drift
- x1.sd - 1
- theta.original - 0.1
- theta.adjustment - 0 if no drift; 0.02 if drift
- risk - 0.25
- n - 10000
- c - 10
- r - 200
- u.prob.original - 0 if no unmeasured confounding; 0.2 if present
- u.prob.adjustment - 0 if no drift; -0.01 if drift
- u.coef - 0.5
- u0.below.prob.original - 0 (this should be kept as 0 for code to work)
- u0.above.prob.original - 1 if there is no U; 0.8 if there is U
- u1.below.prob.original - 0.5
- u1.above.prob.original - 1 (this should be kept as 0 for code to work)

#### Defaults - full results

Here are the default values used for the above parameters:
- prop.y - 0.1; 0.2; 0.3
- b1 - 1
- x1.mean.original - 0
- x1.adjustment - 0 if no drift; 0.1 if drift
- x1.sd - 1
- theta.original - 0.1; 0.2; 0.3
- theta.adjustment - 0 if no drift; 0.02 if drift
- risk - 0.1; 0.2; 0.3; 0.4; 0.5
- n - 10000
- c - 10
- r - 200
- u.prob.original - 0 if no unmeasured confounding; 0.2 if present
- u.prob.adjustment - 0 if no drift; -0.01 if drift
- u.coef - 0.5
- u0.below.prob.original - 0 (this should be kept as 0 for code to work)
- u0.above.prob.original - 1 if there is no U; 0.8 if there is U
- u1.below.prob.original - 0.5
- u1.above.prob.original - 1 (this should be kept as 0 for code to work)

### Results 

Results are stored in the folders with matching code. The name of the results file indicates the update method, presence/absence of U, and if drift was included.
E.g. full.results.refit.U.thetaDrift contains the results when full model refit is used as the update method, unmeasured confounding is included in the simulation, and data drift for the risk threshold for treatment theta is included.

Full results files contains results from all repetitions for each scenario. Summary full results files contain processed results across the repetitions, with means and 95% confidence intervals for the performance measures.

Results are split by model type, where "CPM" indicates the model did not include treatment status as a predictor variable and "CPM with treatment variable" indicates it did.

Results are also split by timing of performance measures. The performance measures reported in the paper are the pre-treatment performance measures, where they are calculated using only the untreated outcomes for individuals. Post-treatment performance measures are also included here, where they are calculated using the outcomes after treatment is allocated to some individuals based on their predicted risk.
