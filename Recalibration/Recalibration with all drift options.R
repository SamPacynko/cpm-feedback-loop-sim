#####################################
## Title: Clinical Prediction Model Feedback Loop Simulation
## Author: Samantha Pacynko
## Date Created: 11/05/2024
## Date Edited: 09/10/2026
#######################################

# Set up ----

## Libraries ----
library(rje) # Expit function for ease
library(tidyverse) # Data summarising and str_replace_all
library(pROC) # AUC / c-statistic

# Create variables ----
## Parameters ----
prop.y <- c(0.1,0.2,0.3) # Prevalence of outcome
b1 <- 1 # Coefficient of X1 in data generating mechanism
x1.mean.original <- 0 # Mean of X1
x1.adjustment <- 0 # Adjustment of x1 mean - adjustment at the beginning of each cycle
x1.sd <- 1 # SD of X1
theta.original <- c(0.1,0.2,0.3) # Risk threshold for treatment
theta.adjustment <- 0 # Adjustment of treatment threshold - adjustment at the beginning of each cycle
risk <- c(0.1,0.2,0.3,0.4,0.5) # Relative risk reduction with treatment
n <- 10000 # Number of patients per cycle
c <- 10 # Number of cycles
r <- 200 # Repetitions

u.prob.original <- 0 # Probability of confounding
u.prob.adjustment <- 0 # Adjustment of probability of confounding - adjustment at the beginning of each cycle
u.coef <- 0.5 # Coefficient of confounding

# Probabilities of individuals receiving treatment, depending on individual U and if estimated risk is above or below threshold
# If you want there to be no effect of U, u0.below/u0.above and u1.below/u1.above should be equal
u0.below.prob.original <- 0 # Should keep as 0 for current code to work
u0.above.prob.original <- 1
u1.below.prob.original <- 0.5
u1.above.prob.original <- 1 # Should keep as 1 for current code to work

u.treatment.drift <- 0 # Include drift for relationship between U and treatment? 1 for yes

# When drift is included, probability of treatment between two groups of U is equal to start (based on U0)
# Then probabilities for U1 change incrementally to set levels
u.below.change <- (u1.below.prob.original-u0.below.prob.original)/c
u.above.change <- (u1.above.prob.original-u0.above.prob.original)/c

## CPM.t vector ----
mod_c <- c(0,0,0) #Set up vector for manual CPM with treatment variable

## Results table ----
# Data frame to store performance measures from each scenario, repetition, and cycle

full.results <- data.frame(
  prop.y = numeric(),
  theta = numeric(),
  risk = numeric(),
  r = integer(),
  c = integer(),
  model = character(),
  timing = character(),
  CITL = numeric(),
  slope = numeric(),
  c.stat = numeric(),
  brier = numeric(),
  actualYcount = numeric(),
  untreatedYcount = numeric(),
  treatedYcount = numeric(),
  treatedn = numeric(),
  TP = numeric(),
  FP = numeric(),
  FN = numeric(),
  TN = numeric(),
  sensitivity = numeric(),
  specificity = numeric(),
  accuracy = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  b0_true = numeric(),
  b1_true = numeric(),
  b2_true = numeric(),
  b0_estimate = numeric(),
  b1_estimate = numeric(),
  b2_estimate = numeric()
)

## CPM.t functions ----

#Linear predictor function
cpm.t.LP <- function(predict_data){
  mod_c[1]+(mod_c[2]*predict_data$x1)+(mod_c[3]*predict_data$t)
}

#Estimated risk
cpm.t.response <- function(predict_data){
  expit(mod_c[1]+(mod_c[2]*predict_data$x1)+(mod_c[3]*predict_data$t))
}

# Simulation ----

for (i in 1:length(prop.y)) {
  for (j in 1:length(theta.original)) {
    for (k in 1:length(risk)) {
      ## Scenario start ----
      
      set.seed(1) # Set seed so that each scenario is based on same seed
      
      for (l in 1:r) {
        ## Repetition start ----
        
        ### Generate data ----
        
        #Reset mod_c
        mod_c <- c(0,0,0)

        #Generate dataset for all cycles in this repetition of this scenario
        theta <- theta.original[j]
        
        x1.mean <- x1.mean.original
        
        x1 <-
          rnorm(n, mean = x1.mean, sd = x1.sd) # X1 from standard normal distribution
        
        u.prob <- u.prob.original
        
        u <- rbinom (n, 1, u.prob)
        
        refuse <- rbinom(n, 1, 0.1)
        data <- data.frame(x1, u)
        b0 <- as.numeric(coef(glm(
          rbinom(n, 1, prop.y[i]) ~ offset((b1 * x1) + (u.coef * u)),
          family = binomial(link = "logit"),
          data = data
        ))[1]) # Based on values of X1, calculate B0 value to be used to give desired prop.y
        
        b2 <- as.numeric(coef(glm(
          rbinom(n, 1, prop.y[i] * (1 - risk[k])) ~ offset(b0 + (b1 * x1) + (u.coef * u)),
          family = binomial(link = "logit"),
          data = data
        ))[1])
        
        b2 <- ifelse(b2>0, 0, b2)
        
        untreated.risk <-
          expit(b0 + (b1 * x1) + (u.coef * u) + (b2 * 0)) # Risk of outcome using data generating mechanism when patient is untreated
        treated.risk <-
          expit(b0 + (b1 * x1) + (u.coef * u) + (b2 * 1)) # Risk of outcome using data generating mechanism when patient is treated
        untreated.y <- rbinom(n, 1, untreated.risk) # Outcome when patient is treated
        treated.y <- ifelse(untreated.y == 0, 0, rbinom(n, 1, (treated.risk /
                                                                 untreated.risk))) # Outcome when patient is untreated
        
        if (u.treatment.drift == 1) {
          u0.below.prob <- u0.below.prob.original
          u0.above.prob <- u0.above.prob.original
          u1.below.prob <- u0.below.prob.original
          u1.above.prob <- u0.above.prob.original
        } else {
          u0.below.prob <- u0.below.prob.original
          u0.above.prob <- u0.above.prob.original
          u1.below.prob <- u1.below.prob.original
          u1.above.prob <- u1.above.prob.original
        }
        
        u0.below.t <- rbinom(n, 1, u0.below.prob)
        u0.above.t <- rbinom(n, 1, u0.above.prob)
        u1.below.t <- rbinom(n, 1, u1.below.prob)
        u1.above.t <- rbinom(n, 1, u1.above.prob)
        
        t <- ifelse(u == 1 &
                      u1.below.t == 1 &
                      refuse == 0, 1, 0) # Change depending if treatment already used based on U or not before model introduced
        
        actual.y <- ifelse(t == 1, treated.y, untreated.y) # Treatment assigned
        
        d <-
          data.frame(
            x1,
            u,
            refuse,
            untreated.risk,
            untreated.y,
            treated.risk,
            treated.y,
            actual.y,
            u0.below.t,
            u0.above.t,
            u1.below.t,
            u1.above.t,
            t
          ) # Combine all of above into one dataset
        
        ### CPM ----
        #### Development ----
        
        #Split baseline group (0) 50/50 into development dataset and test dataset
        sample <- sample.int(
          n = nrow(d),
          size = floor(.5 * nrow(d)),
          replace = F
        )
        d.dev <- d[sample, ]
        d.test  <- d[-sample, ]
        
        #Create copy of d.test for use with CPM.t
        d.test.t <- d.test
        
        #Develop CPM using dev data
        cpm <- glm(actual.y ~ x1, family = binomial, data = d.dev) # Main CPM that doesn't take into account treatment
        
        #### Internal validation ----
        
        #Internally validate main CPM using test data - calculate performance measures
        #Prep
        d.test$predicted.risk <- predict(cpm, newdata = d.test, type = "response") # Calculating predicted risks using CPM
        d.test$LP <- predict(cpm, newdata = d.test) # Calculating linear predictor using coefficients from CPM
        log.cal <- glm(actual.y ~ LP, family = binomial, data = d.test) # Logistic calibration using linear predictor - for slope
        log.cal.offset <- glm(
          actual.y ~ LP,
          family = binomial,
          data = d.test,
          offset = LP
        ) # Logistic calibration using linear predictor, offset by LP - for CITL
        
        #Performance measures
        CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
        slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
        c.stat <- as.numeric(auc(roc(
          d.test$actual.y ~ d.test$predicted.risk, quiet = TRUE
        ))) # C-statistic/area under the curve - greater than 0.5
        brier <- mean((d.test$predicted.risk - d.test$actual.y)^2) # Brier score - 0 is perfect accuracy
        
        #Extra values
        actualYcount <- sum(d.test$actual.y)
        untreatedYcount <- sum(d.test$untreated.y)
        treatedYcount <- sum(d.test$treated.y)
        treatedn <- sum(d.test$t)
        TP <- sum(d.test$untreated.risk >= theta &
                    d.test$predicted.risk >= theta)
        FP <- sum(d.test$untreated.risk < theta &
                    d.test$predicted.risk >= theta)
        FN <- sum(d.test$untreated.risk >= theta &
                    d.test$predicted.risk < theta)
        TN <- sum(d.test$untreated.risk < theta &
                    d.test$predicted.risk < theta)
        sensitivity <- TP / (TP + FN)
        specificity <- TN / (TN + FP)
        accuracy <- (TP + TN) / (TP + TN + FP + FN)
        PPV <- TP / (TP + FP)
        NPV <- TN / (TN + FN)
        b0_true <- b0
        b1_true <- b1
        b2_true <- b2
        b0_estimate <- cpm[["coefficients"]][["(Intercept)"]]
        b1_estimate <- cpm[["coefficients"]][["x1"]]
        b2_estimate <- NA
        
        #Save internal validation baseline results into table - repeated for each performance scenario (pre/post) for ease in analysis/plots
        full.results[nrow(full.results) + 1, ] <-
          list(
            prop.y[i],
            theta.original[j],
            risk[k],
            l,
            0,
            'Main CPM',
            'Pre-treatment',
            CITL,
            slope,
            c.stat,
            brier,
            actualYcount,
            untreatedYcount,
            treatedYcount,
            treatedn,
            TP,
            FP,
            FN,
            TN,
            sensitivity,
            specificity,
            accuracy,
            PPV,
            NPV,
            b0_true,
            b1_true,
            b2_true,
            b0_estimate,
            b1_estimate,
            b2_estimate
          )
        full.results[nrow(full.results) + 1, ] <-
          list(
            prop.y[i],
            theta.original[j],
            risk[k],
            l,
            0,
            'Main CPM',
            'Post-treatment',
            CITL,
            slope,
            c.stat,
            brier,
            actualYcount,
            untreatedYcount,
            treatedYcount,
            treatedn,
            TP,
            FP,
            FN,
            TN,
            sensitivity,
            specificity,
            accuracy,
            PPV,
            NPV,
            b0_true,
            b1_true,
            b2_true,
            b0_estimate,
            b1_estimate,
            b2_estimate
          )
        
        ### CPM.t ----
        #### Development ----
        
        #Update variables in mod_c for CPM with treatment variable (cpm.t)
        mod_c[1] <- cpm[["coefficients"]][["(Intercept)"]]
        mod_c[2] <- cpm[["coefficients"]][["x1"]]
        
        #### Validation ----
        
        #Internally validate CPM.t using test data - calculate performance measures
        #Prep
        d.test.t$predicted.risk <- cpm.t.response(d.test.t) # Calculating predicted risks using CPM.t
        d.test.t$LP <- cpm.t.LP(d.test.t) # Calculating linear predictor using coefficients from CPM.t
        log.cal <- glm(actual.y ~ LP, family = binomial, data = d.test.t) # Logistic calibration using linear predictor - for slope
        log.cal.offset <- glm(
          actual.y ~ LP,
          family = binomial,
          data = d.test.t,
          offset = LP
        ) # Logistic calibration using linear predictor, offset by LP - for CITL
        
        #Performance measures
        CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
        slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
        c.stat <- as.numeric(auc(roc(
          d.test.t$actual.y ~ d.test.t$predicted.risk, quiet = TRUE
        ))) # C-statistic/area under the curve - greater than 0.5
        brier <- mean((d.test.t$predicted.risk - d.test.t$actual.y)^2) # Brier score - 0 is perfect accuracy
        
        #Extra values
        actualYcount <- sum(d.test.t$actual.y)
        untreatedYcount <- sum(d.test.t$untreated.y)
        treatedYcount <- sum(d.test.t$treated.y)
        treatedn <- sum(d.test.t$t)
        TP <- sum(d.test.t$untreated.risk >= theta &
                    d.test.t$predicted.risk >= theta)
        FP <- sum(d.test.t$untreated.risk < theta &
                    d.test.t$predicted.risk >= theta)
        FN <- sum(d.test.t$untreated.risk >= theta &
                    d.test.t$predicted.risk < theta)
        TN <- sum(d.test.t$untreated.risk < theta &
                    d.test.t$predicted.risk < theta)
        sensitivity <- TP / (TP + FN)
        specificity <- TN / (TN + FP)
        accuracy <- (TP + TN) / (TP + TN + FP + FN)
        PPV <- TP / (TP + FP)
        NPV <- TN / (TN + FN)
        b0_true <- b0
        b1_true <- b1
        b2_true <- b2
        b0_estimate <- mod_c[1]
        b1_estimate <- mod_c[2]
        b2_estimate <- mod_c[3]
        
        #Save internal validation baseline results into table - repeated for each performance scenario for ease in analysis/plots
        full.results[nrow(full.results) + 1, ] <-
          list(
            prop.y[i],
            theta.original[j],
            risk[k],
            l,
            0,
            'CPM with treatment variable',
            'Pre-treatment',
            CITL,
            slope,
            c.stat,
            brier,
            actualYcount,
            untreatedYcount,
            treatedYcount,
            treatedn,
            TP,
            FP,
            FN,
            TN,
            sensitivity,
            specificity,
            accuracy,
            PPV,
            NPV,
            b0_true,
            b1_true,
            b2_true,
            b0_estimate,
            b1_estimate,
            b2_estimate
          )
        full.results[nrow(full.results) + 1, ] <-
          list(
            prop.y[i],
            theta.original[j],
            risk[k],
            l,
            0,
            'CPM with treatment variable',
            'Post-treatment',
            CITL,
            slope,
            c.stat,
            brier,
            actualYcount,
            untreatedYcount,
            treatedYcount,
            treatedn,
            TP,
            FP,
            FN,
            TN,
            sensitivity,
            specificity,
            accuracy,
            PPV,
            NPV,
            b0_true,
            b1_true,
            b2_true,
            b0_estimate,
            b1_estimate,
            b2_estimate
          )
        
        for (m in 1:c) {
          # Using data created before, move through cycles, get performance, and update
          ## Cycle start ----
          
          ### Create data for this cycle ----
          
          x1.mean <- x1.mean + x1.adjustment
          
          x1 <- rnorm(n, mean = x1.mean, sd = x1.sd)
          
          u.prob <- u.prob + u.prob.adjustment
          
          u <- rbinom (n, 1, u.prob)
          
          refuse <- rbinom(n, 1, 0.1)
          
          data <- data.frame(x1, u)
          
          untreated.risk <-
            expit(b0 + (b1 * x1) + (u.coef * u) + (b2 * 0)) # Risk of outcome using data generating mechanism when patient is untreated
          treated.risk <-
            expit(b0 + (b1 * x1) + (u.coef * u) + (b2 * 1)) # Risk of outcome using data generating mechanism when patient is treated
          untreated.y <- rbinom(n, 1, untreated.risk) # Outcome when patient is treated
          treated.y <- ifelse(untreated.y == 0, 0, rbinom(n, 1, (treated.risk /
                                                                   untreated.risk))) # Outcome when patient is untreated
          
          
          if (u.treatment.drift == 1) {
            u1.below.prob <- u1.below.prob + u.below.change
            u1.above.prob <- u1.above.prob + u.above.change
          }
          
          if (u1.above.prob > 1) {
            u1.above.prob <- 1
          }
          
          u0.below.t <- rbinom(n, 1, u0.below.prob)
          u0.above.t <- rbinom(n, 1, u0.above.prob)
          u1.below.t <- rbinom(n, 1, u1.below.prob)
          u1.above.t <- rbinom(n, 1, u1.above.prob)
          
          t <- ifelse(u == 1 &
                        u1.below.t == 1 &
                        refuse == 0, 1, 0) # Change depending if treatment already used based on U or not before model introduced
          
          actual.y <- ifelse(t == 1, treated.y, untreated.y) # Treatment assigned
          
          d.cycle <-
            data.frame(
              x1,
              u,
              refuse,
              untreated.risk,
              untreated.y,
              treated.risk,
              treated.y,
              actual.y,
              u0.below.t,
              u0.above.t,
              u1.below.t,
              u1.above.t,
              t
            ) # Combine all of above into one dataset
          
          d.cycle.t <- d.cycle # Create copy of data for use with cpm.t
          
          ### CPM ----
          #### Pre-treatment validation ----
          
          # Validation when all patients are untreated (i.e. before CPM affects treatment decision) & save to table
          #Prep
          d.cycle$predicted.risk <- predict(cpm, newdata = d.cycle, type =
                                              "response") # Calculating predicted risks using CPM
          d.cycle$LP <- predict(cpm, newdata = d.cycle) # Get linear predictor from CPM
          log.cal <- glm(actual.y ~ LP, family = binomial, data = d.cycle) # Logistic calibration using linear predictor - for slope
          log.cal.offset <- glm(
            actual.y ~ LP,
            family = binomial,
            data = d.cycle,
            offset = LP
          ) # Logistic calibration using linear predictor, offset by LP - for CITL
          
          #Performance measures
          CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
          slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
          c.stat <- as.numeric(auc(roc(
            d.cycle$actual.y ~ d.cycle$predicted.risk, quiet = TRUE
          ))) # C-statistic/area under the curve - greater than 0.5
          brier <- mean((d.cycle$predicted.risk - d.cycle$actual.y)^2) # Brier score - 0 is perfect accuracy
          
          #Extra values
          actualYcount <- sum(d.cycle$actual.y)
          untreatedYcount <- sum(d.cycle$untreated.y)
          treatedYcount <- sum(d.cycle$treated.y)
          treatedn <- sum(d.cycle$t)
          TP <- sum(d.cycle$untreated.risk >= theta &
                      d.cycle$predicted.risk >= theta)
          FP <- sum(d.cycle$untreated.risk < theta &
                      d.cycle$predicted.risk >= theta)
          FN <- sum(d.cycle$untreated.risk >= theta &
                      d.cycle$predicted.risk < theta)
          TN <- sum(d.cycle$untreated.risk < theta &
                      d.cycle$predicted.risk < theta)
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          accuracy <- (TP + TN) / (TP + TN + FP + FN)
          PPV <- TP / (TP + FP)
          NPV <- TN / (TN + FN)
          b0_true <- b0
          b1_true <- b1
          b2_true <- b2
          b0_estimate <- cpm[["coefficients"]][["(Intercept)"]]
          b1_estimate <- cpm[["coefficients"]][["x1"]]
          b2_estimate <- NA
          
          # Save results into table
          full.results[nrow(full.results) + 1, ] <-
            list(
              prop.y[i],
              theta.original[j],
              risk[k],
              l,
              m,
              'Main CPM',
              'Pre-treatment',
              CITL,
              slope,
              c.stat,
              brier,
              actualYcount,
              untreatedYcount,
              treatedYcount,
              treatedn,
              TP,
              FP,
              FN,
              TN,
              sensitivity,
              specificity,
              accuracy,
              PPV,
              NPV,
              b0_true,
              b1_true,
              b2_true,
              b0_estimate,
              b1_estimate,
              b2_estimate
            )
          
          #### Allocating treatment ----
          # Using predicted risk from CPM, decide who gets treatment based on risk threshold
          d.cycle$t <- with(
            d.cycle,
            case_when(
              u == 1 & predicted.risk < theta & refuse == 0 ~ u1.below.t,
              u == 1 & predicted.risk >= theta & refuse == 0 ~ u1.above.t,
              u == 0 & predicted.risk < theta & refuse == 0 ~ u0.below.t,
              u == 0 & predicted.risk >= theta & refuse == 0 ~ u0.above.t,
              refuse == 1 ~ 0
            )
          )
          
          d.cycle$actual.y <- ifelse(d.cycle$t == 1,
                                     d.cycle$treated.y,
                                     d.cycle$untreated.y) # Change actual outcome to treated outcome for applicable patients
          
          #### Post-treatment validation ----
          # Validation of main CPM when some patients are treated
          # Prep
          d.cycle$predicted.risk <- predict(cpm, newdata = d.cycle, type =
                                              "response") # Calculating predicted risks using CPM
          d.cycle$LP <- predict(cpm, newdata = d.cycle) # Get linear predictor from CPM
          log.cal <- glm(actual.y ~ LP, family = binomial, data = d.cycle) # Logistic calibration using linear predictor - for slope
          log.intercept <- log.cal[["coefficients"]][["(Intercept)"]] # To be used when updating model
          log.cal.offset <- glm(
            actual.y ~ LP,
            family = binomial,
            data = d.cycle,
            offset = LP
          ) # Logistic calibration using linear predictor, offset by LP - for CITL
          
          #Performance measures
          CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
          slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
          c.stat <- as.numeric(auc(roc(
            d.cycle$actual.y ~ d.cycle$predicted.risk, quiet = TRUE
          ))) # C-statistic/area under the curve - greater than 0.5
          brier <- mean((d.cycle$predicted.risk - d.cycle$actual.y)^2) # Brier score - 0 is perfect accuracy
          
          #Extra values
          actualYcount <- sum(d.cycle$actual.y)
          untreatedYcount <- sum(d.cycle$untreated.y)
          treatedYcount <- sum(d.cycle$treated.y)
          treatedn <- sum(d.cycle$t)
          TP <- sum(d.cycle$untreated.risk >= theta &
                      d.cycle$predicted.risk >= theta)
          FP <- sum(d.cycle$untreated.risk < theta &
                      d.cycle$predicted.risk >= theta)
          FN <- sum(d.cycle$untreated.risk >= theta &
                      d.cycle$predicted.risk < theta)
          TN <- sum(d.cycle$untreated.risk < theta &
                      d.cycle$predicted.risk < theta)
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          accuracy <- (TP + TN) / (TP + TN + FP + FN)
          PPV <- TP / (TP + FP)
          NPV <- TN / (TN + FN)
          b0_true <- b0
          b1_true <- b1
          b2_true <- b2
          b0_estimate <- cpm[["coefficients"]][["(Intercept)"]]
          b1_estimate <- cpm[["coefficients"]][["x1"]]
          b2_estimate <- NA
          
          # Save results into table
          full.results[nrow(full.results) + 1, ] <-
            list(
              prop.y[i],
              theta.original[j],
              risk[k],
              l,
              m,
              'Main CPM',
              'Post-treatment',
              CITL,
              slope,
              c.stat,
              brier,
              actualYcount,
              untreatedYcount,
              treatedYcount,
              treatedn,
              TP,
              FP,
              FN,
              TN,
              sensitivity,
              specificity,
              accuracy,
              PPV,
              NPV,
              b0_true,
              b1_true,
              b2_true,
              b0_estimate,
              b1_estimate,
              b2_estimate
            )
          
          #### Update model ----
          # Update main CPM using treated data using logistic calibration by manually overwriting previous coefficients
          cpm[["coefficients"]][["(Intercept)"]] <- (cpm[["coefficients"]][["(Intercept)"]]*slope) + log.intercept
          cpm[["coefficients"]][["x1"]] <- cpm[["coefficients"]][["x1"]]*slope
          
          ### CPM.t ----
          #### Pre-treatment validation ----
          
          # CPM.t validation when all patients are untreated (i.e. before CPM affects treatment decision)
          #Prep
          d.cycle.t$predicted.risk <- cpm.t.response(d.cycle.t) # Calculating predicted risks using CPM
          d.cycle.t$LP <- cpm.t.LP(d.cycle.t) # Get linear predictor from CPM
          log.cal <- glm(actual.y ~ LP, family = binomial, data = d.cycle.t) # Logistic calibration using linear predictor - for slope
          log.cal.offset <- glm(
            actual.y ~ LP,
            family = binomial,
            data = d.cycle.t,
            offset = LP
          ) # Logistic calibration using linear predictor, offset by LP - for CITL
          
          #Performance measures
          CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
          slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
          c.stat <- as.numeric(auc(roc(
            d.cycle.t$actual.y ~ d.cycle.t$predicted.risk, quiet = TRUE
          ))) # C-statistic/area under the curve - greater than 0.5
          brier <- mean((d.cycle.t$predicted.risk - d.cycle.t$actual.y)^2) # Brier score - 0 is perfect accuracy
          
          #Extra values
          actualYcount <- sum(d.cycle.t$actual.y)
          untreatedYcount <- sum(d.cycle.t$untreated.y)
          treatedYcount <- sum(d.cycle.t$treated.y)
          treatedn <- sum(d.cycle.t$t)
          TP <- sum(d.cycle.t$untreated.risk >= theta &
                      d.cycle.t$predicted.risk >= theta)
          FP <- sum(d.cycle.t$untreated.risk < theta &
                      d.cycle.t$predicted.risk >= theta)
          FN <- sum(d.cycle.t$untreated.risk >= theta &
                      d.cycle.t$predicted.risk < theta)
          TN <- sum(d.cycle.t$untreated.risk < theta &
                      d.cycle.t$predicted.risk < theta)
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          accuracy <- (TP + TN) / (TP + TN + FP + FN)
          PPV <- TP / (TP + FP)
          NPV <- TN / (TN + FN)
          b0_true <- b0
          b1_true <- b1
          b2_true <- b2
          b0_estimate <- mod_c[1]
          b1_estimate <- mod_c[2]
          b2_estimate <- mod_c[3]
          
          # Save results into table
          full.results[nrow(full.results) + 1, ] <-
            list(
              prop.y[i],
              theta.original[j],
              risk[k],
              l,
              m,
              'CPM with treatment variable',
              'Pre-treatment',
              CITL,
              slope,
              c.stat,
              brier,
              actualYcount,
              untreatedYcount,
              treatedYcount,
              treatedn,
              TP,
              FP,
              FN,
              TN,
              sensitivity,
              specificity,
              accuracy,
              PPV,
              NPV,
              b0_true,
              b1_true,
              b2_true,
              b0_estimate,
              b1_estimate,
              b2_estimate
            )
          
          #### Allocating treatment ----
          # Using predicted risk from CPM.t, decide who gets treatment based on risk threshold
          d.cycle.t$t <- with(
            d.cycle.t,
            case_when(
              u == 1 & predicted.risk < theta & refuse == 0 ~ u1.below.t,
              u == 1 & predicted.risk >= theta & refuse == 0 ~ u1.above.t,
              u == 0 & predicted.risk < theta & refuse == 0 ~ u0.below.t,
              u == 0 & predicted.risk >= theta & refuse == 0 ~ u0.above.t,
              refuse == 1 ~ 0
            )
          )
          
          d.cycle.t$actual.y <- ifelse (d.cycle.t$t == 1,
                                        d.cycle.t$treated.y,
                                        d.cycle.t$untreated.y) # Change actual outcome to treated outcome for applicable patients
          
          #### Post-treatment validation ----
          # Validation of CPM.t when some patients are treated
          # Prep
          d.cycle.t$predicted.risk <- cpm.t.response(d.cycle.t) # Calculating predicted risks using CPM
          d.cycle.t$LP <- cpm.t.LP(d.cycle.t) # Get linear predictor from CPM
          log.cal <- glm(actual.y ~ LP, family = binomial, data = d.cycle.t) # Logistic calibration using linear predictor - for slope
          log.intercept <- log.cal[["coefficients"]][["(Intercept)"]] # To be used when updating model
          log.cal.offset <- glm(
            actual.y ~ LP,
            family = binomial,
            data = d.cycle.t,
            offset = LP
          ) # Logistic calibration using linear predictor, offset by LP - for CITL
          
          #Performance measures
          CITL <- log.cal.offset[["coefficients"]][["(Intercept)"]] # Calibration-in-the-large - should be equal to 0
          slope <- log.cal[["coefficients"]][["LP"]] # Calibration slope - should be equal to 1
          c.stat <- as.numeric(auc(roc(
            d.cycle.t$actual.y ~ d.cycle.t$predicted.risk, quiet = TRUE
          ))) # C-statistic/area under the curve - greater than 0.5
          brier <- mean((d.cycle.t$predicted.risk - d.cycle.t$actual.y)^2) # Brier score - 0 is perfect accuracy
          
          #Extra values
          actualYcount <- sum(d.cycle.t$actual.y)
          untreatedYcount <- sum(d.cycle.t$untreated.y)
          treatedYcount <- sum(d.cycle.t$treated.y)
          treatedn <- sum(d.cycle.t$t)
          TP <- sum(d.cycle.t$untreated.risk >= theta &
                      d.cycle.t$predicted.risk >= theta)
          FP <- sum(d.cycle.t$untreated.risk < theta &
                      d.cycle.t$predicted.risk >= theta)
          FN <- sum(d.cycle.t$untreated.risk >= theta &
                      d.cycle.t$predicted.risk < theta)
          TN <- sum(d.cycle.t$untreated.risk < theta &
                      d.cycle.t$predicted.risk < theta)
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          accuracy <- (TP + TN) / (TP + TN + FP + FN)
          PPV <- TP / (TP + FP)
          NPV <- TN / (TN + FN)
          b0_true <- b0
          b1_true <- b1
          b2_true <- b2
          b0_estimate <- mod_c[1]
          b1_estimate <- mod_c[2]
          b2_estimate <- mod_c[3]
          
          # Save results into table
          full.results[nrow(full.results) + 1, ] <-
            list(
              prop.y[i],
              theta.original[j],
              risk[k],
              l,
              m,
              'CPM with treatment variable',
              'Post-treatment',
              CITL,
              slope,
              c.stat,
              brier,
              actualYcount,
              untreatedYcount,
              treatedYcount,
              treatedn,
              TP,
              FP,
              FN,
              TN,
              sensitivity,
              specificity,
              accuracy,
              PPV,
              NPV,
              b0_true,
              b1_true,
              b2_true,
              b0_estimate,
              b1_estimate,
              b2_estimate
            )
          
          #### Update model ----
          # Update CPM.t using treated data
          if (m==1) { # For first cycle only, fit GLM to introduce treatment variable and get first coefficient
            
            log.cal.t <- glm(actual.y~LP+t, family=binomial, data=d.cycle.t)
            mod_c[1] <- log.cal.t[["coefficients"]][["(Intercept)"]] + (log.cal.t[["coefficients"]][["LP"]]*mod_c[1])
            mod_c[2] <- log.cal.t[["coefficients"]][["LP"]]*mod_c[2]
            mod_c[3] <- log.cal.t[["coefficients"]][["t"]]
            
          } else { # For all other cycles, update CPM using same logistic calibration method as used with the main CPM
            
            # Update CPM.t using treated data using logistic calibration by manually overwriting previous coefficients
            mod_c[1] <- (mod_c[1]*slope) + log.intercept
            mod_c[2] <- mod_c[2]*slope
            mod_c[3] <- mod_c[3]*slope
          }
          
          theta <- theta + theta.adjustment
          
        }
      }
    }
  }
}

# Summarise full results ----
summary.full.results <- full.results %>% # Summarise full results with means and CIs for each scenario
  select(
    prop.y,
    theta,
    risk,
    c,
    model,
    timing,
    CITL,
    slope,
    c.stat,
    brier,
    actualYcount,
    untreatedYcount,
    treatedYcount,
    treatedn,
    TP,
    FP,
    FN,
    TN,
    sensitivity,
    specificity,
    accuracy,
    PPV,
    NPV,
    b0_true,
    b1_true,
    b2_true,
    b1_estimate,
    b2_estimate,
    b0_estimate
  ) %>%
  group_by(prop.y, theta, risk, c, model, timing) %>%
  summarise(
    CITL_Mean = mean(CITL),
    CITL_UC = CITL_Mean + (1.96 * (sd(CITL) / sqrt(n(
    )))),
    CITL_LC = CITL_Mean - (1.96 * (sd(CITL) / sqrt(n(
    )))),
    slope_Mean = mean(slope),
    slope_UC = slope_Mean + (1.96 * (sd(slope) / sqrt(n(
    )))),
    slope_LC = slope_Mean - (1.96 * (sd(slope) / sqrt(n(
    )))),
    c.stat_Mean = mean(c.stat),
    c.stat_UC = c.stat_Mean + (1.96 * (sd(c.stat) / sqrt(n(
    )))),
    c.stat_LC = c.stat_Mean - (1.96 * (sd(c.stat) / sqrt(n(
    )))),
    brier_Mean = mean(brier),
    brier_UC = brier_Mean + (1.96 * (sd(brier) / sqrt(n(
    )))),
    brier_LC = brier_Mean - (1.96 * (sd(brier) / sqrt(n(
    )))),
    actualYcount_Mean = mean(actualYcount),
    untreatedYcount_Mean = mean(untreatedYcount),
    treatedYcount_Mean = mean(treatedYcount),
    treatedn_Mean = mean(treatedn),
    TP_Mean = mean(TP),
    FP_Mean = mean(FP),
    FN_Mean = mean(FN),
    TN_Mean = mean(TN),
    sensitivity_Mean = mean(sensitivity),
    specificity_Mean = mean(specificity),
    accuracy_Mean = mean(accuracy),
    PPV_Mean = mean(PPV),
    NPV_Mean = mean(NPV),
    b0_Mean = mean(b0_true),
    b1_Mean = mean(b1_true),
    b2_Mean = mean(b2_true),
    b0_estimate_Mean = mean(b0_estimate),
    b1_estimate_Mean = mean(b1_estimate),
    b2_estimate_Mean = mean(b2_estimate),
    .groups = "keep"
  )
