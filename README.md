---
title: "Local Projections Difference-in-Differences"
author: "Alexander Cardazzi" (Slight help from [Zach Porreca](https://zachporreca.github.io/))
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Local Projections Difference-in-Differences}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Getting Started

<Do we really want to capitalize, the "s" in "started". Feels offensive.>

Please see the working paper on which this package is based [here](https://www.nber.org/papers/w31184).

To install and load this package, use the following code:

```{r include=FALSE}
devtools::load_all(".")
par(mar = c(4.1, 4.1, .1, 4.1))
```

```{r eval=FALSE}
if(!"devtools" %in% installed.packages()) install.packages("devtools")
if(!"lpdid" %in% installed.packages()) devtools::install_github("alexCardazzi/lpdid")
library("lpdid")
```

## A Simple 2x2 DiD Example

To begin, consider a simple 2x2 DiD where treatment begins at $t = 31$.  Only twenty out of fifty units get treated, and the treatment dynamics cause the outcome variable to increase linearly with time.  Below, the outcome is averaged by time and whether the unit receives treatment at some point.

```{r, fig.width=7, fig.height=5, fig.align='center', cache=TRUE}
set.seed(757)
I <- 50; T <- 50; N <- I * T
data.frame(id = rep(1:I, each = T),
           t = rep(1:T, I),
           e = rnorm(N, 0, 4)) -> df
unit_FE <- rnorm(I, 0, 4)
time_FE <- rnorm(T, 0, 4)
df$treated <- ifelse(df$id < 21, 1, 0)
df$time_til <- ifelse(df$treated == 1, df$t - 30, -1)
df$y <- df$e + unit_FE[df$id] + time_FE[df$t] + ifelse(df$id < 21 & df$t > 30, df$t - 30, 0)
aggregate(df$y, list(df$t, df$treated), mean) -> tmp
plot(tmp$Group.1, tmp$x, col = tmp$Group.2 + 1, pch = 19,
     xlab = "Time", ylab = "Outcome")
legend("topleft", legend = c("Treated", "Control"),
       bty = "n", pch = 19, col = c("tomato", "black"))
```

In the most simple of cases, LP-DiD and TWFE preform very similarly.  Of course, TWFE is not biased in this setting, so we observe nearly identical parameter estimates (and standard errors).  The important arguments for ```lpdid()``` are the following:

- ```df```: The analysis dataset.
- ```window```: This is the pre- and post-period lengths in terms of time.  Importantly, the first value must be negative and the second must be positive.
- ```y```: The name of the outcome variable.
- ```unit_index```, ```time_index```: unit and time indices.  For example, if the data
is state-by-year, this is where you would input these names.
- ```rel_time```: This is the relative time variable.  **For untreated units, this must be negative.**

```{r, fig.width=7, fig.height=5, fig.align='center', cache=TRUE}
lpdid(df, window = c(-20, 20), y = "y",
      unit_index = "id", time_index = "t",
      rel_time = "time_til") -> reg
feols(y ~ i(time_til, treated, -1) | id + t, data = df) -> twfe

plot_lpdid(reg, col = "tomato", cex = 1.3)
iplot(twfe, add = TRUE)
legend("topleft", legend = c("TWFE", "LP-DiD"),
       col = c("black", "tomato"), pch = c(20, 19), bty = "n")
```

## Staggered (Exogenous) Treatment Timing

If you are reading this, it is likely that you are aware that TWFE estimates are potentially biased depending on treatment timing and dynamics.  In Dube et al. (2023), the authors simulate two sets of data: one where treatment timing is exogenous and the other where it is endogenous.  First, LP-DiD is applied to the data where treatment is exogenous.

```lpdid``` has a built-in function to replicate (as closely as possible) the simulations in Dube et al. (2023).^[This function is under development.  Currently, the parameters are fixed, so the only options are the exogeneity of the treatment timing and the random number seed.  The intention is to allow for full control over the parameters so users can experiment.]  Using this, we estimate TWFE and LP-DiD models in addition to plotting the "true" average beta values.

```{r, fig.width=7, fig.height=5, fig.align='center', cache=TRUE, message=FALSE, warning=FALSE}
df <- genr_sim_data(exogenous_timing = TRUE, 789)

twfe1 <- feols(y ~ i(rel_time, ever_treat, -1) | i + t, data = df)

lpdid1 <- lpdid(df, window = c(-5, 10), y = "y",
                unit_index = "i", time_index = "t",
                rel_time = "rel_time")

plot_lpdid(lpdid1, col = "tomato", x.shift = 0.1)
iplot(twfe1, add = TRUE, x.shift = -0.1)
points(aggregate(df$beta[df$ever_treat == 1],
                 list(df$rel_time[df$ever_treat == 1]), mean))
legend("topleft", c("True Beta", "TWFE", "LP-DiD"), bty = "n",
       pch = c(1, 19, 19), col = c("black", "black", "tomato"))
```

## Staggered (Endogenous) Treatment Timing

For endogenous treatment timing, the main idea is that treatment is more likely following a large negative shock.  See page 28 of the paper for more details on this.  However, given this, we can include a lag of the outcome variable into the model.  Below, four event studies are plotted: TWFE and LP-DiD, both with and without a lag of the outcome.

```{r, fig.width=7, fig.height=5, fig.align='center', cache=TRUE, message=FALSE, warning=FALSE}
df <- genr_sim_data(exogenous_timing = FALSE, 789)
# The returned "df" is a pdata.frame, so lag() works in a panel setting.
df$lag_y <- lag(df$y, 1)
twfe0 <- feols(y ~ i(rel_time, ever_treat, -1) | i + t, data = df)
twfe1 <- feols(y ~ lag_y + i(rel_time, ever_treat, -1) | i + t, data = df)
lpdid0 <- lpdid(df, window = c(-5, 10), y = "y", outcome_lags = 0,
                unit_index = "i", time_index = "t", rel_time = "rel_time")
lpdid1 <- lpdid(df, window = c(-5, 10), y = "y", outcome_lags = 1,
                unit_index = "i", time_index = "t", rel_time = "rel_time")
plot_lpdid(lpdid0, col = "tomato", x.shift = 0.1)
plot_lpdid(lpdid1, col = "dodgerblue", x.shift = 0.1, add = TRUE)
iplot(twfe0, add = TRUE, x.shift = -0.1)
iplot(twfe1, add = TRUE, x.shift = -0.1, col = "gray")
points(aggregate(df$beta[df$ever_treat == 1], list(df$rel_time[df$ever_treat == 1]), mean))
legend("topleft", bty = "n",
       c("True Beta", "TWFE", "TWFE + Y_lag", "LP-DiD", "LP-DiD + Y_lag"),
       pch = c(1, 19, 19, 19, 19),
       col = c("black", "black", "gray", "tomato", "dodgerblue"))
```

<!-- ## The ```castle``` Dataset -->

<!-- ```{r, fig.width=7, fig.height=5, fig.align='center', cache=TRUE} -->
<!-- library("fixest") -->
<!-- df <- bacondecomp::castle -->

<!-- df$t <- df$year - df$treatment_date -->
<!-- ifelse(is.na(df$t), -1, df$t) -> df$t -->
<!-- df$treat <- ave(df$post, df$state, FUN = max) -->
<!-- twfe1 <- feols(l_homicide ~ i(t, treat, -1) | state + year, data = df) -->

<!-- df$treatment_date1 <- ifelse(is.na(df$treatment_date), 3000, df$treatment_date) -->
<!-- sunab1 <- feols(l_homicide ~ sunab(treatment_date1, year) -->
<!--                | state + year, data = df) -->

<!-- lpdid(df = df, -->
<!--       window = c(-19, 10), -->
<!--       y = "l_homicide", -->
<!--       unit_index = "state", -->
<!--       time_index = "year", -->
<!--       rel_time = "t") -> reg -->

<!-- par(mfrow = c(1,1)) -->
<!-- plot_lpdid(reg, col = "tomato") -->
<!-- iplot(sunab1, add = TRUE, col = "dodgerblue", x.shift = 0.15) -->
<!-- iplot(twfe1, add = TRUE, col = "black", x.shift = -0.15) -->
<!-- legend("bottomright", c("TWFE", "LP-DiD", "Sun & Abraham"), -->
<!--        col = c("black", "tomato", "dodgerblue"), pch = 15, -->
<!--        bty = "n") -->
<!-- ``` -->

## Andrew Baker's Simulated Dataset

Andrew Baker created a simulated data (found [here](https://github.com/scunning1975/mixtape/)) in order to demonstrate the potential pit falls of using TWFE in a setting with staggered treatment timing and temporal treatment dynamics.  Below is an aggregated version of this data where average outcomes are calculated by time and treatment group.

```{r, fig.width=7, fig.height=5, fig.align='center', warning=FALSE, cache=TRUE}
baker <- haven::read_dta('https://github.com/scunning1975/mixtape/raw/master/baker.dta')
baker$treated <- ifelse(baker$treat_date > 0, 1, 0)

tmp <- aggregate(list(y = baker$y),
                 list(year = baker$year,
                      group = baker$group),
                 mean)
plot(tmp$year, tmp$y, pch = 19,
     col = scales::alpha(tmp$group, .4),
     xlab = "Year", ylab = "Outcome")
lines(tmp$year[tmp$group == unique(tmp$group)[1]],
      tmp$y[tmp$group == unique(tmp$group)[1]],
      col = tmp$group[tmp$group == unique(tmp$group)[1]])
lines(tmp$year[tmp$group == unique(tmp$group)[2]],
      tmp$y[tmp$group == unique(tmp$group)[2]],
      col = tmp$group[tmp$group == unique(tmp$group)[2]])
lines(tmp$year[tmp$group == unique(tmp$group)[3]],
      tmp$y[tmp$group == unique(tmp$group)[3]],
      col = tmp$group[tmp$group == unique(tmp$group)[3]])
lines(tmp$year[tmp$group == unique(tmp$group)[4]],
      tmp$y[tmp$group == unique(tmp$group)[4]],
      col = tmp$group[tmp$group == unique(tmp$group)[4]])
abline(v = unique(baker$treat_date)-0.5, col = 1:4)
```

When estimating a simple TWFE model, the following is the estimate generated by OLS:

```{r}
feols(y ~ treat | id + year, data = baker)
```

Of course, the treatment effect is surely non-negative.  Below, we estimate TWFE, Sun & Abraham, and LP-DiD models.  Here, just to exemplify another feature of the package, one can re-weight the regression to obtain an ATT rather than a VWATT.  Here, the results are nearly identical to the estimates via Sun & Abraham.  Also, note that the window (```(-20, 20)```) is wider than what is plotted ```(-17, 17)```.  This is a feature since the FEs and the treatment variable are perfectly colinear outside what is returned.

```{r, fig.width=7, fig.height=5, fig.align='center', warning=FALSE, message=FALSE, cache=TRUE}
res_naive = feols(y ~ i(time_til, treated, ref = -1) | 
                    id + year,
                  baker, cluster = ~state)

res_cohort = feols(y ~ sunab(treat_date, year) | id + year, 
                   baker[baker$year<2004,],
                   cluster = ~state)

lpdid(df = baker, window = c(-20, 20), y = "y",
      unit_index = "id", time_index = "year",
      rel_time = "time_til", reweight = TRUE) -> reg

iplot(res_naive, col = "black", grid = F,
      xlab = "Time to Treatment", main = "")
iplot(res_cohort, col = "dodgerblue", add = TRUE)
plot_lpdid(reg, add = TRUE, col = "tomato", pch = 1, cex = 1.5)
legend("topright", col = c("black", "tomato", "dodgerblue"),
       pch = c(20, 20, 1), bty = "n",
       legend = c("TWFE", "Sun & Abraham", "LP-DiD"))
```

It is also possible to pool together the treatment effects and get standard errors using the ```lpdid``` function.

```{r, fig.width=7, fig.height=5, fig.align='center', warning=FALSE, message = FALSE, cache=TRUE}
lpdid(df = baker, window = c(-20, 20), y = "y",
      unit_index = "id", time_index = "year",
      rel_time = "time_til", reweight = TRUE, pooled = TRUE) -> reg
reg$coeftable
```

## Non-Absorbing Treatment

A final feature of the package (that is still under construction), is handling non-absorbing treatment.  The treatment example used in Dube et al. (2023) is changes in the minimum wage.  It is impossible to find a "clean control" since all states have increased the minimum wage at one point or another.  So, the authors (see also [Cengiz et al., 2019](https://doi.org/10.1093/qje/qjz014)) make an assumption that $L$ time periods after treatment, the effect "stabilizes".  See the text for more details.

Below, we simulate data where there are three groups.  First, there is an untreated group.  Second, there is a group that gets treated only a single time.  Finally, there is a group that is treated twice.  The effect sizes are 15 and 35 for first and second treatments.  The aggregated time series plots can be seen below.

```{r, fig.width=7, fig.height=5, fig.align='center', warning=FALSE, message = FALSE, cache=TRUE}
set.seed(757)
I <- 50; T <- 50; N <- I * T; e_var <- 4
data.frame(id = rep(1:I, each = T),
           t = rep(1:T, I),
           e = rnorm(N, 0, e_var)) -> df
unit_FE <- rnorm(I, 0, e_var)
time_FE <- rnorm(T, 0, e_var)
df$treated <- ifelse(df$id < 21, 1, 0)
df$treat_status <- ifelse(df$id < 21 & df$t %in% 16:25, 1, 0)
df$treat_status <- ifelse(df$id < 11 & df$t %in% 36:45, 1, df$treat_status)
df$y <- df$e + unit_FE[df$id] + time_FE[df$t] +
  ifelse(df$id < 11,
         ifelse(df$t > 35, 50,
                ifelse(df$t > 15, 15, 0)),
         ifelse(df$id < 21,
                ifelse(df$t > 15, 15, 0),
                0))
df$group <- (ave(df$treat_status, df$id)) * 5

agg <- aggregate(df$y, list(df$t, df$group), mean)
plot(agg$Group.1, agg$x, col = agg$Group.2 + 1, pch = 19,
     xlab = "Time", ylab = "Outcome")
```

Notice the differences in this function call versus the previous ones.  First, there is no ```rel_time``` arguement needed.  This is because rel_time does not necessarily make sense for a unit that is treated more than once.  Second, reweighting does not currently work when ```nonabsorbing = TRUE```.  Currently, the reweighting simply gets turned off when ```nonabsorbing = TRUE```.^[As a note, pooling with controls (or outcome lags) is also not working.]  Finally, and perhaps most importantly, notice the ```nonabsorbing_treat_status``` variable.  One must insert the name of a variable that details the treatment status of the unit at every point.  Unfortunately, this portion of the package is still under construction.  To give users an idea of what the variable must look like, consider a unit that is treated in both 1995 and 2010.  Assume that treatment dynamics last only for 5 years. Then, observations of ```nonabsorbing_treat_status``` should be:
- 0 before 1995
- 1 between and including 1995 to 1999
- 0 between and including 2000 to 2009
- 1 between and including 2010 to 2014
- 0 after and including 2015

```{r, fig.width=7, fig.height=5, fig.align='center', warning=FALSE, message = FALSE, cache=TRUE}
lpdid(df = df, window = c(-12, 12), y = "y",
      unit_index = "id", time_index = "t",
      nonabsorbing = TRUE, nonabsorbing_lag = 5,
      nonabsorbing_treat_status = "treat_status") -> reg
plot_lpdid(reg)
abline(h = mean(c(15, 15, 35)), lty = 2)
```


## Under Construction

This is the first iteration of this package.  So, of course, there are some known (and surely many unknown) bugs.  Some of them are as follows:

- Thing 1
- Thing 2
