#' ---
#' title: "Mosquitoes - Analysis"
#' author: "Steffi LaZerte"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: cosmo
#'     toc: true
#'     toc_float:
#'       collapsed: false
#' ---
#'
#' <!--- This style block is to tweak the formatting of the data tables in the report --->
#' <style>
#'
#'   h2 {
#'   font-weight: bold;
#'   }
#'
#'   table {
#'     display: block;
#'     overflow: auto;
#'   }
#'
#'   blockquote {
#'     background-color: #88cabf;
#'  }
#'
#' </style>

#+ echo = FALSE
# RMarkdown -------------------------------------------------------------------

#' R Code to run and save this report:
#+ eval = FALSE
# rmarkdown::render(input = "Scripts/03_analysis.R",
#                   output_dir = "Results",
#                   output_file = paste0("03_analysis_", Sys.Date(), '.html'),
#                   envir = new.env())

#+ echo = FALSE
# Setup -------------------------------------------------------------------
#' # Setup

# renv::restore() # If using these scripts to restore dependencies (do not knit)

#+ setup, include = FALSE
knitr::opts_chunk$set(fig.width = 8, fig.asp = 0.5, out.width = "90%")

#+ echo = FALSE
## Load packages -------------------------------------------
#' ## Load packages
#+ message = FALSE
library(tidyverse)
library(lubridate)
library(here)
library(assertr)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(car)
library(performance)
library(gt)
library(tmaptools)

#+ echo = FALSE
## Load data -------------------------------------------
#' ## Load data
mosq <- read_csv(here("Data/Datasets/mosquitoes.csv"), show_col_types = FALSE) %>%
  filter(after_fogged == FALSE) %>%
  mutate(year = factor(year),
         cdcweek_year = paste0(year, "_", cdcweek)) %>%
  select(date, cdcweek_year, cdcweek, site, sitespecific, species, trapcount, ddm14,
         pt14, rha14, ttmin, year, lat, lon) %>%
  drop_na()

mosq %>%
  filter(site == "Brandon") %>%
  group_by(year) %>%
  summarize(min = min(cdcweek), max = max(cdcweek))


#' Create species specific data sets
mosq_av <- filter(mosq, species == "Ae. vexans")
mosq_cp <- filter(mosq, species == "Cq. perturbans")
mosq_ct <- filter(mosq, species == "Cx. tarsalis")
mosq_od <- filter(mosq, species == "Oc. dorsalis")

#' Outliers removed
mosq_av_o <- filter(mosq, species == "Ae. vexans", site != "Cypress River")
mosq_cp_o <- filter(mosq, species == "Cq. perturbans", site != "Cypress River")
mosq_ct_o <- filter(mosq, species == "Cx. tarsalis", site != "Cypress River")
mosq_od_o <- filter(mosq, species == "Oc. dorsalis", site != "Cypress River")


#+ echo = FALSE
## Figure themes -------------------------------------------
#' ## Figure themes
#'
#' Since we're creating a series of plots based on the same themes, it's easier
#' to bundle them here, and then add the bundle.
#'
#' **Data limits**
#'
#' Make sure that predicted data and figure legends reflect, consistently the
#' limits of the data set (calculated on the whole data set, not individually by
#' species).

limits <- mosq %>%
  select(ddm14, pt14, rha14) %>%
  summarize(across(everything(), list(min = ~floor(min(.)), max = ~ceiling(max(.)))),
            .groups = "drop")

limits

lim_ddm14 <- select(limits, contains("ddm")) %>% unlist()
lim_pt14 <- select(limits, contains("pt14")) %>% unlist()
lim_rha14 <- select(limits, contains("rha14")) %>% unlist()

range_ddm14 <- seq(lim_ddm14[1], lim_ddm14[2],  length.out = 5)
range_pt14 <- seq(lim_pt14[1], lim_pt14[2], length.out = 5)
range_rha14 <- seq(lim_rha14[1], lim_rha14[2],  length.out = 5)


theme_ddm14 <-  function(g) {
  g +
    theme_bw() +
    scale_colour_viridis_c(end = 0.8, option = "inferno",
                           labels = ~paste0(., "\u00B0C"), limits = lim_ddm14) +
    scale_size(breaks = c(1, 5, 10), limits = c(1, 15)) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 5, 10, 50, 100, 500, 1000, 5000)) +
    labs(x = "Week", y = "Trap Count", colour = "Degree Days\n(14-day mean)",
         size = "Trapping\nEvents") +
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2))
}

theme_pt14 <- function(g) {
  g +
    theme_bw() +
    scale_colour_viridis_c(option = "cividis", end = 0.90, direction = -1,
                           labels = ~paste0(., " mm"), limits = lim_pt14) +
    scale_size(breaks = c(1, 5, 10), limits = c(1, 15)) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 5, 10, 50, 100, 500, 1000, 5000)) +
    labs(x = "Week", y = "Trap Count", colour = "Precipitation\n(14-day mean)",
         size = "Trapping\nEvents")+
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2))
}

theme_rha14 <- function(g) {
  g +
    theme_bw() +
    scale_colour_viridis_c(end = 0.8, option = "plasma", direction = -1,
                           labels = ~paste0(., "%"), limits = lim_rha14) +
    scale_size(breaks = c(1, 5, 10, 15), limits = c(1, 20)) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 5, 10, 50, 100, 500, 1000, 5000)) +
    labs(x = "Week", y = "Trap Count",  colour = "Relative Humidity\n(14-day mean)",
         size = "Trapping\nEvents")+
    guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2))
}



#+ echo = FALSE
# Sample Sizes ---------------------------------------------------
#' # Sample Sizes
#'
#' After removing all missing variables among the parameters we include in the
#' models, the final samples are sizes are as follows (note the `n` represents
#' trapping events):
#'
#' ### By Species
mosq %>%
  group_by(species) %>%
  summarize(n = n(),
            n_no_cypress = length(trapcount[site != "Cypress River"]),
            .groups = "drop") %>%
  gt()

#' ### By Species by year
mosq %>%
  group_by(species, year) %>%
  summarize(n = n(),
            n_no_cypress = length(trapcount[site != "Cypress River"]),
            .groups = "drop") %>%
  gt()

#' ### By Species by site
#' This is probably not super relevant for publication, but good for your records
count(mosq, species, site) %>%
  gt()


#' ### By Species number of sites
mosq %>%
  group_by(species) %>%
  summarize(n_sites = n_distinct(site),
            n_sites_no_cypress = n_distinct(site[site != "Cypress River"]),
            .groups = "drop") %>%
  gt()

#' ### By sub site
select(mosq, date, site, sitespecific) %>%
  distinct() %>%
  count(site, sitespecific) %>%
  gt()

select(mosq, year, date, site, sitespecific) %>%
  distinct() %>%
  count(year, site, sitespecific) %>%
  gt()



#+ echo = FALSE
# Analysis Background ---------------------------------------------------
#' # Analysis Background
#'

#+ echo = FALSE
## Details -------------------------------------------------------------------
#' ## Details
#'
#' In the exploration step we discovered several things
#'
#' - Negative binomial is probably the distribution we want to use
#' - There is variation among sites
#' - There are polynomial relationships
#' - There is considerable variation among species (i.e. we'll want to do separate models)
#' - There are some outliers (`trapcounts` > 1000 or 5000; site Cypress River) to consider
#'
#' Therefore, we use negative binomial generalized linear mixed models with
#' site as a random factor. We will test for polynomial relationships among the
#' main parameters and will conduct separate analyses for each species.
#' We will also check for differences in our results due to outliers.
#' 
#' > **Note:** For simplicity, here we present only outlier-removed models. For
#' >  Biological reasons, Cypress River was omitted as a site.
#'
#' **Note:** We cannot use year as a random factor because there are too few levels
#' (i.e. years). But we can include it as a covariate.
#' 
#' > **Following Reviewer Feedback:**  
#' > A helpful and knowledgeable reviewer suggested several changes to the models
#' > to account for leftover temporal variation, sub-site variation and yearly
#' > interactions. They also suggested that we ensure there is no temporal or
#' > spatial autocorrelation that we should be dealing with.
#' 
#' Therefore we will add: 
#' 
#' - `cdcweek_year` as a crossed random factor to account for further temporal 
#'   variation associated with cdcweek by year. 
#' - `sitespecific` nested within `site` to account for potential variation 
#'   associated with sub-sites
#' - an assessment of the `year:cdcweek` interaction. If significant it will be
#'   kept, otherwise it will be omitted.
#' - Checks for temporal and spatial autocorrelation
#' 
#'
#+ echo = FALSE
## Parameters -------------------------------------------------------------------
#' ## Parameters
#'
#' ## Main parameters and polynomials
#'
#' In general we're looking at how time (week, `cdcweek`), temperature
#' (`ddm14`), relative humidity (`rha14`) and precipitation (`pt14`) affect trap
#' counts for each species.
#'
#' We will look for potential polynomial effects for each parameter using the
#' `poly()` function which evaluates *orthogonal polynomials*. Because
#' polynomials are usually added as x + x<sup>2</sup> + x<sup>3</sup> etc. it
#' means that the parameters are mathematically correlated which can cause
#' problems (one assumption of most linear models is that parameters are
#' independent). *Orthogonal polynomials* are created so that they are *not*
#' correlated (hooray!).
#'
#' We will look for potential 2-degree polynomials (x + x<sup>2</sup>) which
#' model curves (can be up or down, think mountains or valleys).
#'
#' We will also include `year` and `ttmin` as a potential covariates. I do not
#' include `ttmax` or `windspd` for several reasons:
#'
#' **`ttmax`**
#'
#' `ttmax` correlates well with `ttmin`, so the conditions we're trying to
#' account for are probably taken care of by including `ttmin`. Omitting `ttmax`
#' also reduces the risk of it explaining variation we want to capture with 
#' `ddm14` or `cdcweek`.
#'
#' **`windspd`**
#'
#' `windspd` is just weird. In the exploration, it seems to mostly reduce trap
#' counts when it is very low or very high (but sometimes the reverse!).
#'
#' It doesn't make sense for very low wind to reduce trap counts, nor for high
#' wind to increase trap counts, so I'm don't think that wind speed, as we 
#' measured it, is actually measuring what we expected.
#'
#' It does correlated with other variables so might be reflecting more of the
#' climatic conditions than wind itself (which is what we wanted to capture).
#'
#' Finally, including `windspd` in the models seemed to make the fit poorer, 
#' even if it did explain some variation.
#'
#+ echo = FALSE
## Building Models ------------------------------------
#' ## Building models
#'
#' How to decide which parameters to put in a model can be tricky. There are many
#' different ways, some are better some are just different. The main thing to do
#' is to pick a method and be consistent.
#'
#' Here, I build several models for each species and compare. The first
#' model is the **basic model** with main effects and covariates but
#' without polynomial terms. The second is the 
#' **full polynomial and interaction model** with main effects, their polynomial
#' terms, interactions and covariates. The third is a reduced interaction 
#' model, the fourth a reduced polynomial model.
#'
#' Here, because there is enough data to have reasonably good power, I keep all
#' main parameters and covariates in the model, even if they're not significant.
#' However, I do omit second order polynomial terms (i.e. the squared component)
#' and interactions if they are not significant (using `summary()` look for the
#' `2` for polynomials, use alpha =
#' 0.05; using `Anova()` look for the `year:cdcweek` for interactions;). 
#' This gives us the **reduced model**, or, in extreme cases, means we would
#' end up back with our basic model.
#'
#' I then confirm that this reduced model is "better" than the other models
#' using the `anova()` function to perform
#' an analysis of deviance. The analysis of deviance tests whether models are
#' significantly different from one another. If there is a significant
#' difference, the model with the lowest likelihood is the 'better' model. If
#' there is no significant difference, you take the simplest, most parsimonious
#' model, because adding terms doesn't make the model significantly better.
#' Sometimes this means we don't take the reduced model but a simpler model even.
#'
#' I always expect the polynomial models to be better than the basic,
#' non-polynomial model, but we check to be sure (it may not alway be true). I also
#' expect the reduced model to be better than the basic, non-polynomial models,
#' but again, we check to be sure.
#'
#+ echo = FALSE
## Visualizing results -------------------------------------------------------------------
#' ## Visualizing the results
#'
#' Polynomial coefficients are difficult if not impossible to interpret, so the
#' best way is to plot the predicted model lines against the data.
#'
#+ echo = FALSE
## How to interpret? ------------------------------------------------------------
#' ## How to interpret?
#'
#' Here I use several tools in these analyses
#'
#' **1. [`glmmTMB`](https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf)**
#'
#'   - A package for quite advanced models, including modelling zero-inflation,
#'   mixed factors (i.e. fixed and random) etc.
#'   - Here I use it because it works well for mixed GLMs
#'   - In this example, the model objects are `m0` and `m1`
#'   - The response variable is `trapcount` the explanatory variables are
#'     polynomials of `ddm14` and `rha14` and a linear covariate, `year`.
#'   - We're including `site` and `sitespecific` as random intercepts
#'   - We're using a classic negative binomial family, `nbinom2`

m0 <- glmmTMB(trapcount ~ poly(ddm14, degree = 2) + year + (1|site),
              family = "nbinom2", data = mosq)

m1 <- glmmTMB(trapcount ~ poly(ddm14, degree = 2) + poly(rha14, degree = 2) +
                year + (1|site),
              family = "nbinom2", data = mosq)

#'
#' **2. Simulated residuals from the [DHARMa package](https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html)**
#'
#'   - These are similar to, but IMPORTANTLY, not the same as QQNorm plots and 
#'   plots of heteroscedacity.
#'   - Simply put, these simulated residuals reflect whether or not data matches 
#'   the kind of model you're running
#'   - This makes them excellent for assessing model fit for ANY type of model.
#'   - In this example we can see issues with the residual variance distribution, 
#'   distribution and dispersion (all the red!)

r <- simulateResiduals(m0, plot = TRUE)
r <- simulateResiduals(m1, plot = TRUE)

#' **3. Interpreting coefficients**
#'
#'   - The "Random effects" section contains variance explained by random factors
#'     (here there is variation captured by site)
#'   - The "Conditional model" section contains fixed effects (i.e what you would
#'     expect from a regular linear model).
#'   - Estimates are currently in log units
summary(m1)

#'   - But we an transform them to incident rate ratios which give us a
#'     multiplicative factor
#'   - Here, if we were to interpret the effect of `year`, we could say that
#'     2021 had, on average, 2.45x greater trap counts than 2020

exp(coefficients(summary(m1))$cond["year2021", "Estimate"])


#' **4. `Anova()` vs. `summary()`**
#'
#' - The `Anova()` function gives you an Anova table which is great for interpreting
#'   values that are less straightforward (i.e. polynomials and categorical variables).
#' - The `summary()` function gives you a summary table which is great for
#' interpreting linear relationships, because you can give meaning to the
#' Estimate.
#' - Therefore, when interpreting the polynomials, I suggest reporting
#' P and Chisq from the Anova table, when interpreting linear effects,
#' I suggest reporting Est, P and z from the summary table.

#' **5. Temporal and Spatial autocorrelation**
#' 
#' - We use the `testSpatialAutocorrelation()` and `testTemporalAutocorrelation()`
#'   functions from the DHARMa package.
#' - We test for overal autocorrelation and for autocorrelation by groups (i.e.
#'   within sites or within years).
#' - Unless there is substantial and consistent autocorrelation we'll leave the
#'   models as is
#' - Where there is autocorrelation, we use the `ar1()` function to model it and
#'   check with `anova()` that it improves the model (cannot use the tests
#'   to check as they won't detect it, see their docs).
#' 

#+ echo = FALSE
# Analysis  ---------------------------------------------------
#' # Analysis
#'
#+ echo = FALSE
# Species  ---------------------------------------------------
#' ## Species {.tabset}

#+ echo = FALSE
## Model --------------------------------------------
#' ### Model
#'
#' Arrange species by numbers to make interpretation easier (no affect on analysis)
mosq_sp <- mosq %>%
  filter(site != "Cypress River") %>%
  mutate(species = factor(species, levels = c(
    "Ae. vexans", "Oc. dorsalis", "Cx. tarsalis", "Cq. perturbans")),
    cdcweek_fct = factor(cdcweek),
    cdcweek_year_fct = factor(cdcweek_year),
    site_year = paste0(site, "_", year))
  
m0 <- glmmTMB(trapcount ~ species * year + ttmin + 
                (1|site/sitespecific) +
                ar1(cdcweek_fct + 0|site_year),
              family = "nbinom2", data = mosq_sp)
r_orig <- simulateResiduals(m0, plot = TRUE)

# Add autocorrelation structure
m2 <- glmmTMB(trapcount ~ species * year + ttmin + 
                (1|site/sitespecific) +
                ar1(cdcweek_fct + 0|site_year),
              family = "nbinom2", data = mosq_sp)
r <- simulateResiduals(m2, plot = TRUE)

Anova(m2, type = "III")
summary(m2)

#+ echo = FALSE
## Post-Hoc --------------------------------------------
#' ### Post-Hoc
#'
#' - Using the False Discovery Rate for P-value adjustment for multiple tests

e <- emmeans(m0, specs = "species", by = "year")
pairs(e, adjust = "fdr", type = "response")

#+ echo = FALSE
## Check Autocorrelation -------------------------------
#' ### Check Autocorrelation
#' 
#' > **NOTE:** This can test for existing autocorrelation *before* account, but
#' > cannot test whether the fix worked (see `?testTemporalAutocorrelation`)
#' 
#' **There is autocorrelation, so we add `ar1()` to the model.**
#' 
#' **Temporal autocorrelation**
data <- mosq_sp

r2 <- recalculateResiduals(
  r_orig, group = data$cdcweek_fct)

t <- testTemporalAutocorrelation(
  r2, time = unique(data$cdcweek))
title("Overall", line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r_orig, group = data$cdcweek,
    sel = data$year == y)
  
  t <- testTemporalAutocorrelation(
    r2, time = unique(data$cdcweek[data$year == y]))
  title(paste0(y, " overall "), line = 3)
  message(y)
  print(t)
  
  for(i in unique(data$site)) {
    if(y %in% data$year[data$site == i]) {
      r2 <- recalculateResiduals(
        r_orig, group = data$cdcweek_fct,
        sel = data$year == y & data$site == i)
      
      message(y, " - ", i)
      t <- testTemporalAutocorrelation(
        r2, time = unique(data$cdcweek_fct[data$year == y & data$site == i]))
      title(paste0(y, " - ", i), line = 3)
      print(t)
    }
  }
}


#' **Spatial Autocorrelation**

r2 <- recalculateResiduals(
  r, group = data$site)

t <- testSpatialAutocorrelation(
  r2,
  x = unique(data$lon),
  y = unique(data$lat))
title(paste0(round(t$p.value, 3)), line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$site,
    sel = data$year == y)
  
  t <- testSpatialAutocorrelation(
    r2,
    x = unique(data$lon[data$year == y]),
    y = unique(data$lat[data$year == y]))
  title(paste0(y, " - ", round(t$p.value, 3)), line = 3)
  print(t)
}


#+ echo = FALSE
# Ae. vexans ------------------------------------------------------------------
#' ## Ae. vexans {.tabset}

#+ echo = FALSE
## Models -----------------------------------
#' ### Models

#' **Basic model**
m0 <- glmmTMB(trapcount ~ year + ttmin +
                cdcweek + ddm14 + pt14 + rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
r <- simulateResiduals(m0, plot = TRUE)

#' **Full model**
m1 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                poly(pt14, degree = 2) * year + 
                poly(rha14, degree = 2) * year +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
r <- simulateResiduals(m1, plot = TRUE)
Anova(m1, type = "III")

#' **Interactions reduced model**
m2 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) + 
                poly(rha14, degree = 2) +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
r <- simulateResiduals(m2, plot = TRUE)
Anova(m2, type = "III")
summary(m2)

#' **Polynomials reduced model**
m3 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) + 
                rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
r <- simulateResiduals(m3, plot = TRUE)
Anova(m3, type = "III")
summary(m3)

#' check interactions back 
m4 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) + 
                rha14 * year +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
r <- simulateResiduals(m4, plot = TRUE)
Anova(m4, type = "III")

anova(m0, m1, m2, m3, m4)




#+ echo = FALSE
## Results -----------------------------------
#' ### Results

check_collinearity(m3)
m3a <- glmmTMB(trapcount ~ 
                 ttmin + year +
                 poly(cdcweek, degree = 2) * year +
                 poly(ddm14, degree = 2) +
                 poly(pt14, degree = 2) + 
                 rha14 +
                 (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_av_o)
check_collinearity(m3a)
r <- simulateResiduals(m3, plot = TRUE)
Anova(m3, type = "III")
summary(m3)

#' Get effect of ttmin as incident rate ratios
exp(coefficients(summary(m3))$cond["ttmin","Estimate"])
exp(coefficients(summary(m3))$cond["rha14","Estimate"])

#+ echo = FALSE
## Check Autocorrelation -------------------------------
#' ### Check Autocorrelation
#' 
#' **Temporal autocorrelation**
data <- mosq_av_o
for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$cdcweek,
    sel = data$year == y)
  
  t <- testTemporalAutocorrelation(
    r2, time = unique(data$cdcweek[data$year == y]))
  title(paste0(y, " overall "), line = 3)
  message(y)
  print(t)
  
  for(i in unique(data$site)) {
    if(y %in% data$year[data$site == i]) {
      r2 <- recalculateResiduals(
        r, group = data$cdcweek,
        sel = data$year == y & data$site == i)
      
      message(y, " - ", i)
      t <- testTemporalAutocorrelation(
        r2, time = unique(data$cdcweek[data$year == y & data$site == i]))
      title(paste0(y, " - ", i), line = 3)
      print(t)
    }
  }
}


#' **Spatial Autocorrelation**

r2 <- recalculateResiduals(
  r, group = data$site)

t <- testSpatialAutocorrelation(
  r2,
  x = unique(data$lon),
  y = unique(data$lat))
title(paste0(round(t$p.value, 3)), line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$site,
    sel = data$year == y)
  
  t <- testSpatialAutocorrelation(
    r2,
    x = unique(data$lon[data$year == y]),
    y = unique(data$lat[data$year == y]))
  title(paste0(y, " - ", round(t$p.value, 3)), line = 3)
  print(t)
}


#+ echo = FALSE
## Figures code -----------------------------------
#' ### Figures code
#'
#' Create figures with models plotted on them for various levels of the parameters

data <- mosq_av_o
model <- m3

new_av <- with(data,
               expand_grid(site = unique(site),
                           cdcweek = unique(cdcweek),
                           ddm14 = range_ddm14,
                           pt14 = range_pt14,
                           rha14 = range_rha14,
                           ttmin = mean(ttmin),
                           year = c(2020, 2021))) %>%
  left_join(select(data, site, sitespecific) %>% distinct(), by = "site") %>%
  mutate(cdcweek_year = paste0(year, "_", cdcweek))


#' **First Figure - ddm14**
t <- new_av %>%
  mutate(pt14 = mean(pt14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, ddm14) %>%
  summarize(y = mean(y), .groups = "drop")

g1 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = ddm14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = ddm14), linewidth = 2) +
  facet_wrap(~ year) 

g1 <- theme_ddm14(g1)

#' **Second Figure - pt14**
t <- new_av %>%
  mutate(ddm14 = mean(ddm14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, pt14) %>%
  summarize(y = mean(y), .groups = "drop")

g2 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = pt14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = pt14), size = 2) +
  facet_wrap(~ year) 

g2 <- theme_pt14(g2)

#' **Third Figure - rha14**
t <- new_av %>%
  mutate(ddm14 = mean(ddm14), pt14 = mean(pt14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, rha14) %>%
  summarize(y = mean(y), .groups = "drop")

g3 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = rha14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = rha14), linewidth = 2) +
  facet_wrap(~ year) 

g3 <- theme_rha14(g3)

#+ echo = FALSE
## Pub Figures -----------------------------------
#' ### Publication Figure
#+ fig.asp = 0.5, width = 4
g1

g2

g3


#+ echo = FALSE
# Cq. perturbans -----------------------------------------------
#' ## Cq. perturbans {.tabset}
#'
#+ echo = FALSE
## Models -----------------------------------
#' ### Models
#' **Basic model**
m0 <- glmmTMB(trapcount ~ ttmin + year +
                cdcweek + ddm14 + pt14 + rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_cp_o)
r <- simulateResiduals(m0, plot = TRUE)

#' **Full model**
m1 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                poly(pt14, degree = 2) * year +  
                poly(rha14, degree = 2) * year +
                + (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_cp_o)
r <- simulateResiduals(m1, plot = TRUE)
Anova(m1, "III")

#' **Interactions reduced model**
m2 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) +  
                poly(rha14, degree = 2) +
                ttmin + (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_cp_o)
r <- simulateResiduals(m2, plot = TRUE)
Anova(m2, "III")
summary(m2)


#' **Polynomials reduced model**
m3 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                ddm14 +
                pt14 +  
                poly(rha14, degree = 2) +
                ttmin + (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_cp_o)
r <- simulateResiduals(m3, plot = TRUE)
summary(m3)
Anova(m3, "III")

testZeroInflation(m3)



#' check interactions again
m4 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                ddm14 * year+
                pt14 * year+  
                poly(rha14, degree = 2) +
                ttmin + (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_cp_o)
r <- simulateResiduals(m4, plot = TRUE)
Anova(m4, "III")

anova(m0, m1, m3, m4)

#+ echo = FALSE
## Results -----------------------------------
#' ### Results

check_collinearity(m3)
r <- simulateResiduals(m3, plot = TRUE)
Anova(m3, type = "III")
summary(m3)

#' Get effect of ttmin as incident rate ratios
exp(coefficients(summary(m2))$cond["ttmin","Estimate"])

#+ echo = FALSE
## Check Autocorrelation -------------------------------
#' ### Check Autocorrelation

data <- mosq_cp_o
 
#' **Temporal autocorrelation**

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$cdcweek,
    sel = data$year == y)
  
  t <- testTemporalAutocorrelation(
    r2, time = unique(data$cdcweek[data$year == y]))
  title(paste0(y, " overall "), line = 3)
  message(y)
  print(t)
  
  for(i in unique(data$site)) {
    if(y %in% data$year[data$site == i]) {
      r2 <- recalculateResiduals(
        r, group = data$cdcweek,
        sel = data$year == y & data$site == i)
      
      message(y, " - ", i)
      t <- testTemporalAutocorrelation(
        r2, time = unique(data$cdcweek[data$year == y & data$site == i]))
      title(paste0(y, " - ", i), line = 3)
      print(t)
    }
  }
}


#' **Spatial Autocorrelation**

r2 <- recalculateResiduals(
  r, group = data$site)

t <- testSpatialAutocorrelation(
  r2,
  x = unique(data$lon),
  y = unique(data$lat))
title(paste0(round(t$p.value, 3)), line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$site,
    sel = data$year == y)
  
  t <- testSpatialAutocorrelation(
    r2,
    x = unique(data$lon[data$year == y]),
    y = unique(data$lat[data$year == y]))
  title(paste0(y, " - ", round(t$p.value, 3)), line = 3)
  print(t)
}


#+ echo = FALSE
## Figures code ---------------------------------
#' ### Figures code
data <- mosq_cp_o
model <- m3

new_cp <- with(data,
               expand_grid(site = unique(site),
                           cdcweek = unique(cdcweek),
                           ddm14 = range_ddm14,
                           pt14 = range_pt14,
                           rha14 = range_rha14,
                           ttmin = mean(ttmin),
                           year = c(2020, 2021))) %>%
  left_join(select(data, site, sitespecific) %>% distinct(), by = "site") %>%
  mutate(cdcweek_year = paste0(year, "_", cdcweek))

#' **First Figure - ddm14**
#' 
#' NON-Significant, therefore only one line showing the seasonal effect
t <- new_cp %>%
  mutate(ddm14 = mean(ddm14), pt14 = mean(pt14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, ddm14) %>%
  summarize(y = mean(y), .groups = "drop")

g1 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = ddm14)) +
  geom_count() +
  geom_line(data = t, aes(y = y), size = 2, colour = "black") +
  facet_wrap(~ year)

g1 <- theme_ddm14(g1)

#' **Second Figure - pt14**
#' 
#' NON-Significant, therefore only one line showing the seasonal effect
t <- new_cp %>%
  mutate(ddm14 = mean(ddm14), pt14 = mean(pt14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, pt14) %>%
  summarize(y = mean(y), .groups = "drop")

g2 <- ggplot(data = mosq_cp_o, aes(x = cdcweek, y = trapcount, colour = pt14)) +
  geom_count() +
  geom_line(data = t, aes(y = y), size = 2, colour = "black") +
  facet_wrap(~ year)

g2 <- theme_pt14(g2)


#' **Third Figure - rha14**
t <- new_cp %>%
  mutate(ddm14 = mean(ddm14), pt14 = mean(pt14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, rha14) %>%
  summarize(y = mean(y), .groups = "drop")

g3 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = rha14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = rha14), size = 2) +
  facet_wrap(~ year)

g3 <- theme_rha14(g3)

#+ echo = FALSE
## Pub Figures -----------------------------------
#' ### Publication Figure
#+ fig.asp = 0.5, width = 4

g1

g2

g3

#+ echo = FALSE
# Cx. tarsalis -----------------------------------------------
#' ## Cx. tarsalis {.tabset}
#'
#+ echo = FALSE
## Models -----------------------------------
#' ### Models
#' **Basic model**
m0 <- glmmTMB(trapcount ~ ttmin + year + 
                cdcweek + ddm14 + pt14 + rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
r <- simulateResiduals(m0, plot = TRUE)

#' **Full model**
m1 <- glmmTMB(trapcount ~ year + ttmin + 
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                poly(pt14, degree = 2) * year + 
                poly(rha14, degree = 2) * year +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
r <- simulateResiduals(m1, plot = TRUE)
Anova(m1, "III")

#' **Reduced interactions model**
m2 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) + 
                poly(rha14, degree = 2) +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
r <- simulateResiduals(m2, plot = TRUE)
Anova(m2, "III")
summary(m2)


#' **Reduced polynomials model**
m3 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                ddm14 +
                pt14 + 
                rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
r <- simulateResiduals(m3, plot = TRUE)
summary(m3)
Anova(m3, "III")

#' Check interactions
m4 <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                ddm14 +
                pt14 * year + 
                rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
r <- simulateResiduals(m4, plot = TRUE)
Anova(m4, "III")

anova(m0, m1, m2, m3, m4)
anova(m3, m4)


#+ echo = FALSE
## Results -----------------------------------
#' ### Results

check_collinearity(m4)
m4a <- glmmTMB(trapcount ~ ttmin + year +
                poly(cdcweek, degree = 2) + year + 
                ddm14 +
                pt14 + year + 
                rha14 +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_ct_o)
check_collinearity(m4a) # Good!

r <- simulateResiduals(m4, plot = TRUE)
Anova(m4, type = "III")
summary(m4)

#' Get effect of year as incident rate ratios
exp(coefficients(summary(m4))$cond["year2021:pt14", "Estimate"])
exp(coefficients(summary(m4))$cond["ddm14", "Estimate"])

#' Check if slope sig in either year...
e <- emmeans(m4, specs = "pt14", by = "year")
regrid(e, type = "response")


#+ echo = FALSE
## Check Autocorrelation -------------------------------
#' ### Check Autocorrelation

data <- mosq_ct_o

#' **Temporal autocorrelation**

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$cdcweek,
    sel = data$year == y)
  
  t <- testTemporalAutocorrelation(
    r2, time = unique(data$cdcweek[data$year == y]))
  title(paste0(y, " overall "), line = 3)
  message(y)
  print(t)
  
  for(i in unique(data$site)) {
    if(y %in% data$year[data$site == i]) {
      r2 <- recalculateResiduals(
        r, group = data$cdcweek,
        sel = data$year == y & data$site == i)
      
      message(y, " - ", i)
      t <- testTemporalAutocorrelation(
        r2, time = unique(data$cdcweek[data$year == y & data$site == i]))
      title(paste0(y, " - ", i), line = 3)
      print(t)
    }
  }
}


#' **Spatial Autocorrelation**

r2 <- recalculateResiduals(
  r, group = data$site)

t <- testSpatialAutocorrelation(
  r2,
  x = unique(data$lon),
  y = unique(data$lat))
title(paste0(round(t$p.value, 3)), line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$site,
    sel = data$year == y)
  
  t <- testSpatialAutocorrelation(
    r2,
    x = unique(data$lon[data$year == y]),
    y = unique(data$lat[data$year == y]))
  title(paste0(y, " - ", round(t$p.value, 3)), line = 3)
  print(t)
}




#+ echo = FALSE
## Figures code -----------------------------------
#' ### Figures code
#'
#' Create figures with models plotted on them for various levels of the parameters
data <- mosq_ct_o
model <- m4

new_ct <- with(data,
               expand_grid(site = unique(site),
                           cdcweek = unique(cdcweek),
                           ddm14 = range_ddm14,
                           pt14 = range_pt14,
                           rha14 = range_rha14,
                           ttmin = mean(ttmin),
                           year = c(2020, 2021))) %>%
  left_join(select(data, site, sitespecific) %>% distinct(), by = "site") %>%
  mutate(cdcweek_year = paste0(year, "_", cdcweek))

#' **First Figure - ddm14**
t <- new_ct %>%
  mutate(pt14 = mean(pt14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, ddm14) %>%
  summarize(y = mean(y), .groups = "drop")

g1 <- ggplot(data = data, aes(x = cdcweek, 
                                   y = trapcount, colour = ddm14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = ddm14), linewidth = 2) +
  facet_wrap(~year)

g1 <- theme_ddm14(g1)


#' **Second Figure - pt14**
#' 
#' Significant interaction
t <- new_ct %>%
  mutate(ddm14 = mean(ddm14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, pt14) %>%
  summarize(y = mean(y), .groups = "drop")

g2 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = pt14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = pt14), linewidth = 2) +
  facet_wrap(~year)

g2 <- theme_pt14(g2)


#' **Third Figure - rha14**
#' 
#' NON-Significant, therefore only one line showing the seasonal effect
t <- new_ct %>%
  mutate(ddm14 = mean(ddm14), rha14 = mean(rha14), pt14 = mean(pt14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, rha14) %>%
  summarize(y = mean(y), .groups = "drop")

g3 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = rha14)) +
  geom_count() +
  geom_line(data = t, aes(y = y), linewidth = 2, colour = "black") +
  facet_wrap(~year)
g3 <- theme_rha14(g3)


#+ echo = FALSE
## Pub Figures -----------------------------------
#' ### Publication Figure
#+ fig.asp = 0.5, width = 4

g1

g2

g3


#+ echo = FALSE
# Oc. dorsalis -----------------------------------------------
#' ## Oc. dorsalis {.tabset}
#'
#+ echo = FALSE
## Models -----------------------------------
#' ### Models
#' **Basic model**
m0 <- glmmTMB(trapcount ~ ttmin + year +
                cdcweek + ddm14 + pt14 + rha14 + 
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m0, plot = TRUE)

#' **Full model**
m1 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                poly(pt14, degree = 2) * year + 
                poly(rha14, degree = 2) * year +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m1, plot = TRUE)
Anova(m1, "III")

#' **Reduced interactions models**
m2 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                poly(pt14, degree = 2) + 
                poly(rha14, degree = 2) +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m2, plot = TRUE)
summary(m2)
Anova(m2, "III")

m2a <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) +
                poly(pt14, degree = 2) + 
                poly(rha14, degree = 2) +
                (1|cdcweek_year) + (1|site/sitespecific),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m2a, plot = TRUE)
summary(m2a)
Anova(m2a, "III")

#' **Reduced polynomial models**
m3 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                pt14 + 
                rha14 +
                (1|cdcweek_year) + (1|site),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m3, plot = TRUE)
summary(m3)
Anova(m3)

m3a <- glmmTMB(trapcount ~ 
                 ttmin + year +
                 poly(cdcweek, degree = 2) * year + 
                 ddm14 +
                 pt14 + 
                 rha14 +
                 (1|cdcweek_year) + (1|site/sitespecific),
               family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m3a, plot = TRUE)
summary(m3a)
Anova(m3a, "III")

#' Check interaction back in
m4 <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                poly(ddm14, degree = 2) * year +
                pt14 + 
                rha14 * year +
                (1|cdcweek_year) + (1|site),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m4, plot = TRUE)
Anova(m4)

m4a <- glmmTMB(trapcount ~ 
                ttmin + year +
                poly(cdcweek, degree = 2) * year + 
                ddm14 * year +
                pt14 + 
                rha14 * year +
                (1|cdcweek_year) + (1|site),
              family = "nbinom2", data = mosq_od_o)
r <- simulateResiduals(m4a, plot = TRUE)
Anova(m4a)

anova(m0, m1, m2, m2a, m3, m3a, m4, m4a)

#'
#+ echo = FALSE
## Results -----------------------------------
#' ### Results

check_collinearity(m3a)

r <- simulateResiduals(m3a, plot = TRUE)
Anova(m3a, type = "III")
summary(m3a)

#' Get effect of ttmin as incident rate ratios
exp(coefficients(summary(m3a))$cond["pt14","Estimate"])

#+ echo = FALSE
## Check Autocorrelation -------------------------------
#' ### Check Autocorrelation

data <- mosq_od_o

#' **Temporal autocorrelation**

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$cdcweek,
    sel = data$year == y)
  
  t <- testTemporalAutocorrelation(
    r2, time = unique(data$cdcweek[data$year == y]))
  title(paste0(y, " overall "), line = 3)
  message(y)
  print(t)
  
  for(i in unique(data$site)) {
    if(y %in% data$year[data$site == i]) {
      r2 <- recalculateResiduals(
        r, group = data$cdcweek,
        sel = data$year == y & data$site == i)
      
      message(y, " - ", i)
      t <- testTemporalAutocorrelation(
        r2, time = unique(data$cdcweek[data$year == y & data$site == i]))
      title(paste0(y, " - ", i), line = 3)
      print(t)
    }
  }
}


#' **Spatial Autocorrelation**

r2 <- recalculateResiduals(
  r, group = data$site)

t <- testSpatialAutocorrelation(
  r2,
  x = unique(data$lon),
  y = unique(data$lat))
title(paste0(round(t$p.value, 3)), line = 3)
print(t)

for(y in 2020:2021) {
  r2 <- recalculateResiduals(
    r, group = data$site,
    sel = data$year == y)
  
  t <- testSpatialAutocorrelation(
    r2,
    x = unique(data$lon[data$year == y]),
    y = unique(data$lat[data$year == y]))
  title(paste0(y, " - ", round(t$p.value, 3)), line = 3)
  print(t)
}



#+ echo = FALSE
## Figures code ---------------------------------
#' ### Figures code
#'
data <- mosq_od_o
model <- m3a

new_od <- with(data,
               expand_grid(site = unique(site),
                           cdcweek = unique(cdcweek),
                           ddm14 = range_ddm14,
                           pt14 = range_pt14,
                           rha14 = range_rha14,
                           ttmin = mean(ttmin),
                           year = c(2020, 2021))) %>%
  left_join(select(data, site, sitespecific) %>% distinct(), by = "site") %>%
  mutate(cdcweek_year = paste0(year, "_", cdcweek))

#' **First Figure - ddm14**
#' 
#' NON-Significant, therefore only one line showing the seasonal effect
t <- new_od %>%
  mutate(ddm14 = mean(ddm14), pt14 = mean(pt14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, ddm14) %>%
  summarize(y = mean(y), .groups = "drop")

g1 <- ggplot(data = data, aes(x = cdcweek, y = trapcount, colour = ddm14)) +
  geom_count() +
  geom_line(data = t, aes(y = y), linewidth = 2, colour = "black") +
  facet_wrap(~year)
g1 <- theme_ddm14(g1)


#' **Second Figure - pt14**
#' 
t <- new_od %>%
  mutate(ddm14 = mean(ddm14), rha14 = mean(rha14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, pt14) %>%
  summarize(y = mean(y), .groups = "drop")

g2 <- ggplot(data = mosq_od_o, aes(x = cdcweek, y = trapcount, colour = pt14)) +
  geom_count() +
  geom_line(data = t, aes(y = y, group = pt14), size = 2) +
  facet_wrap(~year) 
g2 <- theme_pt14(g2)

#' **Third Figure - rha14**
#' NON-Significant, therefore only one line showing the seasonal effect
t <- new_od %>%
  mutate(ddm14 = mean(ddm14), rha14 = mean(rha14), pt14 = mean(pt14)) %>%
  mutate(y = predict(model, newdata = ., type = "response", re.form = NA)) %>%
  select(-site, -sitespecific, -cdcweek_year) %>%
  distinct() %>%
  group_by(year, cdcweek, rha14) %>%
  summarize(y = mean(y), .groups = "drop")

g3 <- ggplot(data = mosq_od_o, aes(x = cdcweek, y = trapcount, colour = rha14)) +
  geom_count() +
  geom_line(data = t, aes(y = y), linewidth = 2, colour = "black") +
  facet_wrap(~year)
g3 <- theme_rha14(g3)

#+ echo = FALSE
## Pub Figures -----------------------------------
#' ### Publication Figure
#+ fig.asp = 0.5, width = 4
g1

g2

g3

#+ echo = FALSE
# Methods ---------------------------------------------------------------
#' # Methods
#'
#' Analyses were conducted with R statistical software
#' (v`r paste(R.version$major, R.version$minor, sep = ".")`; R Core Team 2022).
#'
#'
#' We explored two different sets of models: 1) differences among mosquito
#' species trap counts; and 2) effects of average degree days, relative
#' humidity, precipitation and time (week) on mosquito trap counts.
#' All models were Generalized Linear Mixed Models
#' (GLMMs; glmmTMB package v`r packageVersion("glmmTMB")`; Brooks et al. 2017),
#' using a negative binomial distribution and random intercepts for site and subsite to
#' control for site-level effects.
#' We used a negative binomial distribution because we have count data,
#' but they did not match a Poisson distribution (overdispersed).
#'
#' Differences in trap counts among species were modelled using the complete
#' data set with trap counts as the response and species, year, and their
#' interaction explanatory variables. We also included minimum temperature of
#' the trap-day as a covariate. We conducted Post-Hoc analyses to compare species
#' difference among years (emmeans package v`r packageVersion("emmeans")`;
#' Lenth 2021) with the false discovery rate P-value adjustment. We included an 
#' AR1 model of temporal autocorrelation.
#'
#' Preliminary exploration suggested that trap count differences with local
#' weather variation were species-specific. Therefore we modelled weather variables
#' separately for each species. To account for curvilinear relationships,
#' variables were modelled as 2nd degree orthologonal polynomials. We also included
#' year and minimum temperature of the trap-day as linear covariates. We investigated
#' an interaction with year. 
#'
#' Type III ANOVA tables were computed with the car package
#' (v`r packageVersion("car")`; Fox and Weisberg 2019)
#' and used to assess polynomial terms. Where non-significant, 2nd order
#' polynomial terms were omitted from the models.
#'
#' All model fits and assumptions, and potential temporal and spatial autocorrelation
#' were assessed with the DHARMa package
#' (v`r packageVersion("DHARMa")`; Hartig 2020);
#' multicollinearity was assessed with the performance package
#' (v`r packageVersion("performance")`; Lüdecke et al., 2021);
#' and figures were created with the ggplot2 package
#' (v`r packageVersion("ggplot2")`; Wickham 2016). Note that figure scales
#' are log10 transformed after first adding 1 to better visualize patterns.
#'
#'
#' **References**
#'
#+ results = "asis"
list("base", "car", "glmmTMB", "emmeans", "DHARMa", "performance", "ggplot2") %>%
  map(~citation(.) %>%
        print(style = "text") %>%
        capture.output() %>%
        paste0(collapse = "\n")) %>%
  unlist() %>%
  sort() %>%
  cat(sep = "\n\n")


#+ echo = FALSE
# Reproducible ---------------------------------------------------------------
#' # Reproducible

#' ## Data sets
DT::datatable(mosq, extensions = 'Buttons',
              options = list(dom = 'Bfrtip', buttons = c('csv', 'excel')))

#' ## Session Info
#+ R.options = list(width = 100)
devtools::session_info()
