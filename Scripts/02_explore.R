#' ---
#' title: "Mosquitoes - Data Exploration"
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
# rmarkdown::render(input = "Scripts/02_explore.R",
#                   output_dir = "Results",
#                   output_file = paste0("02_explore_", Sys.Date(), '.html'),
#                   envir = new.env())

#+ echo = FALSE
# Setup -------------------------------------------------------------------
#' # Setup

# renv::restore() # If using these scripts to restore dependencies (do not knit)

#+ setup, include = FALSE
knitr::opts_chunk$set(fig.width = 10, fig.asp = 0.8, out.width = "100%")

#' ## Load packages
#+ message = FALSE
library(tidyverse)
library(lubridate)
library(here)
library(gt)
library(GGally)

#' ## Load data
#+ message = FALSE
mosq <- read_csv(here("Data/Datasets/mosquitoes.csv"), show_col_types = FALSE) %>%
  mutate(year = factor(year))


#+ echo = FALSE
# About the figures -----------------------------------------------------------
#' # About the figures
#'
#' **Log scale**
#'
#' We'll be using a negative binomial distribution (see Distribution below),
#' which uses a log based link function by default.
#'
#' Therefore, most of the figures are plotted on a log scale (after adding one
#' to the data to be able to deal with 0's) for illustration. But note that this
#' is NOT how they'll be analyzed!
#'
#' **Polynomials**
#'
#' Because there does seem to be a polynomial (curved) relationship with many
#' variables, I often add a polynomial (2) line to the data to illustrate
#' relationships. Remember that this is just exploration, and doesn't represent
#' a necessarily significant, AND that we haven't account for other parameters,
#' either.
#'
#'
#+ echo = FALSE
# Distribution -------------------------------------------------------------------
#' # Distribution
#'
#' - Definitely not normal
#' - Poisson / Negative Binomial
#' - Some mosquito papers use Gamma but only if averaged results

ggplot(data = mosq, aes(x = trapcount, fill = site)) +
  theme_bw() +
  geom_histogram(position = "dodge") +
  facet_grid(year ~ species, scales = "free_x") +
  scale_fill_viridis_d(end = 0.8)

#' Negative Binomial seems the best
#'
#'
#' ## Explore fit
#'
#' > **Using `fitdistrplus` package to explore fit**
#' >
#' > Note that this isn't definitive! Depends on the model parameters, etc.
#' > This is a good place to start

m <- mosq %>%
  select(species, trapcount) %>%
  drop_na()


#' Looks like negative binomial is the distribution that matches the data best.

#+ fig.asp = 1
par0 <- par(mfrow = c(2,2))
for(s in unique(mosq$species)) {
  fitdistrplus::descdist(m$trapcount[m$species == s], discrete = TRUE)
  title(sub = s)
}
par(par0)


#+ echo = FALSE
# Outliers -------------------------------------------------------------------
#' # Outliers
#' Let's look for outliers (look where counts > 500)
ggplot(data = filter(mosq, trapcount > 500), aes(x = trapcount, fill = site)) +
  theme_bw() +
  geom_histogram(position = "dodge") +
  facet_grid(year ~ species) +
  scale_fill_viridis_d(end = 0.8) +
  geom_vline(xintercept = c(1000, 5000), linetype = "dotted")


#' Trap counts greater than 1000 by site
filter(mosq, trapcount >= 1000) %>%
  count(site, species) %>%
  gt()


#' - There are a couple of instances where counts are much higher than 1000 or
#'   even 5000, but that's not common
#' - Cypress River seems to be an outlier site, especially with respect to
#'   Cq. perturbans
#'
#' > Consider omitting Cypress River if makes sense biologically speaking

#+ echo = FALSE
# Data parameters --------------------------------------------------------
#' # Data parameters
#'
#+ echo = FALSE
## Trapping events --------------------------------------------------------
#' ## Trapping events
#'
#' Number of trapping events per site
#'
#'  - not all sites present in all years
#'  - West Saint Paul 6 has only one event (NA, so will be omitted)
#'
count(mosq, site, year, species) %>%
  complete(site, year, species) %>%
  pivot_wider(names_from = year, values_from = n) %>%
  gt()

filter(mosq, site == "West Saint Paul 6") %>%
  select(date, site, species, trapcount) %>%
  gt()

events <- mosq %>%
  filter(!is.na(trapcount)) %>%
  count(site, year, species) %>%
  complete(site, year, species)

#+ fig.asp = 0.5
ggplot(events, aes(x = species, y = n, fill = site)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~year) +
  scale_fill_viridis_d() +
  labs(y = "No. Trapping events")

#' - Fewer trapping events in 2021
#' - Fewer trapping events for Cx. tarsalsi in 2020

#+ echo = FALSE
## Species --------------------------------------------------------
#' ## Species
#+ fig.asp = 0.5
ggplot(data = mosq, aes(x = species, y = trapcount, fill = species)) +
  theme_bw() +
  stat_summary(geom = "bar", fun = "mean") +
  scale_fill_viridis_d(end = 0.8) +
  facet_wrap(~ year) +
  labs(y = "Mean trapcounts", title = "Mean species counts")

ggplot(data = filter(mosq, site != "Cypress River"),
       aes(x = species, y = trapcount, fill = species)) +
  theme_bw() +
  stat_summary(geom = "bar", fun = "mean") +
  scale_fill_viridis_d(end = 0.8) +
  facet_wrap(~ year) +
  labs(y = "Mean trapcounts",
       title = "Mean species counts", subtitle = "Cypress River omitted")


#+ echo = FALSE
## Time --------------------------------------------------------
#' ## Time
#'
#' - Definitely a non-linear relationship

ggplot(data = mosq, aes(x = date, y = trapcount, colour = site)) +
  theme_bw() +
  geom_count() +
  facet_grid(species ~ year, scales = "free") +
  scale_colour_viridis_d(end = 0.8) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p")

#' Lets look at cumulative counts and see when species appear
cumulative_counts <- mosq %>%
  select(date, site, species, trapcount, year) %>%
  filter(!is.na(trapcount)) %>%
  group_by(year, species) %>%
  arrange(date) %>%
  mutate(cnts = cumsum(trapcount)) %>%
  ungroup()

cc_perc <- cumulative_counts %>%
  group_by(year, species) %>%
  summarize(p = c(0.1, 0.5, 0.9),
            q = quantile(cnts, prob = p), .groups = "drop") %>%
  left_join(cumulative_counts, by = c("year", "species")) %>%
  filter(cnts >= q) %>%
  group_by(year, species, p) %>%
  slice(1)


#' No really consistent patterns
ggplot(data = cumulative_counts, aes(x = date, y = cnts, colour = species)) +
  theme_bw() +
  geom_line(size = 1) +
  facet_grid(species ~ year, scales = "free_x") +
  scale_colour_viridis_d(end = 0.7) +
  geom_vline(data = cc_perc,
             aes(xintercept = date, linetype = factor(p))) +
  scale_y_continuous(trans = "log1p") +
  labs(linetype = "Percentile", y = "Cumulative trap counts")



#+ echo = FALSE
## Main parameters --------------------------------------------------------
#' ## Main parameters
#'
#+ echo = FALSE
### Temperature --------------------------------------------------------
#' ### Temperature
ggplot(data = mosq, aes(x = ddm14, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_count() +
  facet_grid(species ~ ., scales = "free") +
  scale_colour_viridis_c(end = 0.8)  +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p")

#+ echo = FALSE
### Precipitation --------------------------------------------------------
#' ### Precipitation
ggplot(data = mosq, aes(x = pt14, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_count() +
  facet_grid(species ~ ., scales = "free") +
  scale_colour_viridis_c(end = 0.8)  +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p")

#+ echo = FALSE
### Humidity --------------------------------------------------------
#' ### Humidity
#' - Some really high humidity days at the end of the year... humidity or date?
ggplot(data = mosq, aes(x = rha14, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_point() +
  facet_grid(species ~ ., scales = "free") +
  scale_colour_viridis_c(end = 0.8)  +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p")

#+ echo = FALSE
## Trapping conditions -------------------------------------------
#' ## Trapping conditions (covariates)
#'
#' - There was an option of using ttmin as a cutoff value, rather than a covariate
#' - However, looking here, ttmin seems to have a positive relationship with
#'   trapcount numbers, even above 4C (red line).

ggplot(data = mosq, aes(x = ttmin, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_count() +
  facet_grid(species ~ ., scales = "free") +
  scale_colour_viridis_c(end = 0.8) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  stat_smooth(data = filter(mosq, ttmin > 4),
              method = "lm", formula = y ~ x+ I(x^2), colour = "red", linetype = "dashed") +
  scale_y_continuous(trans = "log1p") +
  geom_vline(xintercept = 4, linetype = "dotted")


ggplot(data = mosq, aes(x = ttmax, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_point() +
  facet_grid(species ~ .) +
  scale_colour_viridis_c(end = 0.8)   +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p") +
  labs(caption = "Mostly no relationship (Except maybe Oc. dorsalis)")

ggplot(data = mosq, aes(x = windspd, y = trapcount, colour = yday(date))) +
  theme_bw() +
  geom_point() +
  facet_grid(species ~ .) +
  scale_colour_viridis_c(end = 0.8) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black") +
  scale_y_continuous(trans = "log1p") +
  labs(caption = "Mostly negative relationship (Except maybe Oc. dorsalis)")

#' Compare among covariates
#'
#' - We see that some are correlated (e.g., ttmin and ttmax and ttmin and windspeed)

cov <- select(mosq, date, ttmin, ttmax, windspd) %>%
  distinct() %>%
  mutate(date = yday(date))

#+ message = FALSE, warning = FALSE
ggpairs(cov, lower = list(continuous = "smooth",
                          mapping = aes(colour = date))) +
  scale_colour_viridis_c(end = 0.8)


#' Compare wind speed to main parameters
cov <- select(mosq, date, ddm14, rha14, pt14, windspd) %>%
  mutate(date = yday(date)) %>%
  distinct()

#+ message = FALSE, warning = FALSE
ggpairs(cov, lower = list(continuous = "smooth",
                          mapping = aes(colour = date))) +
  scale_colour_viridis_c(end = 0.8)



#+ echo = FALSE
# Reproducible ---------------------------------------------------------------
#' # Reproducible

#' ## Data sets
DT::datatable(mosq, extensions = 'Buttons',
              options = list(dom = 'Bfrtip', buttons = c('csv', 'excel')))

#' ## Session Info
#+ R.options = list(width = 100)
devtools::session_info()
